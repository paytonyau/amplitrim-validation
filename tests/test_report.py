"""REPORT.md, the figure, the tables and the screen summary, built from a hand-made result."""
import json
import re
from types import SimpleNamespace

from amplitrim.config import Config
from amplitrim.pipeline import Result, _warnings
from amplitrim.ranking import Candidate, select_recommendation
from amplitrim.report import _share_shorter, screen_summary, write_outputs


def candidate(trim_f, trim_r, yield_est, overlap):
    return Candidate(trim_f=trim_f, trim_r=trim_r, ee_pass_rate=1.0, merge_rate_strict=0.8,
                     merge_rate_tolerant=yield_est, mean_overlap=overlap, mean_overlap_diffs=0.2,
                     max_mergeable_insert=trim_f + trim_r - 12, covers_longest=True, yield_est=yield_est)


def make_result(recommend=True, warnings=(), uniform=True):
    sub = SimpleNamespace(n_pairs=12000, n_input_pairs=12000, stride=1, read_len_f=231, read_len_r=227,
                          uniform_length=uniform, min_len_f=200, min_len_r=200, quality_values=40)
    cands = [candidate(222, 226, 0.995, 18.0), candidate(224, 226, 0.998, 20.0), candidate(229, 225, 1.0, 24.0)]
    res = Result(Config(), sub, [422] * 5 + [430] * 10 + [434] * 20, 1.0, 434, 446, 458,
                 totals=[448, 450, 454], candidates=cands, warnings=list(warnings))
    if recommend:
        res.recommendation = select_recommendation(cands, 0.01)
    return res


def test_screen_summary_is_short_and_gives_the_settings(tmp_path):
    lines = screen_summary(make_result(), tmp_path / "REPORT.md")
    assert lines == ["", "Recommended DADA2 settings",
                     "  truncLenF = 222   truncLenR = 226   maxEE = 2   minOverlap = 12", "",
                     "  99.5% of the 12,000 read pairs analysed should pass filtering and merge.",
                     "  Covers amplicons up to 436 nt (your longest: 434 nt).", "",
                     f"Report: {tmp_path / 'REPORT.md'}"]


def test_screen_summary_lists_warnings_and_handles_no_recommendation(tmp_path):
    text = "\n".join(screen_summary(make_result(warnings=["Check this."]), tmp_path / "R.md"))
    assert "Warning: Check this." in text
    none = "\n".join(screen_summary(make_result(recommend=False), tmp_path / "R.md"))
    assert "No recommendation" in none and "truncLenF" not in none


def test_report_leads_with_settings_snippets_and_what_to_expect(tmp_path):
    text = write_outputs(make_result(), tmp_path).read_text()
    assert text.startswith("# AmpliTrim report\n")
    assert text.index("## Recommended settings") < text.index("## What to expect") < text.index("## Limits")
    assert "| **222** | **226** | 2 | 12 |" in text
    assert "truncLen = c(222, 226)" in text and "--p-trunc-len-f 222 --p-trunc-len-r 226" in text
    assert "minOverlap = 12, maxMismatch = 0" in text
    assert "About **99.5%**" in text and "up to 436 nt" in text
    assert "## Things to check" not in text                         # only shown when there is something to check
    assert "adapters and primers removed" in text and "zero-mismatch" in text
    assert "before DADA2 denoises them" in text and "not as the final retained fraction" in text
    assert not re.search(r"\bnan\b|Stage [0-9]|FIGARO score|16S|V3", text, re.IGNORECASE)


def test_report_lists_warnings_when_there_are_any(tmp_path):
    text = write_outputs(make_result(warnings=["Only 40% of pairs overlap."]), tmp_path).read_text()
    assert "## Things to check" in text and "- Only 40% of pairs overlap." in text


def test_report_without_a_recommendation_says_so(tmp_path):
    text = write_outputs(make_result(recommend=False), tmp_path).read_text()
    assert "## No recommendation" in text and "truncLen = c(" not in text


def test_outputs_are_the_four_files_with_tidy_contents(tmp_path):
    write_outputs(make_result(), tmp_path)
    assert {p.name for p in tmp_path.iterdir()} == {"REPORT.md", "summary.png", "candidates.csv", "summary.json"}
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["recommended"] == {"truncLenF": 222, "truncLenR": 226, "maxEE": 2.0, "minOverlap": 12,
                                      "maxMismatch": 0, "estimated_yield": 0.995, "longest_mergeable_insert": 436}
    assert summary["longest_amplicon"] == 434 and summary["warnings"] == []
    rows = (tmp_path / "candidates.csv").read_text().splitlines()
    assert len(rows) == 4 and rows[0].startswith("truncLenF,truncLenR,total_length")
    assert (tmp_path / "summary.png").stat().st_size > 5000


def test_outputs_are_written_even_without_a_recommendation(tmp_path):
    write_outputs(make_result(recommend=False), tmp_path)
    assert json.loads((tmp_path / "summary.json").read_text())["recommended"] is None


def test_warns_when_figaros_split_passes_clearly_fewer_pairs_than_the_best_split():
    res = make_result()
    res.split_checks = [(448, 0.95, 0.99, (225, 223))]
    assert any("FIGARO's suggested splits of 448 nt pass 95.0%" in w and "99.0% for 225 + 223" in w
               for w in _warnings(res, False))
    res.split_checks = [(448, 0.985, 0.99, (225, 223))]            # within 1 point: nothing to say
    assert not any("FIGARO's suggested" in w for w in _warnings(res, False))


def test_warns_how_many_pairs_are_lost_when_read_lengths_vary():
    res = make_result(uniform=False)
    res.short_lost = 0.10
    assert any("Read lengths vary (forward 200-231 nt, reverse 200-227 nt)" in w and "loses 10%" in w
               for w in _warnings(res, False))
    res.short_lost = 0.001
    assert not any("Read lengths vary" in w for w in _warnings(res, False))


def test_low_merge_rate_warning_names_shorter_reads_among_the_possible_causes():
    res = make_result()
    res.frac_overlapping = 0.75
    text = " ".join(_warnings(res, False))
    assert "Only 75% of read pairs could be merged at full read length" in text
    assert "some reads may be shorter than the rest" in text


def test_report_has_a_plain_your_data_section_and_no_read_through_note_when_amplicons_are_longer(tmp_path):
    text = write_outputs(make_result(), tmp_path).read_text()
    assert text.index("## Your data") < text.index("![summary]")
    assert "- Reads: 12,000 pairs analysed (231 + 227 nt)." in text
    assert "- Amplicons: 422-434 nt (median 434)" in text and "the 100% of pairs that merge" in text
    assert "- Overlap at the recommended settings: about 18 nt on average." in text
    assert "shorter than the reads" not in text


def test_report_notes_when_amplicons_are_shorter_than_the_reads(tmp_path):
    res = make_result()
    res.inserts = [200] * 10 + [434] * 90                          # 10% shorter than the 231 nt reads
    assert _share_shorter(res) == 0.1
    text = write_outputs(res, tmp_path).read_text()
    assert "10% of amplicons are shorter than the reads" in text and "remove them first" in text
    assert "## Things to check" not in text                         # a plain note, not a warning


def test_amplicon_line_shows_where_99_percent_fall_when_a_few_outliers_stretch_the_range(tmp_path):
    res = make_result()
    res.inserts = sorted([58] + [253] * 998 + [369])
    text = write_outputs(res, tmp_path).read_text()
    assert "median 253 nt; 99% are 253-253 nt (all: 58-369)" in text


def test_summary_json_describes_the_dataset_for_comparing_studies(tmp_path):
    write_outputs(make_result(), tmp_path)
    data = json.loads((tmp_path / "summary.json").read_text())["data"]
    assert data == {"pairs_analysed": 12000, "read_length_f": 231, "read_length_r": 227, "uniform_read_length": True,
                    "amplicon_min": 422, "amplicon_median": 434, "amplicon_max": 434,
                    "amplicon_99pct_min": 422, "amplicon_99pct_max": 434,
                    "share_merged_at_full_length": 1.0, "share_amplicons_shorter_than_reads": 0.0,
                    "quality_values": 40}


def make_its_result():
    res = make_result()
    cand = candidate(0, 0, 0.9, 0.0)
    cand.ee_pass_rate, cand.merge_rate_tolerant, cand.yield_est, cand.mean_overlap = 0.9, 1.0, 0.9, 60.0
    cand.max_mergeable_insert, cand.source = 480, "no truncation"
    res.marker, res.min_len, res.candidates, res.recommendation = "its", 50, [cand], cand
    res.cfg, res.totals = Config(), []
    return res


def test_its_report_says_do_not_truncate_and_gives_the_dada2_settings(tmp_path):
    text = write_outputs(make_its_result(), tmp_path).read_text()
    assert "Do not truncate." in text and "| **0 (none)** | **0 (none)** | 2 | 50 | 12 |" in text
    assert "truncQ = 2, minLen = 50" in text and "--p-trunc-len-f 0 --p-trunc-len-r 0" in text
    assert "90.0% pass the filter, and 100.0% of those merge" in text
    assert "## How these settings were chosen" in text and "No truncation is tested" in text
    assert "total length" not in text.lower()                        # nothing about searching truncation lengths
    summary = json.loads((tmp_path / "summary.json").read_text())["recommended"]
    assert (summary["truncLenF"], summary["truncLenR"], summary["minLen"]) == (0, 0, 50)


def test_its_screen_summary_is_short(tmp_path):
    lines = screen_summary(make_its_result(), tmp_path / "R.md")
    assert lines[1] == "Recommended DADA2 settings"
    assert lines[2] == "  no truncation (truncLen = 0)   maxEE = 2   minLen = 50   minOverlap = 12"
    assert "90.0% of the 12,000 read pairs analysed" in lines[4] and len(lines) <= 7


def test_a_few_short_outliers_do_not_trigger_the_shorter_than_reads_note(tmp_path):
    res = make_result()
    res.inserts = [200] * 2 + [434] * 998                           # 0.2% shorter than the 231 nt reads
    assert "shorter than the reads" not in write_outputs(res, tmp_path).read_text()


def test_two_quality_values_give_a_clear_warning_and_binned_scores_a_plain_note(tmp_path):
    flat = make_result()
    flat.sub.quality_values = 2
    assert any("Quality scores take only 2 distinct values" in w and "original FASTQ" in w for w in _warnings(flat, False))
    binned = make_result()
    binned.sub.quality_values = 4
    assert not any("Quality scores" in w for w in _warnings(binned, False))          # a note, not a warning
    assert "Quality scores take only 4 distinct values" in write_outputs(binned, tmp_path).read_text()
    assert not any("Quality scores" in w for w in _warnings(make_result(), False))   # normal data: nothing


def test_the_shorter_than_reads_note_needs_reads_of_one_length(tmp_path):
    """After read-through trimming every read is as long as its own insert, so the note would mislead."""
    res = make_result(uniform=False)
    res.inserts = [200] * 10 + [434] * 90
    assert "shorter than the reads" not in write_outputs(res, tmp_path).read_text()


def test_the_its_report_also_warns_that_the_yield_ignores_denoising(tmp_path):
    assert "before DADA2 denoises them" in write_outputs(make_its_result(), tmp_path).read_text()


def test_two_value_maxEE_appears_in_the_dada2_and_qiime_snippets_and_the_json(tmp_path):
    import json
    res = make_result()
    res.cfg.max_ee = (2.0, 6.0)
    text = write_outputs(res, tmp_path).read_text()
    assert "maxEE = c(2, 6)" in text and "--p-max-ee-f 2" in text and "--p-max-ee-r 6" in text
    summary = json.loads((tmp_path / "summary.json").read_text())
    assert summary["recommended"]["maxEE"] == [2.0, 6.0]
    res.cfg.max_ee = 2.0
    text = write_outputs(res, tmp_path).read_text()
    assert "maxEE = c(2, 2)" in text and "--p-max-ee-r 2" in text


def _rows(*triples):
    return [{"ee": ee, "pass_rate": p, "yield_est": p - 0.01, "score": 0.0, "current": False, "ee_f": 0.4, "ee_r": 1.1} for ee, p in triples]


def test_maxEE_advice_states_dada2_authors_guidance_and_never_asks_for_a_change():
    from amplitrim.report import ee_advice
    healthy = " ".join(ee_advice(_rows(((2, 2), 0.95), ((2, 5), 0.98)), 2.0))
    assert "95.0%" in healthy and "general-purpose value" in healthy and "number of reads kept is not what matters" in healthy
    assert "could make them worse" in healthy and "tutorial" not in healthy and "--max-ee" not in healthy
    low = " ".join(ee_advice(_rows(((2, 2), 0.40), ((2, 5), 0.80)), 2.0))
    assert "40.0%" in low and "not a problem in itself" in low and "maxEE = c(2, 5)" in low and "--max-ee" not in low
    assert ee_advice(_rows(((1, 1), 0.5)), 2.0) == []


def test_the_maxEE_table_is_exact_and_consistent(cli_run):
    import json
    from amplitrim.trimmer import ReadPairs
    summary = json.loads((cli_run.out / "summary.json").read_text())
    rec, table = summary["recommended"], summary["maxEE_table"]
    rows = {tuple(x["maxEE"]): x for x in table}
    assert rows[(2.0, 2.0)]["estimated_yield"] == rec["estimated_yield"]                 # the row in use is the recommendation itself
    pairs = ReadPairs.load(cli_run.work / "subsample" / "subsample_R1.fastq.gz", cli_run.work / "subsample" / "subsample_R2.fastq.gz")
    for ee, x in rows.items():
        assert abs(pairs.pass_table(ee).rate(rec["truncLenF"], rec["truncLenR"]) - x["pass_rate"]) < 1e-3
        assert x["estimated_yield"] <= x["pass_rate"] + 1e-9
    assert rows[(1.0, 1.0)]["pass_rate"] <= rows[(2.0, 2.0)]["pass_rate"] <= rows[(2.0, 3.0)]["pass_rate"] <= rows[(2.0, 5.0)]["pass_rate"]
    report = (cli_run.out / "REPORT.md").read_text()
    assert "## What maxEE keeps" in report and "**2 (yours)**" in report and "FIGARO's score" in report
    load = {k: v["expected_errors_kept"] for k, v in rows.items()}                  # looser limits keep noisier reads
    assert load[(2.0, 2.0)][1] <= load[(2.0, 3.0)][1] <= load[(2.0, 4.0)][1] <= load[(2.0, 5.0)][1] and load[(1.0, 1.0)][1] <= load[(2.0, 2.0)][1]


def test_the_maxEE_in_use_is_marked_when_it_is_a_pair(cli_run_ee25):
    report = (cli_run_ee25.out / "REPORT.md").read_text()
    assert "**2,5 (yours)**" in report and "maxEE = c(2, 5)" in report


def test_the_per_maxEE_search_is_consistent_with_the_main_recommendation(cli_run_explore):
    import json
    summary = json.loads((cli_run_explore.out / "summary.json").read_text())
    expl, rec = summary["maxEE_exploration"], summary["recommended"]
    sug = {"rows": expl["rows"], "estimated_yield": expl["suggested"]["estimated_yield"], "suggested": expl["suggested"]["maxEE"]}
    rows = {tuple(x["maxEE"]): x for x in sug["rows"]}
    mine = rows[(2.0, 2.0)]
    assert (mine["truncLenF"], mine["truncLenR"], mine["estimated_yield"]) == (rec["truncLenF"], rec["truncLenR"], rec["estimated_yield"])
    assert tuple(sug["suggested"]) in rows and tuple(sug["suggested"]) >= (2.0, 2.0)               # never stricter than the default
    assert expl["minRetention"] == 0.5 and (expl["suggested"]["estimated_yield"] >= 0.5 or not expl["suggested"]["reachesFloor"])
    assert all(x["truncLenF"] + x["truncLenR"] >= rec["longest_mergeable_insert"] + 1 for x in sug["rows"] if x["truncLenF"])   # every row covers the amplicons
    for x in sug["rows"]:                                                           # every "beaten by" claim is true on both counts
        for other in x["beatenBy"]:
            o = rows[tuple(other)]
            assert o["pass_rate"] >= x["pass_rate"] - 1e-9
    report = (cli_run_explore.out / "REPORT.md").read_text()
    assert "## maxEE: retention against read quality" in report and "beaten by" in report and "Quality first" in report and "**2 (yours)**" in report and "never changes your maxEE" in report
    assert "general-purpose value" in report and "expected errors per kept read" in report
    assert "## What maxEE keeps" not in report                                     # the richer table replaces the fixed-truncation one


def test_coverage_rows_keep_their_promise_and_the_report_shows_them(cli_run_coverage):
    import json
    summary = json.loads((cli_run_coverage.out / "summary.json").read_text())
    rows, rec = summary["coverage_exploration"]["rows"], summary["recommended"]
    by_level = {x["allowUnmergeable"]: x for x in rows}
    assert 0.002 in by_level and 0.01 in by_level and by_level[0.002]["truncLenF"] == rec["truncLenF"] and by_level[0.002]["estimated_yield"] == rec["estimated_yield"]       # the default, and the old 1% beside it
    for x in rows:
        if x["totalLength"] is not None:
            assert x["leftOut"] <= x["allowUnmergeable"] + 1e-9                 # no row leaves out more amplicons than it says
            assert x["truncLenF"] + x["truncLenR"] == x["totalLength"]
    report = (cli_run_coverage.out / "REPORT.md").read_text()
    assert "## Coverage against quality" in report and "**0.2% (yours)**" in report and "1.0%" in report
    import matplotlib.image as mpimg
    img = mpimg.imread(cli_run_coverage.out / "summary.png")
    assert img.shape[1] / img.shape[0] > 3.5                                    # three panels side by side


def test_allow_unmergeable_is_validated_and_limits_the_share_of_amplicons_left_out(reads, tmp_path):
    import pytest
    from amplitrim.pipeline import PipelineError, Settings, run
    from amplitrim.ranking import frac_unmergeable
    for share in (0.0, 0.05):
        res = run(Settings(input_dir=reads, work_dir=tmp_path / f"w{share}", max_unmergeable=share))
        assert frac_unmergeable(res.inserts, res.recommendation.max_mergeable_insert) <= share + 1e-9
    with pytest.raises(PipelineError):
        run(Settings(input_dir=reads, work_dir=tmp_path / "bad", max_unmergeable=0.7))


def test_the_report_shows_how_noisy_the_extra_reads_are(cli_run):
    import json
    summary = json.loads((cli_run.out / "summary.json").read_text())
    rows = {tuple(x["maxEE"]): x for x in summary["maxEE_table"]}
    assert rows[(2.0, 2.0)]["expected_errors_added"] is None and rows[(1.0, 1.0)]["expected_errors_added"] is None
    for ee in ((2.0, 3.0), (2.0, 4.0), (2.0, 5.0)):
        added, kept = rows[ee]["expected_errors_added"], rows[ee]["expected_errors_kept"]
        if added:
            assert added[1] > 2.0 and added[1] >= kept[1] - 1e-9          # the extra reads are noisier than the average kept read
    assert "expected errors in the reads added beyond maxEE 2" in (cli_run.out / "REPORT.md").read_text()


def test_min_retention_moves_the_quality_first_suggestion_but_never_the_default(reads, tmp_path):
    import pytest
    from amplitrim.pipeline import PipelineError, Settings, run
    from amplitrim.ranking import quality_first_suggestion
    lenient = run(Settings(input_dir=reads, work_dir=tmp_path / "a", explore_max_ee=True, min_retention=0.2))
    strict = run(Settings(input_dir=reads, work_dir=tmp_path / "b", explore_max_ee=True, min_retention=0.999))
    pick_l, ok_l = quality_first_suggestion(lenient.ee_search, 0.2)
    pick_s, ok_s = quality_first_suggestion(strict.ee_search, 0.999)
    assert ok_l and pick_l["ee"] == (2.0, 2.0)                             # a low floor keeps the default
    assert pick_s is not None and min(pick_s["ee"]) >= 2.0                 # an unreachable floor still never goes below 2
    assert lenient.recommendation.max_ee_unchanged if hasattr(lenient.recommendation, "max_ee_unchanged") else lenient.cfg.max_ee == 2.0
    with pytest.raises(PipelineError):
        run(Settings(input_dir=reads, work_dir=tmp_path / "bad", min_retention=0))


def test_the_its_report_gives_the_maxEE_table_with_length_and_noise_and_a_quality_first_line(reads, tmp_path):
    import json
    from amplitrim.pipeline import Settings, run
    from amplitrim.report import write_outputs
    from amplitrim.pipeline import EE_GRID
    res = run(Settings(input_dir=reads, work_dir=tmp_path / "w", marker="its", min_len=60))
    rows = res.ee_search
    assert len(rows) == len({tuple(float(x) for x in (e if isinstance(e, tuple) else (e, e))) for e in EE_GRID} | {(2.0, 2.0)})
    by = {tuple(x["ee"]): x for x in rows}
    assert by[(1.0, 1.0)]["pass_rate"] <= by[(2.0, 2.0)]["pass_rate"] <= by[(2.0, 5.0)]["pass_rate"] <= by[(5.0, 5.0)]["pass_rate"]       # a looser limit never keeps fewer pairs
    assert by[(1.0, 1.0)]["len_f"] <= by[(5.0, 5.0)]["len_f"] + 1e-9                         # and the kept reads are no shorter: long reads are the first to fail
    now = next(x for x in rows if x["current"])
    assert now["ee"] == (2.0, 2.0) and abs(now["yield_est"] - res.recommendation.yield_est) < 1e-9
    assert by[(2.0, 5.0)]["add_r"] > 2.0 and by[(2.0, 2.0)]["add_f"] != by[(2.0, 2.0)]["add_f"]    # the reads added beyond maxEE 2 exceed it; the baseline adds none
    write_outputs(res, tmp_path / "o")
    out = tmp_path / "o"
    text = (out / "REPORT.md").read_text(encoding="utf-8")
    assert "## maxEE: retention against read quality" in text and "mean length kept (F / R)" in text and "Quality first" in text and "upper bound" in text
    assert text.index("## Things to check") < text.index("## maxEE: retention against read quality") < text.index("## Your data") if "## Things to check" in text else True
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    expl = summary["maxEE_exploration"]
    assert expl["suggested"]["maxEE"] == [2.0, 2.0] and expl["minRetention"] == 0.5 and all(r["meanLengthKept"] is not None for r in expl["rows"] if r["estimated_yield"])


def test_the_default_share_of_longest_amplicons_left_out_is_point_two_per_cent_everywhere(tmp_path):
    """The default was 1% until a rule fixed in advance (TEST_REGISTER.md, section 6) was met by the validation; the three places that define it must agree."""
    import inspect
    from amplitrim import pipeline
    from amplitrim.cli import build_parser
    from amplitrim.config import Config
    args = build_parser().parse_args(["--input-dir", str(tmp_path), "--work-dir", str(tmp_path / "w"), "--out-dir", str(tmp_path / "o")])
    assert args.allow_unmergeable == 0.002 and Config().max_unmergeable_frac == 0.002
    assert inspect.signature(pipeline.Settings).parameters["max_unmergeable"].default == 0.002
    assert 0.002 in pipeline.COVERAGE_LEVELS and 0.01 in pipeline.COVERAGE_LEVELS            # the report's coverage table shows the new default and the old one
