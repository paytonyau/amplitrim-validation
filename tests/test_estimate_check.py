"""validation/estimate_check.py: strict and tolerant estimates against DADA2 (part A, tables), and the vsearch cut-off for differences against DADA2's merge step (part B). vsearch is a stand-in."""
import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "validation"))


def load():
    spec = importlib.util.spec_from_file_location("estimate_check", HERE.parent / "validation" / "estimate_check.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


COLS = ["name", "truncF", "truncR", "maxEE", "est_yield", "est_pass", "est_strict", "pass", "merged", "error"]


def write_compare(runs: Path, run: str, dataset: str, rows: list, tracks: dict = None):
    from compare_methods import method_key
    folder = runs / run / dataset / "compare"
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / "compare_methods.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in COLS})
    for name, lines in (tracks or {}).items():
        with open(folder / f"track_{method_key(name)}.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "input", "filtered", "denoised_f", "denoised_r", "merged", "nonchim"])
            w.writerows(lines)
    return folder


def row(name, tf, tr, ee, est, strict, merged, **kw):
    return dict({"name": name, "truncF": tf, "truncR": tr, "maxEE": ee, "est_yield": est, "est_strict": strict, "est_pass": 0.9, "pass": 0.9, "merged": merged}, **kw)


def test_fractions_are_read_whatever_the_format():
    ec = load()
    assert ec.to_fraction("0.469") == pytest.approx(0.469) and ec.to_fraction("46.9%") == pytest.approx(0.469) and ec.to_fraction(46.9) == pytest.approx(0.469) and ec.to_fraction("1") == 1.0
    assert ec.to_fraction("") is None and ec.to_fraction("nan") is None and ec.to_fraction(None) is None
    assert ec.parse_ee("2") == (2.0, 2.0) and ec.parse_ee("2,5") == (2.0, 5.0) and ec.ee_text((2.0, 2.0)) == "2" and ec.ee_text((2.0, 5.0)) == "2,5"


def test_the_tables_on_disk_are_read_and_a_setting_run_on_several_nights_counts_once(tmp_path):
    ec = load()
    write_compare(tmp_path, "night1_2026-10-07", "v4", [row("AmpliTrim (maxEE 2)", 124, 115, "2", 0.90, 0.85, 0.60), row("FIGARO", 118, 121, "2.0", 0.88, 0.80, 0.58, error="boom"),
                                                          row("no truncation", 0, 0, "2", None, None, "", )])
    write_compare(tmp_path, "night3_2026-10-09", "v4", [row("AmpliTrim (maxEE 2)", 124, 115, "2.0", 0.91, 0.86, "62.6%"), row("AmpliTrim (maxEE 2,5)", 124, 115, "2,5", 0.95, 0.90, 0.70)])
    rows = ec.load_rows(tmp_path)
    assert len(rows) == 3                                                                                   # the error row and the row without a DADA2 result are left out
    d = ec.dedupe(rows)
    assert [(r["name"], r["maxEE"]) for r in d] == [("AmpliTrim (maxEE 2)", "2"), ("AmpliTrim (maxEE 2,5)", "2,5")]
    newest = next(r for r in d if r["maxEE"] == "2")
    assert newest["run"] == "night3_2026-10-09" and newest["merged"] == pytest.approx(0.626) and newest["est_yield"] == pytest.approx(0.91)      # the newest run of the same setting


def test_the_rank_correlation_is_the_spearman_of_the_ranks():
    ec = load()
    assert ec.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0) and ec.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    assert ec.spearman([1, 1, 2, 3], [1, 2, 3, 4]) == pytest.approx(4.5 / 22.5 ** 0.5)                      # ties share the mean rank: 4.5 / sqrt(4.5 * 5.0)
    assert ec.spearman([1, 2], [1, 2]) is None and ec.spearman([1, 1, 1], [1, 2, 3]) is None                  # too few, or no variation


def test_the_overshoot_and_the_ranking_are_computed_per_dataset():
    ec = load()
    def r(ds, est, strict, merged): return {"dataset": ds, "est_yield": est, "est_strict": strict, "merged": merged}
    rows = [r("X", 0.70, 0.60, 0.50), r("X", 0.60, 0.55, 0.45), r("X", 0.50, 0.45, 0.40), r("X", 0.40, 0.35, 0.30), r("Y", 0.80, 0.70, 0.75), r("Y", 0.50, 0.60, 0.55), r("Z", None, 0.5, 0.4),
            r("W", 0.70, 0.40, 0.50), r("W", 0.60, 0.55, 0.45), r("W", 0.50, 0.60, 0.40), r("W", 0.40, 0.30, 0.30)]
    s = ec.dataset_stats(rows)
    assert set(s) == {"X", "Y", "W"} and s["X"]["n"] == 4                                                    # Z has no tolerant estimate
    assert s["W"]["rho_tol"] == pytest.approx(1.0) and s["W"]["rho_strict"] == pytest.approx(0.2)          # the two estimates rank W differently: strict ranks [2, 3, 4, 1] against DADA2's [4, 3, 2, 1], 1 - 6 * 8 / 60
    assert s["X"]["tol"] == pytest.approx(13.75) and s["X"]["strict"] == pytest.approx(7.5) and s["X"]["tol_abs"] == pytest.approx(13.75) and s["X"]["strict_abs"] == pytest.approx(7.5)
    assert s["X"]["rho_tol"] == pytest.approx(1.0) and s["X"]["rho_strict"] == pytest.approx(1.0)
    assert s["Y"]["n"] == 2 and s["Y"]["rho_tol"] is None                                                    # two settings do not make a ranking
    assert s["Y"]["tol"] == pytest.approx((5 - 5) / 2) and s["Y"]["tol_abs"] == pytest.approx(5.0) and s["Y"]["strict_abs"] == pytest.approx(5.0)


def test_the_verdicts_of_part_a_follow_the_registered_criteria():
    ec = load()
    def st(tol_abs, strict_abs, rt=None, rs=None): return {"n": 4, "tol": tol_abs, "strict": strict_abs, "tol_abs": tol_abs, "strict_abs": strict_abs, "rho_tol": rt, "rho_strict": rs}
    v = {p: ok for p, _, ok in ec.verdicts_a({"a": st(10, 4, 0.85, 0.81), "b": st(3, 5, 0.85, 0.81)})}
    assert v == {"P36": True, "P37": True}                                                                   # 7.5 against 6.5 points; 0.81 is within 0.05 of 0.85
    v = {p: ok for p, _, ok in ec.verdicts_a({"a": st(5, 5, 0.85, 0.79)})}
    assert v == {"P36": False, "P37": False}                                                                 # equal is not closer; 0.79 is more than 0.05 below 0.85
    assert {p: ok for p, _, ok in ec.verdicts_a({"a": st(5, 4)})}["P37"] is None and all(ok is None for _, _, ok in ec.verdicts_a({}))


def test_dadas_merge_step_is_read_from_its_track(tmp_path):
    ec = load()
    p = tmp_path / "track.tsv"
    p.write_text("sample\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\ns1\t120\t100\t80\t70\t56\t55\ns2\t220\t200\t150\t160\t120\t118\n")
    d = ec.dada2_merge_step(p)
    assert d["step"] == pytest.approx(176 / 220)                                                             # merged over the sum of min(denoised forward, denoised reverse): 70 + 150
    assert d["merged"] == pytest.approx(176 / 340) and d["filtered"] == pytest.approx(300 / 340) and d["denoise_loss"] == pytest.approx(80 / 300)


@pytest.fixture(scope="module")
def ec_sim(tmp_path_factory):
    from simulate import write_dataset
    from amplitrim.subsample import read_fastq
    from amplitrim.trimmer import ReadPairs
    d = tmp_path_factory.mktemp("ec")
    write_dataset(d, n_samples=2, pairs_per_sample=300, seed=5)
    r1, r2 = [], []
    for s in ("sampleA", "sampleB"):
        r1 += list(read_fastq(d / "primer_trimmed" / f"{s}_R1.fastq.gz"))
        r2 += list(read_fastq(d / "primer_trimmed" / f"{s}_R2.fastq.gz"))
    for name, recs in (("subsample_R1.fastq.gz", r1), ("subsample_R2.fastq.gz", r2)):
        with gzip.open(d / name, "wt", newline="\n") as fh:
            fh.write("".join(f"{h}\n{s}\n{p}\n{q}\n" for h, s, p, q in recs))
    return SimpleNamespace(pairs=ReadPairs(r1, r2), dir=d)


def make_merge(calls):
    """vsearch as a function of the allowed differences only: it merges 20%, 40%, 60% and 80% of the pairs with 0, 1, 2 and 3 differences allowed, and all of them with more."""
    share = {0: 0.2, 1: 0.4, 2: 0.6, 3: 0.8}

    def merge(r1, r2, out_dir, tag, min_overlap, max_diffs, setting, keep_merged=True, threads=1):
        n = sum(1 for _ in open(r1)) // 4
        calls.append({"n": n, "min_overlap": min_overlap, "max_diffs": max_diffs, "setting": setting, "keep": keep_merged})
        m = round(n * share.get(max_diffs, 1.0))
        return SimpleNamespace(total_pairs=n, merged_pairs=m, merge_rate=m / n if n else 0.0)
    return merge


def test_vsearch_is_run_at_each_cut_off_on_the_filtered_pairs_only(ec_sim, tmp_path):
    ec = load()
    calls = []
    rates = ec.vsearch_rates(ec_sim.pairs, (225, 220, (2.0, 2.0)), make_merge(calls), tmp_path, diffs=(0, 1, 3, 10), min_overlap=12)
    keep = ec_sim.pairs._masks(225, 220, (2.0, 2.0))[4]
    assert [c["max_diffs"] for c in calls] == [0, 1, 3, 10] and all(c["n"] == int(keep.sum()) for c in calls)    # the filtered pairs, not all of them
    assert all(c["min_overlap"] == 12 and not c["keep"] for c in calls) and [c["setting"] for c in calls] == ["maxdiffs0", "maxdiffs1", "maxdiffs3", "maxdiffs10"]
    assert rates == {0: pytest.approx(0.2, abs=0.01), 1: pytest.approx(0.4, abs=0.01), 3: pytest.approx(0.8, abs=0.01), 10: 1.0} and 0 < int(keep.sum()) < ec_sim.pairs.n
    calls.clear()
    assert ec.vsearch_rates(ec_sim.pairs, (900, 900, (2.0, 2.0)), make_merge(calls), tmp_path / "b", diffs=(0, 10)) == {0: 0.0, 10: 0.0} and not calls          # nothing passes: nothing is merged
    calls.clear()
    ec.vsearch_rates(ec_sim.pairs, (0, 0, (2.0, 2.0)), make_merge(calls), tmp_path / "c", diffs=(0,))                    # ITS: no truncation, the no-truncLen filter
    _, _, short, has_n, ee_fail = ec_sim.pairs.untruncated((2.0, 2.0), 50)
    assert len(calls) == 1 and calls[0]["n"] == int((~(short | has_n | ee_fail)).sum()) > 0


def test_the_cut_off_summary_is_the_mean_difference_to_dadas_merge_step():
    ec = load()
    arms = [{"rates": {0: 0.70, 3: 0.82, 10: 0.99}, "dada2": {"step": 0.80}}, {"rates": {0: 0.88, 3: 0.90, 10: 0.99}, "dada2": {"step": 0.90}}]
    s = ec.cutoff_summary(arms, diffs=(0, 3, 10))
    assert s[0]["signed"] == pytest.approx(-6.0) and s[0]["abs"] == pytest.approx(6.0)
    assert s[3]["signed"] == pytest.approx(1.0) and s[3]["abs"] == pytest.approx(1.0)
    assert s[10]["signed"] == pytest.approx(14.0) and s[10]["abs"] == pytest.approx(14.0) and s[3]["n"] == 2
    def summ(a0, a3, a10): return {0: {"n": 2, "abs": a0, "signed": 0}, 3: {"n": 2, "abs": a3, "signed": 0}, 10: {"n": 2, "abs": a10, "signed": 0}}
    p, txt, ok = ec.verdict_b(summ(6.0, 1.0, 14.0))
    assert ok is True and "the closest is 3" in txt                                                          # 6.0 - 1.0 = 5.0 points closer than the better extreme
    assert ec.verdict_b(summ(0.9, 1.0, 14.0))[2] is False                                                   # the best cut-off is an extreme
    assert ec.verdict_b(summ(6.0, 4.1, 14.0))[2] is False                                                   # 1.9 points closer: under the 2 points set in advance
    assert ec.verdict_b(summ(6.0, 4.0, 14.0))[2] is True                                                    # exactly 2.0
    assert ec.verdict_b({0: {"n": 2, "abs": 1.0, "signed": 0}, 10: {"n": 2, "abs": 2.0, "signed": 0}})[2] is None


def setup_runs(tmp_path, ec_sim):
    runs = tmp_path / "runs"
    reads = tmp_path / "reads"
    reads.mkdir()
    cfg = tmp_path / "paths.json"
    cfg.write_text(json.dumps({"extras": {"mockB": {"ready": str(reads).replace("\\", "/")}}}))
    tr = [("s1", 120, 100, 80, 70, 56, 55), ("s2", 220, 200, 150, 160, 120, 118)]
    rows = [row("AmpliTrim (maxEE 2)", 225, 220, "2", 0.70, 0.62, 0.55), row("AmpliTrim (maxEE 2,5)", 225, 220, "2,5", 0.60, 0.54, 0.47), row("FIGARO (a)", 220, 215, "2", 0.50, 0.44, 0.40),
            row("FIGARO (b)", 215, 210, "3", 0.40, 0.36, 0.30)]
    write_compare(runs, "night1_2026-10-07", "mockB", [row("AmpliTrim (maxEE 2)", 225, 220, "2", 0.90, 0.80, 0.50)])      # an older night with the same setting and other numbers: it must not count
    write_compare(runs, "night3_2026-10-09", "mockB", rows, tracks={"AmpliTrim (maxEE 2)": tr, "AmpliTrim (maxEE 2,5)": tr, "FIGARO (a)": tr})
    build = lambda folder, out: SimpleNamespace(r1=ec_sim.dir / "subsample_R1.fastq.gz", r2=ec_sim.dir / "subsample_R2.fastq.gz")      # noqa: E731
    return runs, cfg, build


def test_part_a_alone_reads_the_tables_and_says_that_part_b_was_not_run(ec_sim, tmp_path):
    ec = load()
    runs, cfg, build = setup_runs(tmp_path, ec_sim)
    assert ec.main(["--runs", str(runs), "--out", str(tmp_path / "o")]) == 0
    rep = (tmp_path / "o" / "ESTIMATE_CHECK_REPORT.md").read_text()
    assert "- **P36** MET" in rep and "- **P37** MET" in rep and "- **P38** not run" in rep and "| mockB | 4 | +12.0 | +6.0 | strict | +1.00 | +1.00 |" in rep      # tolerant 15, 13, 10, 10 over DADA2; strict 7, 7, 4, 6; both rank the four settings as DADA2 does
    lines = [r for r in csv.reader(open(tmp_path / "o" / "estimate_check.tsv"), delimiter="\t")]
    assert len(lines) == 5 and lines[0][0] == "dataset" and not (tmp_path / "o" / "cutoff.tsv").exists()


def test_part_b_runs_vsearch_on_the_settings_that_have_a_track_and_reports_the_cut_offs(ec_sim, tmp_path):
    ec = load()
    runs, cfg, build = setup_runs(tmp_path, ec_sim)
    calls = []
    rc = ec.main(["--runs", str(runs), "--cutoff", "--config", str(cfg), "--diffs", "0,1,2,3,10", "--out", str(tmp_path / "o")], merge=make_merge(calls), build=build)
    assert rc == 0
    rep = (tmp_path / "o" / "ESTIMATE_CHECK_REPORT.md").read_text()
    assert "## B. The number of differences" in rep and "mockB: AmpliTrim (maxEE 2) |" in rep and "mockB: AmpliTrim (maxEE 2,5) |" in rep and "FIGARO" not in rep.split("## B.")[1]      # only the AmpliTrim settings (the name filter)
    assert "- **P38**" in rep and "vsearch minus DADA2, mean absolute (points)" in rep
    assert len(calls) == 10 and {c["max_diffs"] for c in calls} == {0, 1, 2, 3, 10}                              # 2 settings x 5 cut-offs
    rows = [r for r in csv.reader(open(tmp_path / "o" / "cutoff.tsv"), delimiter="\t")]
    assert rows[0][:2] == ["dataset", "name"] and "diffs_3" in rows[0] and len(rows) == 3 and float(rows[1][-2]) == pytest.approx(176 / 220, abs=1e-4)


def test_missing_things_are_said_plainly(ec_sim, tmp_path, capsys, monkeypatch):
    ec = load()
    assert ec.main(["--runs", str(tmp_path / "absent"), "--out", str(tmp_path / "o")]) == 2 and "No such folder" in capsys.readouterr().out
    (tmp_path / "empty").mkdir()
    assert ec.main(["--runs", str(tmp_path / "empty"), "--out", str(tmp_path / "o")]) == 2 and "No <run>/<dataset>/compare/compare_methods.tsv found" in capsys.readouterr().out
    runs, cfg, build = setup_runs(tmp_path, ec_sim)
    monkeypatch.setattr("shutil.which", lambda name: None)
    assert ec.main(["--runs", str(runs), "--cutoff", "--config", str(cfg), "--out", str(tmp_path / "o")]) == 2 and "vsearch was not found" in capsys.readouterr().out
    for t in (runs / "night3_2026-10-09" / "mockB" / "compare").glob("track_*.tsv"):
        t.unlink()
    assert ec.main(["--runs", str(runs), "--cutoff", "--config", str(cfg), "--out", str(tmp_path / "o")], merge=make_merge([]), build=build) == 2
    assert "Part B found no setting with a track and a reads folder" in capsys.readouterr().out
