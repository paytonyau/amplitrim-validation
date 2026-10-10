"""The search: where it starts, how it steps, when it stops, and what it accepts as input. FIGARO is
replaced by a stand-in so only vsearch is needed."""
import gzip

import pytest

from amplitrim import pipeline
from amplitrim.figaro_wrapper import FigaroCandidate
from conftest import need_vsearch

pytestmark = need_vsearch


@pytest.fixture
def settings(reads, tmp_path, monkeypatch):
    """Settings for the simulated reads, in FIGARO cross-check mode with a stand-in FIGARO that
    proposes two fixed splits per total."""
    requested = []

    def fake_figaro(script, py, inp, work, total, min_overlap, top_n=3, compat=True):
        requested.append(total)
        return [FigaroCandidate(trim_f=total // 2, trim_r=total - total // 2),
                FigaroCandidate(trim_f=total // 2 - 3, trim_r=total - total // 2 + 3)]

    monkeypatch.setattr(pipeline, "candidates_for_total", fake_figaro)
    figaro = tmp_path / "figaro.py"
    figaro.touch()
    s = pipeline.Settings(input_dir=reads, work_dir=tmp_path / "work", figaro_script=figaro, splits="both")
    s.requested = requested
    return s


def test_search_starts_at_the_exact_coverage_floor_and_steps_upward(settings):
    res = pipeline.run(settings)
    assert res.totals[0] == res.floor == res.longest_insert + 12
    assert all(b - a == settings.step for a, b in zip(res.totals, res.totals[1:]))
    assert min(settings.requested) >= res.floor                    # nothing shorter is ever tested


def test_search_stops_when_yield_stops_improving(settings):
    stopped = pipeline.run(settings)
    settings.exhaustive = True
    full = pipeline.run(settings)
    assert len(stopped.totals) < len(full.totals)
    assert full.totals[-1] + settings.step > full.max_combined      # exhaustive goes to the read limit
    assert stopped.recommendation.total_length <= full.recommendation.total_length


def test_the_full_reads_are_tested_without_figaro_when_nothing_is_trimmed(settings):
    settings.exhaustive = True
    res = pipeline.run(settings)
    assert res.max_combined not in settings.requested              # FIGARO is not asked for this total
    full = [c for c in res.candidates if c.total_length == res.max_combined]
    assert [(c.trim_f, c.trim_r) for c in full] == [(res.sub.read_len_f, res.sub.read_len_r)]


def test_recommendation_covers_the_longest_amplicon(settings):
    res = pipeline.run(settings)
    assert res.recommendation.covers_longest
    assert res.recommendation.max_mergeable_insert >= res.longest_insert


def test_a_given_longest_insert_sets_the_floor(settings):
    settings.longest_insert = 430
    assert pipeline.run(settings).floor == 442


def test_min_overlap_moves_the_floor_and_the_rule_merging_is_tested_under(settings):
    settings.splits = "quality"
    settings.min_overlap = 20
    res = pipeline.run(settings)
    assert res.floor == res.longest_insert + 20 and res.cfg.dada2_min_overlap == 20
    assert res.recommendation.max_mergeable_insert == res.recommendation.total_length - 20


def test_unestimable_amplicon_lengths_ask_for_longest_insert(settings, monkeypatch):
    monkeypatch.setattr(pipeline, "estimate_insert_sizes", lambda *a, **k: [430] * 10)
    with pytest.raises(pipeline.PipelineError, match="--longest-insert"):
        pipeline.run(settings)


def test_reads_too_short_for_the_longest_amplicon_give_no_recommendation(settings):
    settings.longest_insert = 600
    res = pipeline.run(settings)
    assert res.recommendation is None and not res.totals
    assert any("No setting can merge" in w for w in res.warnings)


def test_a_work_folder_vsearch_for_windows_could_not_open_is_refused_early(tmp_path):
    from pathlib import Path
    pipeline.check_vsearch_path(Path("C:/data/work"), "ascii")                    # plain characters: fine
    pipeline.check_vsearch_path(Path("C:/donn\u00e9es/work"), "cp1252")           # accented letters exist in cp1252: fine
    for path, encoding in ((Path("C:/donn\u00e9es/work"), "ascii"), (Path("C:/\u6578\u64da/work"), "cp1252")):
        with pytest.raises(pipeline.PipelineError, match="work folder .* cannot open"):
            pipeline.check_vsearch_path(path, encoding)


def test_missing_figaro_script_is_a_clear_error(settings):
    settings.figaro_script = settings.work_dir / "missing.py"
    with pytest.raises(pipeline.PipelineError, match="FIGARO script not found"):
        pipeline.run(settings)


def test_input_folder_problems_are_plain_errors(settings, tmp_path):
    settings.input_dir = tmp_path / "does_not_exist"
    with pytest.raises(pipeline.PipelineError, match="Input folder not found"):
        pipeline.run(settings)
    empty = tmp_path / "empty"
    empty.mkdir()
    settings.input_dir = empty
    with pytest.raises(pipeline.PipelineError, match="No paired FASTQ files found"):
        pipeline.run(settings)


# ---- who proposes the forward/reverse splits ---------------------------------------------------
def test_each_total_also_tests_the_best_pass_rate_split_and_labels_who_proposed_it(settings):
    res = pipeline.run(settings)
    assert all(c.source for c in res.candidates)
    assert any("quality table" in c.source for c in res.candidates)
    assert any(c.source == "FIGARO" for c in res.candidates)       # FIGARO's own splits are still tested


def test_quality_mode_does_not_need_figaro(settings):
    settings.splits = "quality"
    res = pipeline.run(settings)
    assert settings.requested == [] and res.recommendation is not None
    assert {c.source for c in res.candidates} <= {"quality table", "full reads"}


def test_figaro_mode_adds_no_quality_table_candidates(settings):
    settings.splits = "figaro"
    res = pipeline.run(settings)
    assert {c.source for c in res.candidates} <= {"FIGARO", "full reads"}


# ---- reads of different lengths ---------------------------------------------------------------
@pytest.fixture
def mixed(reads, tmp_path):
    """The simulated reads with every 10th pair cut short (as after adapter read-through trimming)."""
    inp = tmp_path / "mixed"
    inp.mkdir()
    for mate in ("R1", "R2"):
        with gzip.open(reads / f"sampleA_{mate}.fastq.gz", "rt") as src, \
                gzip.open(inp / f"s_{mate}.fastq.gz", "wt", newline="\n") as out:
            lines = src.read().splitlines()
            for i in range(0, len(lines), 4):
                h, seq, plus, qual = lines[i:i + 4]
                if (i // 4) % 10 == 0:
                    seq, qual = seq[:150], qual[:150]
                out.write(f"{h}\n{seq}\n{plus}\n{qual}\n")
    return inp


def test_reads_of_different_lengths_are_fine_by_default_and_the_loss_is_reported(settings, mixed):
    settings.splits, settings.input_dir = "quality", mixed
    res = pipeline.run(settings)
    assert not res.sub.uniform_length and res.recommendation is not None
    assert res.short_lost == pytest.approx(0.1, abs=0.01)           # the cut-short 10% of pairs
    assert res.recommendation.ee_pass_rate <= 0.9 + 1e-9            # ... can never pass
    assert any("Read lengths vary" in w and "loses 10%" in w for w in res.warnings)


def test_figaro_needs_reads_of_one_length_and_says_so(settings, mixed):
    settings.input_dir = mixed
    with pytest.raises(pipeline.PipelineError, match="FIGARO needs reads of one length"):
        pipeline.run(settings)


# ---- ITS: no truncation ------------------------------------------------------------------------
def test_its_mode_evaluates_no_truncation_on_reads_of_different_lengths(settings, mixed):
    settings.marker, settings.splits, settings.input_dir = "its", "quality", mixed
    res = pipeline.run(settings)
    r = res.recommendation
    assert (r.trim_f, r.trim_r) == (0, 0) and res.totals == [] and len(res.candidates) == 1
    assert r.source == "no truncation" and 0 < r.yield_est <= 1
    assert r.ee_pass_rate * r.merge_rate_tolerant == pytest.approx(r.yield_est)
    assert not any("Read lengths vary" in w for w in res.warnings)   # nothing is truncated, so no length loss


def test_unknown_marker_is_a_plain_error(settings):
    settings.marker = "18s"
    with pytest.raises(pipeline.PipelineError, match="Unknown marker"):
        pipeline.run(settings)
