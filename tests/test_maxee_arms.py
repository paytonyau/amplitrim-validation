"""The measures and rules of the maxEE head-to-head (validation/maxee_arms.py), on reads and tables built here."""
import gzip
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "validation"))
sys.path.insert(0, str(HERE.parent))
import maxee_arms as ma                                    # noqa: E402
from amplitrim.trimmer import ReadPairs                    # noqa: E402


def write_fq(path, reads):
    with gzip.open(path, "wt") as fh:
        for i, (seq, q) in enumerate(reads):
            fh.write(f"@p{i}\n{seq}\n+\n{q}\n")


def pairs_with(tmp_path, n=200, f_good=100, r_good=60, length=150):
    """Forward reads are Q37 for f_good bases and Q10 after; reverse reads are Q37 for r_good bases and Q15 after."""
    q37, q10, q15 = chr(33 + 37), chr(33 + 10), chr(33 + 15)
    r1 = [("A" * length, q37 * f_good + q10 * (length - f_good)) for _ in range(n)]
    r2 = [("C" * length, q37 * r_good + q15 * (length - r_good)) for _ in range(n)]
    write_fq(tmp_path / "a_R1.fq.gz", r1)
    write_fq(tmp_path / "a_R2.fq.gz", r2)
    return ReadPairs.load(tmp_path / "a_R1.fq.gz", tmp_path / "a_R2.fq.gz")


# ------------------------------------------------------------------------------------------------ the rules
def test_the_median_q30_cut_is_where_the_median_quality_first_drops_below_30(tmp_path):
    assert ma.q30_cut(pairs_with(tmp_path)) == (100, 60)
    assert ma.q30_cut(pairs_with(tmp_path, f_good=150, r_good=150)) == (150, 150)         # never below 30: no cut, the read length
    assert ma.q30_cut(pairs_with(tmp_path, f_good=20, r_good=10)) == (50, 50)             # never shorter than the minimum


def test_the_figaro_like_rule_is_the_mean_expected_errors_at_the_cut_plus_one_rounded(tmp_path):
    pairs = pairs_with(tmp_path)
    # forward cut 100: 100 bases of Q37 = 0.02 errors, +1 -> 1; reverse cut 150: 60 x Q37 + 90 x Q15 = 0.012 + 2.846, +1 -> 3.86 -> 4
    assert ma.rule_ee(pairs, 100, 150) == (1, 4)
    assert ma.rule_ee(pairs, 100, 60) == (1, 1)                                           # clean reads: never below 1
    assert [ma.round_ee(x) for x in (0.2, 1.0, 2.49, 2.5, 3.5, 6.5)] == [1, 1, 2, 3, 4, 7]                 # halves go up; never below 1


def fake_result(rows, inserts, trim=(250, 240), longest=None):
    return SimpleNamespace(recommendation=SimpleNamespace(trim_f=trim[0], trim_r=trim[1]), ee_search=rows, inserts=inserts, longest_insert=longest or max(inserts))


def row(ee, y, ef, er, tf=250, tr=240):
    return {"ee": ee, "yield_est": y, "pass_rate": y + 0.01, "ee_f": ef, "ee_r": er, "trim_f": tf, "trim_r": tr, "current": ee == (2, 2)}


def test_explore_gives_the_quality_first_pick_at_each_floor_and_the_amplicon_lengths(tmp_path):
    rows = [row((1, 1), 0.05, 0.1, 0.5), row((2, 2), 0.20, 0.5, 1.0), row((2, 3), 0.35, 0.6, 1.5), row((2, 4), 0.52, 0.7, 2.0),
            row((2, 5), 0.70, 0.8, 2.4), row((3, 5), 0.85, 1.0, 2.5), row((5, 5), 0.90, 1.5, 3.0)]
    info = ma.explore(tmp_path, tmp_path / "w", run_fn=lambda settings: fake_result(rows, [370] * 98 + [550] * 2))
    assert info["picks"][0.5]["ee"] == (2, 4) and info["picks"][0.5]["reached"]          # the cleanest row that keeps half the pairs, never stricter than 2
    assert info["picks"][0.8]["ee"] == (3, 5)
    assert info["trim"] == (250, 240) and info["longest"] == 550 and info["typical"] == 370
    assert ma.has_tail(info["longest"], info["typical"]) and not ma.has_tail(430, 428)
    plan = ma.arm_plan(info, rule=(2, 4))
    assert plan == {"A": "2", "B": "2,4", "B2": "3,5", "B4": "2,4"} and ma.unique_ees(plan) == ["2", "2,4", "3,5"]      # B and B4 are one run
    assert ma.arm_plan(info)["B4"] is None and ma.unique_ees(ma.arm_plan(info)) == ["2", "2,4", "3,5"]


def test_a_floor_nothing_reaches_still_gives_the_row_that_keeps_most(tmp_path):
    rows = [row((2, 2), 0.20, 0.5, 1.0), row((2, 3), 0.30, 0.6, 1.5), row((5, 5), 0.40, 1.5, 3.0)]
    info = ma.explore(tmp_path, tmp_path / "w", run_fn=lambda s: fake_result(rows, [300] * 10))
    assert info["picks"][0.5]["ee"] == (5, 5) and not info["picks"][0.5]["reached"]


# ------------------------------------------------------------------------------------------------ the measures
def test_paired_difference_has_an_interval_only_from_ten_samples():
    a = {f"s{i}": 0.80 + 0.001 * i for i in range(12)}
    b = {f"s{i}": 0.70 for i in range(12)}
    d = ma.paired_diff(a, b)
    assert d["n"] == 12 and abs(d["mean"] - 10.55) < 1e-9 and d["lo"] < d["mean"] < d["hi"] and d["lo"] > 10.0
    assert ma.paired_diff(a, b) == d                                                      # deterministic
    few = ma.paired_diff({k: a[k] for k in list(a)[:5]}, b)
    assert few["n"] == 5 and few["lo"] is None and few["mean"] is not None
    assert ma.paired_diff({}, b) == {"n": 0, "mean": None, "lo": None, "hi": None}
    assert ma.paired_diff({"x": 1.0}, {"y": 1.0})["n"] == 0                               # no common sample


def test_retention_and_noise_classes_use_the_fixed_margins():
    assert [ma.retention_class(x) for x in (2.5, 2.0, 0.0, -2.0, -2.5, None)] == ["more", "same", "same", "same", "fewer", "unknown"]
    assert ma.noise_class((1.0, 2.0), (1.0, 2.0)) == "about the same"
    assert ma.noise_class((1.01, 2.01), (1.0, 2.0)) == "about the same"                   # inside 0.02
    assert ma.noise_class((0.5, 1.9), (1.0, 2.0)) == "cleaner"
    assert ma.noise_class((1.5, 2.1), (1.0, 2.0)) == "noisier"
    assert ma.noise_class((0.5, 2.5), (1.0, 2.0)) == "mixed"                              # one read each way
    assert ma.noise_class(None, (1.0, 2.0)) == ma.noise_class((float("nan"), 1), (1, 1)) == "unknown"


def test_the_verdict_has_a_winner_only_when_one_pick_is_better_without_being_worse():
    v = ma.verdict
    assert v("same", "about the same") == "matches"
    assert v("more", "about the same") == v("same", "cleaner") == v("more", "cleaner") == "better"
    assert v("fewer", "about the same") == v("same", "noisier") == v("fewer", "noisier") == "worse"
    assert v("more", "noisier") == "trade-off: more pairs, noisier"
    assert v("fewer", "cleaner") == "trade-off: fewer pairs, cleaner"
    assert v("more", "mixed") == "mixed noise (more pairs)" and v("unknown", "cleaner") == "no result"


def test_accuracy_against_the_fixed_battery_margins():
    c = {"found": "7", "strains": "3/3", "identical": "0.7828"}
    assert ma.accuracy_verdict({"found": "6", "strains": "3/3", "identical": "0.7526"}, c) == "less accurate"           # 3 points fewer reads in exact ASVs
    assert ma.accuracy_verdict({"found": "7", "strains": "3/3", "identical": "0.78285"}, c) == "as accurate"
    assert ma.accuracy_verdict({"found": "9", "strains": "3/3", "identical": "0.70"}, c) == "mixed"                      # two more variants but fewer exact reads
    assert ma.accuracy_verdict({"found": "8", "strains": "3/3", "identical": "0.79"}, {"found": "6", "strains": "3/3", "identical": "0.78"}) == "more accurate"
    assert ma.accuracy_verdict({"strains": None}, {"strains": None}) == "no truth"                                        # a real dataset: nothing to judge


def test_arm_names_are_stable_for_the_report():
    assert ma.ee_str((2, 2)) == "2" and ma.ee_str((2, 5)) == "2,5" and ma.ee_str((6, 9)) == "6,9"


def test_the_no_truncation_arm_gets_a_minimum_length_and_a_truncated_one_does_not():
    import compare_methods as cm
    assert cm.min_len_for("16s", 0, 0, 50) == 50 and cm.min_len_for("16s", 250, 0, 50) == 50          # any arm with a 0 length
    assert cm.min_len_for("16s", 250, 240, 50) == 0                                                    # truncation discards the shorter reads itself
    assert cm.min_len_for("its", 200, 200, 50) == 50                                                   # the ITS workflow always has one
