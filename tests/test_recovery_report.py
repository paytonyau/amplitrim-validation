"""validation/recovery_report.py: the counting of members found (per sample, in any, in all; at 0.1% and at any read) and the 16S strains."""
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "validation"))
import recovery_report as R  # noqa: E402


def test_counting_per_sample_any_and_all():
    found = {"s1": {"A": 0.5, "B": 0.0005, "C": 0.2}, "s2": {"A": 0.4, "B": 0.01}, "s3": {}}
    expected = {"s1": ["A", "B", "C"], "s2": ["A", "B", "C"], "s3": ["A", "B", "C"]}
    c = R.count(found, expected)
    assert c["samples"] == 3 and c["expected"] == 3
    assert c["found01"] == pytest.approx((2 + 2 + 0) / 3) and c["found1"] == pytest.approx((3 + 2 + 0) / 3)
    assert (c["any01"], c["any1"], c["any_expected"]) == (3, 3, 3) and (c["all01"], c["all1"]) == (0, 0)
    c = R.count({"s1": found["s1"], "s2": found["s2"]}, {"s1": expected["s1"], "s2": expected["s2"]})
    assert (c["all01"], c["all1"]) == (1, 2)                                  # A at 0.1% in both; A and B with any read in both


def test_members_not_expected_in_a_sample_do_not_count_there():
    c = R.count({"MixB-a": {"K": 0.5, "T": 0.5}}, {"MixB-a": ["K", "C"]})     # T is not in mix B
    assert c["found1"] == 1 and c["any_expected"] == 2 and c["any1"] == 1


def test_sixteen_strains_through_their_variants(tmp_path):
    base, run = tmp_path / "data", tmp_path / "run"
    ref = base / R.SIXTEEN["mockB"]
    ref.parent.mkdir(parents=True)
    ref.write_text(">v1|S1|a\nAAAA\n>v2|S1|b\nCCCC\n>v3|S2|a\nGGGG\n")
    d = run / "16s_mockB" / "dada2" / "std"
    d.mkdir(parents=True)
    (d / "samples.tsv").write_text("sample\tinput\tfiltered\tdenoised\tmerged\tnonchim\tasvs\nx\t10\t10\t10\t10\t10\t2\ny\t10\t10\t10\t10\t10\t1\n")
    (d / "seqs.tsv").write_text("sample\tsequence\tabundance\nx\tCCCC\t90\nx\tGGGT\t10\ny\tGGGG\t100\n")
    hits = lambda seqs, reference: {"CCCC": ("v2|S1|b", 100.0), "GGGT": ("v3|S2|a", 99.5), "GGGG": ("v3|S2|a", 100.0)}
    rows = {r["dataset"]: r for r in R.sixteen_rows(run, base, hits)}
    m = rows["mockB"]
    assert m["variants"] == "2 of 3" and m["expected"] == 2                   # S1 through its variant v2; GGGT only near, not exact
    assert m["found1"] == 1 and (m["any1"], m["all1"]) == (2, 0)
    assert "missing" in rows["mockA_high"]
