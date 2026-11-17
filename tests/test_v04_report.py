"""validation/v04_report.py: Pauvert's truth, the safety rule, the best trim (P130), the criteria, H and the bases cut."""
import gzip
import json
import math
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "validation"))
import v04_report as V  # noqa: E402

A = "ACGTTGCAAGCTTCGA" * 15                                               # 240 bases


def test_pauverts_truth_complete_partial_and_expected():
    refs = [("M1", "complete", A[:200]), ("M2", "partial", A[50:230]), ("M3", "complete", "T" * 200)]
    t = V.PartialTruth(refs, ["M1", "M2"])
    assert t.members == ["M1", "M2"]
    assert t.classify(A[:200])[:2] == ("true", "M1") and t.classify(A[:199] + "G")[0] == "true"
    longer = A[50:230] + "GGGGCCCC" * 6                                  # the strain's amplicon runs past its partial reference
    assert t.classify(longer)[:2] == ("true", "M2")
    assert t.classify(A[50:120] + "C" + A[121:230] + "T" * 30)[0] == "true"   # one edit inside the partial reference
    assert t.classify("A" * 300)[0] == "spurious" and t.classify(A[50:150])[0] == "spurious"   # too short to be the strain
    assert t.classify(A[50:170], prefix=True)[:2] == ("true", "M2")      # a forward read shorter than the partial reference
    assert t.classify("T" * 200)[:2] == ("true", "M3")                    # a reference that is not expected still judges true


def R_with(acc, flow=None, cut=None, merge=None):
    return {"acc": acc, "flow": flow or {}, "cut": cut or {}, "merge": merge or {}, "est": {}, "pred": {}}


def a(sp, rec, found=5.0):
    return {"spurious": sp, "recall": rec, "exact": 0.9, "found": found, "unjudged": 0.0}


def test_the_safety_rule_at_its_limits():
    R = R_with({("x", "t"): a(0.02, 0.88), ("x", "std"): a(0.01, 0.90), ("y", "t"): a(0.021, 0.90), ("y", "std"): a(0.01, 0.90),
                ("z", "t"): a(0.0, 0.879), ("z", "std"): a(0.0, 0.90)})
    assert V.safe(R, "x", "t", "std")[0] is True                          # +1.00 and -2.0: at the limits, safe
    assert V.safe(R, "y", "t", "std")[0] is False and V.safe(R, "z", "t", "std")[0] is False
    assert V.safe(R, "w", "t", "std")[0] is None


def trim_R(yields, cuts, unsafe=()):
    acc, flow, cut = {}, {}, {}
    for lab in V.TRUTH:
        acc[(lab, "std")] = a(0.0, 0.9)
        flow[(lab, "std")] = {"yield": 0.5, "pass": 0.5, "merged": 0.5, "input": 100, "asv_median": 10}
        cut[lab] = {}
        for arm in V.TRIMS:
            acc[(lab, arm)] = a(0.05 if arm in unsafe else 0.0, 0.9)
            flow[(lab, arm)] = {"yield": yields.get(arm, 0.6), "pass": 0.9, "merged": 0.8, "input": 100, "asv_median": 10}
            cut[lab][arm] = cuts.get(arm, 30.0)
    return R_with(acc, flow, cut)


def test_the_best_trim_rule():
    R = trim_R({"sw15": 0.80, "w8": 0.803, "w2": 0.90, "q12": 0.70}, {"sw15": 30.0, "w8": 22.0, "w2": 5.0}, unsafe={"w2"})
    b = V.best_trim(R)
    assert b["best"] == "w8" and set(b["tied"]) == {"sw15", "w8"}          # w2 keeps more but is unsafe; w8 ties sw15 and cuts less
    R = trim_R({"sw15": 0.804, "w8": 0.800}, {"sw15": 30.0, "w8": 22.0})
    assert V.best_trim(R)["best"] == "w8"                                 # tied: the gentler trim wins although sw15 keeps 0.4 points more
    R = trim_R({"sw15": 0.80, "w8": 0.81}, {"sw15": 30.0, "w8": 40.0})
    assert V.best_trim(R)["best"] == "w8"                                 # 1 point above: no tie, the highest wins even if it cuts more


def test_a_gate_that_fails_decides_nothing_and_the_tq_that_changes_nothing_is_not_judged():
    R = trim_R({"sw15": 0.8}, {})
    for lab in V.TRUTH:
        for arm in V.TQS + V.MINLENS + ["eer", "fwd", "fwd15"]:
            R["acc"][(lab, arm)] = a(0.0, 0.9)
            R["flow"][(lab, arm)] = dict(R["flow"][(lab, "std")])
    R["h"] = (True, "16S same")
    assert all(ok is None for _, ok, _ in V.verdicts(R, False))
    V_ = {c: (ok, t) for c, ok, t in V.verdicts(R, True)}
    assert V_["P122 tq3"][0] is None and "changes nothing" in V_["P122 tq3"][1] and V_["P127"][0] is True
    assert V_["P117"][0] is True and V_["P130"][0] is True


def summary(path, rec):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"recommended": rec}))


def test_h_the_16s_comparison(tmp_path):
    rec = {"truncLenF": 240, "truncLenR": 200, "maxEE": 2.0, "maxEE_f": 2.0, "maxEE_r": 2.0, "estimated_yield": 0.9, "minOverlap": 12, "maxMismatch": 0}
    summary(tmp_path / "cli" / "v4" / "report" / "summary.json", rec)
    new = dict(rec, truncQ=2)
    summary(tmp_path / "night" / "16s_v4" / "a04" / "summary.json", new)
    summary(tmp_path / "night" / "16s_v4" / "a04_explicit" / "summary.json", new)
    assert V.compare_16s(tmp_path / "night", tmp_path / "cli")[0] is True
    summary(tmp_path / "night" / "16s_v4" / "a04_explicit" / "summary.json", dict(new, truncLenF=239))
    ok, text = V.compare_16s(tmp_path / "night", tmp_path / "cli")
    assert ok is False and "DIFFERENT" in text


def test_bases_cut(tmp_path):
    for name, n in (("reads", 100), ("reads_sw15", 80), ("reads_w8", 90)):
        (tmp_path / name).mkdir()
        with gzip.open(tmp_path / name / "s_R1.fastq.gz", "wt") as fh:
            fh.write(f"@a\n{'A' * n}\n+\n{'I' * n}\n@b\n{'A' * n}\n+\n{'I' * n}\n")
    c = V.bases_cut(tmp_path)
    assert c["sw15"] == 20 and c["w8"] == 10 and "w2" not in c


def test_the_tidy_step_records_the_bases_cut_then_deletes_the_bulk(tmp_path):
    import its_night as NIGHT
    d = tmp_path / "lab"
    for name, n in (("reads", 100), ("reads_sw15", 80), ("reads_w8", 90), ("reads_q12", 95)):
        (d / name).mkdir(parents=True)
        with gzip.open(d / name / "s_R1.fastq.gz", "wt") as fh:
            fh.write(f"@a\n{'A' * n}\n+\n{'I' * n}\n")
    (d / "dada2" / "std" / "filtered").mkdir(parents=True)
    (d / "dada2" / "std" / "filtered" / "x.fastq.gz").write_bytes(b"x")
    (d / "dada2" / "std" / "samples.tsv").write_text("kept")
    NIGHT.Night.step_tidy(None, d)
    assert json.loads((d / "cut.json").read_text())["mean_bases_cut"] == {"sw15": 20.0, "w8": 10.0, "q12": 5.0}
    assert not (d / "reads_sw15").exists() and not (d / "dada2" / "std" / "filtered").exists()
    assert (d / "reads").is_dir() and (d / "dada2" / "std" / "samples.tsv").read_text() == "kept"
    assert V.bases_cut(d) == {"sw15": 20.0, "w8": 10.0, "q12": 5.0}
