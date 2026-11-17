"""validation/pauvert_prep.py: the truth of the fourth mock community from Table S1, and its reads made ready."""
import gzip
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "validation"))
import pauvert_prep as P  # noqa: E402
import bakker_prep as B  # noqa: E402


def levenshtein(a, b):
    prev = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        cur = [i] + [0] * len(b)
        for j in range(1, len(b) + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
        prev = cur
    return prev[-1]


def test_the_banded_distance_agrees_with_the_plain_one():
    rng = random.Random(1)
    for _ in range(400):
        a = "".join(rng.choice("ACGT") for _ in range(rng.randint(20, 60)))
        b = list(a)
        for _ in range(rng.randint(0, 4)):
            k = rng.randint(0, len(b) - 1)
            op = rng.choice("sid")
            if op == "s":
                b[k] = rng.choice("ACGT")
            elif op == "i":
                b.insert(k, rng.choice("ACGT"))
            elif len(b) > 1:
                del b[k]
        b = "".join(b)
        assert P.within(a, b) == (levenshtein(a, b) <= 2), (a, b)


SITE = B.rc(P.ITS2)


def test_the_cut_and_the_prefix_rule():
    its1 = "ACGT" * 20                                                     # 80 bases, then the ITS2 site
    assert P.cut_its1(its1 + SITE + "TTTT") == its1 and P.cut_its1(its1) == ""
    part = its1[:60]
    assert P.prefix_within(part, its1 + SITE) and P.prefix_within(part[:30] + "G" + part[31:], its1)
    assert not P.prefix_within("T" * 60, its1)


def table(path, rows):
    head = "Mock_Strain;Phylum;Present_In_Raw_MiSeq_Data;Unique_ITS1_Sanger_Sequence;Sanger_Sequence\n"
    path.write_text(head + "".join(f"{n};X;{p};Yes;{s}\n" for n, p, s in rows), encoding="utf-8")


def test_the_truth_from_a_small_table(tmp_path):
    rng = random.Random(2)
    seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    a, b, c = seq(200), seq(220), seq(240)
    rows = [("Alpha one", "Yes", a + SITE + "AAAA"),                       # complete
            ("Alpha two", "Yes", a[:100] + "T" + a[101:] + SITE),          # 1 edit from Alpha one: one member
            ("Beta", "Yes", b[:230]),                                      # partial, judged as a prefix
            ("Gamma", "No", c + SITE),                                     # complete, not found by the authors: reference, not expected
            ("Short", "Yes", seq(81)),                                     # partial and too short: left out
            ("Coded", "Yes", "".join(rng.choice("ACGTRY") for _ in range(200)) + SITE)]   # too many variants (very likely)
    table(tmp_path / "t.csv", rows)
    res = P.truth(P.table(tmp_path / "t.csv"), tmp_path / "ref")
    assert res["kept"] == 4 or (res["kept"] == 5 and not res["too_coded"])
    assert res["partial"] == 1 and any("Short" in x for x in res["no_site"]) and res["members"] == res["kept"] - 1
    assert res["expected"] == res["members"] - 1                          # Gamma is not expected
    fa = (tmp_path / "ref" / "references_pauvert_its1.fasta").read_text()
    assert "|partial|" in fa and "|complete|" in fa and "Alpha_one+1" in fa
    exp = (tmp_path / "ref" / "expected_pauvert_its1.tsv").read_text()
    assert "Gamma" not in exp and "Beta" in exp and "PAUVERT" not in exp
    assert "Partial (the Sanger read stopped before the ITS2 site; judged as a prefix): 1" in (tmp_path / "ref" / "PAUVERT_REFERENCES.md").read_text()


def fq(path, recs):
    with gzip.open(path, "wt") as fh:
        fh.write("".join(f"@{n}\n{s}\n+\n{'I' * len(s)}\n" for n, s in recs))


def test_the_reads_are_paired_by_replicate_and_cut(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    amp = "ACGTTGCA" * 15                                                  # 120 bases of amplicon
    r1 = P.ITS1F + amp[:100]
    r2 = P.ITS2 + B.rc(amp)[:100]
    for rep in ("DikA", "DikB"):
        fq(raw / f"Forward_R1_ITSEnz.{rep}.fastq.gz", [("p1", r1), ("p2", "T" * 120)])
        fq(raw / f"Reverse_R2_ITSEnz.{rep}.fastq.gz", [("p1", r2), ("p2", "T" * 120)])
    rows = P.reads(raw, tmp_path / "ready")
    assert [(s, n, k) for s, n, k, _, _ in rows] == [("DikA", 2, 1), ("DikB", 2, 1)]
    with gzip.open(tmp_path / "ready" / "DikA_R1.fastq.gz", "rt") as fh:
        out = fh.read().splitlines()
    assert out[0] == "@p1" and out[1] == amp[:100] and not out[1].startswith(P.ITS1F[:10])
