"""validation/its_night.R through a stand-in for dada2 (tests/its_night_stub.R): the arms, the chimera removal, the reuse of the default arm's denoising, the samples
that keep no read, and the restart. Needs Rscript; skipped without it. The numbers follow from the stand-in: dada gives SEQA 60, SEQB 30 and CHIM1 10 reads per sample;
mergePairs gives MERGEDA 50 (overlap 40), MERGEDB 25 (25), CHIM2 10 (40) and MERGEDC 15 (15); filterAndTrim keeps all reads less max(0, 6 - maxEEF - maxEER)."""
import csv
import os
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "validation" / "its_night.R"
STUB = HERE / "its_night_stub.R"
sys.path.insert(0, str(HERE.parent / "validation"))
import rtools  # noqa: E402  (finds Rscript as the night does: on PATH, under R_HOME, or in the usual Windows install folders)

RSCRIPT = rtools.find_rscript()
pytestmark = pytest.mark.skipif(RSCRIPT is None, reason="Rscript was not found")
COLS = ["arm", "mode", "reads", "truncF", "truncR", "maxEEF", "maxEER", "truncQ", "minOverlap"]


def fastq(path, n):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(n)))


def reads(folder, samples, mates=(1, 2)):
    for s, n in samples.items():
        for m in mates:
            fastq(folder / f"{s}_R{m}.fastq", n)
    return folder


def arms_file(path, rows):
    path.write_text("\t".join(COLS) + "\n" + "".join("\t".join(str(x) for x in r) + "\n" for r in rows))
    return path


def run(arms, out, tmp_path):
    env = dict(os.environ, ITS_NIGHT_STUB=str(STUB), ITS_NIGHT_LOG=str(tmp_path / "calls.log"))
    return subprocess.run([RSCRIPT, str(SCRIPT), str(arms), str(out), "FALSE"], capture_output=True, text=True, env=env)


def table(path):
    with open(path, newline="") as fh:
        return {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}


def num(row, col):
    return float(row[col])


@pytest.fixture
def night(tmp_path):
    paired = reads(tmp_path / "reads", {"s1": 10, "s2": 8, "s3empty": 5})
    merged = reads(tmp_path / "merged", {"m1": 12}, mates=(1,))
    arms = arms_file(tmp_path / "arms.tsv", [
        ("std", "paired", paired, 0, 0, 2, 2, 2, 12),
        ("ee25", "paired", paired, 0, 0, 2, 5, 2, 12),
        ("tq10", "paired", paired, 0, 0, 2, 2, 10, 12),
        ("fwd", "fwd", "std", 0, 0, 2, 2, 2, 12),
        ("mo30", "overlap", "std", 0, 0, 2, 2, 2, 30),
        ("mfirst", "single", merged, 0, 0, 2, 2, 2, 12)])
    p = run(arms, tmp_path / "out", tmp_path)
    assert p.returncode == 0, p.stderr
    return tmp_path, arms, p


def test_the_default_arm_counts_every_step_and_removes_the_chimeras(night):
    tmp, _, _ = night
    t = table(tmp / "out" / "std" / "samples.tsv")
    s1 = t["s1"]
    assert (num(s1, "input"), num(s1, "filtered"), num(s1, "denoised"), num(s1, "merged"), num(s1, "nonchim"), num(s1, "asvs")) == (10, 8, 100, 100, 90, 3)
    seqs = [r for r in csv.DictReader(open(tmp / "out" / "std" / "seqs.tsv"), delimiter="\t") if r["sample"] == "s1"]
    assert sorted(r["sequence"] for r in seqs) == ["MERGEDA", "MERGEDB", "MERGEDC"]               # CHIM2 is gone
    assert (tmp / "out" / "std" / "filtered").is_dir() and (tmp / "out" / "std" / "objects.rds").is_file() and (tmp / "out" / "std" / "DONE").is_file()


def test_a_sample_that_keeps_no_read_gives_zeros_not_an_error(night):
    tmp, _, _ = night
    e = table(tmp / "out" / "std" / "samples.tsv")["s3empty"]
    assert num(e, "input") == 5 and all(num(e, c) == 0 for c in ("filtered", "denoised", "merged", "nonchim", "asvs"))


def test_each_arm_passes_its_own_settings(night):
    tmp, _, _ = night
    assert num(table(tmp / "out" / "ee25" / "samples.tsv")["s1"], "filtered") == 10                 # maxEE 2,5 keeps every read in the stand-in
    calls = (tmp / "calls.log").read_text().splitlines()
    filt = [c.split("\t") for c in calls if c.startswith("filterAndTrim")]
    assert ["filterAndTrim", "paired", "0,0", "2,5", "2", "50"] in filt
    assert ["filterAndTrim", "paired", "0,0", "2,2", "10", "50"] in filt                            # tq10
    assert ["filterAndTrim", "single", "0", "2", "2", "50"] in filt                                 # mfirst
    assert all(c.split("\t")[1] == "consensus" for c in calls if c.startswith("removeBimeraDenovo"))
    assert ["mergePairs", "30", "0"] in [c.split("\t") for c in calls if c.startswith("mergePairs")]


def test_the_overlap_and_forward_arms_reuse_the_denoising_of_the_default_arm(night):
    tmp, _, _ = night
    mo = table(tmp / "out" / "mo30" / "samples.tsv")["s1"]
    assert (num(mo, "filtered"), num(mo, "denoised"), num(mo, "merged"), num(mo, "nonchim")) == (8, 100, 60, 50)      # only the pairs with an overlap of 30 or more
    fw = table(tmp / "out" / "fwd" / "samples.tsv")["s1"]
    assert (num(fw, "denoised"), num(fw, "merged"), num(fw, "nonchim"), num(fw, "asvs")) == (100, 100, 90, 2)        # SEQA and SEQB; CHIM1 is gone
    calls = (tmp / "calls.log").read_text().splitlines()
    assert sum(c.startswith("learnErrors") for c in calls) == 2 * 3 + 1                              # std, ee25 and tq10 (two each) and mfirst; not fwd or mo30


def test_a_single_end_arm_denoises_the_r1_files(night):
    tmp, _, _ = night
    m = table(tmp / "out" / "mfirst" / "samples.tsv")["m1"]
    assert (num(m, "input"), num(m, "filtered"), num(m, "denoised"), num(m, "merged"), num(m, "nonchim")) == (12, 8, 100, 100, 90)


def test_a_finished_arm_is_not_run_again(night):
    tmp, arms, _ = night
    before = (tmp / "out" / "std" / "samples.tsv").stat().st_mtime_ns
    p = run(arms, tmp / "out", tmp)
    assert p.returncode == 0 and p.stdout.count("done before") == 6
    assert (tmp / "out" / "std" / "samples.tsv").stat().st_mtime_ns == before


def test_one_sample_alone_works_although_dada_then_returns_no_list(tmp_path):
    one = reads(tmp_path / "one", {"only": 10})
    arms = arms_file(tmp_path / "arms.tsv", [("std", "paired", one, 0, 0, 2, 2, 2, 12), ("fwd", "fwd", "std", 0, 0, 2, 2, 2, 12)])
    p = run(arms, tmp_path / "out", tmp_path)
    assert p.returncode == 0, p.stderr
    assert num(table(tmp_path / "out" / "std" / "samples.tsv")["only"], "nonchim") == 90


def test_an_arm_that_reuses_another_must_come_after_it(tmp_path):
    paired = reads(tmp_path / "reads", {"s1": 10})
    arms = arms_file(tmp_path / "arms.tsv", [("fwd", "fwd", "std", 0, 0, 2, 2, 2, 12), ("std", "paired", paired, 0, 0, 2, 2, 2, 12)])
    p = run(arms, tmp_path / "out", tmp_path)
    assert p.returncode != 0 and "must run before" in p.stderr
