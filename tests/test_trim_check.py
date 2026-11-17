"""validation/trim_check.py: the comparison and the verdicts, and, where AmpliTrim 0.3 is available (AMPLITRIM_REPO), the whole check against an
independent computation of DADA2's filter on reads with poor ends."""
import gzip
import importlib.util
import os
import random
import re
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "validation"))
spec = importlib.util.spec_from_file_location("trim_check", HERE.parent / "validation" / "trim_check.py")
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)


def samples(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("sample\tinput\tfiltered\tdenoised\tmerged\tnonchim\tasvs\n" + "".join(f"{s}\t{i}\t{f}\t{f}\t{m}\t{m}\t3\n" for s, i, f, m in rows))


def test_dada2s_shares_are_summed_over_the_samples(tmp_path):
    samples(tmp_path / "s.tsv", [("a", 1000, 900, 800), ("b", 3000, 1500, 1200)])
    d = T.dada2_share(tmp_path / "s.tsv")
    assert d["input"] == 4000 and d["pass"] == pytest.approx(0.6) and d["merged_share"] == pytest.approx(0.5)


def test_the_verdicts():
    row = lambda p_trim, d_trim, p_def, d_def: {"pass_trim": p_trim, "dada2_sw15": d_trim, "pass_default": p_def, "dada2_std": d_def}
    good = {lab: row(0.90, 0.905, 0.40, 0.395) for lab in T.DATASETS}
    V = dict((c, ok) for c, ok, _ in T.verdicts(good))
    assert V == {"P115": True, "P116": True}
    bad = dict(good, its2=row(0.90, 0.88, 0.40, 0.40))                                   # 2 points off after the trim
    V = {c: (ok, t) for c, ok, t in T.verdicts(bad)}
    assert V["P115"][0] is False and "+2.00 points (its2)" in V["P115"][1] and V["P116"][0] is True
    short = {lab: good[lab] for lab in T.DATASETS[:5]}
    assert all(ok is None for _, ok, _ in T.verdicts(short))                              # a dataset missing: nothing decided


def q(xs):
    return "".join(chr(33 + x) for x in xs)


def make_night(root, rng):
    """One dataset folder with reads whose reverse ends fade, and nothing else."""
    reads = root / "atcc_its1" / "reads"
    reads.mkdir(parents=True)
    pairs = []
    for i in range(1500):
        f = [37] * 250
        cut = rng.randint(60, 250)
        r = [37] * cut + [rng.choice([8, 10, 12, 20, 30]) for _ in range(250 - cut)]
        if rng.random() < 0.05:
            f[rng.randint(0, 249)] = 2
        pairs.append((f, r))
    for m in (1, 2):
        with gzip.open(reads / f"s_R{m}.fastq.gz", "wt") as fh:
            fh.write("".join(f"@r{i}\n{''.join(rng.choice('ACGT') for _ in p[m - 1])}\n+\n{q(p[m - 1])}\n" for i, p in enumerate(pairs)))
    return pairs


def sliding_cut(qual, window=4, q=15):
    """The cut of the validation (validation/its_night.py), written out again here so that the expectation does not come from the code under test."""
    s = [ord(c) - 33 for c in qual]
    for i in range(0, max(0, len(s) - window + 1)):
        if sum(s[i:i + window]) / window < q:
            return i
    return len(s)


def dada2_filter(quals, trim):
    """An independent DADA2 no-truncLen filter (truncQ 2, minLen 50, maxEE 2, maxN 0: the reads have no N), with or without the trim first."""
    def one(x):
        if trim:
            x = x[:sliding_cut(q(x))]
        k = next((i for i, v in enumerate(x) if v <= 2), len(x))
        x = x[:k]
        return len(x) >= 50 and sum(10 ** (-v / 10) for v in x) <= 2
    return sum(one(f) and one(r) for f, r in quals)


@pytest.mark.skipif(not os.environ.get("AMPLITRIM_REPO"), reason="needs AmpliTrim 0.3 or later: set AMPLITRIM_REPO to its repository folder")
def test_amplitrims_prediction_equals_an_independent_dada2_filter(tmp_path):
    rng = random.Random(3)
    pairs = make_night(tmp_path / "night", rng)
    n = len(pairs)
    for arm, trim in (("std", False), ("sw15", True)):
        kept = dada2_filter(pairs, trim)
        samples(tmp_path / "results" / "atcc_its1" / "dada2" / arm / "samples.tsv", [("s", n, kept, kept)])
    saved_modules = {m: sys.modules[m] for m in list(sys.modules) if m == "amplitrim" or m.startswith("amplitrim.")}
    saved_path, saved_datasets = list(sys.path), list(T.DATASETS)
    try:
        repo = Path(os.environ["AMPLITRIM_REPO"])
        held = re.search(r'__version__ = "(\d+\.\d+)\.', (repo / "amplitrim" / "__init__.py").read_text(encoding="utf-8")).group(1)
        version, pipeline = T.load_amplitrim(repo, want=held + ".")         # the version the repository holds (0.3 for run 30; 0.4 and later since)
        T.DATASETS[:] = ["atcc_its1"]
        rc = T.main(["--night", str(tmp_path / "night"), "--results", str(tmp_path / "results"), "--out", str(tmp_path / "out")], pipeline=pipeline)
    finally:                                                                      # the other tests keep the AmpliTrim they had
        for m in [m for m in sys.modules if m == "amplitrim" or m.startswith("amplitrim.")]:
            del sys.modules[m]
        sys.modules.update(saved_modules)
        sys.path[:] = saved_path
        T.DATASETS[:] = saved_datasets
    s = (tmp_path / "out" / "TRIM_CHECK_SUMMARY.md").read_text()
    assert rc == 0 and "**P115** MET" in s and "the largest +0.00 points" in s and "**P116** MET" in s
