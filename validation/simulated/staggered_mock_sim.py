"""Simulated staggered mock community (15 strains: 5 x 18.02%, 5 x 1.80%, 5 x 0.18%; three pairs are close relatives that differ by 4-7 bases) whose reads degrade
along their length, the reverse read faster: the test of what a looser maxEE keeps (relaxing the reverse limit from 2 to 5 raised the share passing from 78% to 96%).

  python validation/simulated/staggered_mock_sim.py OUT_DIR

Writes OUT_DIR/reads/*.fastq.gz (no primers) and OUT_DIR/truth.fasta.
"""
import os
import sys
from pathlib import Path

OUT = Path(sys.argv[1])
(OUT / "reads").mkdir(parents=True, exist_ok=True)
os.chdir(OUT)
import gzip
import numpy as np
rng = np.random.default_rng(5)
RL1, RL2 = 280, 283
comp = np.array([3, 2, 1, 0])
def posq(n, pts):                                     # mean quality by position, piecewise linear
    xs, ys = zip(*pts); return np.interp(np.arange(n), xs, ys)
Q1 = posq(RL1, [(0, 38), (120, 37), (200, 35), (280, 31)])
Q2 = posq(RL2, [(0, 36), (60, 32), (120, 27), (180, 22), (240, 17), (283, 13)])
# 15 strains (a staggered mock): 5 x 18.02%, 5 x 1.80%, 5 x 0.18%; three pairs are close relatives (a few differences)
taxa = [rng.integers(0, 4, int(rng.integers(455, 496))) for _ in range(15)]
for i, j in ((5, 0), (6, 1), (7, 2)):                 # replace strains 5-7 by variants of 0-2 with 4-7 substitutions
    v = taxa[j].copy(); pos = rng.choice(len(v), int(rng.integers(4, 8)), replace=False); v[pos] = (v[pos] + rng.integers(1, 4, len(pos))) % 4; taxa[i] = v
ab = np.array([18.02] * 5 + [1.80] * 5 + [0.18] * 5); ab /= ab.sum()
B = np.array(list("ACGT"))
open("truth.fasta", "w").write("".join(f">strain{i+1}\n{''.join(B[t])}\n" for i, t in enumerate(taxa)))
def make(seqs, meanq):
    n = len(seqs)
    q = np.clip(rng.normal(meanq + rng.normal(0, 4.0, (n, 1)), 3.0, size=(n, len(meanq))), 3, 40)      # reads differ in overall quality
    dips = rng.random(q.shape) < 0.002; q[dips] = rng.integers(4, 14, dips.sum())          # short bad stretches
    q = np.round(q).astype(int)
    err = rng.random(q.shape) < 10.0 ** (-q / 10)                                         # errors that match the qualities
    arr = np.array(seqs); arr[err] = (arr[err] + rng.integers(1, 4, err.sum())) % 4
    return ["".join(B[r]) for r in arr], ["".join(chr(33 + int(x)) for x in r) for r in q]
for sample in ("s1", "s2", "s3"):
    n = 12000; picks = rng.choice(15, n, p=ab)
    f = np.array([taxa[t][:RL1] for t in picks]); r = np.array([comp[taxa[t][::-1]][:RL2] for t in picks])
    s1, q1 = make(f, Q1); s2, q2 = make(r, Q2)
    with gzip.open(f"reads/{sample}_R1.fastq.gz", "wt", newline="\n") as o1, gzip.open(f"reads/{sample}_R2.fastq.gz", "wt", newline="\n") as o2:
        for i in range(n):
            o1.write(f"@{sample}_{i}\n{s1[i]}\n+\n{q1[i]}\n"); o2.write(f"@{sample}_{i}\n{s2[i]}\n+\n{q2[i]}\n")
print("generated 3 x 20,000 pairs; amplicon lengths", sorted(len(t) for t in taxa)[:3], "...", sorted(len(t) for t in taxa)[-3:])
