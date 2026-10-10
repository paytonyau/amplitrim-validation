"""Simulated variable-length community with planted chimeras: 40 taxa of 120-360 nt (ITS-like, no primers), four samples of 9,000 pairs, 8% of the pairs from 12 recurring
chimeric templates (the start of one abundant taxon joined to the end of another). The test of whether DADA2's chimera removal, as the validation scripts run it, finds
chimeras in reads of variable length and without truncation (it removed 6.5% of the merged reads when 7.8% were planted).

  python validation/simulated/chimera_sim.py OUT_DIR

Chimeras that occur once are unique sequences that DADA2 cannot confirm and drops at denoising; only recurring ones reach chimera removal, which is why these recur.
Writes OUT_DIR/reads/*.fastq.gz.
"""
import gzip
import os
import sys
from pathlib import Path

import numpy as np

OUT = Path(sys.argv[1])
(OUT / "reads").mkdir(parents=True, exist_ok=True)
os.chdir(OUT)
rng = np.random.default_rng(23)
B = np.array(list("ACGT")); comp = np.array([3, 2, 1, 0]); rc = lambda a: comp[a[::-1]]
N = 40; L = np.clip(np.round(rng.normal(230, 60, N)).astype(int), 120, 360)
taxa = [rng.integers(0, 4, l) for l in L]
ab = rng.lognormal(0, 1.0, N); ab /= ab.sum()
top = list(np.argsort(-ab)[:14])                                       # chimeras form between abundant parents, at recurring junctions
types = []
for _ in range(12):
    a, b = rng.choice(top, 2, replace=False)
    types.append(np.concatenate([taxa[a][:int(rng.integers(50, len(taxa[a]) - 50))], taxa[b][int(rng.integers(50, len(taxa[b]) - 50)):]]))


def posq(n, pts):
    xs, ys = zip(*pts)
    return np.interp(np.arange(n), xs, ys)


Q1, Q2 = posq(251, [(0, 37), (150, 34), (251, 27)]), posq(251, [(0, 36), (100, 31), (251, 17)])


def make(seqs, q0):
    S, Q = [], []
    for s in seqs:
        q = np.clip(rng.normal(q0[:len(s)] + rng.normal(0, 3.0), 3.0), 3, 40)
        q = np.round(q).astype(int)
        e = rng.random(len(s)) < 10.0 ** (-q / 10)
        s = s.copy()
        s[e] = (s[e] + rng.integers(1, 4, e.sum())) % 4
        S.append("".join(B[s]))
        Q.append("".join(chr(33 + int(x)) for x in q))
    return S, Q


total = chim = 0
for name in ("S1", "S2", "S3", "S4"):
    n = 9000
    frags = []
    for t in rng.choice(N, n, p=ab):
        if rng.random() < 0.08:
            frags.append(types[int(rng.integers(12))])
            chim += 1
        else:
            frags.append(taxa[t])
        total += 1
    s1, q1 = make([f[:251] for f in frags], Q1)
    s2, q2 = make([rc(f)[:251] for f in frags], Q2)
    with gzip.open(f"reads/{name}_R1.fastq.gz", "wt", newline="\n") as o1, gzip.open(f"reads/{name}_R2.fastq.gz", "wt", newline="\n") as o2:
        for i in range(n):
            o1.write(f"@{name}_{i}\n{s1[i]}\n+\n{q1[i]}\n")
            o2.write(f"@{name}_{i}\n{s2[i]}\n+\n{q2[i]}\n")
print(f"{total} pairs, {chim} ({chim / total:.1%}) from 12 recurring chimera types; true taxa: {N}")
