"""Simulated ITS-like community: 60 taxa whose amplicons differ in length (fragments of 80-380 nt between the ITS1-Fl2 and ITS2 primers, so many are shorter than the
250 nt reads and both primers are present as read-through), five samples of 9,000 pairs and a negative control of 80. Used for the ITS mode: the table of maxEE, the
mean cut length of the kept reads and their expected errors against DADA2's filterAndTrim with no truncLen.

  python validation/simulated/its_sim.py OUT_DIR

Writes OUT_DIR/raw/*.fastq.gz (WITH primers, as a sequencer delivers them).
"""
import os
import sys
from pathlib import Path

OUT = Path(sys.argv[1])
(OUT / "raw").mkdir(parents=True, exist_ok=True)
os.chdir(OUT)
import gzip
import numpy as np
rng = np.random.default_rng(7)
B = np.array(list("ACGT")); comp = np.array([3, 2, 1, 0])
enc = lambda s: np.array(["ACGT".index(c) for c in s])
FWD, REV = enc("CTTGGTCATTTAGAGGAAGTAA"), enc("GCTGCGTTCTTCATCGATGC")
def rc(a): return comp[a[::-1]]
N = 60
L = np.clip(np.round(rng.normal(220, 70, N)).astype(int), 80, 380)
frags = [np.concatenate([FWD, rng.integers(0, 4, l), rc(REV)]) for l in L]
ab = rng.lognormal(0, 1.1, N); ab /= ab.sum()
def posq(n, pts): xs, ys = zip(*pts); return np.interp(np.arange(n), xs, ys)
def make(seqs, q0):
    out_s, out_q = [], []
    for s in seqs:
        q = np.clip(rng.normal(q0[:len(s)] + rng.normal(0, 3.5), 3.0), 3, 40); d = rng.random(len(s)) < 0.002; q[d] = rng.integers(4, 14, d.sum()); q = np.round(q).astype(int)
        err = rng.random(len(s)) < 10.0 ** (-q / 10); s = s.copy(); s[err] = (s[err] + rng.integers(1, 4, err.sum())) % 4
        out_s.append("".join(B[s])); out_q.append("".join(chr(33 + int(x)) for x in q))
    return out_s, out_q
Q1, Q2 = posq(250, [(0, 38), (150, 35), (250, 27)]), posq(250, [(0, 36), (100, 31), (200, 22), (250, 14)])
for name, n in (("S1", 9000), ("S2", 9000), ("S3", 9000), ("S4", 9000), ("S5", 9000), ("neg", 80)):
    picks = rng.choice(N, n, p=ab)
    r1 = [frags[t][:250] for t in picks]; r2 = [rc(frags[t])[:250] for t in picks]      # fragments shorter than 250 give short reads, both primers present
    s1, q1 = make(r1, Q1); s2, q2 = make(r2, Q2)
    with gzip.open(f"raw/{name}_R1.fastq.gz", "wt", newline="\n") as o1, gzip.open(f"raw/{name}_R2.fastq.gz", "wt", newline="\n") as o2:
        for i in range(n):
            o1.write(f"@{name}_{i}\n{s1[i]}\n+\n{q1[i]}\n"); o2.write(f"@{name}_{i}\n{s2[i]}\n+\n{q2[i]}\n")
fl = sorted(len(f) for f in frags); print(f"{N} taxa; fragment lengths with primers {fl[0]}-{fl[-1]} (median {fl[N//2]}); {sum(1 for f in fl if f < 250)} of {N} taxa are shorter than the 250 nt reads")
