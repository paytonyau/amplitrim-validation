"""Simulated community with a spread of amplicon lengths: the tests of the coverage rule, of the "beaten by" rule and of the comparison with FIGARO.

  python validation/simulated/coverage_sim.py OUT_DIR [PROFILE]

PROFILE noisy_reverse (default; seed 41): the forward read stays good (Q38 falling to Q31) and the reverse read degrades (Q36 to Q13). Used for the coverage rule:
  a taxon can merge only if it is at most total - 12 nt long (real DADA2 at totals 511, 505, 480 and 470: the longest taxon recovered was total - 12 to the nucleotide).
PROFILE both_degrade (seed 43): the forward read degrades too (Q38 to Q17). Used for the "beaten by" rule and the comparison with FIGARO at equal maxEE.

120 taxa of 430-525 nt (median about 472), abundances log-normal, two samples of 14,000 pairs of 2 x 280/283 nt reads without primers. Writes OUT_DIR/reads/*.fastq.gz and
OUT_DIR/truth.fasta (every taxon with its length in the label).
"""
import os
import sys
from pathlib import Path

PROFILES = {"noisy_reverse": (41, "[(0, 38), (120, 37), (200, 35), (280, 31)]"), "both_degrade": (43, "[(0, 38), (100, 35), (180, 29), (240, 23), (280, 17)]")}
OUT = Path(sys.argv[1])
PROFILE = sys.argv[2] if len(sys.argv) > 2 else "noisy_reverse"
if PROFILE not in PROFILES:
    raise SystemExit(f"PROFILE must be one of {', '.join(PROFILES)}")
SEED = PROFILES[PROFILE][0]
Q1_POINTS = eval(PROFILES[PROFILE][1])
(OUT / "reads").mkdir(parents=True, exist_ok=True)
os.chdir(OUT)
import gzip, numpy as np
rng = np.random.default_rng(SEED)
RL1, RL2 = 280, 283
comp = np.array([3, 2, 1, 0]); B = np.array(list("ACGT"))
def posq(n, pts):
    xs, ys = zip(*pts); return np.interp(np.arange(n), xs, ys)
Q1 = posq(RL1, Q1_POINTS)
Q2 = posq(RL2, [(0, 36), (60, 32), (120, 27), (180, 22), (240, 17), (283, 13)])
N = 120
lens = np.clip(np.round(rng.normal(472, 16, N)).astype(int), 430, 525)
taxa = [rng.integers(0, 4, L) for L in lens]
ab = rng.lognormal(0, 1.2, N); ab /= ab.sum()
open("truth.fasta", "w").write("".join(f">taxon{i}_len{len(t)}\n{''.join(B[t])}\n" for i, t in enumerate(taxa)))
def make(seqs, meanq):
    n = len(seqs)
    q = np.clip(rng.normal(meanq + rng.normal(0, 4.0, (n, 1)), 3.0, size=(n, len(meanq))), 3, 40)
    d = rng.random(q.shape) < 0.002; q[d] = rng.integers(4, 14, d.sum()); q = np.round(q).astype(int)
    err = rng.random(q.shape) < 10.0 ** (-q / 10)
    arr = np.array(seqs); arr[err] = (arr[err] + rng.integers(1, 4, err.sum())) % 4
    return ["".join(B[r]) for r in arr], ["".join(chr(33 + int(x)) for x in r) for r in q]
for sample in ("s1", "s2"):
    n = 14000; picks = rng.choice(N, n, p=ab)
    f = np.array([taxa[t][:RL1] for t in picks]); r = np.array([comp[taxa[t][::-1]][:RL2] for t in picks])
    s1, q1 = make(f, Q1); s2, q2 = make(r, Q2)
    with gzip.open(f"reads/{sample}_R1.fastq.gz", "wt", newline="\n") as o1, gzip.open(f"reads/{sample}_R2.fastq.gz", "wt", newline="\n") as o2:
        for i in range(n):
            o1.write(f"@{sample}_{i}\n{s1[i]}\n+\n{q1[i]}\n"); o2.write(f"@{sample}_{i}\n{s2[i]}\n+\n{q2[i]}\n")
sh = sorted(lens)
print(f"{N} taxa, lengths {sh[0]}-{sh[-1]} (median {sh[N//2]}); 2 x 14,000 pairs")
