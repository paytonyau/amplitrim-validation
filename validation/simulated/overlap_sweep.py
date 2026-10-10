"""
Overlap-regime sweep. For each condition: simulate substitution-error paired reads from known
templates, apply the tool's DADA2-style maxEE<=2 filter, then run the package's own wrappers
(vsearch tolerant and strict, min_overlap 12) and score merges against truth.
DADA2 is run separately (dada2_check.R) on the same files as the real-pipeline reference.

Simplifications: substitution errors only, no indels/chimeras/repeats, and the
2x300 quality profile is the 2x250 one stretched (real 2x300 tails are worse). Inserts are always
longer than the reads, so read-through into adapters is NOT modelled (short ITS not covered).
"""
import gzip, random, subprocess, sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rtools import find_rscript                                          # noqa: E402
from simulate_hard import qual_profile, mutate, rc                       # noqa: E402
from amplitrim.merge_wrapper import run_vsearch                          # noqa: E402

BASE = Path(sys.argv[2] if len(sys.argv) > 2 else "sweep_out"); PF, PR, N = 17, 21, 6000
# name, read length, insert-length spec, truncation (tf, tr) of primer-trimmed reads or None
CONDS = [
    ("V3V4_2x250_ov12", 250, ("fix", 428), (214, 226)),
    ("V3V4_2x250_ov16", 250, ("fix", 428), (218, 226)),
    ("V3V4_2x250_ov20", 250, ("fix", 428), (222, 226)),
    ("V3V4_2x250_ov34", 250, ("fix", 428), None),
    ("V1V3_2x300_ov62", 300, ("fix", 500), None),
    ("V3V4_2x300_ov134", 300, ("fix", 428), None),
    ("V1V2_2x250_ov162", 250, ("fix", 300), None),
    ("V4_2x250_ov209", 250, ("fix", 253), None),
    ("ITS_2x300_var", 300, ("var", 300, 520), None),
    ("ITS_2x250_var", 250, ("var", 300, 520), None),
]


ee = lambda q: sum(10 ** (-(ord(c) - 33) / 10) for c in q)


def build(name, R, spec, trunc, seed):
    rng = random.Random(seed); d = BASE / name; d.mkdir(parents=True, exist_ok=True)
    pf = "".join(rng.choice("ACGT") for _ in range(PF)); pr = "".join(rng.choice("ACGT") for _ in range(PR))
    T = {}
    for i in range(20):
        L = spec[1] + rng.randint(-6, 6) if spec[0] == "fix" else rng.randint(spec[1], spec[2])
        T[f"taxon{i:02d}"] = "".join(rng.choice("ACGT") for _ in range(L))
    labs = list(T)
    with gzip.open(d / "r1.fq.gz", "wt", newline="\n") as f1, gzip.open(d / "r2.fq.gz", "wt", newline="\n") as f2:
        for n in range(N):
            lab = rng.choice(labs); ins = T[lab]
            qf = qual_profile(R, 38, 31, rng); qr = qual_profile(R, 37, 20, rng)
            f = mutate((pf + ins)[:R], qf, rng, 0.0)[PF:]; r = mutate((pr + rc(ins))[:R], qr, rng, 0.0)[PR:]
            qf, qr = qf[PF:], qr[PR:]
            if trunc:
                f, qf, r, qr = f[:trunc[0]], qf[:trunc[0]], r[:trunc[1]], qr[:trunc[1]]
            qfs, qrs = "".join(chr(q + 33) for q in qf), "".join(chr(q + 33) for q in qr)
            f1.write(f"@sim{n}|{lab} 1\n{f}\n+\n{qfs}\n"); f2.write(f"@sim{n}|{lab} 2\n{r}\n+\n{qrs}\n")
    with open(d / "truth.tsv", "w") as t:
        t.write("label\tkind\tseq\n"); [t.write(f"{k}\tnormal\t{v}\n") for k, v in T.items()]
    return d, T, len(f), len(r)


def fq(p):
    with gzip.open(p, "rt") as fh:
        while True:
            h = fh.readline()
            if not h: return
            s = fh.readline().strip(); fh.readline(); q = fh.readline().strip(); yield h.strip(), s, q


def correct_frac(path, T):
    """fraction of merged reads whose length equals the true insert length"""
    ok = tot = 0
    if str(path).endswith(".fasta"):
        name, seq, recs = None, [], []
        for line in open(path):
            if line.startswith(">"):
                if name: recs.append((name, "".join(seq)))
                name, seq = line[1:].split()[0], []
            else: seq.append(line.strip())
        if name: recs.append((name, "".join(seq)))
    else:
        recs = [(h[1:].split()[0], s) for h, s, _ in fq(path)]
    for nm, s in recs:
        tot += 1; ok += len(s) == len(T[nm.split("|")[1]])
    return ok / tot if tot else float("nan")




def mode_engines():
    rows = []
    for i, (name, R, spec, trunc) in enumerate(CONDS):
        d, T, lf, lr = build(name, R, spec, trunc, seed=100 + i)
        # maxEE filter (the tool's Stage-2 filter) -> files given to the engines
        with gzip.open(d / "f1.fq.gz", "wt", newline="\n") as o1, gzip.open(d / "f2.fq.gz", "wt", newline="\n") as o2:
            for (h1, s1, q1), (h2, s2, q2) in zip(fq(d / "r1.fq.gz"), fq(d / "r2.fq.gz")):
                if ee(q1) <= 2 and ee(q2) <= 2:
                    o1.write(f"{h1}\n{s1}\n+\n{q1}\n"); o2.write(f"{h2}\n{s2}\n+\n{q2}\n")
        L = sorted(len(v) for v in T.values())
        ov = (lf + lr) - (spec[1] if spec[0] == "fix" else (spec[1] + spec[2]) // 2)
        vt = run_vsearch(d / "f1.fq.gz", d / "f2.fq.gz", d / "m", "x", 12, None, "tolerant")
        vs = run_vsearch(d / "f1.fq.gz", d / "f2.fq.gz", d / "m", "x", 12, 0, "strict")
        rows.append(dict(
            condition=name, read_lens=f"{lf}/{lr}", nominal_overlap=ov, insert=f"{L[0]}-{L[-1]}",
            vsearch_tolerant=vt.merged_pairs / N, vsearch_strict=vs.merged_pairs / N,
            vs_correct_len=correct_frac(vt.merged_path, T),
            trunc_f=lf, trunc_r=lr))
        print(f"{name:20s} overlap~{ov:4d}  vsearch {rows[-1]['vsearch_tolerant']:.3f}", flush=True)
    pd.DataFrame(rows).to_csv(BASE / "engines.csv", index=False)


def mode_dada2():
    """Real DADA2 on the same (unfiltered) reads; needs R + dada2. Slow: ~15 s per condition."""
    e = pd.read_csv(BASE / "engines.csv")
    for _, r in e.iterrows():
        lf, lr = r.read_lens.split("/"); c = BASE / r.condition
        subprocess.run([find_rscript() or "Rscript", str(Path(__file__).parent / "dada2_check.R"), str(c / "r1.fq.gz"), str(c / "r2.fq.gz"),
                        str(c / "truth.tsv"), str(c / "dada2.csv"), f"{lf},{lr}"], check=True, capture_output=True)
        print("dada2 done:", r.condition, flush=True)


def mode_report():
    e = pd.read_csv(BASE / "engines.csv")
    d = pd.concat([pd.read_csv(BASE / c / "dada2.csv").assign(condition=c) for c in e.condition])
    m = e.merge(d[["condition", "dada2_merged"]], on="condition")
    m["vsearch - DADA2 (pts)"] = (m.vsearch_tolerant - m.dada2_merged) * 100
    pd.set_option("display.width", 250)
    print(m[["condition", "nominal_overlap", "vsearch_strict", "vsearch_tolerant", "dada2_merged",
             "vsearch - DADA2 (pts)"]].round(3).to_string(index=False))


if __name__ == "__main__":
    {"engines": mode_engines, "dada2": mode_dada2, "report": mode_report}[sys.argv[1]]()
