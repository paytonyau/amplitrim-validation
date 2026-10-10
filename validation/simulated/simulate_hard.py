"""
Simulates paired-end amplicon reads with known ground truth and three confounders. The reads are
written without primers, as AmpliTrim expects them.

  * repeat taxa  : 5 of 25 templates carry a perfect (AT)16 repeat that spans the whole
                   overlap zone -> any even shift aligns perfectly, overlap is ambiguous.
  * chimeras     : 5% of pairs come from 8 recurring chimeric templates (left part of one
                   taxon + right part of another). They merge cleanly but are not real taxa.
  * indels       : 1e-4 per base per read (exaggerated; real Illumina rates are lower).
  * substitution errors drawn from the quality scores, reverse reads worse (as before).

Header: @sim<sample>_<n>|<label>|<kind>   kind in {normal, repeat, chimera}
Truth : ground_truth.json  label -> {seq, kind}
"""
import gzip, json, random, sys
from pathlib import Path

COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]


def qual_profile(n, q0, q1, rng, noise=2.5):
    out = []
    for i in range(n):
        f = i / max(n - 1, 1)
        out.append(max(2, min(40, int(round(rng.gauss(q0 + (q1 - q0) * f ** 2.2, noise))))))
    return out


def mutate(seq, quals, rng, indel_rate):
    out = []
    for b, q in zip(seq, quals):
        r = rng.random()
        if r < indel_rate / 2:            # deletion
            continue
        if r < indel_rate:                # insertion
            out.append(rng.choice("ACGT"))
        if rng.random() < 10 ** (-q / 10):
            b = rng.choice([x for x in "ACGT" if x != b])
        out.append(b)
    return "".join(out)


def make_templates(rng, n_taxa=25, base=428, jitter=6, n_repeat=5, n_chim=8):
    t = {}
    for i in range(n_taxa):
        L = base + rng.randint(-jitter, jitter)
        s = [rng.choice("ACGT") for _ in range(L)]
        kind = "normal"
        if i < n_repeat:                  # perfect dinucleotide repeat across the overlap zone
            s[200:232] = list("AT" * 16)
            kind = "repeat"
        t[f"taxon{i:02d}"] = {"seq": "".join(s), "kind": kind}
    normals = [k for k, v in t.items() if v["kind"] == "normal"]
    for j in range(n_chim):
        a, b = rng.sample(normals, 2)
        bp = rng.randint(60, 360)
        t[f"chim{j}"] = {"seq": t[a]["seq"][:bp] + t[b]["seq"][bp:], "kind": "chimera"}
    return t


def simulate(out_dir, n_samples=3, pairs=4000, read_len=248, chim_frac=0.05, indel=1e-4, seed=11,
             pf_len=17, pr_len=21):
    rng = random.Random(seed)
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    pf = "".join(rng.choice("ACGT") for _ in range(pf_len))
    pr = "".join(rng.choice("ACGT") for _ in range(pr_len))
    T = make_templates(rng)
    chims = [k for k, v in T.items() if v["kind"] == "chimera"]
    reals = [k for k, v in T.items() if v["kind"] != "chimera"]
    for s in range(1, n_samples + 1):
        with gzip.open(out / f"sample{s}_S{s}_L001_R1_001.fastq.gz", "wt", newline="\n") as f1, \
             gzip.open(out / f"sample{s}_S{s}_L001_R2_001.fastq.gz", "wt", newline="\n") as f2:
            for n in range(pairs):
                lab = rng.choice(chims) if rng.random() < chim_frac else rng.choice(reals)
                ins = T[lab]["seq"]
                qf = qual_profile(read_len + 8, 38, 31, rng)
                qr = qual_profile(read_len + 8, 37, 20, rng)
                fwd = mutate((pf + ins)[:read_len + 8], qf, rng, indel)[:read_len]
                rev = mutate((pr + rc(ins))[:read_len + 8], qr, rng, indel)[:read_len]
                qf, qr = qf[:read_len], qr[:read_len]
                h = f"@sim{s}_{n}|{lab}|{T[lab]['kind']}"
                f1.write(f"{h} 1\n{fwd[pf_len:]}\n+\n{''.join(chr(q + 33) for q in qf[pf_len:])}\n")
                f2.write(f"{h} 2\n{rev[pr_len:]}\n+\n{''.join(chr(q + 33) for q in qr[pr_len:])}\n")
    (out / "ground_truth.json").write_text(json.dumps(T))
    L = sorted(len(v["seq"]) for v in T.values())
    print(f"templates {len(T)} (repeat={sum(v['kind']=='repeat' for v in T.values())}, "
          f"chimera={len(chims)}); insert length {L[0]}-{L[-1]}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python simulate_hard.py NEW_OUTPUT_FOLDER")
    target = Path(sys.argv[1])
    if target.exists() and any(target.iterdir()):
        raise SystemExit(f"{target} is not empty. Give a new or empty folder: the simulator writes its own reads "
                         f"there (six FASTQ files and ground_truth.json), and they would mix with whatever is already "
                         f"in it.")
    simulate(sys.argv[1])
