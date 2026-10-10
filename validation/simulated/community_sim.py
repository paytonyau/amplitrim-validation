#!/usr/bin/env python
"""Simulated amplicon communities with known truth, for the richness and the over-merging tests (validation/round4_steps.py).

  python validation/simulated/community_sim.py OUT --taxa 250 --seed 1 [--samples 10 --pairs 5000 --qual-from READY_FOLDER]
  python validation/simulated/community_sim.py OUT --taxa 25 --repeat-len 10 --chimera 0.01 --indel 1e-5 --samples 1 --pairs 12000

What it makes (numpy only; no AmpliTrim code except to read real quality strings with --qual-from):
  * taxa related like real ones: families 5 to 10% apart (substitutions, and short indels in most families) and strains within 1 to 3 bases of their family;
  * a log-normal abundance distribution (--even for an even community), shared by the samples with some noise and absences, as in a study;
  * paired reads whose quality strings are resampled from real reads (--qual-from) or follow a simple profile (a better forward read and a worse reverse read),
    with substitution errors exactly as the quality scores say;
  * optional stress: a perfect AT repeat across the middle of the amplicon in a share of the taxa (--repeat-len), chimeras (--chimera, a share of the pairs),
    and indels in the reads (--indel, per base).
Every read header carries its true template and kind: @sim<sample>_<n>|<template>|<normal|repeat|chimera>, with |indel added when a read of the pair carries an indel.
Writes S<i>_R1.fastq.gz / S<i>_R2.fastq.gz, truth.fasta (>T0001|T0001|kind) and truth.json (lengths, kinds, reads per template).
"""
import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

LETTERS = np.frombuffer(b"ACGT", dtype=np.uint8)
_LUT = np.zeros(256, dtype=np.uint8)
for _i, _ch in enumerate(b"ACGT"):
    _LUT[_ch] = _i


def encode(seq: str) -> np.ndarray:
    return _LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]


def decode(codes: np.ndarray) -> str:
    return LETTERS[codes].tobytes().decode()


def revcomp_codes(codes: np.ndarray) -> np.ndarray:
    return (3 - codes)[::-1]


def _mutate(seq: np.ndarray, rng, rate: float) -> np.ndarray:
    out = seq.copy()
    hit = rng.random(len(out)) < rate
    out[hit] = (out[hit] + rng.integers(1, 4, int(hit.sum()))) % 4
    return out


def _indel(seq: np.ndarray, rng) -> np.ndarray:
    size, pos = int(rng.integers(1, 5)), int(rng.integers(50, len(seq) - 50))
    if rng.random() < 0.5:
        return np.concatenate([seq[:pos], rng.integers(0, 4, size).astype(np.uint8), seq[pos:]])
    return np.concatenate([seq[:pos], seq[pos + size:]])


def make_taxa(n_taxa: int, rng, length: int = 465, family_size: int = 5, family_div: float = 0.075, strain_edits=(1, 3), repeat_len: int = 0,
              repeat_frac: float = 0.2, max_length: int = 478) -> list:
    """n_taxa distinct templates: families derived from one ancestor, strains within 1 to 3 substitutions of their family; optionally a perfect AT repeat in the middle of a share of them."""
    n_fam = max(1, -(-n_taxa // family_size))
    ancestor = rng.integers(0, 4, length).astype(np.uint8)
    fams = []
    for _ in range(n_fam):
        f = _mutate(ancestor, rng, family_div)
        if rng.random() < 0.6:
            f = _indel(f, rng)
        if rng.random() < 0.3:
            f = _indel(f, rng)
        fams.append(f)
    taxa, seen = [], set()
    for i in range(n_taxa):
        fam = i % n_fam
        s = fams[fam].copy()
        edits = int(rng.integers(strain_edits[0], strain_edits[1] + 1))
        pos = rng.choice(len(s), edits, replace=False)
        s[pos] = (s[pos] + rng.integers(1, 4, edits)) % 4
        while s.tobytes() in seen:
            p = int(rng.integers(len(s)))
            s[p] = (s[p] + int(rng.integers(1, 4))) % 4
        seen.add(s.tobytes())
        taxa.append({"id": f"T{i + 1:04d}", "family": fam, "codes": s, "kind": "normal"})
    if repeat_len > 0:
        for i in rng.choice(n_taxa, max(1, int(round(repeat_frac * n_taxa))), replace=False):
            s = taxa[i]["codes"].copy()
            mid, half = len(s) // 2, repeat_len // 2
            s[mid - half:mid - half + repeat_len] = np.array([0, 3] * ((repeat_len + 1) // 2), dtype=np.uint8)[:repeat_len]
            taxa[i]["codes"], taxa[i]["kind"] = s, "repeat"
    for t in taxa:
        t["codes"] = t["codes"][:max_length]
    return taxa


def abundances(n: int, rng, sigma: float = 1.5, even: bool = False) -> np.ndarray:
    w = np.ones(n) if even else rng.lognormal(0.0, sigma, n)
    return w / w.sum()


def sample_profiles(base: np.ndarray, n_samples: int, rng, presence: float = 0.8, noise: float = 0.5) -> list:
    """The abundance of every taxon in each sample: the shared abundance, with noise and some absences, as in a study."""
    out = []
    for _ in range(n_samples):
        w = base * rng.lognormal(0.0, noise, len(base)) * (rng.random(len(base)) < presence)
        out.append((base if w.sum() == 0 else w) / (base if w.sum() == 0 else w).sum())
    return out


def parametric_pool(rng, n: int, read_len: int, drop_f: float = 9.0, drop_r: float = 18.0):
    """Quality strings (as Phred arrays) of a simple profile: the forward read stays better than the reverse read."""
    pos = (np.arange(read_len) / read_len) ** 1.5
    pools = []
    for drop in (drop_f, drop_r):
        q = 37.0 - drop * pos[None, :] + rng.normal(0, 1.5, (n, 1)) + rng.normal(0, 2.0, (n, read_len))
        pools.append(np.clip(np.round(q), 2, 40).astype(np.uint8))
    return pools[0], pools[1]


def load_quality_pool(reads_dir, read_len: int, n: int = 4000, per_sample: int = 500):
    """Quality strings of real reads (the first reads of up to twelve samples), padded or cut to read_len; two matrices (forward, reverse)."""
    from amplitrim.subsample import find_pairs, read_fastq
    pairs = find_pairs(Path(reads_dir))
    chosen = pairs[::max(1, len(pairs) // 12)][:12]
    qf, qr = [], []
    for _name, r1, r2 in chosen:
        for k, (a, b) in enumerate(zip(read_fastq(r1), read_fastq(r2))):
            if k >= per_sample:
                break
            qf.append(a[3])
            qr.append(b[3])
    if not qf:
        raise ValueError(f"no reads found in {reads_dir}")

    def matrix(quals):
        quals = quals[:n]
        m = np.zeros((len(quals), read_len), dtype=np.uint8)
        for i, q in enumerate(quals):
            a = np.frombuffer(q.encode(), dtype=np.uint8) - 33
            a = a[:read_len]
            m[i, :len(a)] = a
            m[i, len(a):] = a[-1] if len(a) else 2
        return m
    return matrix(qf), matrix(qr)


def _edit_read(codes: np.ndarray, qual: np.ndarray, rng):
    """One insertion or deletion in a read (the quality string follows the base)."""
    pos = int(rng.integers(5, len(codes) - 5))
    if rng.random() < 0.5:
        return np.insert(codes, pos, int(rng.integers(0, 4))).astype(np.uint8), np.insert(qual, pos, qual[pos]).astype(np.uint8)
    return np.delete(codes, pos), np.delete(qual, pos)


def simulate_sample(taxa: list, profile: np.ndarray, n_pairs: int, rng, read_len: int, pool_f: np.ndarray, pool_r: np.ndarray, chimera_frac: float, indel_rate: float, tag: str):
    """(headers, forward reads, reverse reads, forward qualities, reverse qualities, template of each pair) of one sample; reads are bytes."""
    k = len(taxa)
    idx = rng.choice(k, n_pairs, p=profile)
    R1 = np.stack([t["codes"][:read_len] for t in taxa])
    R2 = np.stack([revcomp_codes(t["codes"])[:read_len] for t in taxa])
    r1, r2 = R1[idx].copy(), R2[idx].copy()
    labels, kinds = [taxa[i]["id"] for i in idx], [taxa[i]["kind"] for i in idx]
    for n in np.flatnonzero(rng.random(n_pairs) < chimera_frac):
        a, b = rng.choice(k, 2, replace=False)
        ta, tb = taxa[a]["codes"], taxa[b]["codes"]
        bp = int(rng.integers(int(0.25 * min(len(ta), len(tb))), int(0.75 * min(len(ta), len(tb)))))
        chim = np.concatenate([ta[:bp], tb[bp:]])
        r1[n], r2[n] = chim[:read_len], revcomp_codes(chim)[:read_len]
        labels[n], kinds[n] = f"{taxa[a]['id']}x{taxa[b]['id']}", "chimera"
    qf, qr = pool_f[rng.integers(len(pool_f), size=n_pairs)], pool_r[rng.integers(len(pool_r), size=n_pairs)]
    for arr, q in ((r1, qf), (r2, qr)):
        err = rng.random(arr.shape) < 10.0 ** (-q.astype(float) / 10.0)
        arr[err] = (arr[err] + rng.integers(1, 4, int(err.sum()))) % 4
    seqs1, seqs2 = [LETTERS[r].tobytes() for r in r1], [LETTERS[r].tobytes() for r in r2]
    quals1, quals2 = [(q + 33).astype(np.uint8).tobytes() for q in qf], [(q + 33).astype(np.uint8).tobytes() for q in qr]
    indel_pairs = set()
    if indel_rate > 0:
        for which, arr, q, seqs, quals in ((1, r1, qf, seqs1, quals1), (2, r2, qr, seqs2, quals2)):
            for n in np.flatnonzero(rng.random(n_pairs) < indel_rate * read_len):
                c, qq = _edit_read(arr[n], q[n], rng)
                seqs[n], quals[n] = LETTERS[c].tobytes(), (qq + 33).astype(np.uint8).tobytes()
                indel_pairs.add(int(n))
    headers = [f"@sim{tag}_{n}|{labels[n]}|{kinds[n]}" + ("|indel" if n in indel_pairs else "") for n in range(n_pairs)]
    return headers, seqs1, seqs2, quals1, quals2, labels


def _write_fastq(path: Path, headers, seqs, quals) -> None:
    with gzip.open(path, "wb") as fh:
        for h, s, q in zip(headers, seqs, quals):
            fh.write(h.encode() + b"\n" + s + b"\n+\n" + q + b"\n")


def write_community(out_dir, n_taxa: int, seed: int, n_samples: int = 10, pairs_per_sample: int = 5000, read_len: int = 250, qual_from=None, sigma: float = 1.5,
                    even: bool = False, repeat_len: int = 0, chimera_frac: float = 0.0, indel_rate: float = 0.0, max_length: int = 478) -> dict:
    """Make one community and write its reads, truth.fasta and truth.json; returns the information the tests need."""
    out_dir = Path(out_dir)
    if out_dir.exists() and any(out_dir.iterdir()):
        raise SystemExit(f"{out_dir} is not empty: use a new folder")
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    taxa = make_taxa(n_taxa, rng, repeat_len=repeat_len, max_length=max_length)
    base = abundances(n_taxa, rng, sigma, even)
    profiles = sample_profiles(base, n_samples, rng)
    if qual_from:
        pool_f, pool_r = load_quality_pool(qual_from, read_len)
    else:
        pool_f, pool_r = parametric_pool(rng, 2000, read_len)
    counts = {t["id"]: 0 for t in taxa}
    for s, profile in enumerate(profiles, 1):
        headers, s1, s2, q1, q2, labels = simulate_sample(taxa, profile, pairs_per_sample, rng, read_len, pool_f, pool_r, chimera_frac, indel_rate, f"S{s}")
        _write_fastq(out_dir / f"S{s}_R1.fastq.gz", headers, s1, q1)
        _write_fastq(out_dir / f"S{s}_R2.fastq.gz", headers, s2, q2)
        for lab in labels:
            if lab in counts:
                counts[lab] += 1
    (out_dir / "truth.fasta").write_text("".join(f">{t['id']}|{t['id']}|{t['kind']}\n{decode(t['codes'])}\n" for t in taxa), encoding="utf-8")
    info = {"n_taxa": n_taxa, "seed": seed, "read_len": read_len, "lengths": {t["id"]: int(len(t["codes"])) for t in taxa}, "kinds": {t["id"]: t["kind"] for t in taxa},
            "counts": counts, "max_len": int(max(len(t["codes"]) for t in taxa)), "truth_fasta": str(out_dir / "truth.fasta")}
    (out_dir / "truth.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    return info


def abundance_class(reads: int) -> str:
    """The share of the truth DADA2 can be asked to recover: taxa with fewer than two reads cannot be, and are left out."""
    return "abundant" if reads >= 50 else ("medium" if reads >= 10 else ("rare" if reads >= 2 else "unseen"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--taxa", type=int, default=25)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--samples", type=int, default=10)
    ap.add_argument("--pairs", type=int, default=5000, help="pairs per sample")
    ap.add_argument("--read-len", type=int, default=250)
    ap.add_argument("--qual-from", help="a folder of real primer-free reads whose quality strings are resampled")
    ap.add_argument("--even", action="store_true")
    ap.add_argument("--repeat-len", type=int, default=0)
    ap.add_argument("--chimera", type=float, default=0.0)
    ap.add_argument("--indel", type=float, default=0.0)
    args = ap.parse_args(argv)
    info = write_community(args.out, args.taxa, args.seed, args.samples, args.pairs, args.read_len, args.qual_from, even=args.even, repeat_len=args.repeat_len,
                           chimera_frac=args.chimera, indel_rate=args.indel)
    classes = {}
    for c in info["counts"].values():
        classes[abundance_class(c)] = classes.get(abundance_class(c), 0) + 1
    print(f"{info['n_taxa']} taxa (longest {info['max_len']} nt); taxa by reads: {classes}; written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
