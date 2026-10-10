#!/usr/bin/env python
"""Show what the reads of some samples look like, in one or more folders (for example before and after trimming).

For each sample and folder: the most common first bases of R1 and R2 with their shares, the read lengths, and a few example reads. It is
for questions like "did the trimming leave a primer behind?" or "which read starts with which primer in this sample?".

  python validation/peek_reads.py C:\\data\\raw C:\\data\\results\\study\\ready --sample ERR1 ERR2 [--n 5] [--scan 3000] [--k 10]
"""

import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def read_head(path: Path, n: int) -> list:
    """The first n (sequence, quality) pairs of a gzipped FASTQ."""
    out = []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        while len(out) < n:
            rec = [fh.readline() for _ in range(4)]
            if not rec[0]:
                break
            out.append((rec[1].strip(), rec[3].strip()))
    return out


def start_table(seqs: list, k: int = 10, top: int = 4) -> str:
    """The most common first k bases and the share of reads that start with each."""
    counts = Counter(s[:k] for s in seqs if len(s) >= k)
    n = max(sum(counts.values()), 1)
    return " | ".join(f"{kmer} {c / n:.1%}" for kmer, c in counts.most_common(top)) or "(no reads)"


def length_summary(seqs: list) -> str:
    if not seqs:
        return "(no reads)"
    lens = sorted(len(s) for s in seqs)
    mode, count = Counter(lens).most_common(1)[0]
    return f"min {lens[0]}, median {lens[len(lens) // 2]}, max {lens[-1]}; most common {mode} nt ({count / len(lens):.0%})"


def peek(folder: Path, sample: str, scan: int = 3000, n_examples: int = 5, k: int = 10, width: int = 60) -> list:
    from amplitrim.subsample import find_pairs
    pair = next(((r1, r2) for s, r1, r2 in find_pairs(folder) if s == sample), None)
    if pair is None:
        return [f"== {sample} in {folder}: not found"]
    r1, r2 = read_head(pair[0], scan), read_head(pair[1], scan)
    lines = [f"== {sample} in {folder}  (first {len(r1):,} pairs)"]
    for label, reads in (("R1", r1), ("R2", r2)):
        seqs = [s for s, _q in reads]
        lines.append(f"  {label} starts (first {k} bases): {start_table(seqs, k)}")
        lines.append(f"  {label} lengths: {length_summary(seqs)}")
    lines.append(f"  example pairs (first {width} bases):")
    for (s1, _a), (s2, _b) in zip(r1[:n_examples], r2[:n_examples]):
        lines.append(f"    R1 {s1[:width]:<{width}s} ({len(s1)} nt)")
        lines.append(f"    R2 {s2[:width]:<{width}s} ({len(s2)} nt)")
    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folders", nargs="+", type=Path, help="one or more folders of paired FASTQs")
    ap.add_argument("--sample", nargs="+", required=True, help="sample name(s), as in the file names before _R1/_R2")
    ap.add_argument("--n", type=int, default=5, help="example pairs to print (default 5)")
    ap.add_argument("--scan", type=int, default=3000, help="pairs used for the start and length summaries (default 3000)")
    ap.add_argument("--k", type=int, default=10, help="bases at the start counted together (default 10)")
    args = ap.parse_args()
    for sample in args.sample:
        for folder in args.folders:
            print("\n".join(peek(folder, sample, args.scan, args.n, args.k)) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
