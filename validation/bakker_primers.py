#!/usr/bin/env python
"""Which primers start the reads? A census for the Bakker (2018) runs, before any processing (the paper does not print its ITS2 primers).

  python validation\\bakker_primers.py --base C:\\data            (both Bakker folders; the first 2,000 reads of every file)

For every read pair it looks for a known fungal ITS primer at the start of each read (after a spacer of up to --max-offset bases, with up to --max-mismatch
mismatches; IUPAC codes in the primers are honoured), and counts the primer of R1, the primer of R2 and the pair of the two. Unmatched read starts are counted too
(their most common 20 bases), so that a primer that is not in the list shows itself. Writes PRIMER_CENSUS.md to --out (default --base) and prints it.
"""
import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path

PRIMERS = {  # name: 5'->3' as the primer reads (forward and reverse primers alike); the common fungal ITS primers
    "ITS1F": "CTTGGTCATTTAGAGGAAGTAA", "ITS1": "TCCGTAGGTGAACCTGCGG", "ITS1-F_KYO2": "TAGAGGAAGTAAAAGTCGTAA", "ITS5": "GGAAGTAAAAGTCGTAACAAGG",
    "ITS2": "GCTGCGTTCTTCATCGATGC", "ITS2_KYO2": "TTYRCTRCGTTCTTCATC", "ITS3": "GCATCGATGAAGAACGCAGC", "ITS3_KYO2": "GATGAAGAACGYAGYRAA",
    "fITS7": "GTGARTCATCGAATCTTTG", "gITS7": "GTGARTCATCGARTCTTTG", "ITS86F": "GTGAATCATCGAATCTTTGAA", "5.8S-Fun": "AACTTTYRRCAAYGGATCWCT",
    "BITS": "ACCTGCGGARGGATCA", "B58S3": "GAGATCCRTTGYTRAAAGTT",                                    # Bokulich & Mills 2013 (Bakker's first run)
    "ITS4": "TCCTCCGCTTATTGATATGC", "ITS4-Fun": "AGCCTCCGCTTATTGATATGCTTAART", "ITS4_KYO1": "TCCTCCGCTTWTTGWTWTGC", "ITS4_KYO3": "CTBTTVCCKCTTCACTCG",
}
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
FOLDERS = ["22_mock_ITS1_bakker_pcr_PRJNA377530", "23_mock_ITSmix_bakker_polymerase_PRJNA377530"]


def mismatches(read_part: str, primer: str) -> int:
    return sum(1 for b, p in zip(read_part, primer) if b not in IUPAC.get(p, p))


def best_primer(seq: str, max_mm: int = 2, max_offset: int = 10):
    """(name, offset, mismatches) of the best primer at the start of seq (fewest mismatches, then the smallest offset), or None."""
    best = None
    for name, pr in PRIMERS.items():
        for off in range(0, max_offset + 1):
            part = seq[off:off + len(pr)]
            if len(part) < len(pr):
                break
            mm = mismatches(part, pr)
            if mm <= max_mm and (best is None or (mm, off, -len(pr)) < (best[2], best[1], -len(PRIMERS[best[0]]))):
                best = (name, off, mm)
    return best


def reads(path: Path, n: int):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for i in range(n):
            h = fh.readline()
            if not h:
                return
            s = fh.readline().strip()
            fh.readline()
            fh.readline()
            yield s.upper()


def census(folder: Path, n: int, max_mm: int, max_offset: int) -> dict:
    r1c, r2c, pairs, offsets, unmatched, total = Counter(), Counter(), Counter(), Counter(), Counter(), 0
    for r1 in sorted(Path(folder).glob("*_R1*.f*q*")):
        r2 = r1.with_name(r1.name.replace("_R1", "_R2", 1))
        if not r2.is_file():
            continue
        for a, b in zip(reads(r1, n), reads(r2, n)):
            total += 1
            pa, pb = best_primer(a, max_mm, max_offset), best_primer(b, max_mm, max_offset)
            r1c[pa[0] if pa else "-"] += 1
            r2c[pb[0] if pb else "-"] += 1
            pairs[(pa[0] if pa else "-", pb[0] if pb else "-")] += 1
            for p in (pa, pb):
                if p:
                    offsets[(p[0], p[1])] += 1
            if not pa:
                unmatched[("R1", a[:20])] += 1
            if not pb:
                unmatched[("R2", b[:20])] += 1
    return {"total": total, "r1": r1c, "r2": r2c, "pairs": pairs, "offsets": offsets, "unmatched": unmatched}


def report(name: str, c: dict) -> list:
    t = c["total"] or 1
    L = [f"## {name}", "", f"{c['total']:,} read pairs looked at.", "", "| Primer on R1 | Primer on R2 | Pairs | Share |", "|---|---|---|---|"]
    L += [f"| {a} | {b} | {k:,} | {k / t:.1%} |" for (a, b), k in c["pairs"].most_common() if k / t >= 0.005]
    offs = {}
    for (p, off), k in c["offsets"].items():
        offs.setdefault(p, Counter())[off] += k
    L += ["", "Spacer before the primer (bases: reads): " + "; ".join(f"{p} " + ", ".join(f"{o}: {k:,}" for o, k in sorted(v.items())) for p, v in sorted(offs.items()))]
    if c["unmatched"]:
        L += ["", "Most common starts of unmatched reads: " + "; ".join(f"{m} {s} ({k:,})" for (m, s), k in c["unmatched"].most_common(6))]
    return L + [""]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--folders", nargs="*", default=FOLDERS, help="folders under --base, each with a raw\\ folder (default: the two Bakker folders)")
    ap.add_argument("--reads", type=int, default=2000, help="read pairs looked at per sample (default 2,000)")
    ap.add_argument("--max-mismatch", type=int, default=2)
    ap.add_argument("--max-offset", type=int, default=10)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    L = ["# Primer census of the Bakker (2018) runs", "", f"The first {a.reads:,} read pairs of every sample; a primer is found after a spacer of up to {a.max_offset} bases "
         f"with up to {a.max_mismatch} mismatches. '-' means no known primer.", ""]
    found = False
    for f in a.folders:
        raw = a.base / f / "raw"
        if not raw.is_dir():
            L += [f"## {f}", "", "not found", ""]
            continue
        found = True
        L += report(f, census(raw, a.reads, a.max_mismatch, a.max_offset))
    out = a.out or a.base
    out.mkdir(parents=True, exist_ok=True)
    (out / "PRIMER_CENSUS.md").write_text("\n".join(L) + "\nNext: send this file.\n", encoding="utf-8")
    print("\n".join(L))
    print(f"written: {out / 'PRIMER_CENSUS.md'}")
    return 0 if found else 2


if __name__ == "__main__":
    sys.exit(main())
