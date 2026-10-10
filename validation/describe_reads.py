#!/usr/bin/env python3
"""
Describe the paired FASTQ files in a folder before analysing them: reads, read lengths, quality scheme, and which
common 16S/ITS primers the reads start with (and after how many bases, which shows barcodes or spacers), plus the
share of pairs sequenced in the reverse orientation. One line per sample, then a list of anything unusual.

    python validation/describe_reads.py C:\\data\\fastq [--reads 2000] [--csv describe.csv]
"""
import argparse
import csv
import gzip
import re
import statistics
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from amplitrim.subsample import find_pairs  # noqa: E402

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "M": "[AC]", "K": "[GT]", "S": "[CG]",
         "W": "[AT]", "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
# name: (sequence, F or R). Only the first 14 bases are searched for, so degenerate and truncated versions match.
PRIMERS = {
    "27F": ("AGAGTTTGATYMTGGCTCAG", "F"), "338F": ("ACTCCTACGGGAGGCAGCAG", "F"), "341F": ("CCTACGGGNGGCWGCAG", "F"),
    "357F": ("CCTACGGGAGGCAGCAG", "F"), "515F": ("GTGYCAGCMGCCGCGGTAA", "F"), "799F": ("AACMGGATTAGATACCCKG", "F"),
    "926F": ("AAACTYAAAKGAATTGRCGG", "F"), "967F": ("CAACGCGAAGAACCTTACC", "F"),
    "ITS1-Fl2": ("GAACCWGCGGARGGATCA", "F"), "ITS1F-KYO1": ("CTHGGTCATTTAGAGGAASTAA", "F"), "ITS2-KYO2": ("TTYRCTRCGTTCTTCATC", "R"),
    "563F": ("AYTGGGYDTAAAGNG", "F"), "338R": ("GCTGCCTCCCGTAGGAGT", "R"), "802R": ("TACNVGGGTATCTAATCC", "R"), "519R": ("GWATTACCGCGGCKGCTG", "R"), "534R": ("ATTACCGCGGCTGCTGG", "R"),
    "785R/805R": ("GACTACHVGGGTATCTAATCC", "R"), "806R": ("GGACTACNVGGGTWTCTAAT", "R"), "926R": ("CCGYCAATTYMTTTRAGTTT", "R"),
    "1193R": ("ACGTCATCCCCACCTTCC", "R"), "1391R": ("GACGGGCGGTGWGTRCA", "R"), "1492R": ("GGTTACCTTGTTACGACTT", "R"),
    "ITS1F": ("CTTGGTCATTTAGAGGAAGTAA", "F"), "ITS1": ("TCCGTAGGTGAACCTGCGG", "F"), "ITS3": ("GCATCGATGAAGAACGCAGC", "F"),
    "ITS5": ("GGAAGTAAAAGTCGTAACAAGG", "F"), "gITS7": ("GTGARTCATCGARTCTTTG", "F"), "ITS86F": ("GTGAATCATCGAATCTTTGAA", "F"),
    "ITS2": ("GCTGCGTTCTTCATCGATGC", "R"), "ITS4": ("TCCTCCGCTTATTGATATGC", "R"),
}
SEEDS = {name: re.compile("".join(IUPAC[c] for c in seq[:14])) for name, (seq, _d) in PRIMERS.items()}
SEEDS_SHIFTED = {name: [(k, re.compile("".join(IUPAC[c] for c in seq[k:k + 14]))) for k in (1, 2)] for name, (seq, _d) in PRIMERS.items()}
LENGTH_VARIES = 0.05       # a sample's read lengths count as varying when more than this share of its reads is off the most common length
HEAD = 20               # the latest base at which a built-in primer may start (barcodes and spacers); a match further in is the gene itself,
                        # e.g. 563F inside any amplicon that begins at 515F


def read_head(path: Path, n: int):
    opener = gzip.open if str(path).endswith(".gz") else open
    out = []
    with opener(path, "rt") as fh:
        while len(out) < n:
            h = fh.readline()
            if not h:
                break
            seq, _, qual = fh.readline().strip(), fh.readline(), fh.readline().strip()
            out.append((seq, qual))
    return out


def count_reads(path: Path) -> int:
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rb") as fh:
        while chunk := fh.read(1 << 22):
            n += chunk.count(b"\n")
    return n // 4


def _full_mismatches(head: str, primer: str, start: int) -> int:
    """Bases of the whole primer (IUPAC codes allowed) that differ from the read, with the primer starting at `start` (negative: its first bases are missing)."""
    skip = max(-start, 0)
    return sum(1 for i, c in enumerate(primer[skip:], start=max(start, 0)) if i < len(head) and not re.match(IUPAC[c], head[i]))


def find_primer(seq: str):
    """(name, start offset) of the primer whose seed occurs earliest in the first HEAD bases, or None. A read that begins inside a primer
    (its first one or two bases are missing, which is common) gives offset -1 or -2."""
    head = seq[:HEAD + 30]
    found = []
    for name, rx in SEEDS.items():
        m = rx.search(head)                                 # HEAD limits where a match may START, not how much is searched
        if m and m.start() <= HEAD:
            found.append((m.start(), name))
    for name, shifted in SEEDS_SHIFTED.items():
        for k, rx in shifted:
            if rx.match(head):                              # anchored: the read must begin with the rest of the primer
                found.append((-k, name))
    if not found:
        return None
    # the earliest start wins; primers that start together (a degenerate primer contains a plain one, as ITS1F-KYO1 contains ITS1F) are
    # told apart by how well the WHOLE primer fits the read, then by which is more specific
    start, _m, _d, name = min((s, _full_mismatches(head, PRIMERS[n][0], s), sum(c not in "ACGT" for c in PRIMERS[n][0]), n) for s, n in found)
    return name, start


MIN_SHARE = 0.05         # a primer in fewer reads than this is chance (for example a conserved site), not a primer


def top_primer(reads):
    found = [find_primer(s) for s, _q in reads]
    counts = Counter(f[0] for f in found if f)
    if not counts:
        return None, 0.0, None
    name, n = counts.most_common(1)[0]
    if n / len(reads) < MIN_SHARE:
        return None, 0.0, None
    offsets = Counter(f[1] for f in found if f and f[0] == name)
    return name, n / len(reads), offsets.most_common(1)[0][0]


def start_profile(reads, motif_length: int = 20, positions: int = 12):
    """(commonest first `motif_length` bases, its share of reads, mean share of the commonest base over the first
    `positions` positions). Reads that begin at the same place, as after primer removal, are strongly conserved."""
    seqs = [s for s, _q in reads if len(s) >= motif_length]
    if not seqs:
        return "", 0.0, 0.0
    motif, n = Counter(s[:motif_length] for s in seqs).most_common(1)[0]
    conservation = sum(Counter(s[i] for s in seqs).most_common(1)[0][1] / len(seqs) for i in range(positions)) / positions
    return motif, n / len(seqs), conservation


ALLOWED = {c: set(v.strip("[]")) for c, v in IUPAC.items()}
COMPLEMENT = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def mismatches(window: str, primer: str) -> int:
    return sum(base not in ALLOWED.get(p, {p}) for base, p in zip(window, primer))


def best_match(seq: str, primer: str, max_mismatches: int, max_offset: int):
    """(mismatches, offset) of the best occurrence of `primer` starting within max_offset, or None."""
    best = None
    for offset in range(0, min(max_offset, len(seq) - len(primer)) + 1):
        mm = mismatches(seq[offset:offset + len(primer)], primer)
        if mm <= max_mismatches and (best is None or (mm, offset) < best):
            best = (mm, offset)
    return best


def search_primer(reads, primer: str, max_mismatches: int = 2, max_offset: int = 100):
    """Share of reads with the primer starting within the first `max_offset` bases (and its commonest offset), and
    the share containing its reverse complement anywhere (read-through into the other end of the amplicon)."""
    forward, offsets, rc, rc_offsets = 0, Counter(), 0, []
    rc_primer = reverse_complement(primer)
    for seq, _q in reads:
        hit = best_match(seq, primer, max_mismatches, max_offset)
        if hit:
            forward += 1
            offsets[hit[1]] += 1
        hit_rc = best_match(seq, rc_primer, max_mismatches, len(seq))
        if hit_rc:
            rc += 1
            rc_offsets.append(hit_rc[1])
    n = max(len(reads), 1)
    median = sorted(rc_offsets)[len(rc_offsets) // 2] if rc_offsets else None
    # Read-through starts where the amplicon ends, so its matches cluster; chance matches are spread along the read.
    clustered = sum(abs(o - median) <= 10 for o in rc_offsets) / len(rc_offsets) if rc_offsets else 0.0
    return {"share": forward / n, "offset": offsets.most_common(1)[0][0] if offsets else None, "rc_share": rc / n,
            "rc_median_offset": median, "rc_clustered": clustered}


CODES = {"A": "A", "C": "C", "G": "G", "T": "T", "AG": "R", "CT": "Y", "AC": "M", "GT": "K", "CG": "S", "AT": "W",
         "CGT": "B", "AGT": "D", "ACT": "H", "ACG": "V", "ACGT": "N"}


def start_consensus(reads, positions: int = 20, min_share: float = 0.12) -> str:
    """The first `positions` bases as an IUPAC string: at each position every base seen in at least `min_share` of the
    reads. A primer with degenerate positions shows up as a short string with a few ambiguity codes, for example
    AYTGGGYDTAAAGNG; unrelated sequence quickly turns into N."""
    out = []
    for i in range(positions):
        counts = Counter(s[i] for s, _q in reads if len(s) > i and s[i] in "ACGT")
        total = sum(counts.values())
        keep = "".join(sorted(b for b, n in counts.items() if total and n / total >= min_share))
        out.append(CODES.get(keep, "N"))
    return "".join(out)


def position_table(reads, positions: int = 15):
    """Per position, the share of each base among the reads, and the commonest base."""
    rows = []
    for i in range(positions):
        counts = Counter(s[i] for s, _q in reads if len(s) > i)
        total = max(sum(counts.values()), 1)
        rows.append({"position": i + 1, **{b: counts.get(b, 0) / total for b in "ACGTN"}, "top": counts.most_common(1)[0][0] if counts else ""})
    return rows


def pooled_heads(folder: Path, per_sample: int, limit: int = 3000):
    """The first reads of R1 and of R2 pooled over all samples (at most `limit` of each)."""
    r1, r2 = [], []
    for _sample, p1, p2 in find_pairs(folder):
        r1 += read_head(p1, per_sample)
        r2 += read_head(p2, per_sample)
    step = max(1, len(r1) // limit)
    return r1[::step][:limit], r2[::step][:limit]


LENGTH_BINS = [(0, 99, "<100"), (100, 149, "100-149"), (150, 189, "150-189"), (190, 199, "190-199"), (200, 205, "200-205"),
               (206, 209, "206-209"), (210, 219, "210-219"), (220, 229, "220-229"), (230, 10 ** 6, ">=230")]


def pair_length_summary(len1, len2) -> str:
    """How the two reads of a pair compare in length. A constant difference (R2 is 3 nt longer because its primer is 3 nt shorter) is
    normal and is said as such; only a mixture of differences is a sign of trouble."""
    from collections import Counter
    diffs = Counter(b - a for a, b in zip(len1, len2))
    n = max(sum(diffs.values()), 1)
    d, k = diffs.most_common(1)[0] if diffs else (0, 0)
    if d == 0:
        return f"R1 and R2 have the same length in {k / n:.0%} of pairs."
    side = "longer" if d > 0 else "shorter"
    tail = "this is normal when the two primers differ in length" if k / n >= 0.8 else "the rest differ by other amounts"
    return f"R2 is {abs(d)} nt {side} than R1 in {k / n:.0%} of pairs; {tail}."


def length_bins(lengths):
    """[(label, share)] of the lengths in each bin."""
    n = max(len(lengths), 1)
    return [(label, sum(lo <= x <= hi for x in lengths) / n) for lo, hi, label in LENGTH_BINS]


def describe(folder: Path, n_reads: int = 2000):
    rows = []
    for sample, p1, p2 in find_pairs(folder):
        a, b = read_head(p1, n_reads), read_head(p2, n_reads)
        qa, qb = Counter(), Counter()
        for _s, q in a:
            qa.update(q)
        for _s, q in b:
            qb.update(q)
        lens = Counter(len(s) for s, _q in a)
        n1, share1, off1 = top_primer(a)
        n2, share2, off2 = top_primer(b)
        motif, motif_share, conservation = start_profile(a)
        pf = [find_primer(sa) for (sa, _q) in a]
        pr = [find_primer(sb) for (sb, _q) in b]
        reversed_pairs = forward_pairs = 0
        for x, y in zip(pf, pr):
            if x and y:
                dx, dy = PRIMERS[x[0]][1], PRIMERS[y[0]][1]
                forward_pairs += (dx, dy) == ("F", "R")
                reversed_pairs += (dx, dy) == ("R", "F")
        oriented = forward_pairs + reversed_pairs
        rows.append({"sample": sample, "reads": count_reads(p1), "read_length": lens.most_common(1)[0][0],
                     "off_modal_share": round(1 - lens.most_common(1)[0][1] / max(sum(lens.values()), 1), 4),
                     "lengths_differ": 1 - lens.most_common(1)[0][1] / max(sum(lens.values()), 1) > LENGTH_VARIES, "quality_values": len(set(qa) | set(qb)),
                     "r1_primer": n1 or "", "r1_share": round(share1, 3), "r1_offset": off1 if off1 is not None else "",
                     "r2_primer": n2 or "", "r2_share": round(share2, 3), "r2_offset": off2 if off2 is not None else "",
                     "start_motif": motif, "start_motif_share": round(motif_share, 3), "start_conservation": round(conservation, 3),
                     "reversed_share": round(reversed_pairs / oriented, 3) if oriented >= MIN_SHARE * len(a) else ""})
    return rows


def summarise(rows):
    """Plain-language findings, one per line."""
    if not rows:
        return ["No paired FASTQ files were found."]
    notes = []
    flat = [r["sample"] for r in rows if r["quality_values"] <= 2]
    if flat:
        notes.append(f"{len(flat)} of {len(rows)} samples have 2 or fewer distinct quality values ({', '.join(flat[:4])}"
                     f"{', ...' if len(flat) > 4 else ''}): the per-base quality is missing and DADA2 cannot learn an error "
                     f"model from them. Get the original-quality files.")
    binned = [r for r in rows if 3 <= r["quality_values"] <= 5]
    if binned:
        notes.append(f"{len(binned)} sample{'s' if len(binned) != 1 else ''} {'have' if len(binned) != 1 else 'has'} only 3-5 distinct quality values (binned scores).")
    missing = [r["sample"] for r in rows if not r["r1_primer"] or r["r1_share"] < 0.5]
    if missing:
        notes.append(f"{len(missing)} sample{'s' if len(missing) != 1 else ''} {'show' if len(missing) != 1 else 'shows'} no known primer at the start of R1 ({', '.join(missing[:4])}"
                     f"{', ...' if len(missing) > 4 else ''}): the primers may already be removed, or are not in the "
                     f"built-in list.")
    if missing:
        sample = max((r for r in rows if r["sample"] in missing), key=lambda r: r["start_conservation"])
        notes.append(f"Without a primer, R1 of {sample['sample']} begins with {sample['start_motif']} in "
                     f"{sample['start_motif_share']:.0%} of reads, and its first 12 positions are {sample['start_conservation']:.0%} "
                     f"conserved on average. Reads that all begin at the same place, as after primer removal, are strongly "
                     f"conserved; if that is what you see, run AmpliTrim directly on these files. To confirm, run again with "
                     f"--positions 15 --primer YOUR_PRIMER_SEQUENCE.")
    offsets = Counter((r["r1_offset"], r["r2_offset"]) for r in rows if r["r1_primer"])
    if len(offsets) > 1:
        notes.append("The primer starts at different offsets in different samples (R1, R2): "
                     + ", ".join(f"{k} in {v}" for k, v in offsets.most_common(4)) + ".")
    elif offsets:
        (o1, o2), _n = offsets.most_common(1)[0]
        if (o1 or 0) < 0 or (o2 or 0) < 0:
            notes.append(f"The reads start inside the primer: {abs(min(o1 or 0, 0))} base(s) of it are missing on R1 and {abs(min(o2 or 0, 0))} on R2 "
                         f"(common); cutadapt removes what remains.")
        if (o1 or 0) > 0 or (o2 or 0) > 0:
            notes.append(f"Every sample's primer starts after {o1} bases on R1 and {o2} on R2 (barcodes or spacers); "
                         f"cutadapt removes them with the primer.")
    mixed = [r["sample"] for r in rows if r["reversed_share"] != "" and 0.02 < r["reversed_share"]]
    if mixed:
        notes.append(f"{len(mixed)} sample{'s' if len(mixed) != 1 else ''} {'contain' if len(mixed) != 1 else 'contains'} pairs sequenced in the reverse orientation; use --mixed-orientation "
                     f"when trimming.")
    varying = [r for r in rows if r["lengths_differ"]]
    if varying:
        worst = max(varying, key=lambda r: r["off_modal_share"])
        notes.append(f"Read lengths vary: in {len(varying)} sample{'s' if len(varying) != 1 else ''} more than {LENGTH_VARIES:.0%} of the reads are not the most "
                     f"common length (most in {worst['sample']}: {worst['off_modal_share']:.0%}).")
    counts = [r["reads"] for r in rows]
    notes.append(f"Reads per sample: min {min(counts):,}, median {int(statistics.median(counts)):,}, max {max(counts):,}.")
    return notes


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", type=Path)
    ap.add_argument("--reads", type=int, default=2000, help="reads per file to inspect (default 2000)")
    ap.add_argument("--csv", type=Path, help="also write the table to this CSV file")
    ap.add_argument("--primer", action="append", default=[], metavar="SEQ",
                    help="also search for this primer (IUPAC letters allowed), with mismatches and at any start offset; repeat for more")
    ap.add_argument("--mismatches", type=int, default=2, help="mismatches allowed in --primer searches (default 2)")
    ap.add_argument("--lengths", action="store_true", help="also print the distribution of read lengths (useful after trimming)")
    ap.add_argument("--positions", type=int, default=0, metavar="N",
                    help="also print the base composition of the first N positions of R1 and R2 (try 15)")
    args = ap.parse_args()
    if not args.folder.is_dir():
        print(f"Folder not found: {args.folder}")
        return 1
    rows = describe(args.folder, args.reads)
    print(f"{'sample':24s} {'reads':>8s} {'len':>4s} {'Qvals':>5s}  {'R1 primer (share, offset)':28s} {'R2 primer (share, offset)':28s} {'reversed':>8s}  R1 starts with (share)")
    for r in rows:
        p1 = f"{r['r1_primer']} ({r['r1_share']:.0%}, {r['r1_offset']})" if r["r1_primer"] else "none"
        p2 = f"{r['r2_primer']} ({r['r2_share']:.0%}, {r['r2_offset']})" if r["r2_primer"] else "none"
        rev = "" if r["reversed_share"] == "" else f"{r['reversed_share']:.0%}"
        print(f"{r['sample'][:24]:24s} {r['reads']:>8,d} {r['read_length']:>4d} {r['quality_values']:>5d}  {p1:28s} {p2:28s} {rev:>8s}  {r['start_motif']} ({r['start_motif_share']:.0%})")
    print()
    for note in summarise(rows):
        print("-", note)
    r1, r2 = pooled_heads(args.folder, args.reads)
    print(f"\nWhat the reads start with (every base seen in at least 12% of reads at each of the first 20 positions; "
          f"{len(r1):,} reads pooled):\n  R1: {start_consensus(r1)}\n  R2: {start_consensus(r2)}\n  A primer shows up here as a "
          f"conserved string with a few ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N). Compare it with your protocol, or "
          f"search it with --primer.")
    if args.lengths:
        r1_len, r2_len = pooled_heads(args.folder, args.reads, limit=20000)
        print(f"\nRead lengths ({len(r1_len):,} read pairs pooled over all samples):")
        for label, reads in (("R1", r1_len), ("R2", r2_len)):
            print(f"  {label}: " + ", ".join(f"{name} {share:.1%}" for name, share in length_bins([len(s) for s, _q in reads])))
        print("  " + pair_length_summary([len(x[0]) for x in r1_len], [len(y[0]) for y in r2_len]))
    if args.primer or args.positions:
        if args.primer:
            print(f"\nPrimer search ({args.mismatches} mismatches allowed, {len(r1):,} reads pooled over all samples):")
            for primer in args.primer:
                f1, f2 = search_primer(r1, primer, args.mismatches), search_primer(r2, primer, args.mismatches)
                where = lambda r: f"{r['share']:.1%}" + (f" (starts at base {r['offset']})" if r["offset"] is not None else "")
                def rc_text(r):
                    if r["rc_share"] < 0.005:
                        return f"{r['rc_share']:.1%}"
                    return (f"{r['rc_share']:.1%} (starting near base {r['rc_median_offset']}; {r['rc_clustered']:.0%} of the "
                            f"matches within 10 bases of that)")
                print(f"  {primer}: R1 {where(f1)}, R2 {where(f2)}\n      reverse complement anywhere: R1 {rc_text(f1)}; "
                      f"R2 {rc_text(f2)}")
        if args.positions:
            for label, reads in (("R1", r1), ("R2", r2)):
                print(f"\n{label}: base composition of the first {args.positions} positions")
                print("  pos     A     C     G     T     N  top")
                for t in position_table(reads, args.positions):
                    print(f"  {t['position']:>3d} " + " ".join(f"{t[b]:>5.0%}" for b in "ACGTN") + f"  {t['top']}")
    if args.csv and rows:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
