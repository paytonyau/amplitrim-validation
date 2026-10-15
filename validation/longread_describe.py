#!/usr/bin/env python
"""Describe long-read amplicon reads (PacBio HiFi, Nanopore, ...) before any analysis. One line per FASTQ file, then a list of anything unusual.

  python validation/longread_describe.py reads/ --out longread --primer-f AGRGTTYGATYMTGGCTCAG --primer-r RGYTACCTTGTTACGACTT
  python validation/longread_describe.py reads/ --out longread                       (no primers given: every read is used as it is)

For every file (single-end, .fastq or .fastq.gz):
  primers      the share of reads with both primers (in the forward or the reverse orientation), with one only, with none, and how often each primer is found. A primer must match
               within --window bases of the end of the read with at most --error-rate of its bases wrong, where a substitution, an insertion and a deletion each count as one (an
               edit distance: HiFi and Nanopore reads often have an indel inside a primer, which a mismatches-only search cannot see).
  lengths      of the reads between the primers: percentiles, the commonest 10 nt bin, the share within 10% of the median, the shares much
               longer (concatemer-like) or much shorter than it
  errors       the expected errors of those reads (the sum of the error probabilities of the quality scores), per read and per 1,000 nt, and
               the share of reads that stays at each maxEE of --grid
  quality      the number of distinct quality values and the highest one (a copy of the reads that kept only a few values cannot give expected
               errors), and the mean quality of the first 50, the middle half and the last 50 bases of the reads between the primers

Needs only numpy. It makes no decision: it measures what a long-read mode would need to decide, and says where the reads are unusual.
Writes <out>/longread_describe.tsv and <out>/LONGREAD_DESCRIBE.md.
"""
import argparse
import csv
import gzip
import math
import sys
from pathlib import Path

import numpy as np

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
LUT = np.full(256, 4, dtype=np.uint8)                                  # A, C, G, T -> 0..3; anything else (N) -> 4, which never matches
for _i, _c in enumerate("ACGT"):
    LUT[ord(_c)] = LUT[ord(_c.lower())] = _i
ERROR_PROBABILITY = 10.0 ** (-np.arange(94) / 10.0)                     # Phred 0..93 (PacBio HiFi goes up to 93)
DEFAULT_GRID = (1, 2, 3, 5, 8, 10, 15, 20, 30)
SEGMENT = 50                                                            # bases at each end that are compared with the middle
CHUNK = 2000                                                            # reads whose primers are searched for at once
MIN_BETWEEN_PRIMERS = 50


def revcomp(seq: str) -> str:
    return seq.upper().translate(COMPLEMENT)[::-1]


class Primer:
    """A primer (IUPAC codes allowed) as a table: which read bases are accepted at each of its positions."""

    def __init__(self, seq: str):
        self.seq = seq.upper()
        if not self.seq or any(c not in IUPAC for c in self.seq):
            raise ValueError(f"not a primer sequence: {seq!r}")
        self.allowed = np.zeros((len(self.seq), 5), dtype=bool)
        for i, ch in enumerate(self.seq):
            for base in IUPAC[ch]:
                self.allowed[i, "ACGT".index(base)] = True

    def __len__(self):
        return len(self.seq)


def edit_matches(regions: np.ndarray, primer: Primer):
    """For each row of `regions` (read bases as 0..4; 4 never matches) the fewest edits (substitutions, insertions and deletions, each one) that turn the primer into some stretch of the
    row, and the position just after that stretch. The IUPAC codes of the primer accept their bases. Where several stretches are equally good the one that ends last is taken, so that a
    substitution at the end of the primer is not read as a missing base. All rows are done at once."""
    rows, width = regions.shape
    m = len(primer)
    column = np.arange(width + 1, dtype=np.int16)
    previous = np.zeros((rows, width + 1), dtype=np.int16)               # no primer base used yet: the stretch may start anywhere
    for i in range(1, m + 1):
        mismatch = ~primer.allowed[i - 1][regions]                       # the read base does not fit primer base i
        t = np.empty_like(previous)
        t[:, 0] = i
        t[:, 1:] = np.minimum(previous[:, :-1] + mismatch, previous[:, 1:] + 1)    # use a read base for primer base i, fitting or not; or skip primer base i (a deletion)
        previous = np.minimum.accumulate(t - column, axis=1) + column    # or let read bases that belong to no primer base go by (insertions)
    fewest = previous.min(axis=1)
    end = width - previous[:, ::-1].argmin(axis=1)
    return fewest, end


class Locator:
    """Finds the primers at the two ends of reads and says in which orientation each read was sequenced."""

    def __init__(self, primer_f: str, primer_r: str, window: int = 60, error_rate: float = 0.15):
        self.f, self.r = Primer(primer_f), Primer(primer_r)
        self.f_rc, self.r_rc = Primer(revcomp(primer_f)), Primer(revcomp(primer_r))
        self.window, self.error_rate = window, error_rate
        self._from_the_end = {id(p): Primer(p.seq[::-1]) for p in (self.f_rc, self.r_rc)}      # a primer at the 3' end is matched from the end of the read inwards

    def _limit(self, primer: Primer) -> int:
        return int(self.error_rate * len(primer))

    def _five(self, batch: list, primer: Primer):
        """Per read: where the amplicon starts after a primer at the 5' end (-1 when there is none), and the edits."""
        width = len(primer) + self.window
        regions = np.full((len(batch), width), 4, dtype=np.uint8)
        for k, codes in enumerate(batch):
            part = codes[:width]
            regions[k, :len(part)] = part
        edits, end = edit_matches(regions, primer)
        return np.where(edits <= self._limit(primer), end, -1), edits

    def _three(self, batch: list, primer: Primer):
        """Per read: where the amplicon ends before a primer at the 3' end (-1 when there is none), and the edits."""
        width = len(primer) + self.window
        regions = np.full((len(batch), width), 4, dtype=np.uint8)
        lengths = np.empty(len(batch), dtype=np.int64)
        for k, codes in enumerate(batch):
            tail = codes[-width:][::-1]
            regions[k, :len(tail)] = tail
            lengths[k] = len(codes)
        edits, consumed = edit_matches(regions, self._from_the_end[id(primer)])
        return np.where(edits <= self._limit(primer), lengths - consumed, -1), edits

    def locate_many(self, seqs: list) -> list:
        """For every read (orientation '+', '-' or None, category, start, end): the category is 'both', 'five_only', 'three_only' or 'none'; start and end are the positions of the
        amplicon between the primers (a missing primer leaves that end of the read as it is)."""
        batch = [LUT[np.frombuffer(s.encode("ascii"), dtype=np.uint8)] for s in seqs]
        f5, f5e = self._five(batch, self.f)
        r3, r3e = self._three(batch, self.r_rc)
        r5, r5e = self._five(batch, self.r)
        f3, f3e = self._three(batch, self.f_rc)
        out = []
        for k, codes in enumerate(batch):
            n = len(codes)
            plus = ((int(f5[k]), int(f5e[k])) if f5[k] >= 0 else None, (int(r3[k]), int(r3e[k])) if r3[k] >= 0 else None)
            minus = ((int(r5[k]), int(r5e[k])) if r5[k] >= 0 else None, (int(f3[k]), int(f3e[k])) if f3[k] >= 0 else None)
            found = lambda pair: sum(x is not None for x in pair)                                # noqa: E731
            cost = lambda pair: sum(x[1] for x in pair if x is not None)                         # noqa: E731
            if found(plus) == 0 and found(minus) == 0:
                out.append((None, "none", 0, n))
                continue
            orientation, pair = ("+", plus) if (found(plus), -cost(plus)) >= (found(minus), -cost(minus)) else ("-", minus)
            start = pair[0][0] if pair[0] else 0
            end = pair[1][0] if pair[1] else n
            if end - start < MIN_BETWEEN_PRIMERS:
                out.append((None, "none", 0, n))
                continue
            out.append((orientation, "both" if all(pair) else ("five_only" if pair[0] else "three_only"), start, end))
        return out

    def locate(self, seq: str):
        """The same for one read."""
        return self.locate_many([seq])[0]


def open_text(path: Path):
    return gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else open(path, "rt", newline="")


def read_fastq(path: Path):
    """(name, sequence, quality) of every record of a single-end FASTQ file."""
    with open_text(path) as fh:
        while True:
            head = fh.readline()
            if not head:
                return
            seq, _, qual = fh.readline().strip(), fh.readline(), fh.readline().strip()
            yield head.rstrip("\r\n")[1:], seq, qual


def find_files(folder: Path) -> list:
    files = [p for p in sorted(Path(folder).iterdir()) if p.name.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))]
    return files


def sample_name(path: Path) -> str:
    name = path.name
    for suffix in (".gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


def mean_q(error_probability: float) -> float:
    return -10.0 * math.log10(error_probability) if error_probability > 0 else float("inf")


def describe_file(path: Path, locator=None, max_reads: int = 0, grid=DEFAULT_GRID) -> dict:
    """The numbers of one FASTQ file (see the module text)."""
    counts = {"both_forward": 0, "both_reverse": 0, "five_only": 0, "three_only": 0, "none": 0}
    raw_lengths, lengths, ees, per_1000 = [], [], [], []
    seen = np.zeros(94, dtype=bool)
    seg_start, seg_mid, seg_end = [], [], []
    bad = 0
    chunk = []

    def process(chunk: list) -> None:
        if locator is None:
            places = [("+", "both", 0, len(seq)) for seq, _ in chunk]
        else:
            places = locator.locate_many([seq for seq, _ in chunk])
        for (seq, qual), (orientation, category, start, end) in zip(chunk, places):
            q = np.clip(np.frombuffer(qual.encode("ascii"), dtype=np.uint8).astype(np.int16) - 33, 0, 93)
            seen[q] = True
            if category == "both":
                counts["both_forward" if orientation == "+" else "both_reverse"] += 1
            else:
                counts[category] += 1
                continue
            pe = ERROR_PROBABILITY[q[start:end]]
            if orientation == "-":
                pe = pe[::-1]                                          # the amplicon read from its own start
            n = len(pe)
            lengths.append(n)
            ees.append(float(pe.sum()))
            per_1000.append(float(pe.sum()) * 1000.0 / n)
            if n >= 3 * SEGMENT:
                seg_start.append(float(pe[:SEGMENT].mean()))
                seg_end.append(float(pe[-SEGMENT:].mean()))
                seg_mid.append(float(pe[n // 4: 3 * n // 4].mean()))

    for i, (_, seq, qual) in enumerate(read_fastq(path)):
        if max_reads and i >= max_reads:
            break
        if not seq or len(seq) != len(qual):
            bad += 1
            continue
        raw_lengths.append(len(seq))
        chunk.append((seq, qual))
        if len(chunk) >= CHUNK:
            process(chunk)
            chunk = []
    if chunk:
        process(chunk)
    total = len(raw_lengths)
    row = {"sample": sample_name(path), "reads": total, "unreadable": bad,
           "raw_len_median": float(np.median(raw_lengths)) if total else None,
           "raw_len_min": min(raw_lengths) if total else None, "raw_len_max": max(raw_lengths) if total else None}
    for k, v in counts.items():
        row[k] = v
    both = counts["both_forward"] + counts["both_reverse"]
    row["share_both"] = both / total if total else None
    row["share_both_reverse"] = counts["both_reverse"] / both if both else None
    row["share_five"] = (both + counts["five_only"]) / total if total else None              # the primer at the 5' end of the read was found
    row["share_three"] = (both + counts["three_only"]) / total if total else None            # the primer at the 3' end of the read was found
    row["between_primers"] = len(lengths)
    if lengths:
        L, E = np.array(lengths), np.array(ees)
        median = float(np.median(L))
        bins = np.bincount(L // 10)
        row.update({"len_min": int(L.min()), "len_p1": float(np.percentile(L, 1)), "len_p5": float(np.percentile(L, 5)), "len_median": median,
                    "len_p95": float(np.percentile(L, 95)), "len_p99": float(np.percentile(L, 99)), "len_max": int(L.max()),
                    "len_mode_bin": int(bins.argmax()) * 10,
                    "share_within_10pct": float(np.mean(np.abs(L - median) <= 0.1 * median)),
                    "share_over_1p5x": float(np.mean(L > 1.5 * median)), "share_under_half": float(np.mean(L < 0.5 * median)),
                    "ee_p25": float(np.percentile(E, 25)), "ee_median": float(np.median(E)), "ee_p75": float(np.percentile(E, 75)), "ee_p95": float(np.percentile(E, 95)),
                    "ee_per_1000_median": float(np.median(per_1000))})
        for g in grid:
            row[f"keep_maxEE_{g:g}"] = float(np.mean(E <= g))
    q_values = np.flatnonzero(seen)
    row["distinct_q"] = int(len(q_values))
    row["max_q"] = int(q_values.max()) if len(q_values) else None
    if seg_mid:
        for name, values in (("start", seg_start), ("mid", seg_mid), ("end", seg_end)):
            row[f"p_{name}"] = float(np.mean(values))
            row[f"q_{name}"] = mean_q(float(np.mean(values)))
    return row


def flags(row: dict) -> list:
    """What is unusual about one file."""
    out = []
    if row["reads"] == 0:
        return ["no reads"]
    if row.get("share_both") is not None and row["share_both"] < 0.9:
        out.append(f"only {100 * row['share_both']:.0f}% of the reads have both primers")
    if row.get("share_both_reverse") and row["share_both_reverse"] > 0.05:
        out.append(f"{100 * row['share_both_reverse']:.0f}% of the reads with both primers are in the reverse orientation")
    if row["distinct_q"] < 20:
        out.append(f"only {row['distinct_q']} distinct quality values (highest {row['max_q']}): the quality scores look simplified, so expected errors may not be reliable")
    if row.get("share_over_1p5x", 0) > 0.01:
        out.append(f"{100 * row['share_over_1p5x']:.1f}% of the reads are more than 1.5 times the median length (concatemers?)")
    if row.get("share_under_half", 0) > 0.01:
        out.append(f"{100 * row['share_under_half']:.1f}% of the reads are under half the median length")
    if "q_mid" in row and min(row["q_start"], row["q_end"]) < row["q_mid"] - 3:
        out.append(f"the quality at an end (start Q{row['q_start']:.0f}, end Q{row['q_end']:.0f}) is more than 3 below the middle (Q{row['q_mid']:.0f})")
    return out


def pooled(rows: list) -> dict:
    """Totals over files: reads, the share with both primers, the weighted quality of the segments."""
    reads = sum(r["reads"] for r in rows)
    both = sum(r["both_forward"] + r["both_reverse"] for r in rows)
    n = sum(r["between_primers"] for r in rows)
    out = {"files": len(rows), "reads": reads, "share_both": both / reads if reads else None, "between_primers": n,
           "share_five": sum(r["share_five"] * r["reads"] for r in rows if r.get("share_five") is not None) / reads if reads else None,
           "share_three": sum(r["share_three"] * r["reads"] for r in rows if r.get("share_three") is not None) / reads if reads else None}
    with_q = [r for r in rows if "p_mid" in r]
    for name in ("start", "mid", "end"):
        weights = sum(r["between_primers"] for r in with_q)
        out[f"q_{name}"] = mean_q(sum(r[f"p_{name}"] * r["between_primers"] for r in with_q) / weights) if weights else None
    medians = [r["len_median"] for r in rows if "len_median" in r]
    out["len_median_of_medians"] = float(np.median(medians)) if medians else None
    out["share_within_10pct"] = (sum(r["share_within_10pct"] * r["between_primers"] for r in rows if "len_median" in r) / n) if n else None
    out["distinct_q_min"] = min((r["distinct_q"] for r in rows), default=None)
    out["max_q_min"] = min((r["max_q"] for r in rows if r["max_q"] is not None), default=None)
    return out


def table_lines(rows: list, grid) -> list:
    head = ["sample", "reads", "both%", "rev%", "median", "mode", "in10%", ">1.5x%", "EE med", "EE/kb"] + [f"<=EE{g:g}%" for g in (2, 5, 10) if g in grid] + ["Qs", "Qmax", "Q 5'", "Q mid", "Q 3'"]
    lines = ["| " + " | ".join(head) + " |", "|" + "|".join("---" for _ in head) + "|"]
    f = lambda r, k, scale=1.0, fmt="%.0f": "-" if r.get(k) is None else fmt % (r[k] * scale)       # noqa: E731
    for r in rows:
        cells = [r["sample"], f"{r['reads']:,}", f(r, "share_both", 100), f(r, "share_both_reverse", 100), f(r, "len_median"), f(r, "len_mode_bin"),
                 f(r, "share_within_10pct", 100), f(r, "share_over_1p5x", 100, "%.1f"), f(r, "ee_median", 1, "%.2f"), f(r, "ee_per_1000_median", 1, "%.2f")]
        cells += [f(r, f"keep_maxEE_{g:g}", 100) for g in (2, 5, 10) if g in grid]
        cells += [str(r["distinct_q"]), str(r["max_q"]), f(r, "q_start"), f(r, "q_mid"), f(r, "q_end")]
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def write_outputs(rows: list, out: Path, folder, with_primers: bool, grid) -> None:
    out.mkdir(parents=True, exist_ok=True)
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(out / "longread_describe.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval="")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    p = pooled(rows)
    L = [f"# Long-read description of {folder}", "",
         f"{p['files']} files, {p['reads']:,} reads" + (f"; {100 * p['share_both']:.1f}% have both primers" if with_primers and p["share_both"] is not None else "; no primers were given, so every read is used as it is") + ".",
         (f"The primer at the 5' end of the read is found in {100 * p['share_five']:.1f}% of the reads and the one at the 3' end in {100 * p['share_three']:.1f}% (insertions and deletions in a primer count as errors, up to the allowed share)." if with_primers and p["share_five"] is not None else ""),
         f"Reads between the primers: {p['between_primers']:,}; median of the file medians {p['len_median_of_medians'] or 0:.0f} nt; share within 10% of the median {100 * (p['share_within_10pct'] or 0):.1f}%.",
         f"Quality (mean of the error probabilities, as Phred): first 50 bases Q{p['q_start'] or 0:.1f}, middle half Q{p['q_mid'] or 0:.1f}, last 50 bases Q{p['q_end'] or 0:.1f}; fewest distinct quality values in a file: {p['distinct_q_min']}, lowest maximum: Q{p['max_q_min']}.",
         "", "Columns: both% reads with both primers; rev% of those, in the reverse orientation; median, mode: length between the primers (mode = commonest 10 nt bin); in10%: share within 10% of the median; "
         ">1.5x%: share longer than 1.5 times the median; EE med, EE/kb: expected errors per read and per 1,000 nt (median); <=EEn%: share of reads at or under maxEE n; Qs, Qmax: distinct quality values and the highest; "
         "Q 5', Q mid, Q 3': mean quality of the first 50 bases, the middle half and the last 50 bases.", ""] + table_lines(rows, grid) + ["", "## What is unusual", ""]
    unusual = [(r["sample"], f) for r in rows for f in flags(r)]
    L += [f"- {s}: {f}" for s, f in unusual] if unusual else ["- nothing"]
    (out / "LONGREAD_DESCRIBE.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def describe_folder(folder: Path, out: Path, primer_f=None, primer_r=None, max_reads: int = 0, window: int = 60, error_rate: float = 0.15,
                    grid=DEFAULT_GRID, log=print) -> list:
    files = find_files(folder)
    if not files:
        raise SystemExit(f"No FASTQ files (.fastq, .fastq.gz, .fq, .fq.gz) in {folder}")
    if bool(primer_f) != bool(primer_r):
        raise SystemExit("Give both --primer-f and --primer-r, or neither.")
    locator = Locator(primer_f, primer_r, window, error_rate) if primer_f else None
    rows = []
    for p in files:
        row = describe_file(p, locator, max_reads, grid)
        rows.append(row)
        log(f"{row['sample']}: {row['reads']:,} reads" + (f", {100 * row['share_both']:.0f}% with both primers" if locator and row["reads"] else "") + (f", median {row['len_median']:.0f} nt" if "len_median" in row else ""))
    write_outputs(rows, out, folder, locator is not None, grid)
    log(f"Saved: {out / 'longread_describe.tsv'} and {out / 'LONGREAD_DESCRIBE.md'}")
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("folder", type=Path, help="folder of single-end FASTQ files")
    ap.add_argument("--out", type=Path, default=Path("longread_describe"))
    ap.add_argument("--primer-f", help="forward primer (IUPAC codes allowed), 5'->3'")
    ap.add_argument("--primer-r", help="reverse primer, 5'->3' as ordered (it is reverse-complemented here)")
    ap.add_argument("--reads", type=int, default=0, metavar="N", help="use only the first N reads of each file (default: all)")
    ap.add_argument("--window", type=int, default=60, help="how far from the end of a read a primer may start (default 60)")
    ap.add_argument("--error-rate", type=float, default=0.15, help="share of a primer's bases that may be wrong, an insertion or a deletion counting as one (default 0.15, so 3 of 20)")
    ap.add_argument("--grid", type=float, nargs="+", default=list(DEFAULT_GRID), help="maxEE values for the retention columns")
    a = ap.parse_args(argv)
    describe_folder(a.folder, a.out, a.primer_f, a.primer_r, a.reads, a.window, a.error_rate, tuple(a.grid))
    return 0


if __name__ == "__main__":
    sys.exit(main())
