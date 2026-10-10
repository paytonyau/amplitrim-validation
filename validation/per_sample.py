#!/usr/bin/env python3
"""
Per-sample check of the DADA2 settings AmpliTrim recommended: for every sample, the share of read pairs that pass
filtering and merge. AmpliTrim's own estimate pools all samples, so a few weak samples can hide in it.

    python validation/per_sample.py C:\\data\\results\\study1             (the folder real_data.py created)
    python validation/per_sample.py --reads C:\\data\\ready --trunc-f 220 --trunc-r 227
    python validation/per_sample.py C:\\data\\results\\study1 --runinfo SraRunInfo.csv     (adds sample names and groups)

Only 16S-style truncation is handled (not --marker its).
"""
import argparse
import csv
import json
import math
import re
import statistics
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from amplitrim.config import Config, parse_ee  # noqa: E402
from amplitrim.pipeline import _evaluate  # noqa: E402
from amplitrim.subsample import count_pairs, find_pairs, read_fastq  # noqa: E402
from amplitrim.trimmer import ReadPairs  # noqa: E402


def evenly_spaced_pairs(p1: Path, p2: Path, n: int) -> ReadPairs:
    stride = max(1, math.ceil(count_pairs(p1) / n))
    keep = [(a, b) for i, (a, b) in enumerate(zip(read_fastq(p1), read_fastq(p2))) if i % stride == 0]
    return ReadPairs([a for a, _ in keep], [b for _, b in keep])


def per_sample(reads: Path, trunc_f: int, trunc_r: int, pairs_per_sample: int = 4000, min_overlap: int = 12,
               max_ee=2.0):
    cfg, rows = Config(dada2_min_overlap=min_overlap, max_ee=max_ee), []
    with tempfile.TemporaryDirectory() as tmp:
        for sample, p1, p2 in find_pairs(reads):
            pairs = evenly_spaced_pairs(p1, p2, pairs_per_sample)
            if pairs.n == 0:
                continue
            c = _evaluate(pairs, trunc_f, trunc_r, cfg, 0, Path(tmp) / sample, "", (pairs.width_f, pairs.width_r))
            rows.append({"sample": sample, "pairs_used": pairs.n, "pass_filter": round(c.ee_pass_rate, 4),
                         "merge_rate": round(c.merge_rate_tolerant, 4), "estimated_yield": round(c.yield_est, 4),
                         # DADA2's own rule (no mismatch allowed in the overlap) applied to the raw reads, so sequencing errors count against it
                         "strict_merge_rate": round(c.merge_rate_strict, 4), "strict_yield": round(c.ee_pass_rate * c.merge_rate_strict, 4),
                         "mean_overlap": round(c.mean_overlap, 1) if c.mean_overlap == c.mean_overlap else ""})
    return rows


def load_design(csv_path: Path) -> dict:
    """{run accession: (sample name, read pairs in the archive)} from an SRA Run Selector CSV (Run, SampleName, spots) or an ENA report
    TSV (run_accession, sample_alias, read_count). ENA counts the reads of both mates, so its pairs are read_count / 2."""
    import io
    text = csv_path.read_text(encoding="utf-8-sig")
    delimiter = "\t" if "\t" in text.splitlines()[0] else ","
    design = {}
    for r in csv.DictReader(io.StringIO(text), delimiter=delimiter):
        run = r.get("Run") or r.get("run_accession")
        name = r.get("SampleName") or r.get("sample_alias") or r.get("LibraryName") or r.get("library_name") or ""
        spots = int(r["spots"]) if r.get("spots") else (int(r["read_count"]) // 2 if r.get("read_count") else 0)
        if run:
            design[run] = (name, spots)
    return design


def group_of(sample_name: str, pattern: str = r"^(.*?)L\d+$") -> str:
    """The sample name without its replicate marker (G1''L2 -> G1''); the whole name when the pattern does not match."""
    m = re.match(pattern, sample_name)
    return m.group(1) if m and m.group(1) else sample_name


def add_design(rows, design: dict, pattern: str = r"^(.*?)L\d+$"):
    for r in rows:
        name, spots = design.get(r["sample"], ("", 0))
        r["sample_name"], r["archive_reads"], r["group"] = name, spots, group_of(name, pattern) if name else ""
    return rows


def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    out, i = [0.0] * len(values), 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            out[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return out


def spearman(x, y):
    """Rank correlation (ties get average ranks); None when it cannot be computed."""
    if len(x) < 3 or len(set(x)) < 2 or len(set(y)) < 2:
        return None
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else None


def add_dada2(rows, track_path: Path):
    """Add DADA2's per-sample merged share (of the input pairs it was given) from the table real_data.py wrote."""
    with open(track_path, encoding="utf-8", newline="") as fh:
        table = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for r in rows:
        t = table.get(r["sample"])
        r["dada2_merged"] = float(t["merged"]) / float(t["input"]) if t and float(t["input"]) else None
        r["dada2_input"] = int(float(t["input"])) if t else None                 # the pairs DADA2 was given (it may be thinned)
    return rows


def summarise_dada2(rows):
    have = [r for r in rows if r.get("dada2_merged") is not None]
    if len(have) < 2:
        return ["No DADA2 per-sample results matched these samples."]
    d = sorted(r["dada2_merged"] for r in have)
    lines = [f"DADA2 merged {d[0]:.1%} to {d[-1]:.1%} of the pairs it was given across {len(have)} samples (median {statistics.median(d):.1%}); "
             f"AmpliTrim's estimates span {min(r['estimated_yield'] for r in have):.1%} to {max(r['estimated_yield'] for r in have):.1%}."]
    if len(have) < 10:
        lines.append("Too few samples for a meaningful rank correlation (it needs at least 10).")
        return lines
    rho = spearman([r["estimated_yield"] for r in have], [r["dada2_merged"] for r in have])
    if rho is not None:
        lines.append(f"Rank correlation between AmpliTrim's estimate and DADA2's merged share, sample by sample (n = {len(have)}): {rho:+.2f}.")
    if all("strict_yield" in r for r in have):
        rho_strict = spearman([r["strict_yield"] for r in have], [r["dada2_merged"] for r in have])
        if rho_strict is not None:
            lines.append(f"Rank correlation between the strict estimate (no mismatch allowed in the overlap, on raw reads) and DADA2's merged share: {rho_strict:+.2f}.")
    rho_depth = spearman([r["dada2_input"] for r in have], [r["dada2_merged"] for r in have])
    if rho_depth is not None:
        lines.append(f"Rank correlation between the number of pairs DADA2 was given and its merged share: {rho_depth:+.2f}.")
    return lines


def summarise_groups(rows):
    """Estimated yield by group: samples, median, lowest, highest."""
    groups = {}
    for r in rows:
        if r.get("group"):
            groups.setdefault(r["group"], []).append(r["estimated_yield"])
    return [f"{g}: {len(v)} samples, median {statistics.median(v):.1%}, range {min(v):.1%}-{max(v):.1%}"
            for g, v in sorted(groups.items())]


def summarise(rows, pooled=None):
    if not rows:
        return ["No samples."]
    ys = sorted(r["estimated_yield"] for r in rows)
    med = statistics.median(ys)
    weighted = sum(r["estimated_yield"] * r["pairs_used"] for r in rows) / sum(r["pairs_used"] for r in rows)
    p10 = ys[int(0.1 * (len(ys) - 1))]
    lines = [f"{len(rows)} samples: estimated yield median {med:.1%}, 10th percentile {p10:.1%}, "
             f"lowest {ys[0]:.1%}, highest {ys[-1]:.1%}; pair-weighted mean {weighted:.1%}"
             + (f" (AmpliTrim's pooled estimate: {pooled:.1%})" if pooled is not None else "") + "."]
    weak = [r for r in rows if r["estimated_yield"] < med - 0.10]
    if weak:
        lines.append(f"{len(weak)} sample{'s' if len(weak) != 1 else ''} more than 10 points below the median: "
                     + ", ".join(f"{r['sample']} ({r['estimated_yield']:.0%})" for r in sorted(weak, key=lambda r: r['estimated_yield'])[:8])
                     + ".")
    else:
        lines.append("No sample is more than 10 points below the median.")
    return lines


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("workdir", nargs="?", type=Path, help="folder created by real_data.py (has ready/ and report/)")
    ap.add_argument("--reads", type=Path, help="folder of primer-free paired FASTQs (instead of workdir)")
    ap.add_argument("--trunc-f", type=int)
    ap.add_argument("--trunc-r", type=int)
    ap.add_argument("--pairs", type=int, default=4000, help="pairs per sample to evaluate (default 4000)")
    ap.add_argument("--min-overlap", type=int, default=12)
    ap.add_argument("--max-ee", type=parse_ee, default=2.0, help="maxEE: one number, or F,R")
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--runinfo", type=Path, help="SRA Run Selector CSV (SraRunInfo.csv): adds sample names, archive read counts and groups")
    ap.add_argument("--dada2-track", type=Path, help="per-sample DADA2 table written by real_data.py --dada2 (dada2_per_sample_*.tsv)")
    ap.add_argument("--group-regex", default=r"^(.*?)L\d+$", help="how to cut the replicate off a sample name to get its group (default: before L<number>)")
    args = ap.parse_args()
    pooled = None
    if args.workdir:
        if not (args.workdir / "report" / "summary.json").exists():
            print(f"{args.workdir} has no AmpliTrim result yet. Run real_data.py first; it creates this folder "
                  f"(--out/--name) and the report in it.")
            return 1
        reads = args.workdir / "ready"
        summary = json.loads((args.workdir / "report" / "summary.json").read_text(encoding="utf-8"))
        rec = summary["recommended"]
        if rec is None:
            print("AmpliTrim made no recommendation for this dataset.")
            return 1
        tf, tr, pooled = rec["truncLenF"], rec["truncLenR"], rec["estimated_yield"]
    elif args.reads and args.trunc_f and args.trunc_r:
        reads, tf, tr = args.reads, args.trunc_f, args.trunc_r
    else:
        ap.error("give the folder real_data.py created, or --reads with --trunc-f and --trunc-r")
    rows = per_sample(reads, tf, tr, args.pairs, args.min_overlap, args.max_ee)
    if args.runinfo:
        add_design(rows, load_design(args.runinfo), args.group_regex)
    if args.dada2_track:
        add_dada2(rows, args.dada2_track)
    print(f"settings tested: truncLenF {tf}, truncLenR {tr}\n")
    named = bool(args.runinfo)
    print((f"{'sample':14s} {'name':10s} " if named else f"{'sample':28s} ")
          + f"{'pairs':>6s} {'pass':>7s} {'merge':>7s} {'yield':>7s} {'strict':>7s} {'overlap':>8s}" + f" {'archive reads':>13s}" * named
          + f" {'DADA2 merged':>12s}" * bool(args.dada2_track))
    for r in rows:
        print((f"{r['sample'][:14]:14s} {r['sample_name'][:10]:10s} " if named else f"{r['sample'][:28]:28s} ")
              + f"{r['pairs_used']:>6d} {r['pass_filter']:>7.1%} {r['merge_rate']:>7.1%} {r['estimated_yield']:>7.1%} {r['strict_yield']:>7.1%} {r['mean_overlap']:>8}"
              + (f" {r['archive_reads']:>13,d}" if named else "")
              + (f" {r['dada2_merged']:>12.1%}" if args.dada2_track and r.get("dada2_merged") is not None else ""))
    print()
    for line in summarise(rows, pooled):
        print("-", line)
    if args.dada2_track:
        print()
        for line in summarise_dada2(rows):
            print("-", line)
    if named:
        print("\nBy group:")
        for line in summarise_groups(rows):
            print("  " + line)
    if args.csv and rows:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
