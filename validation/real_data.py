#!/usr/bin/env python3
"""
Run AmpliTrim (and optionally DADA2) on one dataset and append a row to <out>/results.tsv, so that many
public datasets can be compared in one table.

A dataset is either a list of ENA/SRA run accessions (one per line, or a Run Selector CSV) or a folder of local
paired FASTQs:

  python validation/real_data.py runs.txt --name study1 --out results \
         --primer-f GTGYCAGCMGCCGCGGTAA --primer-r GGACTACNVGGGTWTCTAAT --dada2 --compare 240,160
  python validation/real_data.py --local reads/ --name study2 --out results --dada2
  python validation/real_data.py --local reads/ --name study4 --out results --adapter-f AGATCGGAAGAGC \
         --adapter-r AGATCGGAAGAGC          (adapters only, primers already removed)
  python validation/real_data.py --local reads/ --name study3 --out results --mixed-orientation \
         --primer-f GTGYCAGCMGCCGCGGTAA --primer-r GGACTACNVGGGTWTCTAAT       (pairs in both orientations)

Primers and adapters are removed with cutadapt; the reads are NOT quality-trimmed. Without --primer-f and
--primer-r the reads are used as they are (they must already be free of primers and adapters).
Needs: amplitrim (pip install .), vsearch on PATH, cutadapt (only with primers), R + dada2 (only with --dada2).
The download step needs access to www.ebi.ac.uk.
"""
import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from amplitrim.subsample import find_pairs  # noqa: E402
from download_ena import download, read_accessions  # noqa: E402
from rtools import NOT_FOUND, find_rscript  # noqa: E402

COMPLEMENT = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
HERE = Path(__file__).parent


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


PRIMER_MIN_OVERLAP = 8      # a primer must be present for at least this many bases to count: cutadapt counts a read as trimmed when ANY adapter is found, and
                            # a match of 3 bases happens by chance often enough (about 1 in 64 per end) to keep untrimmed, wrongly oriented pairs


def _cutadapt(in1: Path, in2: Path, out1: Path, out2: Path, primer_f, primer_r, adapter_f=None, adapter_r=None,
              error_rate: float = 0.1) -> None:
    """Remove primers (and anything before them), the other primer's read-through and 3' adapters. With primers, pairs without them are
    discarded. Never quality-trims. With primers this is two passes, because they are different jobs: the first FINDS the primers at the
    starts of the reads (at least PRIMER_MIN_OVERLAP bases, so that chance matches cannot pass a pair) and keeps only pairs that have both;
    the second removes the read-through at the 3' ends, as little as 3 bases of it, and discards nothing."""
    import tempfile
    base = [sys.executable, "-m", "cutadapt", "-e", str(error_rate)]
    three = (["-a", revcomp(primer_r), "-A", revcomp(primer_f)] if primer_f else []) + (["-a", adapter_f, "-A", adapter_r] if adapter_f else [])
    times = ["--times", "2"] if primer_f and adapter_f else []          # a read-through primer and an adapter behind it: two removals
    def run(cmd):
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise SystemExit(f"cutadapt failed on {in1.name}:\n{res.stderr[-500:]}")
    if not primer_f:
        run(base + ["--minimum-length", "50"] + three + ["-o", str(out1), "-p", str(out2), str(in1), str(in2)])
        return
    with tempfile.TemporaryDirectory() as tmp:
        m1, m2 = Path(tmp) / "m1.fastq.gz", Path(tmp) / "m2.fastq.gz"
        run(base + ["-O", str(PRIMER_MIN_OVERLAP), "-g", primer_f, "-G", primer_r, "--discard-untrimmed", "-o", str(m1), "-p", str(m2), str(in1), str(in2)])
        run(base + ["--minimum-length", "50"] + times + three + ["-o", str(out1), "-p", str(out2), str(m1), str(m2)])


def trim_reads(raw_dir: Path, out_dir: Path, primer_f=None, primer_r=None, adapter_f=None, adapter_r=None,
               mixed: bool = False, error_rate: float = 0.1) -> None:
    """Write <sample>_R1/_R2.fastq.gz with primers and/or adapters removed. With `mixed`, pairs sequenced in
    either orientation are found in two passes (the second with the files swapped) and put the same way round."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for sample, r1, r2 in find_pairs(raw_dir):
        o1, o2 = out_dir / f"{sample}_R1.fastq.gz", out_dir / f"{sample}_R2.fastq.gz"
        if not mixed:
            _cutadapt(r1, r2, o1, o2, primer_f, primer_r, adapter_f, adapter_r, error_rate)
            continue
        part = out_dir / "_parts"
        part.mkdir(exist_ok=True)
        _cutadapt(r1, r2, part / "a1.fastq.gz", part / "a2.fastq.gz", primer_f, primer_r, adapter_f, adapter_r, error_rate)
        _cutadapt(r2, r1, part / "b1.fastq.gz", part / "b2.fastq.gz", primer_f, primer_r, adapter_f, adapter_r, error_rate)
        for out, x, y in ((o1, "a1", "b1"), (o2, "a2", "b2")):
            with open(out, "wb") as fh:                      # concatenated gzip files are valid gzip
                for name in (x, y):
                    fh.write((part / f"{name}.fastq.gz").read_bytes())
        shutil.rmtree(part)


def thin_reads(reads: Path, out_dir: Path, per_sample: int) -> Path:
    """Copy at most `per_sample` evenly spaced pairs of every sample, so that DADA2 stays tractable."""
    import gzip
    import math
    from amplitrim.subsample import count_pairs, read_fastq
    out_dir.mkdir(parents=True, exist_ok=True)
    for sample, r1, r2 in find_pairs(reads):
        stride = max(1, math.ceil(count_pairs(r1) / per_sample))
        with gzip.open(out_dir / f"{sample}_R1.fastq.gz", "wt", newline="\n") as o1, \
                gzip.open(out_dir / f"{sample}_R2.fastq.gz", "wt", newline="\n") as o2:
            for i, ((h1, s1, p1, q1), (h2, s2, p2, q2)) in enumerate(zip(read_fastq(r1), read_fastq(r2))):
                if i % stride == 0:
                    o1.write(f"{h1}\n{s1}\n{p1}\n{q1}\n")
                    o2.write(f"{h2}\n{s2}\n{p2}\n{q2}\n")
    return out_dir


def ee_text(value: str) -> str:
    """Validate a --max-ee value ('2' or '2,5') and return it in a normal form."""
    import argparse
    from amplitrim.config import ee_pair, parse_ee
    try:
        f, r = ee_pair(parse_ee(value))
    except argparse.ArgumentTypeError as e:
        raise argparse.ArgumentTypeError(str(e)) from None
    return f"{f:g}" if f == r else f"{f:g},{r:g}"


IUPAC_DNA = set("ACGTURYSWKMBDHVN")


def check_sequence(label: str, seq) -> None:
    """Stop early, in plain words, if a primer or adapter is not a DNA sequence (a placeholder such as FWD, say)."""
    if seq and not set(seq.upper()) <= IUPAC_DNA:
        raise SystemExit(f"{label} '{seq}' is not a DNA sequence. Give the letters of the primer itself (for example "
                         f"GTGYCAGCMGCCGCGGTAA), not a placeholder. Allowed letters: A C G T and the ambiguity codes R Y S W K M B D H V N.")


def dada2_blocker(quality_values):
    """Why DADA2 cannot be run on these reads, or None."""
    if quality_values is not None and quality_values <= 2:
        return (f"the reads have only {quality_values} distinct quality values, so DADA2 cannot learn an error model from them. "
                "AmpliTrim's result above rests on the same flattened quality; get the original-quality files to test DADA2.")
    return None


def keep_only(reads_dir: Path, spec: str) -> list:
    """Delete from the working copy every sample not named in `spec`: comma-separated sample names, or a text file with one name per
    line. A name that is not among the samples stops the run, so that a typo cannot silently drop a sample. Returns the names dropped."""
    p = Path(spec)
    names = [ln.strip() for ln in p.read_text(encoding="utf-8-sig").splitlines()] if p.is_file() else spec.split(",")
    names = [n.strip() for n in names if n.strip() and not n.strip().startswith("#")]
    have = {s: (r1, r2) for s, r1, r2 in find_pairs(reads_dir)}
    missing = [n for n in names if n not in have]
    if missing:
        raise SystemExit(f"--only names samples that are not in the reads: {', '.join(missing)}. "
                         f"Samples found: {', '.join(sorted(have)[:6])}{' ...' if len(have) > 6 else ''}")
    dropped = []
    for sample, (r1, r2) in have.items():
        if sample not in names:
            r1.unlink()
            r2.unlink()
            dropped.append(sample)
    return dropped


def drop_small_samples(reads_dir: Path, minimum: int) -> list:
    """Delete the pairs of every sample with fewer than `minimum` read pairs (controls, failed runs) from the working copy.
    Returns [(sample, pairs)] for the dropped samples."""
    from amplitrim.subsample import count_pairs
    dropped = []
    for sample, r1, r2 in find_pairs(reads_dir):
        n = count_pairs(r1)
        if n < minimum:
            r1.unlink()
            r2.unlink()
            dropped.append((sample, n))
    return dropped


def summary_lines(row: dict, table: Path) -> list:
    """What the runner prints at the end: AmpliTrim's recommendation, and DADA2's result when it was run."""
    if row["truncLenF"] is None:
        lines = [f"{row['dataset']}: AmpliTrim made no recommendation."]
    else:
        lines = [f"{row['dataset']}: AmpliTrim {row['truncLenF']}/{row['truncLenR']}, estimated yield {row['estimated_yield']}"]
    if "dada2_merged" in row:
        lines.append(f"  DADA2 at {row['truncLenF']}/{row['truncLenR']}: {row['dada2_merged']:.1%} of read pairs merged "
                     f"({row['dada2_nonchim']:.1%} after chimera removal), {row['dada2_asvs']} ASVs")
        if "dada2_denoised_f" in row:
            lines.append(f"    where the reads went: {row['dada2_filtered']:.1%} passed the filter, {row['dada2_denoised_f']:.1%} are in "
                         f"denoised forward sequences ({row['dada2_denoised_r']:.1%} reverse), {row['dada2_merged']:.1%} merged; "
                         f"{row['dada2_singletons']:.1%} of the filtered forward reads are singletons before denoising")
    if "compare_dada2_merged" in row:
        lines.append(f"  DADA2 at {row['compare_trunc'].replace(',', '/')}: {row['compare_dada2_merged']:.1%} of read pairs merged, "
                     f"{row['compare_asvs']} ASVs")
    if row.get("suggested_max_ee"):
        lines.append(f"  maxEE, quality first (keeping at least {row['min_retention']:.0%} of pairs): {row['suggested_max_ee']} with truncLen "
                     f"{row['suggested_trunc']}, estimated yield {row['suggested_yield']:.1%}" + ("" if row["suggested_reached"] else " (the floor was not reached)")
                     + f"; yours {row['estimated_yield']:.1%}")
    if row.get("dada2_skipped"):
        lines.append(f"  DADA2 not run: {row['dada2_skipped']}")
    if row.get("dada2_per_sample_table"):
        lines.append(f"  DADA2 per sample: {row['dada2_per_sample_table']}")
    for key in ("dada2_asv_fasta", "compare_asv_fasta"):
        if row.get(key):
            lines.append(f"  DADA2 ASVs (for mock_check.py): {row[key]}")
    lines.append(f"  row added to {table}")
    return lines


def normalise_names(raw_dir: Path, out_dir: Path) -> None:
    """Copy the pairs as <sample>_R1/_R2.fastq.gz, whatever the original file names were."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for sample, r1, r2 in find_pairs(raw_dir):
        for src, mate in ((r1, "R1"), (r2, "R2")):
            shutil.copy(src, out_dir / f"{sample}_{mate}.fastq.gz")


class Dada2Error(Exception):
    """DADA2 ran and failed; the message says why in plain words."""


def explain_dada2_failure(stderr: str) -> str:
    if "Error matrix is NULL" in stderr:
        return ("DADA2 could not learn an error model from these quality scores (it stopped with 'Error matrix is NULL'). "
                "That happens when the quality scores take too few distinct values; AmpliTrim's own result above is unaffected.")
    return f"DADA2 failed: {stderr.strip()[-300:]}"


def dada2_run(reads: Path, trunc_f: int, trunc_r: int, max_ee, min_overlap: int, min_len: int = 0, track: Path = None, asvs: Path = None) -> dict:
    rscript = find_rscript()
    if rscript is None:
        raise SystemExit(NOT_FOUND)
    try:
        res = subprocess.run([rscript, str(HERE / "dada2_real.R"), str(reads), str(trunc_f), str(trunc_r),
                              str(max_ee), str(min_overlap), str(min_len)]
                             + ([str(track) if track else ""] if (track or asvs) else []) + ([str(asvs)] if asvs else []),
                             capture_output=True, text=True)
    except FileNotFoundError:
        raise SystemExit("Rscript was not found. Install R (https://cran.r-project.org), put its bin folder on PATH, and install DADA2 (see the Windows steps in validation/README.md); python validation/check_environment.py shows what is missing.") from None
    if res.returncode != 0:
        raise Dada2Error(explain_dada2_failure(res.stderr))
    f = res.stdout.strip().splitlines()[-1].split("\t")
    out = {"filtered": float(f[0]), "merged": float(f[1]), "nonchim": float(f[2]), "asvs": int(f[3])}
    if len(f) >= 7:                                  # where the reads went
        out.update({"denoised_f": float(f[4]), "denoised_r": float(f[5]), "singletons": float(f[6])})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="?", type=Path, help="file with one ENA/SRA run accession per line")
    ap.add_argument("--local", type=Path, help="folder of paired FASTQs instead of accessions")
    ap.add_argument("--name", help="dataset label (default: the file or folder name)")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--primer-f"), ap.add_argument("--primer-r")
    ap.add_argument("--error-rate", type=float, default=0.1,
                    help="cutadapt error rate for primers and adapters (default 0.1 allows 1 error in a 15-18 nt primer; 0.15 allows 2)")
    ap.add_argument("--adapter-f", help="3' adapter to remove from R1 (e.g. AGATCGGAAGAGC for TruSeq)")
    ap.add_argument("--adapter-r", help="3' adapter to remove from R2")
    ap.add_argument("--mixed-orientation", action="store_true",
                    help="pairs were sequenced in both orientations (needs the primers); they are put the same way round")
    ap.add_argument("--min-overlap", type=int, default=12)
    ap.add_argument("--max-ee", type=ee_text, default="2",
                    help="maxEE: one number, or F,R for different forward and reverse limits (default 2)")
    ap.add_argument("--explore-max-ee", action="store_true", help="also find the best truncation for several maxEE values and show what each keeps and costs")
    ap.add_argument("--explore-coverage", action="store_true", help="also show what covering fewer of the longest amplicons would allow")
    ap.add_argument("--allow-unmergeable", type=float, default=None, metavar="FRAC", help="share of the longest amplicons you accept not being able to merge (default 0.01)")
    ap.add_argument("--only", metavar="NAMES", help="analyse only these samples: comma-separated names, or a text file with one name per line (leave out controls and low-input samples)")
    ap.add_argument("--min-retention", type=float, default=None, metavar="FRAC", help="with --explore-max-ee: the share of pairs the quality-first maxEE suggestion must keep (default 0.5)")
    ap.add_argument("--min-reads", type=int, default=0, metavar="N",
                    help="leave out samples with fewer than N read pairs after trimming (negative controls, failed runs)")
    ap.add_argument("--marker", choices=["16s", "its"], default="16s")
    ap.add_argument("--min-len", type=int, default=50, help="minLen for --marker its")
    ap.add_argument("--dada2", action="store_true", help="also run DADA2 at the recommended settings")
    ap.add_argument("--dada2-pairs-per-sample", type=int, default=0, metavar="N",
                    help="give DADA2 at most N evenly spaced pairs per sample (default: all; use ~3000 for studies with many samples)")
    ap.add_argument("--compare", metavar="F,R", help="also run DADA2 at this truncLenF,truncLenR")
    args = ap.parse_args()
    for label, seq in (("--primer-f", args.primer_f), ("--primer-r", args.primer_r),
                       ("--adapter-f", args.adapter_f), ("--adapter-r", args.adapter_r)):
        check_sequence(label, seq)
    if bool(args.runs) == bool(args.local):
        ap.error("give either a file of run accessions or --local DIR")
    if bool(args.primer_f) != bool(args.primer_r):
        ap.error("give both --primer-f and --primer-r, or neither")
    if bool(args.adapter_f) != bool(args.adapter_r):
        ap.error("give both --adapter-f and --adapter-r, or neither")
    name = args.name or (args.runs.stem if args.runs else args.local.name)
    work = args.out / name
    work.mkdir(parents=True, exist_ok=True)

    if args.runs:
        raw = work / "raw"
        raw.mkdir(exist_ok=True)
        for acc in read_accessions(str(args.runs)):
            download(acc, raw)
    else:
        raw = args.local
    if args.mixed_orientation and not args.primer_f:
        ap.error("--mixed-orientation needs --primer-f and --primer-r")
    reads = work / "ready"
    if args.primer_f or args.adapter_f:
        trim_reads(raw, reads, args.primer_f, args.primer_r, args.adapter_f, args.adapter_r, args.mixed_orientation,
                    args.error_rate)
    else:
        normalise_names(raw, reads)

    if args.only:
        left_out = keep_only(reads, args.only)
        print(f"Analysing only the {len(list(find_pairs(reads)))} sample(s) named in --only; left out {len(left_out)}.")
    dropped = drop_small_samples(reads, args.min_reads) if args.min_reads else []
    if dropped:
        print(f"Left out {len(dropped)} sample(s) with fewer than {args.min_reads} read pairs: "
              + ", ".join(f"{s} ({n})" for s, n in dropped))
        if not list(find_pairs(reads)):
            raise SystemExit(f"Every sample has fewer than {args.min_reads} read pairs; lower --min-reads.")
    res = subprocess.run([sys.executable, "-m", "amplitrim.cli", "--input-dir", str(reads), "--work-dir", str(work / "work"),
                          "--out-dir", str(work / "report"), "--min-overlap", str(args.min_overlap),
                          "--max-ee", args.max_ee, "--marker", args.marker,
                          "--min-len", str(args.min_len)] + (["--explore-max-ee"] if args.explore_max_ee else [])
                         + (["--min-retention", str(args.min_retention)] if args.min_retention is not None else [])
                         + (["--explore-coverage"] if args.explore_coverage else [])
                         + (["--allow-unmergeable", str(args.allow_unmergeable)] if args.allow_unmergeable is not None else []),
                         capture_output=True, text=True)
    if res.returncode not in (0, 2):
        raise SystemExit(f"AmpliTrim failed:\n{res.stderr[-600:]}")
    summary = json.loads((work / "report" / "summary.json").read_text())
    rec, data = summary["recommended"], summary["data"]
    sug = summary.get("maxEE_exploration")
    row = {"dataset": name, "pairs": data["pairs_analysed"], "read_f": data["read_length_f"],
           "read_r": data["read_length_r"], "amplicon_min": data["amplicon_min"], "amplicon_median": data["amplicon_median"],
           "amplicon_max": data["amplicon_max"], "merged_at_full_length": data["share_merged_at_full_length"],
           "truncLenF": rec and rec["truncLenF"], "truncLenR": rec and rec["truncLenR"],
           "estimated_yield": rec and rec["estimated_yield"]}
    if sug:
        row.update({"suggested_max_ee": ",".join(f"{x:g}" for x in sug["suggested"]["maxEE"]), "suggested_trunc": f"{sug['suggested']['truncLenF']}/{sug['suggested']['truncLenR']}",
                    "suggested_yield": sug["suggested"]["estimated_yield"], "suggested_reached": sug["suggested"]["reachesFloor"],
                    "min_retention": sug["minRetention"]}) if sug.get("suggested") else None
    dada_reads = thin_reads(reads, work / "dada2_reads", args.dada2_pairs_per_sample) if args.dada2_pairs_per_sample else reads
    blocked = dada2_blocker(data.get("quality_values")) if args.dada2 else None
    if blocked:
        row["dada2_skipped"] = blocked
    if args.dada2 and rec and not blocked:
        track = work / f"dada2_per_sample_{rec['truncLenF']}x{rec['truncLenR']}.tsv"
        asvs = work / f"asvs_{rec['truncLenF']}x{rec['truncLenR']}.fasta"
        try:
            d = dada2_run(dada_reads, rec["truncLenF"], rec["truncLenR"], args.max_ee, args.min_overlap,
                          args.min_len if args.marker == "its" else 0, track, asvs)
        except Dada2Error as e:
            row["dada2_skipped"], d = str(e), None
        if d:
            row["dada2_per_sample_table"], row["dada2_asv_fasta"] = str(track), str(asvs)
            row.update({"dada2_filtered": d["filtered"], "dada2_merged": d["merged"], "dada2_nonchim": d["nonchim"], "dada2_asvs": d["asvs"]})
            if "singletons" in d:
                row.update({"dada2_denoised_f": d["denoised_f"], "dada2_denoised_r": d["denoised_r"], "dada2_singletons": d["singletons"]})
    if args.compare and not row.get("dada2_skipped"):
        f, r = (int(x) for x in args.compare.split(","))
        asvs_c = work / f"asvs_{f}x{r}.fasta"
        try:
            d = dada2_run(dada_reads, f, r, args.max_ee, args.min_overlap, asvs=asvs_c)
            row.update({"compare_trunc": args.compare, "compare_dada2_merged": d["merged"], "compare_asvs": d["asvs"], "compare_asv_fasta": str(asvs_c)})
        except Dada2Error as e:
            row["dada2_skipped"] = str(e)
    table = args.out / "results.tsv"
    new = not table.exists()
    with open(table, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row), delimiter="\t")
        if new:
            w.writeheader()
        w.writerow(row)
    print("\n".join(summary_lines(row, table)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
