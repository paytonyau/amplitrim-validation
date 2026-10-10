"""Command line interface."""

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .config import parse_ee
from .pipeline import PipelineError, Settings, run
from .report import screen_summary, write_outputs
from .subsample import DEFAULT_MAX_PAIRS

EXHAUSTIVE_MAX_PAIRS = 10_000      # --exhaustive tests every total length, so it gets a smaller default subsample

DESCRIPTION = """Recommend DADA2 truncLenF / truncLenR for paired-end amplicon reads by testing candidate
settings on your own reads: how many pairs pass a DADA2-style quality filter and merge, and whether the
settings can merge your longest amplicons. The reads should have adapters and primers removed and should
not be quality-trimmed; their lengths may differ."""


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="amplitrim", description=DESCRIPTION)
    ap.add_argument("--version", action="version", version=f"amplitrim {__version__}")
    req = ap.add_argument_group("required")
    req.add_argument("--input-dir", required=True, type=Path,
                     help="folder of paired FASTQs named sample_R1/_R2.fastq.gz or "
                          "sample_S1_L001_R1_001.fastq.gz")
    req.add_argument("--work-dir", required=True, type=Path, help="folder for intermediate files")
    req.add_argument("--out-dir", required=True, type=Path, help="folder for the report")
    opt = ap.add_argument_group("optional")
    opt.add_argument("--marker", choices=["16s", "its"], default="16s",
                     help="16s: search truncation lengths (default). its: amplicon lengths vary, so evaluate DADA2's "
                          "ITS workflow with no truncation")
    opt.add_argument("--min-len", type=int, default=50, metavar="NT",
                     help="minLen for --marker its (default 50)")
    opt.add_argument("--min-overlap", type=int, default=12, metavar="NT",
                     help="minimum overlap to merge, as in your DADA2 mergePairs call (default 12)")
    opt.add_argument("--max-ee", type=parse_ee, default=2.0, metavar="X",
                     help="maxEE, as in your DADA2 filterAndTrim call (default 2); two values, for example 2,5, set the "
                          "forward and the reverse read separately, as maxEE = c(2, 5) does")
    opt.add_argument("--allow-unmergeable", type=float, default=0.002, metavar="FRAC",
                     help="share of the longest amplicons you accept not being able to merge (default 0.002 = 0.2%%, which cost at most 0.4 points of merged pairs on six test datasets and recovered one more strain of a mock community than 0.01; see the validation). Lower keeps longer "
                          "reads so long 16S regions merge; higher allows shorter, cleaner reads")
    opt.add_argument("--explore-coverage", action="store_true",
                     help="also find the best setting for several coverage levels and show what covering fewer of the longest "
                          "amplicons would allow (slower)")
    opt.add_argument("--min-retention", type=float, default=0.5, metavar="FRAC",
                     help="with --explore-max-ee: the share of pairs the quality-first maxEE suggestion must keep (default 0.5). "
                          "The suggestion is the lowest-noise setting that reaches it and never goes below DADA2's default of 2")
    opt.add_argument("--explore-max-ee", action="store_true",
                     help="also find the best truncation for several maxEE values and show what each keeps and the expected "
                          "errors of the kept reads, and where retention levels off (slower)")
    opt.add_argument("--max-pairs", type=int, default=None, metavar="N",
                     help=f"read pairs to analyse, taken evenly across all samples (default {DEFAULT_MAX_PAIRS:,}, or "
                          f"{EXHAUSTIVE_MAX_PAIRS:,} with --exhaustive; 0 = all)")
    opt.add_argument("--longest-insert", type=int, metavar="NT",
                     help="longest amplicon to cover, in nt (default: measured from the reads)")
    opt.add_argument("--yield-tolerance", type=float, metavar="X",
                     help="treat settings within X of the best yield as tied and pick the one with the "
                          "shortest overlap (default 0.01)")
    opt.add_argument("--step", type=int, default=2, metavar="NT",
                     help="step between total lengths tested (default 2)")
    opt.add_argument("--exhaustive", action="store_true",
                     help="test every total length up to the read limit instead of stopping when yield "
                          "levels off")
    opt.add_argument("--verbose", action="store_true", help="show progress and detailed messages")
    fig = ap.add_argument_group("FIGARO cross-check (optional; needs reads of one length)")
    fig.add_argument("--figaro-script", type=Path,
                     help="path to figaro.py; if given, FIGARO's suggested trim pairs are tested too")
    fig.add_argument("--figaro-python", default=sys.executable, metavar="PYTHON",
                     help="Python used to run FIGARO (default: the current one)")
    fig.add_argument("--top-n", type=int, default=3, metavar="N",
                     help="FIGARO trim pairs tested per total length (default 3)")
    return ap


def resolve_max_pairs(args) -> int:
    if args.max_pairs is not None:
        return args.max_pairs
    return EXHAUSTIVE_MAX_PAIRS if args.exhaustive else DEFAULT_MAX_PAIRS


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING, format="%(message)s")
    settings = Settings(
        input_dir=args.input_dir, work_dir=args.work_dir, min_overlap=args.min_overlap, max_ee=args.max_ee, explore_max_ee=args.explore_max_ee, min_retention=args.min_retention, max_unmergeable=args.allow_unmergeable, explore_coverage=args.explore_coverage,
        max_pairs=resolve_max_pairs(args), longest_insert=args.longest_insert, yield_tolerance=args.yield_tolerance,
        step=args.step, exhaustive=args.exhaustive, figaro_script=args.figaro_script,
        figaro_python=args.figaro_python, top_n=args.top_n,
        splits="both" if args.figaro_script else "quality", marker=args.marker, min_len=args.min_len)
    print(f"AmpliTrim {__version__}\nAnalysing reads...", flush=True)
    try:
        result = run(settings, say=(lambda msg: print(msg, flush=True)) if args.verbose else (lambda msg: None))
    except (PipelineError, RuntimeError) as err:
        print(f"Error: {err}", file=sys.stderr)
        return 3
    report = write_outputs(result, args.out_dir)
    print("\n".join(screen_summary(result, report)))
    return 0 if result.recommendation else 2


if __name__ == "__main__":
    sys.exit(main())
