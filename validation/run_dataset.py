#!/usr/bin/env python3
"""Run the whole AmpliTrim test protocol on ONE dataset and collect everything into one text file to send back.
Works on Windows, Linux and macOS (no Bash needed).

    python validation/run_dataset.py NAME RAW_DIR [options]

    python validation\\run_dataset.py mockA C:\\data\\mockA\\raw --primer-f AGAGTTTGATYMTGGCTCAG --primer-r ATTACCGCGGCTGCTGG ^
           --error-rate 0.15 --compare 279,263 --pairs 100000 --out C:\\data\\results

Output: <out>/NAME_outcome.txt (send that one file). Needs amplitrim installed, vsearch, cutadapt (with primers) and R + DADA2.
"""
import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from rtools import find_rscript  # noqa: E402


def build_id() -> str:
    """A short fingerprint of the code that is running (the package and the validation scripts), so that an outcome file shows
    exactly which build produced it. It depends only on file contents, not on dates, paths or line endings."""
    h = hashlib.sha256()
    # a fixed order that does not depend on the operating system: sorting the paths themselves gives a different order on Windows (case-insensitive)
    files = sorted(list((HERE.parent / "amplitrim").glob("*.py")) + list(HERE.glob("*.py")) + list(HERE.glob("*.R")), key=lambda p: (p.parent.name, p.name))
    for p in files:
        h.update(p.name.encode())
        h.update(p.read_bytes().replace(b"\r\n", b"\n"))
    return h.hexdigest()[:10]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    from amplitrim import __version__
    ap.add_argument("--version", action="version", version=f"AmpliTrim {__version__} build {build_id()}")
    ap.add_argument("name")
    ap.add_argument("raw", type=Path, help="folder of paired FASTQs")
    ap.add_argument("--primer-f"), ap.add_argument("--primer-r")
    ap.add_argument("--error-rate", help="cutadapt error rate (0.15 allows 2 errors in a 15-18 nt primer)")
    ap.add_argument("--mixed-orientation", action="store_true")
    ap.add_argument("--max-ee", default="2", help="maxEE: one number, or F,R for different forward and reverse limits (default 2)")
    ap.add_argument("--explore-max-ee", action="store_true", help="also find the best truncation for several maxEE values and show what each keeps and costs")
    ap.add_argument("--explore-coverage", action="store_true", help="also show what covering fewer of the longest amplicons would allow")
    ap.add_argument("--allow-unmergeable", default=None, metavar="FRAC", help="share of the longest amplicons you accept not being able to merge (default 0.01)")
    ap.add_argument("--only", metavar="NAMES", help="analyse only these samples: comma-separated names, or a text file with one name per line")
    ap.add_argument("--group-regex", default=None, help="how to group sample names in the per-sample table: a regex whose first group is the group name")
    ap.add_argument("--min-retention", default=None, metavar="FRAC", help="with --explore-max-ee: the share of pairs the quality-first maxEE suggestion must keep (default 0.5)")
    ap.add_argument("--min-reads", type=int, default=0, metavar="N", help="leave out samples with fewer than N read pairs (controls, failed runs)")
    ap.add_argument("--marker", choices=["16s", "its"], help="its for variable-length amplicons (then leave out --compare)")
    ap.add_argument("--compare", metavar="F,R", help="also run DADA2 at this truncLenF,truncLenR (16S), e.g. 240,160")
    ap.add_argument("--pairs", type=int, default=5000, help="pairs per sample given to DADA2 when there are 6 or more samples (default 5000)")
    ap.add_argument("--runinfo", type=Path, help="SRA RunInfo CSV (sample names and archive read counts)")
    ap.add_argument("--reference", type=Path, help="mock community: FASTA of the expected sequences")
    ap.add_argument("--merge-check", action="store_true", help="also classify why DADA2 fails to merge pairs (3 samples; slower)")
    ap.add_argument("--out", type=Path, default=Path("results"))
    args = ap.parse_args()
    if not args.raw.is_dir():
        raise SystemExit(f"RAW_DIR '{args.raw}' is not a folder")
    try:
        sys.stdout.reconfigure(errors="replace")
    except AttributeError:
        pass
    args.out.mkdir(parents=True, exist_ok=True)
    log_path = args.out / f"{args.name}_outcome.txt"
    log = open(log_path, "w", encoding="utf-8")

    def emit(text: str = "") -> None:
        print(text)
        log.write(text + "\n")
        log.flush()

    def section(title: str) -> None:
        emit(f"\n==================== {title}")

    def run(cmd, tail: int = 0) -> int:
        emit("$ " + " ".join(str(c) for c in cmd))
        try:
            proc = subprocess.run([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")
        except FileNotFoundError as e:
            emit(f"could not start {cmd[0]}: {e}")
            return 127
        lines = proc.stdout.splitlines()
        for line in (lines[-tail:] if tail else lines):
            emit(line)
        return proc.returncode

    py = sys.executable
    primers = [p for p in (args.primer_f, args.primer_r) if p]
    primer_args = [x for p in primers for x in ("--primer", p)]
    from amplitrim.subsample import find_pairs
    n_samples = len(list(find_pairs(args.raw)))

    from amplitrim import __version__
    section(f"dataset: {args.name}  ({time.strftime('%Y-%m-%d %H:%M')})  |  AmpliTrim {__version__} build {build_id()}")
    emit(f"raw folder: {args.raw} | samples found: {n_samples} | primers '{args.primer_f or 'none'}' / '{args.primer_r or 'none'}', "
         f"error-rate '{args.error_rate or 'default'}', marker '{args.marker or '16s'}', compare '{args.compare or 'none'}'")

    section("1. environment")
    run([py, HERE / "check_environment.py"], tail=12)

    section("2. the raw reads (describe_reads)")
    run([py, HERE / "describe_reads.py", args.raw, "--lengths", "--positions", "15"] + primer_args)

    section("3. trim, AmpliTrim and DADA2 (real_data)")
    cmd = [py, HERE / "real_data.py", "--local", args.raw, "--name", args.name, "--out", args.out, "--dada2"]
    if args.primer_f:
        cmd += ["--primer-f", args.primer_f]
    if args.primer_r:
        cmd += ["--primer-r", args.primer_r]
    if args.error_rate:
        cmd += ["--error-rate", args.error_rate]
    if args.mixed_orientation:
        cmd += ["--mixed-orientation"]
    if args.min_reads:
        cmd += ["--min-reads", str(args.min_reads)]
    if args.only:
        cmd += ["--only", args.only]
    cmd += ["--max-ee", args.max_ee]
    if args.explore_max_ee:
        cmd += ["--explore-max-ee"]
    if args.min_retention is not None:
        cmd += ["--min-retention", args.min_retention]
    if args.explore_coverage:
        cmd += ["--explore-coverage"]
    if args.allow_unmergeable is not None:
        cmd += ["--allow-unmergeable", args.allow_unmergeable]
    if args.marker:
        cmd += ["--marker", args.marker]
    if args.compare:
        cmd += ["--compare", args.compare]
    if n_samples >= 6:
        cmd += ["--dada2-pairs-per-sample", str(args.pairs)]
    run(cmd)
    res = args.out / args.name
    summary_file = res / "report" / "summary.json"
    if not summary_file.exists():
        emit("No AmpliTrim result was written: see the error above.")
        emit(f"Send this file: {log_path}")
        log.close()
        return 1

    section("4. the reads after trimming (primers and reverse complements should now read about 0%)")
    run([py, HERE / "describe_reads.py", res / "ready", "--lengths"] + primer_args)

    section("5. AmpliTrim's summary.json")
    summary = json.loads(summary_file.read_text(encoding="utf-8"))
    emit(json.dumps(summary, indent=1))
    rec = summary.get("recommended")
    tf, tr = (rec["truncLenF"], rec["truncLenR"]) if rec else (None, None)

    if n_samples >= 6 and rec and args.marker != "its":
        section("6. per sample")
        cmd = [py, HERE / "per_sample.py", res, "--max-ee", args.max_ee]
        if args.runinfo:
            cmd += ["--runinfo", args.runinfo]
        if args.group_regex:
            cmd += ["--group-regex", args.group_regex]
        track = res / f"dada2_per_sample_{tf}x{tr}.tsv"
        if track.exists():
            cmd += ["--dada2-track", track]
        run(cmd)

    if args.merge_check and rec and args.marker != "its":
        section("7. why pairs fail to merge in DADA2 (3 samples)")
        reads = res / "dada2_reads" if (res / "dada2_reads").is_dir() else res / "ready"
        if find_rscript():
            run([find_rscript(), HERE / "dada2_merge_check.R", reads, tf, tr, 3, args.max_ee])
        else:
            emit("Rscript was not found: skipped.")

    if args.reference:
        section("8. mock community: ASVs against the expected sequences")
        for fasta in sorted(res.glob("asvs_*.fasta")):
            run([py, HERE / "mock_check.py", fasta, args.reference])

    section("done")
    emit(f"Send this ONE file back: {log_path}")
    log.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
