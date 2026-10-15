"""Does AmpliTrim's estimate of what a maxEE filter keeps agree with what DADA2's own filter keeps, on long reads? (test L4)

For every single-end FASTQ file in a folder, on the same reads:
  1. the estimate: longread_describe.py's share of the reads with both primers that are at or under each maxEE, counted between the primers;
  2. the primers are removed with cutadapt (both primers, the orientation corrected, insertions and deletions allowed as in the estimate) and
     DADA2's filterAndTrim is run with ONLY maxEE acting (truncQ 0, minQ 0, maxN 0, no length limits): the share it keeps;
  3. the gate: the two agree within the tolerance (default 2 points) for every maxEE in every file.
It also records, as a description, how many reads survive four ways of removing the primers (the estimate's own search; cutadapt allowing
insertions and deletions; cutadapt with none and at most 2 mismatches; DADA2's removePrimers with its defaults), which shows how much of the
loss at the primer step comes from insertions and deletions. Needs: Python 3.9+, numpy, cutadapt (pip install cutadapt) and R with DADA2.

Usage:  python longread_filter_check.py FOLDER --out OUT --primer-f AGRGTTYGATYMTGGCTCAG --primer-r RGYTACCTTGTTACGACTT
        [--reads N] [--grid 1 2 3 5 10] [--tolerance 0.02] [--keep]
Writes OUT/filter_check.tsv, OUT/primer_arms.tsv and OUT/FILTER_CHECK.md.
"""

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from longread_describe import Locator, describe_file, find_files, open_text, revcomp, sample_name  # noqa: E402
from rtools import NOT_FOUND, find_rscript  # noqa: E402

DEFAULT_GRID = (1, 2, 3, 5, 10)
TOLERANCE = 0.02
NO_INDEL_MISMATCHES = 2                  # DADA2's removePrimers allows 2 mismatches by default
R_SCRIPT = HERE / "filter_check.R"


def count_reads(path: Path) -> int:
    with open_text(path) as fh:
        return sum(1 for _ in fh) // 4


def subset(path: Path, n: int, dest: Path) -> Path:
    """The first n reads of a file (all of them when n is 0), written to dest, so that every step sees the same reads."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path) as src, gzip.open(dest, "wt", newline="") as out:
        for i, line in enumerate(src):
            if n and i >= 4 * n:
                break
            out.write(line)
    return dest


def run_cutadapt(src: Path, dst: Path, primer_f: str, primer_r: str, errors, indels: bool) -> None:
    """Both primers as one linked adapter (the reverse primer as it appears on the read), reads in either orientation, reads without both dropped."""
    cmd = [sys.executable, "-m", "cutadapt", "-g", f"{primer_f}...{revcomp(primer_r)}", "--revcomp", "--discard-untrimmed",
           "-e", str(errors), "-j", "1", "--quiet", "-o", str(dst)]
    if not indels:
        cmd.append("--no-indels")
    cmd.append(str(src))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"cutadapt failed on {Path(src).name}: " + (res.stderr or res.stdout).strip()[-500:])


def run_r(args: list) -> None:
    rscript = find_rscript()
    if not rscript:
        raise SystemExit(NOT_FOUND)
    res = subprocess.run([rscript, str(R_SCRIPT), *[str(a) for a in args]], capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError("The DADA2 step failed: " + (res.stderr or res.stdout).strip()[-800:])


def read_table(path: Path) -> list:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def dada2_filter(reads: Path, grid, table: Path) -> dict:
    """{maxEE: (reads in, reads out)} of DADA2's filterAndTrim with only maxEE acting."""
    run_r(["filter", reads, ",".join(f"{g:g}" for g in grid), table])
    return {float(r["maxEE"]): (int(float(r["reads_in"])), int(float(r["reads_out"]))) for r in read_table(table)}


def dada2_remove_primers(reads: Path, primer_f: str, primer_r: str, table: Path) -> tuple:
    run_r(["primers", reads, primer_f, primer_r, table])
    row = read_table(table)[0]
    return int(float(row["reads_in"])), int(float(row["reads_out"]))


def compare(estimate: dict, kept: dict, grid, tolerance: float) -> list:
    """One row per maxEE: the estimate, the share DADA2 keeps, their difference and whether it is within the tolerance."""
    rows = []
    for g in grid:
        n_in, n_out = kept.get(float(g), (0, 0))
        share = n_out / n_in if n_in else None
        est = estimate.get(f"keep_maxEE_{g:g}")
        diff = None if share is None or est is None else share - est
        rows.append({"maxEE": g, "estimate": est, "dada2": share, "difference": diff, "reads_in": n_in, "reads_out": n_out,
                     "within": diff is not None and abs(diff) <= tolerance})
    return rows


def verdict(results: list, tolerance: float) -> tuple:
    """(met, the row with the largest difference): met only when every row of every file is within the tolerance."""
    rows = [dict(r, sample=name) for name, file_rows in results for r in file_rows]
    if not rows:
        return False, None
    measurable = [r for r in rows if r["difference"] is not None]
    worst = max(measurable, key=lambda r: abs(r["difference"])) if measurable else None
    return all(r["within"] for r in rows), worst


def check_file(path: Path, work: Path, locator: Locator, primer_f: str, primer_r: str, grid, max_reads: int, error_rate: float,
               tolerance: float, log=print) -> dict:
    name = sample_name(path)
    wd = work / name
    reads = subset(path, max_reads, wd / "reads.fastq.gz")
    estimate = describe_file(reads, locator, 0, tuple(grid))
    arms = {"reads": estimate["reads"], "describer": estimate["between_primers"]}
    trimmed = wd / "cutadapt_indels.fastq.gz"
    run_cutadapt(reads, trimmed, primer_f, primer_r, error_rate, True)
    arms["cutadapt_indels"] = count_reads(trimmed)
    other = wd / "cutadapt_no_indels.fastq.gz"
    run_cutadapt(reads, other, primer_f, primer_r, NO_INDEL_MISMATCHES, False)
    arms["cutadapt_no_indels"] = count_reads(other)
    arms["dada2_remove_primers"] = dada2_remove_primers(reads, primer_f, primer_r, wd / "primers.tsv")[1]
    if arms["cutadapt_indels"] == 0:
        log(f"{name}: no read has both primers, so there is nothing to filter")
        kept = {float(g): (0, 0) for g in grid}
    else:
        kept = dada2_filter(trimmed, grid, wd / "filter.tsv")
    rows = compare(estimate, kept, grid, tolerance)
    log(f"{name}: {arms['reads']:,} reads, both primers {arms['cutadapt_indels']:,} (cutadapt), largest difference "
        + ("none" if all(r["difference"] is None for r in rows) else f"{100 * max(abs(r['difference']) for r in rows if r['difference'] is not None):.2f} points"))
    return {"name": name, "arms": arms, "rows": rows}


def pct(x, digits: int = 1) -> str:
    return "-" if x is None else f"{100 * x:.{digits}f}"


def points(x) -> str:
    return "-" if x is None else f"{100 * x:+.2f}"


def write_outputs(results: list, out: Path, folder, primer_f: str, primer_r: str, grid, tolerance: float, max_reads: int) -> None:
    out.mkdir(parents=True, exist_ok=True)
    pairs = [(r["name"], r["rows"]) for r in results]
    met, worst = verdict(pairs, tolerance)
    with open(out / "filter_check.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "maxEE", "estimate", "dada2", "difference", "within_tolerance", "reads_in", "reads_out"])
        for name, rows in pairs:
            for r in rows:
                w.writerow([name, f"{r['maxEE']:g}", r["estimate"], r["dada2"], r["difference"], r["within"], r["reads_in"], r["reads_out"]])
    with open(out / "primer_arms.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "reads", "describer", "cutadapt_indels", "cutadapt_no_indels", "dada2_remove_primers"])
        for r in results:
            a = r["arms"]
            w.writerow([r["name"], a["reads"], a["describer"], a["cutadapt_indels"], a["cutadapt_no_indels"], a["dada2_remove_primers"]])
    total = sum(r["arms"]["reads"] for r in results)
    lines = ["# Filter agreement: AmpliTrim's estimate against DADA2's filterAndTrim", "",
             f"Folder {folder}; forward primer {primer_f}, reverse primer {primer_r}; {len(results)} file(s), {total:,} reads"
             + (f" (the first {max_reads:,} of each file)" if max_reads else "") + f"; tolerance {100 * tolerance:g} points.", ""]
    if worst is None:
        lines += ["**Gate (L4): NOT MET.** No file gave a comparison (no read had both primers).", ""]
    else:
        lines += [f"**Gate (L4): {'MET' if met else 'NOT MET'}.** The largest difference is {100 * abs(worst['difference']):.2f} points "
                  f"({worst['sample']}, maxEE {worst['maxEE']:g}); the gate is {100 * tolerance:g} points for every maxEE in every file.", ""]
    lines += ["The estimate is the share of the reads with both primers that are at or under maxEE, counted between the primers (longread_describe.py). "
              "DADA2's figure is the share of the reads left after cutadapt removed both primers (insertions and deletions allowed) that "
              "filterAndTrim keeps with only maxEE acting.", "",
              "| sample | maxEE | estimate % | DADA2 % | difference (points) | within |", "|---|---|---|---|---|---|"]
    for name, rows in pairs:
        for r in rows:
            lines.append(f"| {name} | {r['maxEE']:g} | {pct(r['estimate'])} | {pct(r['dada2'])} | {points(r['difference'])} | {'yes' if r['within'] else 'NO'} |")
    lines += ["", "## Reads that survive the primer step (description only, no gate)", "",
              "| sample | reads | estimate's search % | cutadapt, indels allowed % | cutadapt, no indels, 2 mismatches % | DADA2 removePrimers % |", "|---|---|---|---|---|---|"]
    for r in results:
        a = r["arms"]
        n = a["reads"]
        lines.append(f"| {r['name']} | {n:,} | {pct(a['describer'] / n if n else None)} | {pct(a['cutadapt_indels'] / n if n else None)} | "
                     f"{pct(a['cutadapt_no_indels'] / n if n else None)} | {pct(a['dada2_remove_primers'] / n if n else None)} |")
    lines += ["", "If the gate is not met, look first at the primer table: a gap between the estimate's search and cutadapt means the two do not find the same "
              "reads, and a gap in the filter alone means the primer boundaries differ. Nothing else is changed until the cause is known."]
    (out / "FILTER_CHECK.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def check_folder(folder: Path, out: Path, primer_f: str, primer_r: str, grid=DEFAULT_GRID, max_reads: int = 0, error_rate: float = 0.15,
                 tolerance: float = TOLERANCE, keep: bool = False, log=print) -> list:
    files = find_files(folder)
    if not files:
        raise SystemExit(f"No FASTQ files (.fastq, .fastq.gz, .fq, .fq.gz) in {folder}")
    locator = Locator(primer_f, primer_r, 60, error_rate)
    work = out / "work"
    results = []
    for p in files:
        results.append(check_file(p, work, locator, primer_f, primer_r, grid, max_reads, error_rate, tolerance, log))
    write_outputs(results, out, folder, primer_f, primer_r, grid, tolerance, max_reads)
    if not keep:
        shutil.rmtree(work, ignore_errors=True)
    met, worst = verdict([(r["name"], r["rows"]) for r in results], tolerance)
    log(f"Saved: {out / 'FILTER_CHECK.md'}")
    log("Gate L4: " + ("MET" if met else "NOT MET") + ("" if worst is None else f" (largest difference {100 * abs(worst['difference']):.2f} points)"))
    return results


def check_tools() -> None:
    res = subprocess.run([sys.executable, "-m", "cutadapt", "--version"], capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit("cutadapt is not installed in this Python: pip install cutadapt")
    if not find_rscript():
        raise SystemExit(NOT_FOUND)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("folder", type=Path, help="folder of single-end FASTQ files")
    ap.add_argument("--out", type=Path, default=Path("filter_check"))
    ap.add_argument("--primer-f", required=True, help="forward primer (IUPAC codes allowed), 5'->3'")
    ap.add_argument("--primer-r", required=True, help="reverse primer, 5'->3' as ordered")
    ap.add_argument("--reads", type=int, default=0, metavar="N", help="use only the first N reads of each file (default: all)")
    ap.add_argument("--grid", type=float, nargs="+", default=list(DEFAULT_GRID), help="maxEE values to compare (default 1 2 3 5 10)")
    ap.add_argument("--error-rate", type=float, default=0.15, help="share of a primer's bases that may be wrong when indels are allowed (default 0.15)")
    ap.add_argument("--tolerance", type=float, default=TOLERANCE, help="the gate: the largest allowed difference as a share (default 0.02 = 2 points)")
    ap.add_argument("--keep", action="store_true", help="keep the intermediate files")
    a = ap.parse_args(argv)
    check_tools()
    check_folder(a.folder, a.out, a.primer_f, a.primer_r, tuple(a.grid), a.reads, a.error_rate, a.tolerance, a.keep)
    return 0


if __name__ == "__main__":
    sys.exit(main())
