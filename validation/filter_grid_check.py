#!/usr/bin/env python
"""Development test D1: AmpliTrim's pass table against DADA2's own filterAndTrim on the DADA2 example reads, over a grid of truncation lengths and maxEE.

  python validation/filter_grid_check.py [--quick] [--out DIR]

DADA2 filters each of the two example samples at every combination of truncLenF, truncLenR and maxEE (maxN 0, truncQ 2, minLen 50) and the number of pairs it keeps
is compared with the one AmpliTrim's quality table gives. They must be identical: the filter replica is exact, not approximate. Needs R with DADA2.
"""

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

FULL = {"tf": [120, 150, 180, 210, 240], "tr": [100, 130, 160, 190, 220], "ee": [(1, 1), (2, 2), (2, 5), (5, 5)]}
QUICK = {"tf": [150, 240], "tr": [130, 220], "ee": [(2, 2), (5, 5)]}
R_CODE = """
suppressMessages(library(dada2))
d <- system.file("extdata", package = "dada2")
sams <- list(s1 = c("sam1F.fastq.gz", "sam1R.fastq.gz"), s2 = c("sam2F.fastq.gz", "sam2R.fastq.gz"))
res <- NULL
for (s in names(sams)) for (tf in c(%(tf)s)) for (tr in c(%(tr)s)) for (ee in list(%(ee)s)) {
  o <- suppressWarnings(filterAndTrim(file.path(d, sams[[s]][1]), tempfile(fileext = ".fastq.gz"), file.path(d, sams[[s]][2]), tempfile(fileext = ".fastq.gz"),
                                      truncLen = c(tf, tr), maxN = 0, maxEE = ee, truncQ = 2, minLen = 50, rm.phix = FALSE, compress = TRUE))
  res <- rbind(res, data.frame(sample = s, tf = tf, tr = tr, ee_f = ee[1], ee_r = ee[2], n_in = o[1, 1], n_out = o[1, 2]))
}
write.csv(res, "%(csv)s", row.names = FALSE)
cat(system.file("extdata", package = "dada2"))
"""


def compare(rows: list, pairs_of) -> list:
    """[(sample, tf, tr, ee_f, ee_r, n_in, n_dada2, n_amplitrim)] from the DADA2 counts and a function sample -> ReadPairs."""
    out, tables = [], {}
    for r in rows:
        ee = (float(r["ee_f"]), float(r["ee_r"]))
        key = (r["sample"], ee)
        if key not in tables:
            tables[key] = pairs_of(r["sample"]).pass_table(ee)
        n_in = int(r["n_in"])
        out.append((r["sample"], int(r["tf"]), int(r["tr"]), ee[0], ee[1], n_in, int(r["n_out"]), round(tables[key].rate(int(r["tf"]), int(r["tr"])) * n_in)))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="a small grid (8 settings per sample) instead of the full 100")
    ap.add_argument("--out", type=Path, default=Path("filter_grid_check"))
    args = ap.parse_args(argv)
    from rtools import find_rscript
    rscript = find_rscript()
    if rscript is None:
        raise SystemExit("Rscript was not found (R with DADA2 is needed).")
    grid = QUICK if args.quick else FULL
    args.out.mkdir(parents=True, exist_ok=True)
    csv_path = args.out / "dada2_counts.csv"
    code = R_CODE % {"tf": ", ".join(map(str, grid["tf"])), "tr": ", ".join(map(str, grid["tr"])), "ee": ", ".join(f"c({a}, {b})" for a, b in grid["ee"]), "csv": csv_path.as_posix()}
    with tempfile.TemporaryDirectory() as tmp:
        script = Path(tmp) / "grid.R"
        script.write_text(code, encoding="utf-8")
        res = subprocess.run([str(rscript), str(script)], capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit(f"R failed: {(res.stderr or res.stdout).strip()[-400:]}")
    example = Path(res.stdout.strip().splitlines()[-1])
    from amplitrim.trimmer import ReadPairs
    names = {"s1": ("sam1F.fastq.gz", "sam1R.fastq.gz"), "s2": ("sam2F.fastq.gz", "sam2R.fastq.gz")}
    cache = {}
    def pairs_of(sample):
        if sample not in cache:
            cache[sample] = ReadPairs.load(example / names[sample][0], example / names[sample][1])
        return cache[sample]
    rows = list(csv.DictReader(open(csv_path)))
    out = compare(rows, pairs_of)
    worst = max(abs(r[6] - r[7]) for r in out)
    rates = [r[6] / r[5] for r in out]
    with open(args.out / "grid_compare.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "truncF", "truncR", "maxEE_F", "maxEE_R", "pairs_in", "kept_dada2", "kept_amplitrim"])
        w.writerows(out)
    print(f"{len(out)} comparisons ({len(names)} samples x {len(grid['tf']) * len(grid['tr'])} truncation pairs x {len(grid['ee'])} maxEE settings); pass rates {min(rates):.1%} to {max(rates):.1%}")
    print(f"largest difference in the number of pairs kept: {worst}  ->  {'IDENTICAL' if worst == 0 else 'DIFFERENT'}")
    return 0 if worst == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
