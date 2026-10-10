#!/usr/bin/env python3
"""
Copy the real MiSeq V4 example reads that ship with the DADA2 R package (2 samples of 1,500 read pairs, 250 nt,
primers already removed) into OUTDIR as sam1_R1.fastq.gz, sam1_R2.fastq.gz, ... so they can be used as a small
real-data test with no download. Needs R with DADA2.

    python validation/example_data.py OUTDIR
"""
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rtools import NOT_FOUND, find_rscript  # noqa: E402


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    rscript = find_rscript()
    if rscript is None:
        print(NOT_FOUND)
        return 1
    try:
        res = subprocess.run([rscript, "-e", "cat(system.file('extdata', package='dada2'))"],
                             capture_output=True, text=True)
    except FileNotFoundError:
        print("Rscript was not found. Install R (https://cran.r-project.org), put its bin folder on PATH, and install DADA2 (see the Windows steps in validation/README.md); python validation/check_environment.py shows what is missing.")
        return 1
    src = Path(res.stdout.strip())
    if res.returncode != 0 or not src.is_dir():
        print("Could not find the DADA2 package. Install R and DADA2 first (see check_environment.py).")
        return 1
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    for sample in ("sam1", "sam2"):
        for mate, tag in (("F", "R1"), ("R", "R2")):
            shutil.copy(src / f"{sample}{mate}.fastq.gz", out / f"{sample}_{tag}.fastq.gz")
    print(f"copied 4 files to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
