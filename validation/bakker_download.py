#!/usr/bin/env python
"""The fungal mock communities of Bakker (2018): download from ENA, by the SRA run table (no accession is written in this file).

Bakker MG (2018). A fungal mock community control for amplicon sequencing experiments. Molecular Ecology Resources 18(3):541-556. doi:10.1111/1755-0998.12760.
BioProject PRJNA377530; 19 fungal isolates in three mock communities (even, staggered A, staggered B), in triplicate, paired MiSeq reads.

  python validation\\bakker_download.py --base C:\\data --runinfo C:\\data\\bakker\\SraRunInfo.csv              (download, resumable)
  python validation\\bakker_download.py --base C:\\data --runinfo C:\\data\\bakker\\SraRunInfo.csv --check      (every file there, the pair counts of the run table)
  --dry-run lists the runs; --only pcr or --only polymerase; --with-negatives adds the negative controls.

Two sequencing runs (Table S2 of the paper; Appendix S1):
  pcr         (SRX2614..., library names FMock*): ITS1 only; the standard workflow, a high PCR cycle number (HC) and a long extension (LE)
              -> 22_mock_ITS1_bakker_pcr_PRJNA377530\\raw
  polymerase  (SRX3014..., library names Mock*_GoTaq / _Phusion / _PreAmp): ITS1 AND ITS2 amplicons in the same files; they are split by primer later
              (bakker_primers.py) -> 23_mock_ITSmix_bakker_polymerase_PRJNA377530\\raw
Files are named <library name>_R1.fastq.gz and _R2.fastq.gz. The selected rows of the run table are saved beside the data (bakker_runs.tsv) as the record.
The machinery (the resumable download, the gzip check, the pair count) is that of arfken_download.py; the ENA address is built here for any accession length.
"""
import argparse
import csv
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import arfken_download as ad  # noqa: E402

FOLDERS = {"pcr": "22_mock_ITS1_bakker_pcr_PRJNA377530", "polymerase": "23_mock_ITSmix_bakker_polymerase_PRJNA377530"}


def ena_url(run: str, mate: int) -> str:
    """ENA's FASTQ address. The sub-folder depends on the number of digits of the accession: 6 digits none; 7 digits '00' and the last digit; 8 digits '0' and the
    last two; 9 digits the last three. (arfken_download.ena_url knows only the 8-digit form of Arfken's runs; Bakker's runs have 7 digits.)"""
    digits = len(run) - 3
    sub = "" if digits <= 6 else f"00{run[-1]}/" if digits == 7 else f"0{run[-2:]}/" if digits == 8 else f"{run[-3:]}/"
    return f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{run[:6]}/{sub}{run}/{run}_{mate}.fastq.gz"


def batch_of(row: dict) -> str:
    exp = row.get("Experiment", "")
    if exp.startswith("SRX2614"):
        return "pcr"
    if exp.startswith("SRX3014"):
        return "polymerase"
    raise ValueError(f"run {row.get('Run')} is in neither sequencing run of Bakker (2018): experiment {exp}")


def read_runinfo(path: Path, with_negatives: bool = False) -> dict:
    """{'pcr': [(run, library, pairs, bases)], 'polymerase': [...]} from SraRunInfo.csv (as NCBI's 'Send to: File, RunInfo' writes it)."""
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8-sig")))
    if not rows or not {"Run", "LibraryName", "spots", "bases", "Experiment", "LibraryLayout", "BioProject"} <= set(rows[0]):
        raise SystemExit(f"{path} is not an SRA run table with the columns Run, LibraryName, spots, bases, Experiment, LibraryLayout and BioProject")
    out = {b: [] for b in FOLDERS}
    for r in rows:
        if r["BioProject"] != "PRJNA377530":
            raise SystemExit(f"run {r['Run']} belongs to {r['BioProject']}, not PRJNA377530")
        if r["LibraryLayout"] != "PAIRED":
            raise SystemExit(f"run {r['Run']} is not paired")
        if r["LibraryName"].lower().startswith("neg") and not with_negatives:
            continue
        out[batch_of(r)].append((r["Run"], r["LibraryName"], int(r["spots"]), int(r["bases"])))
    for b in out:
        out[b].sort(key=lambda x: x[1])
        names = [x[1] for x in out[b]]
        if len(names) != len(set(names)):
            raise SystemExit(f"two runs share a library name in the {b} run")
    return out


def save_table(base: Path, runs: dict) -> Path:
    p = Path(base) / "bakker_runs.tsv"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("run_set\tfolder\trun\tlibrary\tread_pairs\tbases\n" + "".join(f"{b}\t{FOLDERS[b]}\t{r}\t{l}\t{n}\t{x}\n" for b in runs for r, l, n, x in runs[b]), encoding="utf-8")
    return p


def main(argv=None, opener=None, sleep=time.sleep, verify=ad.verify_gzip) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--runinfo", type=Path, required=True, help="SraRunInfo.csv of PRJNA377530")
    ap.add_argument("--only", choices=sorted(FOLDERS))
    ap.add_argument("--with-negatives", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    runs = read_runinfo(a.runinfo, a.with_negatives)
    if a.only:
        runs = {a.only: runs[a.only]}
    for b, rs in runs.items():
        print(f"{b:10s} {len(rs):2d} runs, {sum(r[2] for r in rs):>9,} read pairs, mean read length {min(ad.mean_read_length(r[2], r[3]) for r in rs):.0f} to "
              f"{max(ad.mean_read_length(r[2], r[3]) for r in rs):.0f} nt -> {Path(a.base) / FOLDERS[b] / 'raw'}")
    if a.dry_run:
        return 0
    if a.check:
        bad = []
        for b, rs in runs.items():
            for run, lib, pairs, _ in rs:
                for mate in (1, 2):
                    p = Path(a.base) / FOLDERS[b] / "raw" / f"{lib}_R{mate}.fastq.gz"
                    n = ad.count_pairs(p) if p.is_file() else None
                    if n != pairs:
                        bad.append(f"{lib} mate {mate}: {'missing' if n is None else f'{n} reads'}, the run table says {pairs}")
        print("\n".join(bad) if bad else "every file is there and its read count matches the run table")
        return 2 if bad else 0
    saved = save_table(a.base, runs)
    done = repaired = kept = 0
    for b, rs in runs.items():
        for run, lib, pairs, _ in rs:
            for mate in (1, 2):
                dest = Path(a.base) / FOLDERS[b] / "raw" / f"{lib}_R{mate}.fastq.gz"
                state = ad.download(ena_url(run, mate), dest, opener=opener, sleep=sleep, verify=verify)
                print(f"{state:11s} {lib}_R{mate}")
                done += state == "downloaded"
                repaired += state == "repaired"
                kept += state not in ("downloaded", "repaired")
    print(f"\n{done} downloaded, {repaired} repaired, {kept} already there and complete. The run table used: {saved}")
    print(f"Next: python validation\\bakker_download.py --base {a.base} --runinfo {a.runinfo} --check")
    return 0


if __name__ == "__main__":
    sys.exit(main())
