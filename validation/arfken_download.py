#!/usr/bin/env python
"""Download the single-isolate ITS1 and ITS2 reads of Arfken et al. 2023 (BioProject PRJNA693350) and check them against the ENA run table.

  python validation\\arfken_download.py --dry-run                  (shows the plan, downloads nothing)
  python validation\\arfken_download.py --base C:\\data             (downloads, about 2 to 3 GB; a file that is already there is kept, so it can be run again)
  python validation\\arfken_download.py --base C:\\data --check     (counts the read pairs of every file against the run table)

Two folders are made under --base, one per marker, with the files named <library>_R1.fastq.gz and <library>_R2.fastq.gz (for example C-ITS1a_R1.fastq.gz), which is the form
validation\\its_followup.py and the R script read:
  16_mock_ITS1_arfken_isolates_PRJNA693350\\raw      and      17_mock_ITS2_arfken_isolates_PRJNA693350\\raw
The isolates are C (Candida albicans), K (Kazachstania slooffiae), M (the Mucor / Lichtheimia corymbifera isolate), P (Pichia fermentans) and T (Trichosporon asahii), each in
triplicate (a, b, c). Two more folders hold the mixes MixA (KCTPM), MixB (KC) and MixC (KT) in triplicate, which are the mixes of known composition whose ITS1 and ITS2 files differ in the run table (for MixD
to MixM and ATCC2 the two markers have identical read and base counts, which looks like duplicated deposits, so they are not used):
  18_mock_ITS1_arfken_mixes_PRJNA693350\\raw      and      19_mock_ITS2_arfken_mixes_PRJNA693350\\raw       (--only mixes)
The ATCC mycobiome standard MSA-1010 (ten strains, Table 3 of the paper) in triplicate, libraries ATCC-ITS1a to c and ATCC-ITS2a to c (not ATCC2, whose ITS1 and ITS2 runs
have identical counts), goes to two more folders (--only atcc):
  20_mock_ITS1_arfken_atcc_PRJNA693350\\raw      and      21_mock_ITS2_arfken_atcc_PRJNA693350\\raw
Their accessions are not written here: first
  python validation\\arfken_download.py --base C:\\data --find-atcc
asks ENA for the run table of PRJNA693350, takes the libraries named exactly ATCC-ITS1a to c and ATCC-ITS2a to c, refuses unless there are three of each and the ITS1
and ITS2 counts differ, and saves them to <base>\\arfken_atcc_runs.tsv (read it); --only atcc then downloads from that file, and --check counts against it.
The other mixes, the negative controls and the separate K-Feb sample are not downloaded.
"""
import argparse
import csv
import gzip
import io
import re
import shutil
import sys
import time
import urllib.request
from pathlib import Path

DATASETS = {"its1": "16_mock_ITS1_arfken_isolates_PRJNA693350", "its2": "17_mock_ITS2_arfken_isolates_PRJNA693350",
            "mix_its1": "18_mock_ITS1_arfken_mixes_PRJNA693350", "mix_its2": "19_mock_ITS2_arfken_mixes_PRJNA693350",
            "atcc_its1": "20_mock_ITS1_arfken_atcc_PRJNA693350", "atcc_its2": "21_mock_ITS2_arfken_atcc_PRJNA693350"}
GROUPS = {"isolates": ["its1", "its2"], "mixes": ["mix_its1", "mix_its2"], "atcc": ["atcc_its1", "atcc_its2"]}
DEFAULT_KEYS = ["its1", "its2", "mix_its1", "mix_its2"]                  # without --only; the ATCC standard is asked for by name
ENA_REPORT = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA693350&result=read_run"
              "&fields=run_accession,library_name,sample_alias,sample_title,read_count,base_count&format=tsv&limit=0")
ATCC_NAME = re.compile(r"^ATCC[-_ ]?ITS([12])([abc])$", re.I)             # ATCC-ITS1a ... ATCC-ITS2c; ATCC2-ITS1a does not match
ATCC_TABLE = "arfken_atcc_runs.tsv"
# (run accession, library, read pairs, total bases), copied from the ENA run table of PRJNA693350
RUNS = {
    "its1": [("SRR13489861", "C-ITS1c", 366019, 159304823), ("SRR13489862", "C-ITS1b", 328661, 142933045), ("SRR13489866", "C-ITS1a", 351425, 152993050),
             ("SRR13489860", "K-ITS1a", 249182, 138887880), ("SRR13489920", "K-ITS1b", 223001, 124272130), ("SRR13489919", "K-ITS1c", 302891, 168855125),
             ("SRR13489918", "M-ITS1a", 119058, 66104057), ("SRR13489917", "M-ITS1b", 141893, 78705368), ("SRR13489916", "M-ITS1c", 87197, 48370552),
             ("SRR13489915", "P-ITS1a", 307764, 100101951), ("SRR13489914", "P-ITS1b", 491298, 159713252), ("SRR13489913", "P-ITS1c", 320053, 104064356),
             ("SRR13489912", "T-ITS1a", 218061, 88134781), ("SRR13489911", "T-ITS1b", 245483, 98956355), ("SRR13489909", "T-ITS1c", 277311, 111744456)],
    "its2": [("SRR13496185", "C-ITS2a", 121574, 68059694), ("SRR13496181", "C-ITS2b", 113071, 63293136), ("SRR13496180", "C-ITS2c", 147723, 82689973),
             ("SRR13496179", "K-ITS2a", 111535, 62287457), ("SRR13496239", "K-ITS2b", 164916, 92097911), ("SRR13496238", "K-ITS2c", 65371, 36508176),
             ("SRR13496237", "M-ITS2a", 54302, 30420765), ("SRR13496236", "M-ITS2b", 73561, 41212295), ("SRR13496235", "M-ITS2c", 88060, 49338236),
             ("SRR13496234", "P-ITS2a", 246146, 132925017), ("SRR13496233", "P-ITS2b", 103445, 55851629), ("SRR13496232", "P-ITS2c", 158518, 85596216),
             ("SRR13496231", "T-ITS2a", 171953, 95988570), ("SRR13496230", "T-ITS2b", 200431, 111885449), ("SRR13496228", "T-ITS2c", 151206, 84408708)],
    # the mixes with a known composition and with different ITS1 and ITS2 files: MixA (KCTPM), MixB (KC), MixC (KT)
    "mix_its1": [("SRR13489908", "MixA-ITS1a", 156787, 66385305), ("SRR13489907", "MixA-ITS1b", 328607, 139163008), ("SRR13489906", "MixA-ITS1c", 238299, 101032661),
                 ("SRR13489905", "MixB-ITS1a", 151534, 69193359), ("SRR13489904", "MixB-ITS1b", 154362, 70335520), ("SRR13489903", "MixB-ITS1c", 169808, 77432938),
                 ("SRR13489902", "MixC-ITS1a", 280033, 121203945), ("SRR13489901", "MixC-ITS1b", 263801, 113570459), ("SRR13489900", "MixC-ITS1c", 262200, 113066767)],
    "mix_its2": [("SRR13496227", "MixA-ITS2a", 174743, 96691911), ("SRR13496226", "MixA-ITS2b", 225831, 124953972), ("SRR13496225", "MixA-ITS2c", 147104, 81397637),
                 ("SRR13496224", "MixB-ITS2a", 164203, 91844473), ("SRR13496223", "MixB-ITS2b", 175466, 98140964), ("SRR13496222", "MixB-ITS2c", 117401, 65662421),
                 ("SRR13496221", "MixC-ITS2a", 140463, 78423364), ("SRR13496220", "MixC-ITS2b", 158848, 88684696), ("SRR13496219", "MixC-ITS2c", 191514, 106927032)],
}


def ena_url(run: str, mate: int) -> str:
    """ENA's folder rule for an accession of 11 characters: the first six characters, then '0' and the last two digits."""
    return f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{run[:6]}/0{run[-2:]}/{run}/{run}_{mate}.fastq.gz"


def target(base: Path, key: str, library: str, mate: int) -> Path:
    return Path(base) / DATASETS[key] / "raw" / f"{library}_R{mate}.fastq.gz"


def find_atcc(opener=None) -> dict:
    """{'atcc_its1': [(run, library, pairs, bases), ...], 'atcc_its2': [...]} from ENA's run table of PRJNA693350. A library counts when its library name, sample alias or
    sample title is exactly ATCC-ITS1a to c or ATCC-ITS2a to c. Refuses unless there are three runs per marker, one per replicate, and the two markers differ in counts."""
    opener = opener or (lambda u: urllib.request.urlopen(u, timeout=60))
    with opener(ENA_REPORT) as r:
        text = r.read().decode("utf-8")
    found = {"atcc_its1": {}, "atcc_its2": {}}
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        for field in ("library_name", "sample_alias", "sample_title"):
            m = ATCC_NAME.match((row.get(field) or "").strip())
            if m:
                marker, rep_ = m.groups()
                found[f"atcc_its{marker}"][rep_.lower()] = (row["run_accession"], f"ATCC-ITS{marker}{rep_.lower()}", int(row["read_count"]), int(row["base_count"]))
                break
    for key, reps in found.items():
        if sorted(reps) != ["a", "b", "c"]:
            raise RuntimeError(f"ENA's run table has {len(reps)} ATCC libraries for {key} ({', '.join(sorted(reps)) or 'none'}), not a, b and c")
    one, two = ({v[2:] for v in found[k].values()} for k in ("atcc_its1", "atcc_its2"))
    if one == two:
        raise RuntimeError("the ATCC ITS1 and ITS2 runs have identical read and base counts, like the duplicated deposits: not used")
    return {k: [v[r] for r in sorted(v)] for k, v in found.items()}


def save_atcc(base: Path, runs: dict) -> Path:
    p = Path(base) / ATCC_TABLE
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = ["folder\trun\tlibrary\tread_pairs\tbases"] + [f"{k}\t{r}\t{lib}\t{n}\t{b}" for k in sorted(runs) for r, lib, n, b in runs[k]]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def runs_for(key: str, base: Path) -> list:
    """The runs of a folder: written in this file for the isolates and the mixes; for the ATCC standard, read from the table that --find-atcc saved."""
    if key in RUNS:
        return RUNS[key]
    p = Path(base) / ATCC_TABLE
    if not p.is_file():
        raise SystemExit(f"{p} is missing: run  python validation\\arfken_download.py --base {base} --find-atcc  first, and read the table it saves")
    rows = [line.split("\t") for line in p.read_text(encoding="utf-8").splitlines()[1:] if line.strip()]
    return [(r[1], r[2], int(r[3]), int(r[4])) for r in rows if r[0] == key]


def mean_read_length(pairs: int, bases: int) -> float:
    return bases / pairs / 2


def verify_gzip(path: Path) -> bool:
    """True when the whole gzip stream can be read to its end marker (a file cut short by a dropped connection fails)."""
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 22):
                pass
        return True
    except (OSError, EOFError):
        return False


def download(url: str, dest: Path, opener=None, retries: int = 3, sleep=time.sleep, verify=None) -> str:
    """Returns 'kept' when the file is already there and complete, 'repaired' when a file that failed the check was replaced, and 'downloaded' otherwise.
    The file is written under a temporary name; it is renamed only when its size is the one the server promised (when it says) and, if verify is given, when verify(file) is true.
    A file that is not complete is downloaded again, up to retries times."""
    dest = Path(dest)
    status = "downloaded"
    if dest.is_file() and dest.stat().st_size > 0:
        if verify is None or verify(dest):
            return "kept"
        dest.unlink()
        status = "repaired"
    opener = opener or (lambda u: urllib.request.urlopen(u, timeout=60))
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    last = None
    for attempt in range(1, retries + 1):
        try:
            with opener(url) as r, open(part, "wb") as fh:
                shutil.copyfileobj(r, fh, 1 << 20)
                headers = getattr(r, "headers", None)
            promised = headers.get("Content-Length") if headers is not None else None
            got = part.stat().st_size
            if promised is not None and int(promised) != got:
                raise OSError(f"incomplete: {got:,} of {int(promised):,} bytes")
            if verify is not None and not verify(part):
                raise OSError("the file is not a complete gzip file")
            part.replace(dest)
            return status
        except Exception as e:                                           # a dropped connection is the usual cause; try again
            last = e
            part.unlink(missing_ok=True)
            if attempt < retries:
                sleep(5 * attempt)
    raise RuntimeError(f"could not download {url}: {last}")


def count_pairs(path: Path) -> int:
    with gzip.open(path, "rt") as fh:
        return sum(1 for _ in fh) // 4


def check(base: Path, key: str, runs=None) -> list:
    """Problems found: a missing file, or a read count that differs from the run table (for either mate)."""
    problems = []
    for run, library, pairs, _ in (runs if runs is not None else runs_for(key, base)):
        for mate in (1, 2):
            p = target(base, key, library, mate)
            if not p.is_file():
                problems.append(f"{library} mate {mate}: file missing")
                continue
            try:
                n = count_pairs(p)
            except (OSError, EOFError) as e:
                problems.append(f"{library} mate {mate}: cannot be read ({e})")
                continue
            if n != pairs:
                problems.append(f"{library} mate {mate}: {n:,} reads, the run table says {pairs:,}")
    return problems


def main(argv=None, opener=None, sleep=time.sleep, verify=verify_gzip) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", type=Path, default=Path(r"C:\data"))
    ap.add_argument("--only", choices=sorted(DATASETS) + sorted(GROUPS), help="one folder (its1, its2, mix_its1, mix_its2, atcc_its1, atcc_its2) or a group (isolates, mixes, atcc); "
                                                                             "the default is the isolates and the mixes")
    ap.add_argument("--find-atcc", action="store_true", help="ask ENA for the ATCC standard's runs and save them to <base>\\arfken_atcc_runs.tsv; download nothing")
    ap.add_argument("--dry-run", action="store_true", help="show the plan and download nothing")
    ap.add_argument("--check", action="store_true", help="count the read pairs of the files already downloaded")
    a = ap.parse_args(argv)
    if a.find_atcc:
        try:
            found = find_atcc(opener)
        except RuntimeError as e:
            print(e)
            return 2
        for k in sorted(found):
            for r, lib, n, b in found[k]:
                print(f"{k:10s} {r:12s} {lib:12s} {n:>9,} read pairs  {b:>12,} bases  (mean read length {mean_read_length(n, b):.0f} nt)")
        print(f"\nsaved: {save_atcc(a.base, found)}. Read it; then  python validation\\arfken_download.py --base {a.base} --only atcc")
        return 0
    keys = GROUPS[a.only] if a.only in GROUPS else [a.only] if a.only else DEFAULT_KEYS
    for key in keys:
        runs = runs_for(key, a.base)
        lengths = [mean_read_length(n, b) for _, _, n, b in runs]
        print(f"{key.upper()}: {len(runs)} runs, {sum(r[2] for r in runs):,} read pairs, mean read length {min(lengths):.0f} to {max(lengths):.0f} nt (run table: bases / pairs / 2) -> {a.base / DATASETS[key] / 'raw'}")
    if a.dry_run:
        return 0
    if a.check:
        bad = [p for key in keys for p in check(a.base, key)]
        print("\n".join(bad) if bad else "every file is there and its read count matches the run table")
        return 2 if bad else 0
    done = kept = repaired = 0
    for key in keys:
        for run, library, _, _ in runs_for(key, a.base):
            for mate in (1, 2):
                try:
                    status = download(ena_url(run, mate), target(a.base, key, library, mate), opener=opener, sleep=sleep, verify=verify)
                except RuntimeError as e:
                    print(e)
                    return 2
                done += status == "downloaded"
                kept += status == "kept"
                repaired += status == "repaired"
                print(f"{status:11s} {library}_R{mate}")
    print(f"\n{done} downloaded, {repaired} repaired, {kept} already there and complete. Next: python validation\\arfken_download.py --base {a.base} --check")
    return 0


if __name__ == "__main__":
    sys.exit(main())
