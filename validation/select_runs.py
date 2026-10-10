#!/usr/bin/env python
"""See what a sequencing project contains, and pick a small, evenly spaced set of runs from one part of it.

A project can hold several amplicons (a V3-V4 library, an ITS library, ...) or several instruments. Give this script the project accession
(it asks ENA, which mirrors NCBI and DDBJ) or a table you downloaded (an ENA read_run report as TSV, or an SRA Run Selector RunInfo CSV).
It groups the runs by a label column, with numbers collapsed so that 'sample 1', 'sample 2' ... fall into one group, and shows for each group
the runs, read length and size. Then --pick chooses the group and --n how many runs to take, evenly spaced over the sorted accessions, and
--out writes the accession list that download_ena.py reads.

  python validation/select_runs.py PRJEB58189                          # what is in the project
  python validation/select_runs.py PRJEB58189 --pick "V3-V4" --n 24 --out list_v3v4.csv
  python validation/select_runs.py runs.tsv --by library_name          # a table you downloaded
"""

import argparse
import csv
import io
import re
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import download_ena as dl  # noqa: E402

API = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields={fields}&format=tsv"
FIELDS = ("run_accession,sample_alias,library_name,library_strategy,library_layout,instrument_model,read_count,base_count,fastq_bytes,"
          "experiment_title,study_title,sample_title")
FALLBACK = "run_accession,sample_alias,library_name,library_layout,read_count,base_count,fastq_bytes"
ALIASES = {"run": ("run_accession", "run"), "layout": ("library_layout", "librarylayout"), "instrument": ("instrument_model", "model"),
           "reads": ("read_count", "spots"), "bases": ("base_count", "bases"), "bytes": ("fastq_bytes",), "mb": ("size_mb",)}
LABELS = ("experiment_title", "library_name", "sample_title", "librarylayout", "sample_alias")


def _get(url: str, opener=None) -> str:
    """The text at `url`; a server error or lost connection is tried again after growing pauses, like the downloader does."""
    opener = opener or urllib.request.urlopen
    for attempt in range(1, len(dl.WAITS) + 2):
        try:
            with opener(url, timeout=60) as r:
                return r.read().decode()
        except OSError as err:
            code = getattr(err, "code", None)
            if (code is not None and code < 500 and code != 429) or attempt > len(dl.WAITS):
                raise
            wait = dl.WAITS[attempt - 1] * dl.PAUSE
            print(f"{dl.ena_failure('ENA', err)} Trying again in {wait:g} s (attempt {attempt + 1} of {len(dl.WAITS) + 1}).")
            time.sleep(wait)
    raise AssertionError("unreachable")


def fetch_table(accession: str, opener=None) -> str:
    """ENA's run table for a project, as TSV text. If ENA rejects a field, a smaller set is requested."""
    for fields in (FIELDS, FALLBACK):
        try:
            return _get(API.format(acc=accession, fields=fields), opener)
        except urllib.error.HTTPError as err:
            if err.code != 400 or fields == FALLBACK:
                raise SystemExit(dl.ena_failure(accession, err))
        except OSError as err:
            raise SystemExit(dl.ena_failure(accession, err))
    raise AssertionError("unreachable")


def read_table(text: str) -> list:
    """Rows of an ENA report (TSV) or SRA RunInfo (CSV), with lower-case column names."""
    text = text.lstrip("\ufeff")
    delimiter = "\t" if "\t" in text.splitlines()[0] else ","
    return [{k.strip().lower(): (v or "").strip() for k, v in r.items() if k} for r in csv.DictReader(io.StringIO(text), delimiter=delimiter)]


def field(row: dict, kind: str) -> str:
    return next((row[k] for k in ALIASES[kind] if row.get(k)), "")


def label_of(row: dict, by: str = None) -> str:
    for col in ([by.lower()] if by else LABELS):
        if row.get(col):
            return re.sub(r"\d+", "#", row[col])
    return "(no label)"


def to_float(x: str) -> float:
    try:
        return float(x)
    except ValueError:
        return 0.0


def summarise(rows: list, by: str = None) -> list:
    """One entry per group: label, runs, paired runs, instruments, reads (pairs), average read length where known, megabytes."""
    groups = defaultdict(list)
    for r in rows:
        groups[label_of(r, by)].append(r)
    out = []
    for label, rs in groups.items():
        paired = [r for r in rs if field(r, "layout").upper() == "PAIRED"]
        reads = sum(to_float(field(r, "reads")) for r in rs)
        bases = sum(to_float(field(r, "bases")) for r in rs)
        mb = sum(to_float(field(r, "mb")) or to_float(field(r, "bytes")) / 1e6 for r in rs)
        out.append({"label": label, "runs": len(rs), "paired": len(paired), "instruments": sorted({field(r, "instrument") for r in rs} - {""}),
                    "reads": reads, "avg_len": bases / reads if reads and bases else None, "mb": mb, "rows": rs})
    return sorted(out, key=lambda g: -g["runs"])


def pick_runs(rows: list, text: str, n: int, by: str = None) -> list:
    """`n` paired runs from the groups whose label contains `text`, evenly spaced over the sorted accessions (all of them if fewer)."""
    chosen = sorted(field(r, "run") for r in rows if text.lower() in label_of(r, by).lower() and field(r, "layout").upper() == "PAIRED")
    if n <= 0 or n >= len(chosen):
        return chosen
    step = (len(chosen) - 1) / (n - 1) if n > 1 else 0
    return [chosen[round(i * step)] for i in range(n)]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", help="a project accession (PRJEB..., PRJNA..., PRJDB...) or a table file")
    ap.add_argument("--by", metavar="COLUMN", help="the column that names a run's group (default: the first of experiment_title, library_name, sample_title)")
    ap.add_argument("--pick", metavar="TEXT", help="take runs whose group label contains this text (case does not matter)")
    ap.add_argument("--n", type=int, default=24, metavar="N", help="with --pick: how many runs, evenly spaced (default 24; 0 = all)")
    ap.add_argument("--out", type=Path, help="with --pick: write the accession list here")
    args = ap.parse_args()
    text = Path(args.source).read_text(encoding="utf-8") if Path(args.source).is_file() else fetch_table(args.source)
    rows = read_table(text)
    if not rows:
        raise SystemExit("The table has no runs.")
    groups = summarise(rows, args.by)
    print(f"{len(rows)} run(s) in {len(groups)} group(s):\n")
    print(f"{'group':<58s} {'runs':>5s} {'paired':>6s} {'read length':>11s} {'pairs':>11s} {'MB':>8s}  instrument")
    for g in groups:
        length = f"{g['avg_len']:.0f}" if g["avg_len"] else "?"
        print(f"{g['label'][:58]:<58s} {g['runs']:>5d} {g['paired']:>6d} {length:>11s} {g['reads'] / 2:>11,.0f} {g['mb']:>8,.0f}  {', '.join(g['instruments'])}")
    print("\n(read length is the average over both mates when the table gives bases; pairs are the table's reads divided by two)")
    if args.pick:
        runs = pick_runs(rows, args.pick, args.n, args.by)
        if not runs:
            raise SystemExit(f"No paired run has a group label containing {args.pick!r}. The groups are listed above.")
        size = sum(to_float(field(r, "mb")) or to_float(field(r, "bytes")) / 1e6 for r in rows if field(r, "run") in set(runs))
        print(f"\nPicked {len(runs)} paired run(s) matching {args.pick!r}, about {size:,.0f} MB: {', '.join(runs[:4])}{' ...' if len(runs) > 4 else ''}")
        if args.out:
            args.out.write_text("acc\n" + "\n".join(runs) + "\n", encoding="utf-8")
            print(f"Accession list written to {args.out}; download with: python validation\\download_ena.py {args.out} --out <folder>")
    return 0


if __name__ == "__main__":
    sys.exit(main())
