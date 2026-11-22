"""The names of the datasets and of the settings, for the pages of the website. One list (docs/DATASET_NAMES.tsv) gives every dataset its ID and its name; the archive's
scripts and results keep their own codes, which this list links to. The settings have readable names here (the codes are the register's)."""
import csv
import pathlib
import re

COLUMNS = ["id", "full_name", "short_name", "archive_code", "data_folder", "accession", "marker_region", "platform", "kind", "source"]


def load(root):
    path = pathlib.Path(root) / "docs" / "DATASET_NAMES.tsv"
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows and list(rows[0]) == COLUMNS, f"{path}: the columns are not {COLUMNS}"
    return rows


def slug(row):
    s = row["short_name"].lower().replace("\u00fc", "u")
    s = re.sub(r"[^a-z0-9]+", "-", s).strip("-")
    return f"{row['id'].lower()}-{s}"


def is_its(row):
    return "ITS" in row["marker_region"]


def by_code(rows):
    return {r["archive_code"]: r for r in rows}


SETTINGS = {
    "std": "DADA2 Default (maxEE 2,2, No Truncation)",
    "ee23": "maxEE 2,3", "ee24": "maxEE 2,4", "ee25": "maxEE 2,5", "ee33": "maxEE 3,3", "ee55": "maxEE 5,5",
    "tq3": "truncQ 3", "tq4": "truncQ 4", "tq5": "truncQ 5", "tq6": "truncQ 6", "tq7": "truncQ 7", "tq8": "truncQ 8", "tq10": "truncQ 10", "tq20": "truncQ 20",
    "sw10": "Trim: Window 4, Quality 10", "sw15": "Trim: Window 4, Quality 15", "sw20": "Trim: Window 4, Quality 20", "sw25": "Trim: Window 4, Quality 25",
    "q12": "Trim: Window 4, Quality 12", "q18": "Trim: Window 4, Quality 18",
    "fwd": "Forward Reads Only", "fwd15": "Forward Reads Only, Trimmed",
    "mfirst": "Merge First (vsearch), Then DADA2", "mo20": "mergePairs minOverlap 20", "mo30": "mergePairs minOverlap 30",
    "eer": "maxEE From The Reads", "sweer": "Trim, Then maxEE From The Trimmed Reads",
    "ml20": "Trim, minLen 20", "ml100": "Trim, minLen 100", "ml150": "Trim, minLen 150",
}
for w in (2, 3, 5, 6, 8, 10, 12, 16, 20):
    SETTINGS[f"w{w}"] = f"Trim: Window {w}, Quality 15"
ORDER = list(SETTINGS)            # the order the settings are listed in


def setting(code):
    return SETTINGS[code]


def setting_rank(code):
    return ORDER.index(code)
