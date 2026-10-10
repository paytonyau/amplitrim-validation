#!/usr/bin/env python
"""Collect the results of the runs into the validation archive: the reports, protocols, logs and result tables, and no reads.

  python validation\\collect_results.py --runs C:\\data\\runs --out C:\\amplitrim_release\\amplitrim-validation\\results [--dry-run]

From every run folder in --runs (night1_..., night2_..., night3_..., night4_..., its_check, its_merge_check, release_check, ...) it copies the text outputs (.md, .tsv, .json,
.log, .txt) of at most 2 MB each, with their folder structure, and leaves out the reads and everything DADA2 and vsearch wrote to work folders ('reads', 'work', 'trimmed',
'merge', 'filt'), and the ASV FASTA files. In the copies, a Windows user folder (C:\\Users\\name) is replaced by C:\\Users\\<user>. It writes MANIFEST.tsv (path, size, sha256) to
--out, so that the record can be checked, and prints what it copied and what it left out. Nothing in --runs is changed.
"""
import argparse
import hashlib
import re
import sys
from pathlib import Path

KEEP_SUFFIXES = {".md", ".tsv", ".json", ".log", ".txt"}
SKIP_DIRS = {"reads", "work", "trimmed", "merge", "filt", "__pycache__"}
SKIP_PREFIXES = ("asvs_",)
MAX_BYTES = 2_000_000
USER_DIR = re.compile(r"([A-Za-z]:[\\/]+Users[\\/]+)[^\\/\s\"']+", re.IGNORECASE)


def redact(text: str) -> str:
    return USER_DIR.sub(lambda m: m.group(1) + "<user>", text)


def wanted(path: Path, root: Path):
    """(True, '') for a file to copy, else (False, the reason it is left out)."""
    rel = path.relative_to(root)
    if any(part in SKIP_DIRS for part in rel.parts[:-1]):
        return False, "work folder"
    if path.name.startswith(SKIP_PREFIXES):
        return False, "ASV FASTA"
    if path.suffix.lower() not in KEEP_SUFFIXES:
        return False, "not a text output"
    if path.stat().st_size > MAX_BYTES:
        return False, f"over {MAX_BYTES // 1_000_000} MB"
    return True, ""


def collect(runs: Path, out: Path, dry_run: bool = False) -> dict:
    runs, out = runs.resolve(), out.resolve()
    if out == runs or runs in out.parents or out in runs.parents:
        raise ValueError("--out must be a folder of its own, outside --runs (and not holding it)")
    copied, left = [], {}
    for run in sorted(p for p in runs.iterdir() if p.is_dir() and p.name not in SKIP_DIRS):
        for f in sorted(p for p in run.rglob("*") if p.is_file()):
            ok, why = wanted(f, runs)
            if not ok:
                left.setdefault(why, []).append(str(f.relative_to(runs)))
                continue
            dest = out / f.relative_to(runs)
            data = redact(f.read_text(encoding="utf-8", errors="replace")).encode("utf-8")
            if not dry_run:
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
            copied.append((str(f.relative_to(runs)).replace("\\", "/"), len(data), hashlib.sha256(data).hexdigest()))
    if copied and not dry_run:
        (out / "MANIFEST.tsv").write_text("path\tbytes\tsha256\n" + "\n".join(f"{p}\t{n}\t{h}" for p, n, h in copied) + "\n", encoding="utf-8", newline="\n")
    return {"copied": copied, "left_out": left}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=Path, required=True, help="the folder that holds the run folders (for example C:\\data\\runs)")
    ap.add_argument("--out", type=Path, required=True, help="the folder to write to (for example amplitrim-validation\\results)")
    ap.add_argument("--dry-run", action="store_true", help="show what would be copied, write nothing")
    a = ap.parse_args(argv)
    if not a.runs.is_dir():
        print(f"No such folder: {a.runs}")
        return 2
    try:
        r = collect(a.runs, a.out, a.dry_run)
    except ValueError as err:
        print(err)
        return 2
    total = sum(n for _, n, _ in r["copied"])
    print(f"{'Would copy' if a.dry_run else 'Copied'} {len(r['copied'])} files ({total / 1e6:.1f} MB) from {len({p.split('/')[0] for p, _, _ in r['copied']})} run folders to {a.out}")
    for why, files in sorted(r["left_out"].items()):
        print(f"  left out, {why}: {len(files)} files" + (f" (for example {files[0]})" if files else ""))
    if not r["copied"]:
        print("nothing to copy: is --runs the folder that holds night1_..., night3_... and so on?")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
