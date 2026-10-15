#!/usr/bin/env python
"""Put the data drive in order: number the runs and the datasets, and find the files that are stored twice. It plans first and changes nothing unless you say --apply.

  python validation/organise_data.py --table                       the numbered runs and datasets (the order of TEST_REGISTER.md)
  python validation/organise_data.py --base C:\\data                 the plan: which folders would get which number (nothing is changed)
  python validation/organise_data.py --base C:\\data --apply         number them, and rewrite the folder names in amplitrim_paths.json (a .bak copy is kept)
  python validation/organise_data.py --base C:\\data --undo C:\\data\\organise_data_<stamp>.csv     put everything back
  python validation/organise_data.py --base C:\\data --duplicates                   the files that are identical in several places, and the space they take (nothing is changed)
  python validation/organise_data.py --base C:\\data --duplicates --apply           set the extra copies aside in <base>\\_duplicates\\<stamp>\\ (never deleted)
  python validation/organise_data.py --base C:\\data --duplicates --restore C:\\data\\duplicates_moved_<stamp>.csv     put them back
  python validation/organise_data.py --base C:\\data --next                         the number the next run should get
  python validation/organise_data.py --base C:\\data --check-paths                 every path in amplitrim_paths.json, ok or MISSING (run it before and after --apply)

Runs: every folder of <base>\\runs gets the number of its run as a prefix (night3_2026-10-09 becomes 05_night3_2026-10-09). A folder whose name is not
recognised is left alone and listed; give it a number with  --number NAME=NN.  Datasets: every dataset folder of <base> (mock_V1-V3_mockB_PRJEB67827, ...)
gets its dataset number. Folders that already start with two digits and an underscore are skipped, so running it twice does nothing the second time.

It refuses to rename while a run seems to be going (a log or a state file under runs\\ was written in the last 30 minutes; --force overrides), and it
rolls back if a rename fails halfway. Names inside amplitrim_paths.json (--config, repeatable; by default the one next to <base> or in C:\\amplitrim_v2\\amplitrim)
are rewritten in the same step. Duplicates are found by size, then a hash of the start, then a hash of the whole file; hard links are not duplicates. A copy that amplitrim_paths.json
names (a dataset's ready folder, a reference) is never set aside, and --protect TEXT protects every copy whose path contains TEXT.
"""
import argparse
import csv
import datetime
import hashlib
import json
import os
import re
import shutil
import sys
import time
from collections import defaultdict
from pathlib import Path

# number, name, sections of TEST_REGISTER.md, scripts, patterns of the folder names under runs\ (any one matches; lower case; checked in the order of the numbers)
RUNS = [
    (1, "development rounds", "1-8", "filter_grid_check.py (test D1), real_data.py", [r"^(dev|development)", r"^round_?[123](?!\d)"]),
    (2, "round 4", "9-10", "round4.py, round4_steps.py", [r"^round_?4(?!\d)"]),
    (3, "night 1", "11", "night1_v2.py (not in this archive), describe_reads.py, compare_methods.py", [r"^night_?1(?!\d)"]),
    (4, "night 2", "12-13", "night2.py, maxee_arms.py", [r"^night_?2(?!\d)"]),
    (5, "night 3", "14-15", "night3.py", [r"^night_?3(?!\d)"]),
    (6, "night 4 (offline)", "16-17", "maxee_rules.py, maxee_frontier.py", [r"^night_?4_?(offline|stage|rules|frontier)"]),
    (7, "night 4b", "18", "night4.py", [r"^night_?4(?!\d)"]),
    (8, "ITS check", "19", "its_check.py, its_filter_check.R", [r"^its_check"]),
    (9, "ITS merge check", "20", "its_merge_check.py, its_merge_check.R", [r"^its_merge"]),
    (10, "release check", "21", "release_check.py", [r"^release_check"]),
    (11, "soil V4", "22", "overmerge_real.py, overnight.py", [r"soil", r"overmerge", r"v4_over"]),
    (12, "merge estimates against DADA2", "23", "estimate_check.py", [r"^estimate_check", r"merge_estimate"]),
    (13, "night 5", "24", "nocap_check.py", [r"^night_?5(?!\d)", r"^nocap_check"]),
    (14, "night 6", "25-26", "night6.py, longread_describe.py", [r"^night_?6(?!\d)"]),
]
# number, id of the results site, every pattern must match the folder name (lower case): the accession, and the region or the name where one accession holds two datasets
DATASETS = [
    (1, "mock-a", "Mock community A", [r"prjeb67797"]),
    (2, "mock-b", "Mock community B", [r"prjeb67827"]),
    (3, "interlab-v1v2", "Inter-laboratory mock, V1-V2", [r"prjeb34919", r"v1-v2"]),
    (4, "interlab-v4v6", "Inter-laboratory mock, V4-V6", [r"prjeb34919", r"v4-v6"]),
    (12, "swabs", "Meat-processing swabs, V3-V4", [r"prjna1261676", r"swab"]),
    (5, "d6300-v3v4", "ZymoBIOMICS D6300, V3-V4", [r"prjna1261676"]),
    (6, "d6300-v4v5", "ZymoBIOMICS D6300, V4-V5", [r"srr1696158[123]|d6300.*v4-v5|v4-v5.*d6300"]),
    (7, "v1v2", "Real samples, V1-V2", [r"prjeb42528"]),
    (8, "soil-v4", "Soil, V4", [r"prjna574030"]),
    (9, "its", "Spring wheat, ITS", [r"prjeb58189"]),
    (10, "wastewater", "Wastewater, V3-V4", [r"prjeb75155"]),
    (11, "spike-in", "Mouse gut with spike-ins, V3-V4", [r"prjna1243306"]),
    (13, "buetas-2024", "Periodontitis patients, V3-V4 and full-length (Buetas et al.)", [r"prjna933120"]),
]
NUMBERED = re.compile(r"^\d\d_")
ACTIVE_MINUTES = 30


def run_number(name: str):
    n = name.lower()
    for number, _, _, _, patterns in sorted(RUNS):
        if any(re.search(p, n) for p in patterns):
            return number
    return None


def dataset_number(name: str):
    n = name.lower()
    for number, _, _, patterns in DATASETS:
        if all(re.search(p, n) for p in patterns):
            return number
    return None


def numbered(number: int, name: str) -> str:
    return f"{number:02d}_{name}"


def next_run_number(runs_dir: Path) -> int:
    highest = max(r[0] for r in RUNS)
    if runs_dir.is_dir():
        for d in runs_dir.iterdir():
            if d.is_dir() and NUMBERED.match(d.name):
                highest = max(highest, int(d.name[:2]))
    return highest + 1


def table_lines() -> list:
    L = ["Runs (the order of TEST_REGISTER.md):", f"  {'no.':>3}  {'run':32s} {'register':9s} scripts"]
    L += [f"  {n:02d}   {name:32s} {sec:9s} {scripts}" for n, name, sec, scripts, _ in RUNS]
    L += ["", "Datasets (the order of the results site; 13 is the dataset of night 6):", f"  {'no.':>3}  {'id':16s} name"]
    L += [f"  {n:02d}   {i:16s} {name}" for n, i, name, _ in sorted(DATASETS)]
    return L


# ----------------------------------------------------------------------------- numbering

def recent_activity(runs_dir: Path, minutes: int = ACTIVE_MINUTES):
    """A log or state file under runs\\ written in the last `minutes`, or None."""
    if not runs_dir.is_dir():
        return None
    limit = time.time() - minutes * 60
    for d in runs_dir.iterdir():
        if d.is_dir():
            for pattern in ("*.log", "*_state.json"):
                for f in d.glob(pattern):
                    if f.stat().st_mtime > limit:
                        return f
    return None


def scan(base: Path, what: str = "all", manual: dict = None) -> dict:
    """{'plan': [(kind, number, old, new)], 'skipped': [(path, why)], 'left': [path]} for the folders under `base`."""
    manual = {k.lower(): v for k, v in (manual or {}).items()}
    plan, skipped, left = [], [], []
    runs_dir = base / "runs"
    if what in ("all", "runs") and runs_dir.is_dir():
        for d in sorted(runs_dir.iterdir()):
            if not d.is_dir():
                continue
            if NUMBERED.match(d.name):
                skipped.append((d, "already numbered"))
                continue
            n = manual.get(d.name.lower()) or run_number(d.name)
            if n is None:
                skipped.append((d, "not recognised: give it a number with --number " + d.name + "=NN"))
            else:
                plan.append(("run", n, d, d.with_name(numbered(n, d.name))))
    if what in ("all", "datasets") and base.is_dir():
        for d in sorted(base.iterdir()):
            if not d.is_dir() or d.name == "runs" or d.name.startswith(("_", ".")):
                continue
            if NUMBERED.match(d.name):
                skipped.append((d, "already numbered"))
                continue
            n = manual.get(d.name.lower()) or dataset_number(d.name)
            if n is not None:
                plan.append(("dataset", n, d, d.with_name(numbered(n, d.name))))
            elif re.match(r"^(mock|real)_", d.name, re.I):
                skipped.append((d, "looks like a dataset but no accession of the table matches: give it a number with --number " + d.name + "=NN"))
            else:
                left.append(d)
    return {"plan": plan, "skipped": skipped, "left": left}


def _norm(s: str) -> str:
    return s.replace("\\", "/").rstrip("/")


def rewrite_value(value: str, mapping: list):
    """`value` with the folder name of any (old, new) in `mapping` replaced when it is the start of the path (not just the start of a name); (text, changed)."""
    v = _norm(value)
    for old, new in mapping:
        o = _norm(old)
        if v.lower() == o.lower() or v.lower().startswith(o.lower() + "/"):
            sep = "\\" if "\\" in value else "/"
            rest = v[len(o):]
            return _norm(new).replace("/", sep) + rest.replace("/", sep), True
    return value, False


def rewrite_json(obj, mapping: list):
    """(new object, number of strings changed), walking every value of a parsed JSON document."""
    if isinstance(obj, str):
        new, changed = rewrite_value(obj, mapping)
        return new, int(changed)
    if isinstance(obj, list):
        out, n = [], 0
        for x in obj:
            y, k = rewrite_json(x, mapping)
            out.append(y)
            n += k
        return out, n
    if isinstance(obj, dict):
        out, n = {}, 0
        for key, x in obj.items():
            y, k = rewrite_json(x, mapping)
            out[key] = y
            n += k
        return out, n
    return obj, 0


def default_configs(base: Path) -> list:
    found = [base / "amplitrim_paths.json", base.parent / "amplitrim_paths.json", Path(r"C:\amplitrim_v2\amplitrim\amplitrim_paths.json")]
    seen, out = set(), []
    for p in found:
        if p.is_file() and str(p.resolve()).lower() not in seen:
            seen.add(str(p.resolve()).lower())
            out.append(p)
    return out


def update_configs(configs: list, mapping: list, stamp: str, write: bool) -> list:
    """[(config path, strings changed)]; with `write`, the file is rewritten after a .bak copy."""
    results = []
    for path in configs:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        new, n = rewrite_json(data, mapping)
        if n and write:
            shutil.copy2(path, str(path) + f".bak_{stamp}")
            Path(path).write_text(json.dumps(new, indent=2) + "\n", encoding="utf-8")
        results.append((Path(path), n))
    return results


def apply_plan(plan: list, base: Path, configs: list, stamp: str, log=print) -> Path:
    """Rename the folders of `plan`, then rewrite the configs; returns the mapping file. A failed rename undoes the ones already done."""
    for kind, n, old, new in plan:
        if new.exists():
            raise SystemExit(f"{new} already exists: nothing was renamed.")
    mapping = [(str(old), str(new)) for _, _, old, new in plan]
    update_configs(configs, mapping, stamp, write=False)                  # a config that does not parse stops us before anything is renamed
    done = []
    try:
        for kind, n, old, new in plan:
            old.rename(new)
            done.append((old, new))
            log(f"  {old.name}  ->  {new.name}")
    except OSError as err:
        failed = plan[len(done)][2]
        for was, now in reversed(done):
            try:
                now.rename(was)
            except OSError:
                pass
        raise SystemExit(f"Could not rename {failed.name} ({err}); the {len(done)} folder(s) already renamed were put back. A folder open in Explorer or a program that is using it is the usual cause.")
    changed = update_configs(configs, mapping, stamp, write=True)
    record = base / f"organise_data_{stamp}.csv"
    with open(record, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "number", "old", "new"])
        for kind, n, old, new in plan:
            w.writerow([kind, n, str(old), str(new)])
    for path, n in changed:
        log(f"  {path}: {n} path(s) rewritten" + (f" (copy kept as {path.name}.bak_{stamp})" if n else ""))
    return record


def undo(record: Path, configs: list, stamp: str, log=print) -> None:
    rows = list(csv.DictReader(open(record, newline="", encoding="utf-8")))
    for r in rows:
        if Path(r["old"]).exists():
            raise SystemExit(f"{r['old']} exists again: nothing was undone.")
        if not Path(r["new"]).exists():
            raise SystemExit(f"{r['new']} is missing: nothing was undone.")
    for r in reversed(rows):
        Path(r["new"]).rename(Path(r["old"]))
        log(f"  {Path(r['new']).name}  ->  {Path(r['old']).name}")
    mapping = [(r["new"], r["old"]) for r in rows]
    for path, n in update_configs(configs, mapping, stamp, write=True):
        log(f"  {path}: {n} path(s) rewritten back")


# ----------------------------------------------------------------------------- duplicates

def _hash(path: Path, limit: int = 0) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        remaining = limit or float("inf")
        while remaining > 0:
            chunk = fh.read(min(1 << 20, remaining) if limit else 1 << 20)
            if not chunk:
                break
            h.update(chunk)
            remaining -= len(chunk)
    return h.hexdigest()


def keeper_key(path: Path):
    return (len(path.parts), str(path).lower())


def _key(p) -> str:
    return str(p).replace("\\", "/").lower()


def is_protected(path, protect) -> bool:
    """True when the path contains any of the texts of `protect` (slash style and case do not matter)."""
    k = _key(path)
    return any(_key(t) in k for t in protect if t)


def config_paths(configs: list) -> list:
    """The strings of the config files that are paths that exist: the folders and files the pipeline is told to use, which must never be set aside."""
    found = []

    def walk(o):
        if isinstance(o, str):
            if ("/" in o or "\\" in o) and os.path.exists(o):
                found.append(o)
        elif isinstance(o, list):
            for x in o:
                walk(x)
        elif isinstance(o, dict):
            for x in o.values():
                walk(x)
    for c in configs:
        try:
            walk(json.loads(Path(c).read_text(encoding="utf-8")))
        except (OSError, ValueError):
            pass
    return found


def check_paths(configs: list) -> list:
    """[(config file, path, exists)] for every string of the config files that is a path (it starts with a drive letter or a slash)."""
    rows = []

    def walk(c, o):
        if isinstance(o, str):
            if re.match(r"^([A-Za-z]:[\\/]|/)", o):
                rows.append((Path(c).name, o, os.path.exists(o)))
        elif isinstance(o, list):
            for x in o:
                walk(c, x)
        elif isinstance(o, dict):
            for x in o.values():
                walk(c, x)
    for c in configs:
        try:
            walk(c, json.loads(Path(c).read_text(encoding="utf-8")))
        except (OSError, ValueError) as err:
            rows.append((Path(c).name, f"cannot be read: {err}", False))
    return rows


def find_duplicates(base: Path, min_size: int = 1_000_000, protect=(), log=print) -> list:
    """[(size, hash, [paths])] of files identical in content and stored in more than one place, the copy to keep first (the shortest path, or a protected one)."""
    by_size = defaultdict(list)
    for root, dirs, files in os.walk(base):
        dirs[:] = [d for d in dirs if d not in ("_duplicates", ".git")]
        for f in files:
            p = Path(root) / f
            try:
                st = p.stat()
            except OSError:
                continue
            if st.st_size >= min_size:
                by_size[st.st_size].append((p, st))
    groups, hashed = [], 0
    for size, items in sorted(by_size.items()):
        if len(items) < 2:
            continue
        files_of = defaultdict(list)
        for p, st in items:
            files_of[(st.st_dev, st.st_ino) if st.st_ino else (str(p),)].append(p)         # hard links of one file are one file
        if len(files_of) < 2:
            continue
        reps = [min(paths, key=keeper_key) for paths in files_of.values()]          # of the names of one file, the shortest path stands for it
        start = defaultdict(list)
        for p in reps:
            start[_hash(p, 1 << 20)].append(p)
        for bucket in start.values():
            if len(bucket) < 2:
                continue
            whole = defaultdict(list)
            for p in bucket:
                whole[_hash(p)].append(p)
                hashed += 1
            for h, paths in whole.items():
                if len(paths) >= 2:
                    keep = sorted(paths, key=lambda p: (not is_protected(p, protect), *keeper_key(p)))
                    groups.append((size, h, keep))
    log(f"{hashed} candidate files hashed")
    return sorted(groups, key=lambda g: -g[0] * (len(g[2]) - 1))


def human(n: float) -> str:
    return f"{n / 1e9:.2f} GB" if n >= 1e9 else f"{n / 1e6:.1f} MB"


def duplicate_lines(groups: list, protect=()) -> list:
    reclaim = sum(size * (len(paths) - 1) for size, _, paths in groups)
    L = [f"{len(groups)} group(s) of identical files; {human(reclaim)} would be freed by keeping one copy of each.", ""]
    for size, h, paths in groups:
        L.append(f"{size / 1e6:,.1f} MB x {len(paths)}  (sha256 {h[:12]})")
        L.append(f"   keep  {paths[0]}")
        L += [f"   extra {p}" for p in paths[1:]]
    return L


def set_aside(groups: list, base: Path, stamp: str, protect=(), log=print) -> Path:
    """Move every copy but the first of each group to <base>\\_duplicates\\<stamp>\\<its path under base>; returns the record to restore them from."""
    target = base / "_duplicates" / stamp
    record = base / f"duplicates_moved_{stamp}.csv"
    rows = []
    for size, h, paths in groups:
        for p in paths[1:]:
            if is_protected(p, protect):
                continue                                                  # protected: never moved
            try:
                rel = p.relative_to(base)
            except ValueError:
                rel = Path(p.name)
            dest = target / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(p), str(dest))
            rows.append((str(p), str(dest), size, h, str(paths[0])))
            log(f"  set aside {p}")
    with open(record, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["original", "set_aside", "bytes", "sha256", "kept_copy"])
        w.writerows(rows)
    return record


def restore(record: Path, log=print) -> int:
    n = 0
    for r in csv.DictReader(open(record, newline="", encoding="utf-8")):
        dest, orig = Path(r["set_aside"]), Path(r["original"])
        if not dest.exists():
            log(f"  missing, skipped: {dest}")
        elif orig.exists():
            log(f"  already there, skipped: {orig}")
        else:
            orig.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(dest), str(orig))
            n += 1
    return n


# ----------------------------------------------------------------------------- command line

def parse_numbers(items) -> dict:
    out = {}
    for item in items or []:
        name, _, number = item.rpartition("=")
        if not name or not number.isdigit() or not 0 < int(number) < 100:
            raise SystemExit(f"--number needs NAME=NN with NN from 1 to 99, got {item!r}")
        out[name] = int(number)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", type=Path, default=Path(r"C:\data"), help=r"the data folder (default C:\data)")
    ap.add_argument("--what", choices=["all", "runs", "datasets"], default="all", help="number the runs, the datasets, or both (default both)")
    ap.add_argument("--number", action="append", metavar="NAME=NN", help="give a folder under runs\\ (or a dataset folder) a number by hand; repeat for several")
    ap.add_argument("--config", action="append", type=Path, help="an amplitrim_paths.json whose folder names are rewritten (repeatable)")
    ap.add_argument("--apply", action="store_true", help="do it (without this nothing is changed)")
    ap.add_argument("--force", action="store_true", help="rename even if a run seems to be going")
    ap.add_argument("--undo", type=Path, metavar="CSV", help="put back the folders of an earlier --apply (its organise_data_<stamp>.csv)")
    ap.add_argument("--duplicates", action="store_true", help="find files that are identical in several places")
    ap.add_argument("--min-size", type=float, default=1.0, metavar="MB", help="ignore files smaller than this (default 1)")
    ap.add_argument("--protect", action="append", default=[], metavar="TEXT", help="a copy whose path contains TEXT is never moved (repeatable)")
    ap.add_argument("--restore", type=Path, metavar="CSV", help="with --duplicates: put back the files of an earlier --apply (its duplicates_moved_<stamp>.csv)")
    ap.add_argument("--next", action="store_true", help="print the number the next run should get")
    ap.add_argument("--table", action="store_true", help="print the numbered runs and datasets")
    ap.add_argument("--check-paths", action="store_true", help="list every path in the config files as ok or MISSING (exit code 1 if any is missing); changes nothing")
    a = ap.parse_args(argv)
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    if a.table:
        print("\n".join(table_lines()))
        return 0
    if a.next:
        print(f"{next_run_number(a.base / 'runs'):02d}")
        return 0
    if a.check_paths:
        configs_c = a.config or default_configs(a.base)
        if not configs_c:
            print("No amplitrim_paths.json found: give --config.")
            return 2
        rows = check_paths(configs_c)
        for name, path, ok in rows:
            print(f"  {'ok     ' if ok else 'MISSING'}  {path}   ({name})")
        missing = [r for r in rows if not r[2]]
        print(f"\n{len(rows) - len(missing)} of {len(rows)} path(s) exist" + (f"; {len(missing)} MISSING: compare with the same check made before the numbering." if missing else "."))
        return 1 if missing else 0
    if not a.base.is_dir():
        print(f"The folder {a.base} does not exist: give --base.")
        return 2
    if a.duplicates:
        if a.restore:
            print(f"{restore(a.restore)} file(s) put back")
            return 0
        configs_d = a.config or default_configs(a.base)
        protect = list(a.protect) + config_paths(configs_d)
        if configs_d:
            print(f"Protected, because {', '.join(c.name for c in configs_d)} names them: {len(protect) - len(a.protect)} path(s); these copies are never set aside.")
        groups = find_duplicates(a.base, int(a.min_size * 1_000_000), protect)
        print("\n".join(duplicate_lines(groups, protect)))
        if not groups:
            return 0
        if not a.apply:
            print("\nNothing was changed. Add --apply to set the extra copies aside (they are moved, not deleted).")
            return 0
        record = set_aside(groups, a.base, stamp, protect)
        print(f"\nThe extra copies are in {a.base / '_duplicates' / stamp}. When you are sure nothing needs them, delete that folder. To put them back: --duplicates --restore {record}")
        return 0
    configs = a.config or default_configs(a.base)
    if a.undo:
        undo(a.undo, configs, stamp)
        print("Done.")
        return 0
    s = scan(a.base, a.what, parse_numbers(a.number))
    for kind, n, old, new in s["plan"]:
        print(f"  {kind:8s} {n:02d}  {old.name}  ->  {new.name}")
    for p, why in s["skipped"]:
        if why != "already numbered":
            print(f"  skipped  {p.name}: {why}")
    if s["left"]:
        print("  left alone (not a run or a dataset of the tables): " + ", ".join(p.name for p in s["left"]))
    print(f"  config files that will be rewritten: {', '.join(str(c) for c in configs) or 'none found (give --config)'}")
    if not s["plan"]:
        print("Nothing to number.")
        return 0
    if not a.apply:
        print(f"\n{len(s['plan'])} folder(s) would be renamed. Nothing was changed. Add --apply to do it.")
        return 0
    active = recent_activity(a.base / "runs")
    if active and not a.force:
        print(f"\nA run seems to be going: {active} was written in the last {ACTIVE_MINUTES} minutes. Wait for it to finish (or give --force). Nothing was renamed.")
        return 1
    try:
        record = apply_plan(s["plan"], a.base, configs, stamp)
    except SystemExit as err:
        print(str(err))
        return 1
    print(f"\nDone. To put everything back: --undo {record}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
