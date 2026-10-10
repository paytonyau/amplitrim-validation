#!/usr/bin/env python
"""Release check: does the released package recommend the settings that nights 3 and 4b ran through DADA2?

DADA2's result depends only on the settings it is given (truncLenF, truncLenR, maxEE, minOverlap, minLen). If the released package recommends the same settings on every
dataset, the DADA2 results of those nights are the released package's results, and nothing needs to be re-run. This runs the released package (no DADA2, no R) on the 12
datasets of the battery and compares what it recommends with what was tested: the arm H of night 4b where the rule relaxed maxEE, otherwise the arm A of night 3 (the
default maxEE, the best truncation at it; for ITS, no truncation).

  python validation\\release_check.py --package C:\\amplitrim_release\\amplitrim --night3 C:\\data\\runs\\night3_2026-10-09 --night4 C:\\data\\runs\\night4_2026-10-10
  python validation\\release_check.py --package ... --only mockB v4        some datasets

About 2 to 5 minutes per dataset (vsearch on the same subsample as the nights). The criterion (P30) is in TEST_REGISTER.md, section 21. Writes RELEASE_CHECK_REPORT.md and
release_check.tsv to --out.
"""
import argparse
import csv
import json
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT)]

DATASETS = ["mockA_high", "v1v2", "its", "v4", "mcm_v12", "mockB", "mcm_v46", "d6300_V3-V4", "d6300_V4-V5", "wastewater_V3-V4", "spikein_V3-V4", "swabs_V3-V4"]


def parse_ee(text) -> tuple:
    """'2' -> (2.0, 2.0); '2,3' -> (2.0, 3.0); a number or a pair of numbers is accepted as well."""
    if isinstance(text, (int, float)):
        return float(text), float(text)
    if isinstance(text, (tuple, list)):
        return float(text[0]), float(text[1])
    parts = [float(x) for x in str(text).replace(" ", "").split(",") if x != ""]
    return (parts[0], parts[0]) if len(parts) == 1 else (parts[0], parts[1])


def ee_text(ee) -> str:
    return f"{ee[0]:g}" if ee[0] == ee[1] else f"{ee[0]:g},{ee[1]:g}"


def _arm(folder: Path, tag: str, name: str, arm: str):
    path = folder / name / f"{tag}_result.json"
    if not path.is_file():
        return None
    r = json.loads(path.read_text(encoding="utf-8"))
    return (r.get("arms") or {}).get(arm)


def tested_settings(night3, night4, name: str):
    """What was run through DADA2 for this dataset: (label, truncF, truncR, maxEE pair), or None. Arm H of night 4b where it exists, otherwise arm A of night 3."""
    for folder, tag, arm, label in ((night4, "night4", "H", "night 4b, arm H"), (night3, "night3", "A", "night 3, arm A")):
        if folder is None:
            continue
        a = _arm(Path(folder), tag, name, arm)
        if a and a.get("truncF") is not None and a.get("maxEE") is not None:
            return label, int(float(a["truncF"])), int(float(a["truncR"])), parse_ee(a["maxEE"])
    return None


def compare(name: str, released: tuple, tested) -> dict:
    """released = (truncF, truncR, maxEE pair) or None when the package recommended nothing."""
    row = {"dataset": name, "released": released, "tested": tested}
    if tested is None:
        row["same"] = None
    elif released is None:
        row["same"] = False
    else:
        row["same"] = (released[0], released[1]) == (tested[1], tested[2]) and released[2] == tested[3]
    return row


def verdict(rows: list):
    decided = [r for r in rows if r["same"] is not None]
    if len(decided) < len(DATASETS):
        return None
    return all(r["same"] for r in decided)


def settings_text(x) -> str:
    return "-" if x is None else f"{x[0]}/{x[1]}, maxEE {ee_text(x[2])}"


def report_lines(rows: list) -> list:
    v = verdict(rows)
    same, diff = [r for r in rows if r["same"]], [r for r in rows if r["same"] is False]
    L = ["# Release check", "", "The criterion (P30) is in TEST_REGISTER.md, section 21, and was fixed before this ran.", "", "## Verdict", "",
         f"- **P30** {'MET' if v else 'NOT MET' if v is not None else 'not decided (not all 12 datasets were run, or a tested result was not found)'}: the released package recommends exactly the settings "
         f"that were run through DADA2 on every dataset: {len(same)} of {len(rows)} the same" + (f", {len(diff)} different" if diff else "") + ".", "",
         "| dataset | tested | what was tested | released package | the same |", "|---|---|---|---|---|"]
    for r in rows:
        t = r["tested"]
        L.append(f"| {r['dataset']} | {t[0] if t else 'not found'} | {settings_text(None if not t else (t[1], t[2], t[3]))} | {settings_text(r['released'])} | "
                 f"{'-' if r['same'] is None else 'yes' if r['same'] else '**NO**'} |")
    if diff:
        L += ["", "Where the settings differ, the DADA2 results of the nights are not claimed for the released package: DADA2 has to be run at the released settings."]
    elif v:
        L += ["", "On every dataset the released package recommends the settings that DADA2 ran at, so the DADA2 results of nights 3 and 4b are the released package's results."]
    return L + ["", "Next: send this file."]


def main(argv=None, pipeline_run=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--package", type=Path, help="the folder of the clean repository (the package under test)")
    ap.add_argument("--config", type=Path, help="amplitrim_paths.json (default: the one next to the package)")
    ap.add_argument("--night3", type=Path, required=False, help="the folder of night 3")
    ap.add_argument("--night4", type=Path, required=False, help="the folder of night 4b")
    ap.add_argument("--only", nargs="+", help="some datasets only (the verdict is then not decided)")
    ap.add_argument("--out", type=Path, default=Path("release_check_out"))
    a = ap.parse_args(argv)
    from its_check import resolve_package
    from maxee_rules import reads_folder
    pl = None
    if a.package and pipeline_run is None:
        root, tried = resolve_package(a.package)
        if root is None:
            print(f"The folder given with --package ({a.package}) has no 'amplitrim' package in it. Looked for: " + "; ".join(str(t) for t in tried))
            return 2
        sys.path.insert(0, str(root.resolve()))
        for m in [m for m in sys.modules if m == "amplitrim" or m.startswith("amplitrim.")]:
            del sys.modules[m]
    import importlib
    pl = importlib.import_module("amplitrim.pipeline")
    if pipeline_run is None:
        print(f"Using the amplitrim package at {Path(pl.__file__).parent}")
        if "rule" not in getattr(pl.Result, "__dataclass_fields__", {}):
            print("That package has no percentile rule (Result.rule): it is the development build. Give --package the folder of the clean repository.")
            return 2
    cfg_path = a.config or ROOT / "amplitrim_paths.json"
    cfg = json.loads(cfg_path.read_text(encoding="utf-8")) if cfg_path.exists() else {}
    if not (a.night3 or a.night4):
        print("Give --night3 and --night4: the folders where those nights wrote their results.")
        return 2
    rows = []
    for name in DATASETS:
        if a.only and name not in a.only:
            continue
        folder = reads_folder(cfg, name)
        if folder is None or not folder.is_dir():
            print(f"{name}: no reads folder in amplitrim_paths.json, skipped")
            continue
        with tempfile.TemporaryDirectory(prefix="release_check_") as tmp:
            res = (pipeline_run or pl.run)(pl.Settings(input_dir=folder, work_dir=Path(tmp) / "work", marker="its" if name == "its" else "16s"))
        r = res.recommendation
        released = None if r is None else (int(r.trim_f), int(r.trim_r), parse_ee(res.cfg.max_ee))
        rows.append(compare(name, released, tested_settings(a.night3, a.night4, name)))
        print(f"{name}: released {settings_text(released)} | tested {settings_text(None if not rows[-1]['tested'] else rows[-1]['tested'][1:])} | same: {rows[-1]['same']}", flush=True)
    if not rows:
        print("nothing was run: check amplitrim_paths.json")
        return 2
    a.out.mkdir(parents=True, exist_ok=True)
    L = report_lines(rows)
    (a.out / "RELEASE_CHECK_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    with open(a.out / "release_check.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dataset", "tested", "tested_truncF", "tested_truncR", "tested_maxEE", "released_truncF", "released_truncR", "released_maxEE", "same"])
        for r in rows:
            t, x = r["tested"], r["released"]
            w.writerow([r["dataset"], t[0] if t else "", *(["", "", ""] if not t else [t[1], t[2], ee_text(t[3])]), *(["", "", ""] if not x else [x[0], x[1], ee_text(x[2])]), r["same"]])
    print("\n".join(L[:8]))
    print(f"\nwritten: {a.out / 'RELEASE_CHECK_REPORT.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
