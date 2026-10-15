#!/usr/bin/env python
"""ITS merge check: where does DADA2 lose the pairs after the filter, exactly, and does an overhang explain part of the merging loss?

  python validation\\its_merge_check.py                                   the ITS dataset of amplitrim_paths.json, 5,000 pairs per sample as in night 3
  python validation\\its_merge_check.py --night3 C:\\data\\runs\\night3_2026-10-09   night 3's folder (its DADA2 input and track are reused when found)
  python validation\\its_merge_check.py --out C:\\data\\runs\\its_merge_check

Needs R and DADA2; it runs DADA2 on the ITS reads (learnErrors and dada: a few minutes). One run, then mergePairs twice on the same denoised reads, with the
default settings and with trimOverhang = TRUE; every rejected pair is counted by reason. Part of TEST_REGISTER.md, section 20 (criteria P26 to P29, fixed before the run).
Writes ITS_MERGE_REPORT.md and its_merge.tsv to --out.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT)]

KEY = "its"
MIN_LEN, MAX_EE, MIN_OVERLAP = 50, 2.0, 12
R_SCRIPT = HERE / "its_merge_check.R"
NIGHT3_POINTS = 0.3                   # P26: the default merge reproduces night 3's merged share within this many points
OVERHANG_POINTS = 2.0                 # P28: trimOverhang = TRUE merges at least this many points more of the input pairs
ESTIMATE_GAP_POINTS = 5.0             # P29: the exact merging loss exceeds AmpliTrim's estimate by more than this
ESTIMATE_POINTS = 1.89                # AmpliTrim's 'pass the filter but do not merge' on this dataset (ITS_CHECK_REPORT.md, register section 19)
COLUMNS = ["sample", "input", "filtered", "reach", "merged", "rej_short", "rej_mismatch", "rej_other", "merged_trim_overhang"]


def read_table(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as fh:
        rows = [{k: (float(v) if k != "sample" else v) for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]
    missing = [c for c in COLUMNS if rows and c not in rows[0]]
    if missing:
        raise RuntimeError(f"the table from R lacks the columns {missing}")
    return rows


def flow(rows: list) -> dict:
    """The exact flow of the input pairs, as shares of the input: through the filter, into the merge step, merged; and what the rejected pairs were."""
    inp = sum(r["input"] for r in rows)
    if inp <= 0:
        return {}
    filt, reach, merged = (sum(r[k] for r in rows) for k in ("filtered", "reach", "merged"))
    rej = {k: sum(r[k] for r in rows) / inp for k in ("rej_short", "rej_mismatch", "rej_other")}
    return {"samples": len(rows), "input": inp, "filtered": filt / inp, "reach": reach / inp, "merged": merged / inp, "lost_denoising": (filt - reach) / inp,
            "lost_merging": (reach - merged) / inp, **rej, "merged_trim_overhang": sum(r["merged_trim_overhang"] for r in rows) / inp,
            "overhang_gain": (sum(r["merged_trim_overhang"] for r in rows) - merged) / inp}


def verdicts(f: dict, night3_merged) -> list:
    if not f:
        return [("P26", "not decided: no result", None), ("P27", "not decided: no result", None), ("P28", "not decided: no result", None), ("P29", "not decided: no result", None)]
    gap = None if night3_merged is None else (f["merged"] - night3_merged) * 100
    return [("P26", f"The default merge reproduces night 3's merged share within {NIGHT3_POINTS} points: {f['merged']:.2%} against "
                    + ("night 3 not found" if gap is None else f"{night3_merged:.2%} ({gap:+.2f} points)"), None if gap is None else abs(gap) <= NIGHT3_POINTS),
            ("P27", f"Denoising accounts for more of the loss between the filter and the merged pairs than merging does, exactly: {f['lost_denoising']:.1%} against {f['lost_merging']:.1%} of the input pairs",
             f["lost_denoising"] > f["lost_merging"]),
            ("P28", f"mergePairs with trimOverhang = TRUE merges at least {OVERHANG_POINTS:g} points more of the input pairs than the default: {f['overhang_gain'] * 100:+.2f} points",
             f["overhang_gain"] * 100 >= OVERHANG_POINTS),
            ("P29", f"The exact merging loss exceeds AmpliTrim's estimate ({ESTIMATE_POINTS:g} points) by more than {ESTIMATE_GAP_POINTS:g} points: {f['lost_merging'] * 100:.1f} points, "
                    f"{f['lost_merging'] * 100 - ESTIMATE_POINTS:+.1f} over", f["lost_merging"] * 100 - ESTIMATE_POINTS > ESTIMATE_GAP_POINTS)]


def report_lines(f: dict, night3_merged, night3) -> list:
    L = ["# ITS merge check", "", "The criteria (P26 to P29) are in TEST_REGISTER.md, section 20, and were fixed before this ran.", "", "## Verdicts", ""]
    for p, txt, ok in verdicts(f, night3_merged):
        L.append(f"- **{p}** {'MET' if ok else 'NOT MET' if ok is not None else 'not decided'}: {txt}")
    if not f:
        return L + ["", "No pairs were given to DADA2."]
    L += ["", f"## Where the input pairs go ({f['samples']} samples, {f['input']:,.0f} input pairs; the exact split, every pair that enters the merge step is counted)", "",
          "| step | share of the input pairs |", "|---|---|", f"| pass the filter | {f['filtered']:.2%} |",
          f"| lost in denoising (reads not in a denoised forward and reverse sequence) | {f['lost_denoising']:.2%} |", f"| enter the merge step | {f['reach']:.2%} |",
          f"| rejected: overlap shorter than {MIN_OVERLAP} nt | {f['rej_short']:.2%} |", f"| rejected: mismatches or indels in the overlap | {f['rej_mismatch']:.2%} |",
          f"| rejected: a clean overlap that is still refused (an overhang) | {f['rej_other']:.2%} |", f"| **merged (default)** | **{f['merged']:.2%}** |",
          f"| merged with `trimOverhang = TRUE` | {f['merged_trim_overhang']:.2%} ({f['overhang_gain'] * 100:+.2f} points) |", "",
          f"- Night 3 merged {night3_merged:.2%} on the same reads." if night3_merged is not None else "- Night 3's track was not found, so P26 is not decided.",
          f"- Between the filter and the merged pairs: {f['lost_denoising'] * 100:.1f} points to denoising, {f['lost_merging'] * 100:.1f} points to merging (rejected pairs).",
          f"- AmpliTrim's estimate of 'pass the filter but do not merge' was {ESTIMATE_POINTS:g} points (vsearch, raw reads).", ""]
    return L + ["Next: send this file."]


def main(argv=None, run=subprocess.run, find_rscript=None, thin=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, help="amplitrim_paths.json (default: the one next to the package)")
    ap.add_argument("--night3", type=Path, help="the folder of night 3 (default: the newest night3_* or NN_night3_* in C:\\data\\runs)")
    ap.add_argument("--pairs", type=int, help="pairs per sample given to DADA2 (default: as in night 3)")
    ap.add_argument("--out", type=Path, default=Path("its_merge_check_out"))
    a = ap.parse_args(argv)
    import rtools
    from night2 import DEFAULT_PAIRS, PAIRS
    from its_check import find_night3_track, newest_night3, post_filter_loss, read_track
    cfg_path = a.config or ROOT / "amplitrim_paths.json"
    cfg = json.loads(cfg_path.read_text(encoding="utf-8")) if cfg_path.exists() else {}
    folder = Path(os.path.expandvars(cfg[KEY])) if cfg.get(KEY) else None
    if folder is None or not folder.is_dir():
        print(f"No reads folder for '{KEY}' in {cfg_path}.")
        return 2
    rscript = (find_rscript or rtools.find_rscript)()
    if not rscript:
        print(rtools.NOT_FOUND)
        return 2
    a.out.mkdir(parents=True, exist_ok=True)
    night3 = a.night3 or newest_night3(Path(os.environ.get("AMPLITRIM_RUNS", r"C:\data\runs")))
    track = find_night3_track(night3) if night3 else None
    night3_merged = post_filter_loss(read_track(track)).get("merged") if track else None
    if thin is None:
        from real_data import thin_reads as thin
    given = a.pairs or PAIRS.get(KEY, DEFAULT_PAIRS)
    reads = thin(folder, a.out / "reads", given) if given else folder
    table = a.out / "its_merge.tsv"
    p = run([rscript, str(R_SCRIPT), str(reads), str(MIN_LEN), f"{MAX_EE:g}", str(MIN_OVERLAP), str(table)], capture_output=True, text=True)
    if p.returncode != 0 or not table.is_file():
        print(f"DADA2 could not run: {(p.stderr or '').strip()[-500:]}")
        return 2
    f = flow(read_table(table))
    L = report_lines(f, night3_merged, night3)
    (a.out / "ITS_MERGE_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L[:11]))
    print(f"\nwritten: {a.out / 'ITS_MERGE_REPORT.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
