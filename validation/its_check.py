#!/usr/bin/env python
"""ITS check: are the numbers of AmpliTrim's ITS report (where the pairs go, and what the forward reads alone would keep) what DADA2 itself gives?

  python validation\\its_check.py --package C:\\path\\to\\amplitrim     the folder of the clean repository: the package with the ITS report additions
  python validation\\its_check.py --night3 C:\\data\\runs\\night3_2026-10-09    where the night 3 DADA2 track tables are (default: the newest night3_* or NN_night3_* in C:\\data\\runs)
  python validation\\its_check.py --out C:\\data\\runs\\its_check

Part A (needs R and DADA2, about a minute): AmpliTrim's own ITS run (the code behind `amplitrim --its`, with vsearch) on the ITS reads, and DADA2's filterAndTrim with no
truncLen on the very same pooled reads, as pairs and with the forward reads alone. No denoising is run.
Part B (reads files that night 3 wrote; optional): DADA2's per-sample track for AmpliTrim's setting (maxEE 2) on the ITS dataset, to split DADA2's loss between the filter
and the merged pairs into reads that denoising does not keep and pairs that do not merge, and to put AmpliTrim's 'pass the filter but do not merge' estimate beside it.
This repository's own package is the build that nights 1 to 3 tested and has no 'where the pairs go': --package gives the package that has it (the clean repository).
The criteria (P23 to P25) are in TEST_REGISTER.md, section 19, and were fixed before this was run. It writes ITS_CHECK_REPORT.md and its_check.tsv to --out.
"""
import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT)]

KEY = "its"                                                       # the dataset's key in amplitrim_paths.json
MIN_LEN, MAX_EE = 50, 2.0                                         # DADA2's ITS workflow, and AmpliTrim's ITS defaults
POINTS = 0.2                                                      # P23 and P24: DADA2 and AmpliTrim agree within this many points of the pairs
ARM_A = "AmpliTrim (maxEE 2)"                                     # the arm of night 3 whose track is read in part B
R_SCRIPT = HERE / "its_filter_check.R"


def dada2_filter_counts(rscript: str, r1: Path, r2: Path, min_len: int = MIN_LEN, max_ee: float = MAX_EE, run=subprocess.run) -> dict:
    """DADA2's filterAndTrim on the pooled reads, as a pair and with the forward reads alone."""
    p = run([rscript, str(R_SCRIPT), str(r1), str(r2), str(min_len), f"{max_ee:g}"], capture_output=True, text=True)
    lines = [ln for ln in (p.stdout or "").splitlines() if ln.count("\t") == 3]
    if p.returncode != 0 or not lines:
        raise RuntimeError(f"Rscript ended with code {p.returncode} and gave no result: {(p.stderr or '').strip()[-400:]}")
    a, b, c, d = (int(x) for x in lines[-1].split("\t"))
    return {"pairs_in": a, "pairs_out": b, "forward_in": c, "forward_out": d}


def compare_filters(its_loss: dict, counts: dict) -> dict:
    """AmpliTrim's shares against DADA2's, on the same reads."""
    dada2_pass = counts["pairs_out"] / counts["pairs_in"] if counts["pairs_in"] else float("nan")
    dada2_fwd = counts["forward_out"] / counts["forward_in"] if counts["forward_in"] else float("nan")
    return {"same_reads": its_loss["pairs"] == counts["pairs_in"] == counts["forward_in"], "pairs_amplitrim": its_loss["pairs"], "pairs_dada2": counts["pairs_in"],
            "pass_amplitrim": its_loss["pass"], "pass_dada2": dada2_pass, "pass_gap": (its_loss["pass"] - dada2_pass) * 100,
            "forward_amplitrim": its_loss["forward_only_pass"], "forward_dada2": dada2_fwd, "forward_gap": (its_loss["forward_only_pass"] - dada2_fwd) * 100}


def read_track(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as fh:
        return [{k: (float(v) if k != "sample" else v) for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]


def post_filter_loss(rows: list) -> dict:
    """DADA2's loss between the filter and the merged pairs, in shares of the input pairs. A pair can only merge if both of its reads were denoised, so
    min(denoised forward, denoised reverse) bounds the pairs that reach merging: filtered - min is at least what denoising removed, and min - merged is at most what
    merging removed."""
    inp = sum(r["input"] for r in rows)
    if inp <= 0:
        return {}
    filtered = sum(r["filtered"] for r in rows)
    reach = sum(min(r["denoised_f"], r["denoised_r"]) for r in rows)
    merged = sum(r["merged"] for r in rows)
    return {"samples": len(rows), "input": inp, "filtered": filtered / inp, "denoised_f": sum(r["denoised_f"] for r in rows) / inp, "denoised_r": sum(r["denoised_r"] for r in rows) / inp,
            "merged": merged / inp, "nonchim": sum(r["nonchim"] for r in rows) / inp,
            "lost_denoising": max(0.0, filtered - reach) / inp, "lost_merging": max(0.0, reach - merged) / inp}


def resolve_package(given: Path):
    """The folder that contains the `amplitrim` package, from what the user gave: that folder itself, the package folder inside it, or the folder above a package
    folder (a repository unzipped one level deeper). Returns (root, None) or (None, the places that were looked at)."""
    given = Path(given)
    tried = [given / "amplitrim" / "__init__.py", given / "__init__.py", given / "amplitrim" / "amplitrim" / "__init__.py"]
    if (given / "amplitrim" / "__init__.py").is_file():
        return given, None
    if (given / "__init__.py").is_file() and given.name == "amplitrim":                 # the package folder itself was given
        return given.parent, None
    if (given / "amplitrim" / "amplitrim" / "__init__.py").is_file():                   # the folder above the repository
        return given / "amplitrim", None
    return None, tried


def find_night3_track(night3: Path, name: str = ARM_A):
    """The per-sample track of an arm of night 3 on the ITS dataset, or None."""
    from compare_methods import method_key
    folder = night3 / KEY / "compare"
    track = folder / f"track_{method_key(name)}.tsv"
    return track if track.is_file() else None


def newest_night3(runs: Path):
    """The newest night 3 folder of the runs folder: night3_<date>, or 05_night3_<date> once organise_data.py has numbered it. The date decides, not the number."""
    found = [d for d in runs.iterdir() if d.is_dir() and re.fullmatch(r"(\d\d_)?night3_.+", d.name)] if runs.is_dir() else []
    found.sort(key=lambda d: re.sub(r"^\d\d_", "", d.name))
    return found[-1] if found else None


def verdicts(cmp_: dict, loss: dict) -> list:
    ok23 = cmp_["same_reads"] and abs(cmp_["pass_gap"]) <= POINTS
    ok24 = cmp_["same_reads"] and abs(cmp_["forward_gap"]) <= POINTS
    out = [("P23", f"AmpliTrim's share of pairs passing the filter is within {POINTS} points of DADA2's, on the same reads: {cmp_['pass_amplitrim']:.2%} against {cmp_['pass_dada2']:.2%} ({cmp_['pass_gap']:+.2f} points)"
                   + ("" if cmp_["same_reads"] else "; THE READS DIFFER, so this is not decided"), ok23 if cmp_["same_reads"] else None),
           ("P24", f"AmpliTrim's share of forward reads passing alone is within {POINTS} points of DADA2's single-end filter, on the same reads: {cmp_['forward_amplitrim']:.2%} against "
                   f"{cmp_['forward_dada2']:.2%} ({cmp_['forward_gap']:+.2f} points)" + ("" if cmp_["same_reads"] else "; THE READS DIFFER, so this is not decided"), ok24 if cmp_["same_reads"] else None)]
    if loss:
        out.append(("P25", f"(exploratory) Of DADA2's loss between the filter and the merged pairs, denoising accounts for more than merging: {loss['lost_denoising']:.1%} of the input pairs "
                           f"lost to denoising (at least) against {loss['lost_merging']:.1%} lost to merging (at most)", loss["lost_denoising"] > loss["lost_merging"]))
    else:
        out.append(("P25", "(exploratory) not decided: the night 3 track table was not found", None))
    return out


def report_lines(cmp_: dict, its_loss: dict, loss: dict, night3) -> list:
    L = ["# ITS check", "", "The criteria are in TEST_REGISTER.md, section 19, and were fixed before this ran.", "", "## Verdicts", ""]
    for p, txt, ok in verdicts(cmp_, loss):
        L.append(f"- **{p}** {'MET' if ok else 'NOT MET' if ok is not None else 'not decided'}: {txt}")
    L += ["", "## A. The filter, on the same reads (AmpliTrim's own ITS run, and DADA2's filterAndTrim with no truncLen)", "",
          f"Pairs analysed: AmpliTrim {cmp_['pairs_amplitrim']:,}, DADA2 {cmp_['pairs_dada2']:,}" + ("" if cmp_["same_reads"] else " (**different: the comparison below is not valid**)") + ".", "",
          "| quantity | AmpliTrim | DADA2 | difference (points) |", "|---|---|---|---|",
          f"| pairs passing the filter | {cmp_['pass_amplitrim']:.2%} | {cmp_['pass_dada2']:.2%} | {cmp_['pass_gap']:+.2f} |",
          f"| forward reads passing alone | {cmp_['forward_amplitrim']:.2%} | {cmp_['forward_dada2']:.2%} | {cmp_['forward_gap']:+.2f} |", "",
          "AmpliTrim's rows of 'where the pairs go' for the same run (shares of the pairs analysed):", "",
          "| step | share |", "|---|---|"]
    for key, label in (("too_short", "too short after the cut at the first base of quality 2 or below"), ("has_n", "contain an N"), ("over_maxEE", "more expected errors than maxEE 2"), ("pass", "pass the filter"),
                       ("pass_not_merged", "pass the filter but do not merge (vsearch's estimate)"), ("pass_and_merged", "pass the filter and merge (vsearch's estimate)")):
        L.append(f"| {label} | {its_loss[key]:.2%} |")
    L += ["", "## B. Where DADA2 loses the pairs after the filter (night 3, AmpliTrim's own setting, maxEE 2)", ""]
    if not loss:
        L += [f"Not available: the track table of '{ARM_A}' was not found" + ("" if night3 is None else f" under {night3}") + ". Give the folder of night 3 with --night3.", ""]
    else:
        L += [f"From {night3}, {loss['samples']} samples, {loss['input']:,.0f} input pairs. Shares of the input pairs:", "", "| step | DADA2 |", "|---|---|",
              f"| pass the filter | {loss['filtered']:.1%} |", f"| reads in denoised forward sequences | {loss['denoised_f']:.1%} |", f"| reads in denoised reverse sequences | {loss['denoised_r']:.1%} |",
              f"| merged | {loss['merged']:.1%} |", f"| after chimera removal | {loss['nonchim']:.1%} |", "",
              f"- Between the filter ({loss['filtered']:.1%}) and the merged pairs ({loss['merged']:.1%}) DADA2 loses {(loss['filtered'] - loss['merged']) * 100:.1f} points of the input pairs: "
              f"at least {loss['lost_denoising'] * 100:.1f} points to denoising (a pair needs both reads in a denoised sequence) and at most {loss['lost_merging'] * 100:.1f} points to merging.",
              f"- AmpliTrim's estimate of the pairs that pass the filter but do not merge is {its_loss['pass_not_merged'] * 100:.1f} points (vsearch, on the raw reads: it cannot see denoising).", ""]
    return L + ["Next: send this file."]


def main(argv=None, run=subprocess.run, find_rscript=None, pipeline_run=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, help="amplitrim_paths.json (default: the one next to the package)")
    ap.add_argument("--package", type=Path, help="the folder that contains the amplitrim package with the ITS report additions (the clean repository); default: this repository")
    ap.add_argument("--night3", type=Path, help="the folder of night 3 (default: the newest night3_* or NN_night3_* in C:\\data\\runs)")
    ap.add_argument("--out", type=Path, default=Path("its_check_out"))
    a = ap.parse_args(argv)
    import rtools
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
    if a.package:                                                  # the package under test comes first, whatever else is on the path
        root, tried = resolve_package(a.package)
        if root is None:
            print(f"The folder given with --package ({a.package}) has no 'amplitrim' package in it. Looked for: " + "; ".join(str(t) for t in tried)
                  + ". Give the folder of the clean repository: the one that contains README.md, pyproject.toml and a folder named amplitrim.")
            return 2
        sys.path.insert(0, str(root.resolve()))
        for name in [m for m in sys.modules if m == "amplitrim" or m.startswith("amplitrim.")]:      # nothing imported earlier may shadow it
            del sys.modules[name]
    import importlib
    pl = importlib.import_module("amplitrim.pipeline")
    print(f"Using the amplitrim package at {Path(pl.__file__).parent}")
    if pipeline_run is None and "its_loss" not in getattr(pl.Result, "__dataclass_fields__", {}):
        print(f"That package has no ITS 'where the pairs go' (Result.its_loss): it is the build that nights 1 to 3 tested. Give --package the folder of the clean repository"
              f"{'' if a.package else ' (it was not given)'}.")
        return 2
    a.out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="its_check_") as tmp:
        res = (pipeline_run or pl.run)(pl.Settings(input_dir=folder, work_dir=Path(tmp) / "work", marker="its"))
        if res.its_loss is None:
            print("AmpliTrim gave no ITS result.")
            return 2
        try:
            counts = dada2_filter_counts(rscript, res.sub.r1, res.sub.r2, run=run)
        except RuntimeError as err:
            print(f"DADA2 could not run: {err}")
            return 2
    cmp_ = compare_filters(res.its_loss, counts)
    night3 = a.night3 or newest_night3(Path(os.environ.get("AMPLITRIM_RUNS", r"C:\data\runs")))
    track = find_night3_track(night3) if night3 else None
    loss = post_filter_loss(read_track(track)) if track else {}
    L = report_lines(cmp_, res.its_loss, loss, night3)
    (a.out / "ITS_CHECK_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    with open(a.out / "its_check.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["quantity", "amplitrim", "dada2", "gap_points"])
        w.writerow(["pass", f"{cmp_['pass_amplitrim']:.5f}", f"{cmp_['pass_dada2']:.5f}", f"{cmp_['pass_gap']:.3f}"])
        w.writerow(["forward_only_pass", f"{cmp_['forward_amplitrim']:.5f}", f"{cmp_['forward_dada2']:.5f}", f"{cmp_['forward_gap']:.3f}"])
    print("\n".join(L[:12]))
    print(f"\nwritten: {a.out / 'ITS_CHECK_REPORT.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
