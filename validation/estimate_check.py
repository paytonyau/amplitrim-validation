#!/usr/bin/env python
"""Estimate check: how close are AmpliTrim's merge estimates to DADA2, and which vsearch cut-off for differences in the overlap is closest?

  python validation\\estimate_check.py --runs C:\\data\\runs --out C:\\data\\runs\\estimate_check               Part A (tables already on disk, seconds)
  python validation\\estimate_check.py --runs C:\\data\\runs --cutoff --out ...                              Parts A and B (vsearch on the reads, tens of minutes)

Part A reads every <run>/<dataset>/compare/compare_methods.tsv of the earlier runs: for each setting AmpliTrim's tolerant estimate (the share passing the filter times the merge
rate with vsearch's default of 10 differences allowed), its strict estimate (no difference allowed, on the raw reads) and DADA2's merged share. It asks which estimate is closer
to DADA2 and which ranks the settings better. Part B asks what number of differences vsearch should allow (0, 1, 2, 3, 4, 5, 7 or 10) to merge as DADA2 does: for each AmpliTrim
setting it merges the filtered pairs at each cut-off and compares the rate with DADA2's merge-step success, taken from the per-sample track of that setting. DADA2 merges
denoised reads, so the step success is merged / min(reads in denoised forward sequences, reads in denoised reverse sequences) summed over the samples: a lower bound of the
true success (on the ITS dataset the exact figure was 80.3% where this gives 77.9%). The criteria (P36 to P38) are in TEST_REGISTER.md, section 23, and were fixed before the run.
Writes ESTIMATE_CHECK_REPORT.md, estimate_check.tsv and, with --cutoff, cutoff.tsv.
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

DIFFS = (0, 1, 2, 3, 4, 5, 7, 10)
RHO_SLACK = 0.05                      # P37: the strict estimate ranks at least as well, within this
CUTOFF_POINTS = 2.0                   # P38: the best intermediate cut-off is closer than both extremes by at least this many points
MIN_ARMS_FOR_RANK = 4


def to_fraction(x):
    if x is None:
        return None
    s = str(x).strip()
    if s in ("", "nan", "None", "NA"):
        return None
    pct = s.endswith("%")
    v = float(s.rstrip("%"))
    return v / 100.0 if (pct or v > 1.5) else v


def parse_ee(text) -> tuple:
    parts = [float(x) for x in str(text).replace(" ", "").split(",") if x != ""]
    return (parts[0], parts[0]) if len(parts) == 1 else (parts[0], parts[1])


def ee_text(ee) -> str:
    return f"{ee[0]:g}" if ee[0] == ee[1] else f"{ee[0]:g},{ee[1]:g}"


def load_rows(runs: Path) -> list:
    """Every setting of every comparison on disk that DADA2 could run: dataset, run, name, cut, maxEE, the estimates and DADA2's merged share (as fractions)."""
    out = []
    for tsv in sorted(Path(runs).glob("*/*/compare/compare_methods.tsv")):
        run, dataset = tsv.parent.parent.parent.name, tsv.parent.parent.name
        with open(tsv, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                merged = to_fraction(r.get("merged"))
                if r.get("error") or merged is None:
                    continue
                out.append({"run": run, "dataset": dataset, "name": r["name"], "truncF": int(float(r.get("truncF") or 0)), "truncR": int(float(r.get("truncR") or 0)),
                            "maxEE": ee_text(parse_ee(r.get("maxEE") or "2")), "est_yield": to_fraction(r.get("est_yield")), "est_strict": to_fraction(r.get("est_strict")),
                            "pass": to_fraction(r.get("pass")), "merged": merged, "folder": tsv.parent})
    return out


def dedupe(rows: list) -> list:
    """The same setting run on several nights is one setting: keep the newest run (the results are identical)."""
    best = {}
    for r in sorted(rows, key=lambda r: r["run"]):
        best[(r["dataset"], r["truncF"], r["truncR"], r["maxEE"])] = r
    return sorted(best.values(), key=lambda r: (r["dataset"], r["name"]))


def _ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r, i = [0.0] * len(v), 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(a, b):
    if len(a) < 3:
        return None
    ra, rb = _ranks(a), _ranks(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    den = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
    return None if den == 0 else sum((x - ma) * (y - mb) for x, y in zip(ra, rb)) / den


def dataset_stats(rows: list) -> dict:
    """Per dataset, over the settings that have both estimates: the mean overshoot (estimate minus DADA2, in points), its mean absolute value, and the rank correlation with DADA2."""
    by = {}
    for r in rows:
        if r["est_yield"] is not None and r["est_strict"] is not None:
            by.setdefault(r["dataset"], []).append(r)
    out = {}
    for ds, rs in by.items():
        tol = [100 * (r["est_yield"] - r["merged"]) for r in rs]
        strict = [100 * (r["est_strict"] - r["merged"]) for r in rs]
        ok = len(rs) >= MIN_ARMS_FOR_RANK
        out[ds] = {"n": len(rs), "tol": sum(tol) / len(rs), "strict": sum(strict) / len(rs), "tol_abs": sum(abs(x) for x in tol) / len(rs), "strict_abs": sum(abs(x) for x in strict) / len(rs),
                   "rho_tol": spearman([r["est_yield"] for r in rs], [r["merged"] for r in rs]) if ok else None,
                   "rho_strict": spearman([r["est_strict"] for r in rs], [r["merged"] for r in rs]) if ok else None}
    return out


def verdicts_a(stats: dict) -> list:
    if not stats:
        return [("P36", "not decided: no settings with both estimates", None), ("P37", "not decided: no settings with both estimates", None)]
    tol_abs = sum(s["tol_abs"] for s in stats.values()) / len(stats)
    strict_abs = sum(s["strict_abs"] for s in stats.values()) / len(stats)
    closer = sum(1 for s in stats.values() if s["strict_abs"] < s["tol_abs"])
    out = [("P36", f"The strict estimate is closer to DADA2's merged share than the tolerant one, as the mean over the datasets of the mean absolute difference: {strict_abs:.1f} against {tol_abs:.1f} points "
                   f"({closer} of {len(stats)} datasets closer)", strict_abs < tol_abs)]
    ranked = [s for s in stats.values() if s["rho_tol"] is not None and s["rho_strict"] is not None]
    if ranked:
        rt, rs = sum(s["rho_tol"] for s in ranked) / len(ranked), sum(s["rho_strict"] for s in ranked) / len(ranked)
        out.append(("P37", f"The strict estimate ranks the settings at least as well as the tolerant one (within {RHO_SLACK}), mean rank correlation over the {len(ranked)} datasets with at least "
                           f"{MIN_ARMS_FOR_RANK} settings: {rs:+.2f} against {rt:+.2f}", rs >= rt - RHO_SLACK))
    else:
        out.append(("P37", f"not decided: no dataset has {MIN_ARMS_FOR_RANK} settings", None))
    return out


def read_track(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as fh:
        return [{k: (float(v) if k != "sample" else v) for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]


def dada2_merge_step(track: Path) -> dict:
    """DADA2's merging success from its per-sample track: merged over the pairs that can reach merging (at most min(denoised forward, denoised reverse))."""
    rows = read_track(track)
    inp = sum(r["input"] for r in rows)
    filtered = sum(r["filtered"] for r in rows)
    reach = sum(min(r["denoised_f"], r["denoised_r"]) for r in rows)
    merged = sum(r["merged"] for r in rows)
    nan = float("nan")
    return {"step": merged / reach if reach else nan, "merged": merged / inp if inp else nan, "filtered": filtered / inp if inp else nan,
            "denoise_loss": (filtered - reach) / filtered if filtered else nan}


def vsearch_rates(pairs, arm: tuple, merge, work: Path, diffs=DIFFS, min_overlap: int = 12, min_len: int = 50) -> dict:
    """The share of the filtered pairs that vsearch merges when it may find at most k differences in the overlap, for each k."""
    tf, tr, ee = arm
    o1, o2 = work / "f_R1.fastq", work / "f_R2.fastq"
    st = pairs.filter(tf, tr, o1, o2, ee) if tf else pairs.filter_untruncated(o1, o2, ee, min_len)
    if st.n_out == 0:
        return {k: 0.0 for k in diffs}
    return {k: merge(o1, o2, work, f"d{k}", min_overlap, k, f"maxdiffs{k}", keep_merged=False).merge_rate for k in diffs}


def cutoff_summary(arms: list, diffs=DIFFS) -> dict:
    """Per cut-off: the mean signed and mean absolute difference, in points, between vsearch's rate and DADA2's merge-step success, over the settings."""
    out = {}
    for k in diffs:
        e = [100 * (a["rates"][k] - a["dada2"]["step"]) for a in arms if k in a["rates"] and a["dada2"]["step"] == a["dada2"]["step"]]
        out[k] = {"n": len(e), "signed": sum(e) / len(e) if e else float("nan"), "abs": sum(abs(x) for x in e) / len(e) if e else float("nan")}
    return out


def verdict_b(summary: dict):
    ok = {k: v for k, v in summary.items() if v["n"]}
    if len(ok) < 3 or 0 not in ok or max(ok) not in ok:
        return ("P38", "not decided: too few cut-offs or settings", None)
    best = min(ok, key=lambda k: ok[k]["abs"])
    ends = min(ok[0]["abs"], ok[max(ok)]["abs"])
    return ("P38", f"An intermediate cut-off is closer to DADA2's merge step than both 0 and {max(ok)} differences, by at least {CUTOFF_POINTS:g} points: the closest is {best} "
                   f"({ok[best]['abs']:.1f} points away on average) against {ok[0]['abs']:.1f} at 0 and {ok[max(ok)]['abs']:.1f} at {max(ok)}", ends - ok[best]["abs"] >= CUTOFF_POINTS)      # an extreme as the best leaves a gain of 0


def report_lines(stats: dict, rows: list, arms: list, summary: dict) -> list:
    L = ["# Estimate check", "", "The criteria (P36 to P38) are in TEST_REGISTER.md, section 23, and were fixed before this ran.", "", "## Verdicts", ""]
    for p, txt, ok in verdicts_a(stats) + ([verdict_b(summary)] if summary else []):
        L.append(f"- **{p}** {'MET' if ok else 'NOT MET' if ok is not None else 'not decided'}: {txt}")
    if not summary:
        L.append("- **P38** not run (give --cutoff)")
    L += ["", f"## A. The two estimates against DADA2 ({len(rows)} settings on {len({r['dataset'] for r in rows})} datasets; points of the input pairs)", "",
          "| dataset | settings | tolerant minus DADA2 | strict minus DADA2 | closer | rank correlation, tolerant | rank correlation, strict |", "|---|---|---|---|---|---|---|"]
    f = lambda x: "-" if x is None else f"{x:+.2f}"          # noqa: E731
    for ds, s in sorted(stats.items()):
        L.append(f"| {ds} | {s['n']} | {s['tol']:+.1f} | {s['strict']:+.1f} | {'strict' if s['strict_abs'] < s['tol_abs'] else 'tolerant' if s['tol_abs'] < s['strict_abs'] else 'same'} | {f(s['rho_tol'])} | {f(s['rho_strict'])} |")
    if summary:
        L += ["", "## B. The number of differences vsearch may allow, against DADA2's merging", "", "| setting | " + " | ".join(f"{k}" for k in summary) + " | DADA2 merge step | denoising loss (at least) |", "|---|" + "---|" * (len(summary) + 2)]
        for a in arms:
            L.append(f"| {a['dataset']}: {a['name']} | " + " | ".join(f"{100 * a['rates'][k]:.1f}%" for k in summary) + f" | {100 * a['dada2']['step']:.1f}% | {100 * a['dada2']['denoise_loss']:.1f}% |")
        L += ["", "| allowed differences | " + " | ".join(f"{k}" for k in summary) + " |", "|---|" + "---|" * len(summary),
              "| vsearch minus DADA2, mean (points) | " + " | ".join(f"{v['signed']:+.1f}" for v in summary.values()) + " |",
              "| vsearch minus DADA2, mean absolute (points) | " + " | ".join(f"{v['abs']:.1f}" for v in summary.values()) + " |", ""]
    return L + ["Next: send this file."]


def main(argv=None, merge=None, build=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=Path, required=True, help="the folder that holds the run folders (for example C:\\data\\runs)")
    ap.add_argument("--cutoff", action="store_true", help="also run part B (vsearch on the reads)")
    ap.add_argument("--config", type=Path, help="amplitrim_paths.json (default: the one next to the package)")
    ap.add_argument("--arms", default="AmpliTrim", help="part B: settings whose name contains this (default AmpliTrim)")
    ap.add_argument("--diffs", default=",".join(map(str, DIFFS)), help="part B: the numbers of differences to try")
    ap.add_argument("--min-overlap", type=int, default=12)
    ap.add_argument("--out", type=Path, default=Path("estimate_check_out"))
    a = ap.parse_args(argv)
    if not a.runs.is_dir():
        print(f"No such folder: {a.runs}")
        return 2
    rows = dedupe(load_rows(a.runs))
    if not rows:
        print("No <run>/<dataset>/compare/compare_methods.tsv found: is --runs the folder that holds night1_..., night3_... and so on?")
        return 2
    stats = dataset_stats(rows)
    arms, summary = [], {}
    if a.cutoff:
        diffs = tuple(int(x) for x in a.diffs.split(","))
        from maxee_rules import reads_folder
        cfg_path = a.config or ROOT / "amplitrim_paths.json"
        cfg = json.loads(cfg_path.read_text(encoding="utf-8")) if cfg_path.exists() else {}
        if merge is None:
            import shutil
            if shutil.which("vsearch") is None:
                print("vsearch was not found on the PATH.")
                return 2
            from amplitrim.merge_wrapper import run_vsearch as merge
        if build is None:
            from amplitrim.subsample import build_subsample as build
        from amplitrim.trimmer import ReadPairs
        from compare_methods import method_key
        with tempfile.TemporaryDirectory(prefix="estimate_check_") as tmp:
            loaded = {}
            for r in rows:
                if a.arms not in r["name"]:
                    continue
                track = r["folder"] / f"track_{method_key(r['name'])}.tsv"
                folder = reads_folder(cfg, r["dataset"])
                if not track.is_file() or folder is None or not folder.is_dir():
                    print(f"  {r['dataset']}: {r['name']}: no track or no reads folder, skipped")
                    continue
                if r["dataset"] not in loaded:
                    info = build(folder, Path(tmp) / r["dataset"] / "sub")
                    loaded[r["dataset"]] = ReadPairs.load(info.r1, info.r2)
                work = Path(tmp) / r["dataset"] / f"{r['truncF']}x{r['truncR']}_{r['maxEE'].replace(',', '_')}"
                work.mkdir(parents=True, exist_ok=True)
                rates = vsearch_rates(loaded[r["dataset"]], (r["truncF"], r["truncR"], parse_ee(r["maxEE"])), merge, work, diffs, a.min_overlap)
                arms.append({"dataset": r["dataset"], "name": r["name"], "rates": rates, "dada2": dada2_merge_step(track)})
                print(f"  {r['dataset']}: {r['name']}: " + ", ".join(f"{k}: {100 * v:.1f}%" for k, v in rates.items()), flush=True)
        if not arms:
            print("Part B found no setting with a track and a reads folder.")
            return 2
        summary = cutoff_summary(arms, diffs)
    L = report_lines(stats, rows, arms, summary)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "ESTIMATE_CHECK_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    with open(a.out / "estimate_check.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dataset", "name", "truncF", "truncR", "maxEE", "est_yield", "est_strict", "merged"])
        for r in rows:
            w.writerow([r["dataset"], r["name"], r["truncF"], r["truncR"], r["maxEE"], "" if r["est_yield"] is None else f"{r['est_yield']:.5f}", "" if r["est_strict"] is None else f"{r['est_strict']:.5f}", f"{r['merged']:.5f}"])
    if arms:
        with open(a.out / "cutoff.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["dataset", "name"] + [f"diffs_{k}" for k in summary] + ["dada2_step", "dada2_denoise_loss"])
            for x in arms:
                w.writerow([x["dataset"], x["name"]] + [f"{x['rates'][k]:.5f}" for k in summary] + [f"{x['dada2']['step']:.5f}", f"{x['dada2']['denoise_loss']:.5f}"])
    print("\n".join(L[:9]))
    print(f"\nwritten: {a.out / 'ESTIMATE_CHECK_REPORT.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
