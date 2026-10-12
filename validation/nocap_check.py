#!/usr/bin/env python
"""Night 5: the rule without a cap, whether a simpler cut would do, and how many pairs wastewater loses to length (TEST_REGISTER.md, section 24).

Part 1, offline (no DADA2, no R; about 2 to 5 minutes per dataset, vsearch on the same subsample as the nights). It runs the released package on the datasets of the
battery and writes NOCAP_CHECK_REPORT.md and nocap_check.tsv:

  P39  the rule: the released package recommends what nights 3 and 4b ran through DADA2 on every dataset except mock B, where the cap is gone and the
       rule is predicted to set maxEE 4,8
  P43  the floor-only cut (the best split of the shortest total that covers the amplicons, at the released maxEE) is the released cut on at least 8 of the 11 datasets
       that have a cut
  P45  descriptive: the share of read pairs shorter than the recommended cut on each dataset (it tells length loss from quality loss)

It also prints the DADA2 commands for the two runs that cannot be done offline: mock B at the released rule's setting (P40 to P42) and, for each dataset where the
floor-only cut differs from the released one, DADA2 at that cut (P44).

  python validation\\nocap_check.py --package C:\\amplitrim_release\\amplitrim --night3 C:\\data\\runs\\night3_2026-10-09 --night4 C:\\data\\runs\\night4_2026-10-10
  python validation\\nocap_check.py --package ... --only mockB wastewater_V3-V4      some datasets (P39 and P43 are then not decided)
  python validation\\nocap_check.py --package ... --like-the-nights                     P39b: the released package on the pairs the nights gave AmpliTrim

Part 2: the table of the mock B run (compare_methods.py with --amplitrim-max-ee 2 2,5 4,5 4,8) is read and P40 to P41 are decided with the margins fixed for the
battery (maxee_arms.py):

  python validation\\nocap_check.py --dada2-table C:\\data\\runs\\night5_mockB\\compare_methods.tsv
"""
import argparse
import csv
import json
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT)]

from release_check import DATASETS, compare, ee_text, parse_ee, settings_text, tested_settings   # noqa: E402
import maxee_arms as ma                                                                          # noqa: E402

EXPECTED_CHANGE = {"mockB": (4.0, 8.0)}      # P39: the one dataset where the cap changed the rule (the offline estimate of the register: q = 0.75 gives (4,8), the cap made it (4,5))
CAPPED, RELEASED = "4,5", "4,8"              # the two settings of P40 and P41 on mock B
PREDICTED_MERGED = (60.0, 76.0)              # P40: where the released rule's merged share on mock B is predicted to fall (the capped rule: 55.0%)
FLOOR_ONLY_SAME = 8                          # P43: at least this many datasets with a cut
PASSING_ACCURACY = ("as accurate", "more accurate")
DADA2_PAIRS = 5000                           # pairs per sample given to DADA2, as in nights 1 to 4b (mock A has all its reads)
PAIRS_PER_SAMPLE = {"mockA_high": 0}          # as night2.PAIRS: the pairs per sample that compare_methods.py thinned each dataset to before AmpliTrim and DADA2 (0: not thinned)
DEFAULT_PAIRS_PER_SAMPLE = 5000


def fraction(x):
    """A value of compare_methods.tsv as a fraction (the table may hold percents)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v / 100.0 if v > 1.5 else v


def inspect(name: str, res, load_pairs) -> dict:
    """What the released package did on one dataset, and what the floor-only cut would be."""
    r = res.recommendation
    released = None if r is None else (int(r.trim_f), int(r.trim_r), parse_ee(res.cfg.max_ee))
    row = {"dataset": name, "released": released, "short_lost": getattr(res, "short_lost", None), "floor": getattr(res, "floor", None), "floor_only": None, "noise": None}
    if released is None or name == "its":                                   # ITS has no cut
        return row
    pairs = load_pairs(res.sub.r1, res.sub.r2)
    best = pairs.pass_table(released[2]).best_split(res.floor)
    if best:
        row["floor_only"] = (int(best[0]), int(best[1]), released[2])
    load = getattr(pairs, "error_load", None)
    if load is not None:
        row["noise"] = load(released[0], released[1], released[2])
    return row


def p39(rows: list, tested_by) -> tuple:
    """(verdict, per-dataset status). The verdict needs every dataset; an unexpected difference, or a different change on mock B, decides it as NOT MET."""
    status, undecided, bad = {}, False, False
    for r in rows:
        name = r["dataset"]
        c = compare(name, r["released"], tested_by(name))
        if c["same"] is None:
            status[name], undecided = "tested result not found", True
        elif name in EXPECTED_CHANGE:                                     # the one dataset that is predicted to change
            if c["same"]:
                status[name], bad = "**unchanged: is the cap still in this package?**", True
            elif r["released"] is not None and r["released"][2] == EXPECTED_CHANGE[name]:
                status[name] = "expected change: the cap is gone"
            else:
                status[name], bad = "**changed, but not to the predicted maxEE " + ee_text(EXPECTED_CHANGE[name]) + "**", True
        elif c["same"]:
            status[name] = "same"
        else:
            status[name], bad = "**NO**", True
    complete = {r["dataset"] for r in rows} >= set(DATASETS)
    return (False if bad else (None if undecided or not complete else True)), status


def p43(rows: list) -> tuple:
    """(verdict, same, total): the floor-only cut against the released cut, on the datasets that have a cut."""
    have = [r for r in rows if r["floor_only"] is not None and r["released"] is not None]
    same = [r for r in have if r["floor_only"][:2] == r["released"][:2]]
    complete = {r["dataset"] for r in rows} >= set(DATASETS)
    return (None if not have or not complete else len(same) >= FLOOR_ONLY_SAME), len(same), len(have)


def read_table(path: Path) -> dict:
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8"), delimiter="\t"))
    return {r["name"]: r for r in rows if not r.get("error") and fraction(r.get("merged")) is not None}


def evaluate_table(rows: dict, capped: str = CAPPED, released: str = RELEASED) -> dict:
    """P40 and P41 from the table of the mock B run: the released rule's setting against the capped one."""
    h, n = rows.get(f"AmpliTrim (maxEE {capped})"), rows.get(f"AmpliTrim (maxEE {released})")
    if h is None or n is None:
        return {"error": f"the table has no row 'AmpliTrim (maxEE {capped})' or 'AmpliTrim (maxEE {released})': run compare_methods.py with --amplitrim-max-ee 2 2,5 {capped} {released}"}
    mh, mn = 100 * fraction(h["merged"]), 100 * fraction(n["merged"])
    row = lambda x: {"found": x.get("found"), "strains": x.get("strains"), "identical": fraction(x.get("identical"))}      # noqa: E731
    acc = ma.accuracy_verdict(row(n), row(h))
    return {"capped": (h["truncF"], h["truncR"], capped, mh), "released": (n["truncF"], n["truncR"], released, mn), "gain": mn - mh, "accuracy": acc,
            "p40": mn - mh > ma.MARGIN_MERGED, "p40_in_range": PREDICTED_MERGED[0] <= mn <= PREDICTED_MERGED[1], "p41": acc in PASSING_ACCURACY,
            "found": (h.get("found"), n.get("found")), "strains": (h.get("strains"), n.get("strains")), "identical": (fraction(h.get("identical")), fraction(n.get("identical")))}


def table_lines(e: dict) -> list:
    if "error" in e:
        return [e["error"]]
    yn = lambda ok: "MET" if ok else "NOT MET"                                                    # noqa: E731
    pc = lambda x: "-" if x is None else f"{100 * x:.1f}%"                                          # noqa: E731
    return ["# Night 5: mock B at the released rule's setting", "",
            f"- **P40** {yn(e['p40'])}: the released rule's setting ({e['released'][0]}/{e['released'][1]}, maxEE {e['released'][2]}) merged {e['released'][3]:.1f}% of the pairs against "
            f"{e['capped'][3]:.1f}% for the capped rule ({e['capped'][0]}/{e['capped'][1]}, maxEE {e['capped'][2]}): {e['gain']:+.1f} points (the margin is {ma.MARGIN_MERGED:g}). "
            f"Predicted range {PREDICTED_MERGED[0]:g}% to {PREDICTED_MERGED[1]:g}%: {'inside' if e['p40_in_range'] else 'OUTSIDE'}.",
            f"- **P41** {yn(e['p41'])}: accuracy of the released setting against the capped one: {e['accuracy']} (variants found {e['found'][1]} against {e['found'][0]}; strains "
            f"{e['strains'][1]} against {e['strains'][0]}; reads in identical ASVs {pc(e['identical'][1])} against {pc(e['identical'][0])}). Only 'as accurate' and 'more accurate' pass.",
            "", "Next: send this file and compare_methods.tsv."]


def command(reads, out: str, ees: str, extra: str = "", dataset: str = "") -> str:
    pairs = PAIRS_PER_SAMPLE.get(dataset, DEFAULT_PAIRS_PER_SAMPLE)                    # DADA2 gets what it got in the nights: 5,000 pairs per sample, all the reads of mock A
    return f'python validation\\compare_methods.py "{reads}" --out {out} --amplitrim-max-ee {ees}' + (f" --pairs {pairs}" if pairs else "") + extra


def report_lines(rows: list, tested_by, reads_of, like_the_nights: bool = False) -> list:
    v40, st = p39(rows, tested_by)
    v44, n_same, n_have = p43(rows)
    yn = lambda v: "MET" if v else "NOT MET" if v is not None else "not decided (not all 12 datasets were run, or a tested result was not found)"       # noqa: E731
    L = ["# Night 5: the rule without a cap, the floor-only cut, and reads lost to length", "",
         "The criteria are in TEST_REGISTER.md, section 24, and were fixed before this ran.", "",
         ("**The pairs:** the ones the nights gave AmpliTrim (every dataset thinned to 5,000 pairs per sample, mock A not thinned, as compare_methods.py does); this is P39b." if like_the_nights
          else "**The pairs:** the package's own sample of the full reads (every Nth pair, up to 50,000 in all), not the pairs the nights gave AmpliTrim; this is P39."), "", "## Verdicts", "",
         f"- **P39** {yn(v40)}: the released package recommends what was tested on every dataset except mock B, where it is predicted to set maxEE {ee_text(EXPECTED_CHANGE['mockB'])}.",
         f"- **P43** {yn(v44)}: the floor-only cut is the released cut on {n_same} of {n_have} datasets with a cut (at least {FLOOR_ONLY_SAME} predicted).", "",
         "## P39: the released rule against what DADA2 was shown", "", "| dataset | tested | what was tested | released package | status |", "|---|---|---|---|---|"]
    for r in rows:
        t = tested_by(r["dataset"])
        L.append(f"| {r['dataset']} | {t[0] if t else 'not found'} | {settings_text(None if not t else (t[1], t[2], t[3]))} | {settings_text(r['released'])} | {st[r['dataset']]} |")
    L += ["", "## P43: the floor-only cut", "", "| dataset | floor (nt) | released cut | total | floor-only cut | the same cut |", "|---|---|---|---|---|---|"]
    for r in rows:
        x, f = r["released"], r["floor_only"]
        L.append(f"| {r['dataset']} | {r['floor'] if r['floor'] is not None else '-'} | {'-' if x is None else f'{x[0]}/{x[1]}'} | {'-' if x is None else x[0] + x[1]} | "
                 f"{'-' if f is None else f'{f[0]}/{f[1]}'} | {'-' if f is None or x is None else ('yes' if f[:2] == x[:2] else '**no**')} |")
    L += ["", "## P45: read pairs shorter than the recommended cut (descriptive)", "", "| dataset | pairs lost because a read is shorter than the cut |", "|---|---|"]
    for r in rows:
        s = r["short_lost"]
        L.append(f"| {r['dataset']} | {'-' if s is None else f'{100 * s:.1f}%'}{'  (the package warns from 1%)' if s is not None and s >= 0.01 else ''} |")
    L += ["", "## DADA2 runs that cannot be done offline", ""]
    mb = next((r for r in rows if r["dataset"] == "mockB"), None)
    if mb is not None:
        L += ["P40 to P42, mock B at the released rule's setting (the expected sequences of the mock are given as in night 4b; the table is read with part 2):", "", "```",
              command(reads_of("mockB"), "C:\\data\\runs\\night5_mockB", f"2 2,5 {CAPPED} {RELEASED}", ' --reference "<the expected sequences of mock B, as in night 4b>"', "mockB"), "",
              "python validation\\nocap_check.py --dada2-table C:\\data\\runs\\night5_mockB\\compare_methods.tsv", "```", ""]
        if mb["noise"]:
            L += [f"P42 (descriptive): expected errors per kept read at the released cut and maxEE: {mb['noise'][0]:.2f} forward, {mb['noise'][1]:.2f} reverse (the offline estimate of the register: 1.02 and 3.31).", ""]
    differ = [r for r in rows if r["floor_only"] is not None and r["released"] is not None and r["floor_only"][:2] != r["released"][:2]]
    if differ:
        L += ["P44, DADA2 at the floor-only cut where it differs from the released cut. The table lists 'released' and 'floor_only' and compares every row with 'released' (--baseline), so the line for "
              "'floor_only' is floor-only minus released: at least -2 points on every dataset is MET.", ""]
        if not like_the_nights:
            L += ["**These cuts come from the package's own sample of the full reads, and DADA2 sees the nights' pairs.** Run this check again with `--like-the-nights` and take the commands from that report.", ""]
        L += ["```"]
        for r in differ:
            x, f = r["released"], r["floor_only"]
            ee = ee_text(x[2])
            L.append(command(reads_of(r["dataset"]), f"C:\\data\\runs\\night5_{r['dataset']}", ee, f" --setting released={x[0]},{x[1]},{ee} --setting floor_only={f[0]},{f[1]},{ee} --baseline released", r["dataset"]))
        L += ["```", ""]
    else:
        L += ["P44: the floor-only cut is the released cut on every dataset, so DADA2 would give the same result at both: nothing to run.", ""]
    return L + ["Next: send this file."]


def main(argv=None, pipeline_run=None, load_pairs=None, thin_reads=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--package", type=Path, help="the folder of the clean repository (the package under test)")
    ap.add_argument("--config", type=Path, help="amplitrim_paths.json (default: the one next to the package)")
    ap.add_argument("--night3", type=Path, help="the folder of night 3")
    ap.add_argument("--night4", type=Path, help="the folder of night 4b")
    ap.add_argument("--only", nargs="+", help="some datasets only (P39 and P43 are then not decided)")
    ap.add_argument("--like-the-nights", action="store_true", help="P39b: give the package the pairs the nights gave AmpliTrim (every dataset thinned to 5,000 pairs per sample, mock A not thinned, as compare_methods.py does) and not the full reads")
    ap.add_argument("--out", type=Path, default=Path("nocap_check_out"))
    ap.add_argument("--dada2-table", type=Path, help="part 2: the compare_methods.tsv of the mock B run")
    a = ap.parse_args(argv)
    if a.dada2_table:
        L = table_lines(evaluate_table(read_table(a.dada2_table)))
        print("\n".join(L))
        a.out.mkdir(parents=True, exist_ok=True)
        (a.out / "NOCAP_DADA2_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
        return 0 if not L[0].startswith("the table has no row") else 2
    import importlib
    from maxee_rules import reads_folder
    if a.package and pipeline_run is None:
        from its_check import resolve_package
        root, tried = resolve_package(a.package)
        if root is None:
            print(f"The folder given with --package ({a.package}) has no 'amplitrim' package in it. Looked for: " + "; ".join(str(t) for t in tried))
            return 2
        sys.path.insert(0, str(root.resolve()))
        for m in [m for m in sys.modules if m == "amplitrim" or m.startswith("amplitrim.")]:
            del sys.modules[m]
    pl = importlib.import_module("amplitrim.pipeline")
    if pipeline_run is None:
        print(f"Using the amplitrim package at {Path(pl.__file__).parent}")
        if "rule" not in getattr(pl.Result, "__dataclass_fields__", {}) or not hasattr(pl, "WARN_ABOVE"):
            print("That package is not the release without a cap (it has no WARN_ABOVE): give --package the folder of the clean repository, after the change.")
            return 2
    if load_pairs is None:
        from amplitrim.trimmer import ReadPairs
        load_pairs = ReadPairs.load
    if a.like_the_nights and thin_reads is None:
        from real_data import thin_reads                                   # the function compare_methods.py used; it reads with the package on the path
    cfg_path = a.config or ROOT / "amplitrim_paths.json"
    if not cfg_path.is_file():
        print(f"No amplitrim_paths.json was found at {cfg_path}.\nGive --config with the path of the file the nights used (it is not in the GitHub clone; search for it with:  where /r C:\\ amplitrim_paths.json ).")
        return 2
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
    print(f"Using the datasets in {cfg_path}")
    tested_by = lambda name: tested_settings(a.night3, a.night4, name) if (a.night3 or a.night4) else None          # noqa: E731
    reads_of = lambda name: reads_folder(cfg, name) or f"<the reads folder of {name}>"                              # noqa: E731
    rows = []
    for name in DATASETS:
        if a.only and name not in a.only:
            continue
        folder = reads_folder(cfg, name)
        if folder is None or not folder.is_dir():
            print(f"{name}: no reads folder for it in {cfg_path} (the key is 'extras' > '{name}' > 'ready'), skipped")
            continue
        with tempfile.TemporaryDirectory(prefix="nocap_check_") as tmp:
            reads = folder
            if a.like_the_nights and PAIRS_PER_SAMPLE.get(name, DEFAULT_PAIRS_PER_SAMPLE):
                reads = thin_reads(folder, Path(tmp) / "thin", PAIRS_PER_SAMPLE.get(name, DEFAULT_PAIRS_PER_SAMPLE))
            res = (pipeline_run or pl.run)(pl.Settings(input_dir=reads, work_dir=Path(tmp) / "work", marker="its" if name == "its" else "16s"))
            rows.append(inspect(name, res, load_pairs))
        x, f, lost = rows[-1]["released"], rows[-1]["floor_only"], rows[-1]["short_lost"]
        floor_text = "-" if f is None else "%d/%d" % (f[0], f[1])
        lost_text = "-" if lost is None else "%.1f%%" % (100 * lost)
        print(f"{name}: released {settings_text(x)} | floor-only {floor_text} | pairs shorter than the cut: {lost_text}", flush=True)
    if not rows:
        print("nothing was run: check amplitrim_paths.json")
        return 2
    a.out.mkdir(parents=True, exist_ok=True)
    L = report_lines(rows, tested_by, reads_of, a.like_the_nights)
    (a.out / "NOCAP_CHECK_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    with open(a.out / "nocap_check.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dataset", "released_truncF", "released_truncR", "released_maxEE", "floor", "floor_only_F", "floor_only_R", "short_lost", "noise_F", "noise_R"])
        for r in rows:
            x, f, nz = r["released"], r["floor_only"], r["noise"]
            w.writerow([r["dataset"], *(["", "", ""] if not x else [x[0], x[1], ee_text(x[2])]), "" if r["floor"] is None else r["floor"], *(["", ""] if not f else [f[0], f[1]]),
                        "" if r["short_lost"] is None else round(r["short_lost"], 4), *(["", ""] if not nz else [round(nz[0], 3), round(nz[1], 3)])])
    print("\n".join(L[:9]))
    print(f"\nwritten: {a.out / 'NOCAP_CHECK_REPORT.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
