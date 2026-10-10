#!/usr/bin/env python
"""Can a maxEE rule do better than FIGARO's, and what would it cost? Answered from the reads, without DADA2 or vsearch.

A looser maxEE always keeps more pairs and always keeps noisier ones, so "better than FIGARO" is not one number. This script puts FIGARO's own
choice (its truncation and its maxEE, from the table that compare_methods.py saved) on the two axes that move together, using AmpliTrim's exact
quality table: the share of pairs that pass DADA2's filter, and the mean expected errors per kept read (forward / reverse). Then, at AmpliTrim's
truncation, it shows

  * the loosest-needed maxEE that reaches FIGARO's retention, and whether its kept reads are cleaner or noisier than FIGARO's;
  * what each candidate rule would pick: the default (2,2), the quality-first rule of --explore-max-ee, a retention target (60, 70, 80%), and
    the rule of FIGARO's score (retention in percent minus (F-1)^2 minus (R-1)^2) applied to exact pass rates.

  python validation\\maxee_frontier.py --reads C:\\data\\new\\mockB\\prepared\\mockB\\ready --compare C:\\data\\overnight\\extra_datasets\\mockB\\compare

It prints a table and writes maxee_frontier.tsv next to the comparison. It runs in under a minute per dataset.
"""
import argparse
import csv
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from amplitrim.config import ee_pair, parse_ee                       # noqa: E402
from amplitrim.subsample import build_subsample                     # noqa: E402
from amplitrim.trimmer import ReadPairs                              # noqa: E402

LADDER = [(1, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 3), (3, 5), (5, 5)]          # what --explore-max-ee offers
F_RANGE, R_RANGE = range(1, 7), range(1, 10)                                       # the wider grid: FIGARO has chosen 6,9 on a real mock


def fig_score(pass_rate: float, ee) -> float:
    """FIGARO's weighting as AmpliTrim's report describes it: retention in percent minus (F-1)^2 minus (R-1)^2."""
    return 100 * pass_rate - ((ee[0] - 1) ** 2 + (ee[1] - 1) ** 2)


def looseness(ee) -> tuple:
    return (ee[0] + ee[1], ee[0], ee[1])


def evaluate(pairs: ReadPairs, trim_f: int, trim_r: int, ee) -> dict:
    ee = ee_pair(ee)
    ef, er = pairs.error_load(trim_f, trim_r, ee)
    rate = pairs.pass_table(ee).rate(trim_f, trim_r)
    return {"ee": ee, "pass": rate, "ee_f": ef, "ee_r": er, "score": fig_score(rate, ee)}


def grid_at(pairs: ReadPairs, trim_f: int, trim_r: int, f_range=F_RANGE, r_range=R_RANGE) -> list:
    """Every maxEE of the wider grid at one truncation, strictest first."""
    return sorted((evaluate(pairs, trim_f, trim_r, (f, r)) for f in f_range for r in r_range), key=lambda x: looseness(x["ee"]))


def strictest_reaching(grid: list, target: float):
    """The strictest maxEE whose pass rate reaches `target`; None when no maxEE of the grid does."""
    for x in grid:
        if x["pass"] >= target - 1e-12:
            return x
    return None


def quality_first(grid: list, min_retention: float = 0.5):
    """The rule of --explore-max-ee: the strictest setting of the ladder that keeps at least `min_retention` of the pairs (never below 2,2); else the best of the ladder."""
    ladder = [x for x in grid if x["ee"] in LADDER and x["ee"] != (1, 1)]
    return strictest_reaching(ladder, min_retention) or max(ladder, key=lambda x: x["pass"])


def best_by_score(grid: list):
    return max(grid, key=lambda x: (x["score"], -looseness(x["ee"])[0]))


def cleaner_or_noisier(mine: dict, fig: dict, tol: float = 0.02) -> str:
    """Compare the kept reads of two settings: 'cleaner' when both reads have no more expected errors and one has fewer, 'noisier' for the reverse, else 'mixed'."""
    d_f, d_r = mine["ee_f"] - fig["ee_f"], mine["ee_r"] - fig["ee_r"]
    if d_f <= tol and d_r <= tol and (d_f < -tol or d_r < -tol):
        return "cleaner"
    if d_f >= -tol and d_r >= -tol and (d_f > tol or d_r > tol):
        return "noisier"
    return "about the same" if abs(d_f) <= tol and abs(d_r) <= tol else "mixed"


def read_methods(compare_folder: Path) -> list:
    table = Path(compare_folder) / "compare_methods.tsv"
    if not table.is_file():
        raise SystemExit(f"{table} was not found: give the output folder of a compare_methods.py run.")
    rows = []
    with open(table, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("error") or not r.get("truncF") or not r.get("truncR"):
                continue
            try:
                tf, tr = int(float(r["truncF"])), int(float(r["truncR"]))
                ee = parse_ee(r["maxEE"])
            except (ValueError, TypeError):
                continue
            if tf > 0 and tr > 0:
                rows.append({"name": r["name"], "tf": tf, "tr": tr, "ee": ee_pair(ee)})
    return rows


def analyse(pairs: ReadPairs, methods: list) -> dict:
    """The whole analysis for one dataset, as data (printed by main)."""
    out = {"methods": [], "figaro": None, "amplitrim": None, "reach": None, "rules": []}
    for m in methods:
        out["methods"].append({**m, **evaluate(pairs, m["tf"], m["tr"], m["ee"])})
    fig = next((x for x in out["methods"] if x["name"].upper().startswith("FIGARO")), None)
    mine = next((x for x in out["methods"] if x["name"].startswith("AmpliTrim (maxEE")), None) or next((x for x in out["methods"] if x["name"].startswith("AmpliTrim")), None)
    out["figaro"], out["amplitrim"] = fig, mine
    if mine is None:
        return out
    grid = grid_at(pairs, mine["tf"], mine["tr"])
    if fig is not None:
        reach = strictest_reaching(grid, fig["pass"])
        out["reach"] = {**reach, "against": cleaner_or_noisier(reach, fig)} if reach else None
    rules = [("default (2,2)", next(x for x in grid if x["ee"] == (2, 2))), ("quality-first (the --explore-max-ee rule, 50%)", quality_first(grid))]
    for t in (0.6, 0.7, 0.8):
        pick = strictest_reaching(grid, t)
        rules.append((f"retention target {t:.0%}", pick if pick else max(grid, key=lambda x: x["pass"])))
    rules.append(("FIGARO's score on exact pass rates", best_by_score(grid)))
    out["rules"] = [{"rule": name, **x, "against": cleaner_or_noisier(x, fig) if fig else ""} for name, x in rules]
    return out


def ee_text(ee) -> str:
    return f"{ee[0]:g},{ee[1]:g}"


def report(name: str, res: dict) -> list:
    L = [f"=== {name}", "", "| setting | truncation | maxEE | pairs passing | expected errors per kept read (F / R) | FIGARO's score |", "|---|---|---|---|---|---|"]
    for x in res["methods"]:
        L.append(f"| {x['name']} | {x['tf']}/{x['tr']} | {ee_text(x['ee'])} | {x['pass']:.1%} | {x['ee_f']:.2f} / {x['ee_r']:.2f} | {x['score']:.1f} |")
    fig, mine, reach = res["figaro"], res["amplitrim"], res["reach"]
    if fig is None or mine is None:
        return L + ["", "No FIGARO row or no AmpliTrim row in the comparison, so nothing to put against FIGARO's choice.", ""]
    L += ["", f"FIGARO keeps {fig['pass']:.1%} of the pairs at {fig['tf']}/{fig['tr']} with maxEE {ee_text(fig['ee'])}; its kept reads carry {fig['ee_f']:.2f} / {fig['ee_r']:.2f} expected errors."]
    if reach is None:
        L.append(f"AmpliTrim's truncation ({mine['tf']}/{mine['tr']}) does not reach that retention at any maxEE up to 6,9.")
    else:
        L.append(f"AmpliTrim's truncation ({mine['tf']}/{mine['tr']}) reaches it at maxEE {ee_text(reach['ee'])} ({reach['pass']:.1%} of the pairs), where the kept reads carry "
                 f"{reach['ee_f']:.2f} / {reach['ee_r']:.2f} expected errors: {reach['against']} than FIGARO's.")
    L += ["", f"What each rule would pick at AmpliTrim's truncation ({mine['tf']}/{mine['tr']}):", "", "| rule | maxEE | pairs passing | expected errors per kept read (F / R) | against FIGARO's kept reads |", "|---|---|---|---|---|"]
    for x in res["rules"]:
        L.append(f"| {x['rule']} | {ee_text(x['ee'])} | {x['pass']:.1%} | {x['ee_f']:.2f} / {x['ee_r']:.2f} | {x['against']} |")
    return L + [""]


def write_tsv(path: Path, name: str, res: dict) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dataset", "what", "name", "truncF", "truncR", "maxEE", "pass", "ee_f", "ee_r", "score", "against_figaro"])
        for x in res["methods"]:
            w.writerow([name, "setting", x["name"], x["tf"], x["tr"], ee_text(x["ee"]), f"{x['pass']:.4f}", f"{x['ee_f']:.3f}", f"{x['ee_r']:.3f}", f"{x['score']:.2f}", ""])
        if res["reach"]:
            x = res["reach"]
            w.writerow([name, "reaches FIGARO's retention", "AmpliTrim truncation", res["amplitrim"]["tf"], res["amplitrim"]["tr"], ee_text(x["ee"]), f"{x['pass']:.4f}", f"{x['ee_f']:.3f}", f"{x['ee_r']:.3f}", f"{x['score']:.2f}", x["against"]])
        for x in res["rules"]:
            w.writerow([name, "rule", x["rule"], res["amplitrim"]["tf"], res["amplitrim"]["tr"], ee_text(x["ee"]), f"{x['pass']:.4f}", f"{x['ee_f']:.3f}", f"{x['ee_r']:.3f}", f"{x['score']:.2f}", x["against"]])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads", type=Path, required=True, help="the folder of reads the comparison was run on (the dataset's ready folder)")
    ap.add_argument("--compare", type=Path, required=True, help="the output folder of compare_methods.py (holds compare_methods.tsv)")
    ap.add_argument("--max-pairs", type=int, default=60000, help="pairs pooled over all samples (default 60,000)")
    ap.add_argument("--name", help="a label for the dataset (default: the name of the folder that holds the comparison's parent)")
    args = ap.parse_args(argv)
    methods = read_methods(args.compare)
    with tempfile.TemporaryDirectory() as tmp:
        info = build_subsample(args.reads, Path(tmp), max_pairs=args.max_pairs)
        pairs = ReadPairs.load(info.r1, info.r2)
        print(f"{pairs.n:,} pairs pooled from {args.reads}\n")
        res = analyse(pairs, methods)
    name = args.name or args.compare.resolve().parent.name
    print("\n".join(report(name, res)))
    write_tsv(args.compare / "maxee_frontier.tsv", name, res)
    print(f"written: {args.compare / 'maxee_frontier.tsv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
