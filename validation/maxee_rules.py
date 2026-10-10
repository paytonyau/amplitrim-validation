#!/usr/bin/env python
"""Night 4 (offline): can one simple formula choose maxEE the way FIGARO does where FIGARO adapts, and leave DADA2's default alone elsewhere?

FIGARO's maxEE is not a black box. At a cut it is  int(expected errors of a typical read at that cut) + 1  (padMaxExpectedError: truncation, then +1),
where "a typical read" comes from a curve FIGARO fits to the reads. This script asks whether the same formula on the reads themselves, with no curve
and no FIGARO, reproduces FIGARO's picks:

  primary rule    per read: int(MEAN expected errors at the cut) + 1
  secondary rule  per read: int(MEDIAN expected errors at the cut) + 1
  floored         each of them never below 2 (DADA2's default): FIGARO's 1 on clean reads is only what the formula gives when a typical read has
                  under one expected error, and it merged 2 to 6 points fewer pairs than the default there

It uses the reads only (AmpliTrim's exact quality table): no DADA2, no vsearch, no FIGARO; a minute or two per dataset. The criteria (P13 to P16) are
in TEST_REGISTER.md, section 16, and were fixed before this script was run on any dataset.

  python validation\\maxee_rules.py                       every registered dataset that has a FIGARO row
  python validation\\maxee_rules.py --only mockB v4        some of them
  python validation\\maxee_rules.py --out C:\\data\\runs\\night4_offline

It writes MAXEE_RULES_REPORT.md and maxee_rules.tsv to --out.
"""
import argparse
import csv
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from amplitrim.config import ee_pair                                   # noqa: E402
from amplitrim.subsample import build_subsample                       # noqa: E402
from amplitrim.trimmer import ReadPairs, _cut_errors                  # noqa: E402

# FIGARO's cut, its maxEE, and the share of pairs DADA2's own filter passed there; AmpliTrim's cut at maxEE 2. All from NIGHT3_REPORT.md (night 3,
# arms C and A), written here before this script was run. A single maxEE means the same limit on both reads.
ROWS = {
    "mockA_high":       dict(fig=(273, 247, (3, 4)), dada2_pass=0.803, amp=(277, 243)),
    "v1v2":             dict(fig=(225, 138, (2, 1)), dada2_pass=0.882, amp=(201, 157)),
    "v4":               dict(fig=(120, 119, (1, 1)), dada2_pass=0.965, amp=(124, 115)),
    "mcm_v12":          dict(fig=(156, 169, (1, 1)), dada2_pass=0.941, amp=(155, 170)),
    "mockB":            dict(fig=(279, 248, (6, 9)), dada2_pass=0.759, amp=(279, 248)),
    "mcm_v46":          dict(fig=(270, 268, (4, 4)), dada2_pass=0.704, amp=(273, 265)),
    "d6300_V3-V4":      dict(fig=(224, 218, (1, 1)), dada2_pass=0.882, amp=(224, 218)),
    "wastewater_V3-V4": dict(fig=(247, 193, (5, 5)), dada2_pass=0.831, amp=(254, 186)),
    "spikein_V3-V4":    dict(fig=(232, 207, (2, 1)), dada2_pass=0.816, amp=(222, 218)),
    "swabs_V3-V4":      dict(fig=(224, 218, (1, 1)), dada2_pass=0.786, amp=(224, 218)),
}
STATS = {"mean": np.mean, "median": np.median, "p75": lambda x: np.percentile(x, 75), "p90": lambda x: np.percentile(x, 90)}
VARIANTS = ((0.75, None), (0.75, 5), (0.80, None), (0.80, 5))              # EXPLORATORY: q and cap (a cap of 5 is DADA2's tutorial example, c(2, 5)); added after the percentile run
QGRID = (0.75, 0.80, 0.85, 0.90)                                         # the percentile rule's q, fixed before it was run (section 17 of the register)
POOR = [k for k, v in ROWS.items() if max(v["fig"][2]) >= 3]            # FIGARO went looser than 2 on at least one read
CLEAN = [k for k, v in ROWS.items() if max(v["fig"][2]) <= 2]           # FIGARO's maxEE is 2 or less on both reads
MAX_PAIRS = 60000
SANITY_POINTS = 0.02                                                     # AmpliTrim's filter replica was within 0.8 points of DADA2's: more than 2 means the wrong reads


def figaro_pad(x: float) -> int:
    """FIGARO's padMaxExpectedError, verbatim: -(int(-x)) + 1. For x >= 0 this is int(x) + 1 (truncation, not a ceiling)."""
    return -(int(-x)) + 1


def read_errors(pairs: ReadPairs, trim_f: int, trim_r: int):
    """Expected errors of every read over its first trim bases, forward and reverse, for the reads that are at least that long."""
    ef = _cut_errors(pairs._phred1, np.full(pairs.n, trim_f), pairs.len1 >= trim_f)
    er = _cut_errors(pairs._phred2, np.full(pairs.n, trim_r), pairs.len2 >= trim_r)
    return ef, er


def rules(ef: np.ndarray, er: np.ndarray) -> dict:
    """The rules' maxEE (forward, reverse) from the expected errors of the reads at a cut."""
    out = {}
    for name, stat in (("mean", np.mean), ("median", np.median)):
        raw = (figaro_pad(float(stat(ef))), figaro_pad(float(stat(er))))
        out[name] = raw
        out[name + "_floored"] = (max(2, raw[0]), max(2, raw[1]))
    return out


def q_rule(ef: np.ndarray, er: np.ndarray, q: float, cap=None) -> tuple:
    """The percentile rule: per read, int(q-th percentile of the expected errors at the cut) + 1. Returns (raw, floored): the same never below 2 and, when a cap is
    given, never above it (the cap does not touch `raw`)."""
    raw = (figaro_pad(float(np.percentile(ef, 100 * q))), figaro_pad(float(np.percentile(er, 100 * q))))
    lim = (lambda v: min(cap, v)) if cap else (lambda v: v)
    return raw, (lim(max(2, raw[0])), lim(max(2, raw[1])))


def share_below(x: np.ndarray, t: float) -> float:
    return float((x < t).mean())


def within(a, b, tol: int = 1) -> bool:
    a, b = ee_pair(a), ee_pair(b)
    return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol


def point(pairs: ReadPairs, tf: int, tr: int, ee) -> dict:
    ee = ee_pair(ee)
    ef, er = pairs.error_load(tf, tr, ee)
    return {"ee": ee, "pass": pairs.pass_table(ee).rate(tf, tr), "ee_f": ef, "ee_r": er}


def analyse(pairs: ReadPairs, row: dict) -> dict:
    """Everything for one dataset, as data."""
    ftf, ftr, fee = row["fig"]
    atf, atr = row["amp"]
    f_f, f_r = read_errors(pairs, ftf, ftr)
    a_f, a_r = read_errors(pairs, atf, atr)
    for label, cut, arrays in (("FIGARO's", (ftf, ftr), (f_f, f_r)), ("AmpliTrim's", (atf, atr), (a_f, a_r))):
        if any(len(x) == 0 for x in arrays):
            raise ValueError(f"no reads are long enough for {label} cut {cut[0]}/{cut[1]} (reads are up to {pairs.width_f}/{pairs.width_r} nt): "
                             f"is this the folder of reads that night 3 ran on?")
    at_fig, at_amp = rules(f_f, f_r), rules(a_f, a_r)
    res = {
        "figaro": {"cut": (ftf, ftr), "pick": ee_pair(fee), "replica": point(pairs, ftf, ftr, fee), "dada2_pass": row["dada2_pass"],
                   "mean_ee": (float(f_f.mean()), float(f_r.mean())), "median_ee": (float(np.median(f_f)), float(np.median(f_r))), "rules": at_fig,
                   "stats": {name: (float(fn(f_f)), float(fn(f_r))) for name, fn in STATS.items()},
                   "q": {q: {"raw": q_rule(f_f, f_r, q)[0], "match": within(q_rule(f_f, f_r, q)[0], fee)} for q in QGRID},
                   "cdf": ((share_below(f_f, ee_pair(fee)[0] - 1), share_below(f_f, ee_pair(fee)[0])), (share_below(f_r, ee_pair(fee)[1] - 1), share_below(f_r, ee_pair(fee)[1]))),
                   "match": {name: within(v, fee) for name, v in at_fig.items()}},
        "amplitrim": {"cut": (atf, atr), "mean_ee": (float(a_f.mean()), float(a_r.mean())), "median_ee": (float(np.median(a_f)), float(np.median(a_r))),
                      "rules": at_amp, "match": {name: within(v, fee) for name, v in at_amp.items()},
                      "points": {"default (2,2)": point(pairs, atf, atr, (2, 2)), "mean rule, floored": point(pairs, atf, atr, at_amp["mean_floored"]),
                                 "FIGARO's pick": point(pairs, atf, atr, fee)},
                      "q": {q: {**point(pairs, atf, atr, q_rule(a_f, a_r, q)[1]), "raw": q_rule(a_f, a_r, q)[0], "match": within(q_rule(a_f, a_r, q)[1], fee)} for q in QGRID},
                      "variants": {f"{q:g}|{cap or 0}": {**point(pairs, atf, atr, q_rule(a_f, a_r, q, cap)[1]), "q": q, "cap": cap} for q, cap in VARIANTS}},
    }
    res["sanity_gap"] = abs(res["figaro"]["replica"]["pass"] - row["dada2_pass"])
    return res


def verdicts(results: dict) -> list:
    """P13 to P16 of TEST_REGISTER.md section 16, from the results of the datasets that were run."""
    have = [k for k in ROWS if k in results]
    n = len(have)
    m13 = sum(results[k]["figaro"]["match"]["mean"] for k in have)
    m16 = sum(results[k]["figaro"]["match"]["median"] for k in have)
    poor, clean = [k for k in POOR if k in results], [k for k in CLEAN if k in results]
    m14 = sum(results[k]["amplitrim"]["match"]["mean_floored"] for k in poor)
    m15 = sum(ee_pair(results[k]["amplitrim"]["rules"]["mean_floored"]) == (2.0, 2.0) for k in clean)
    out = [("P13", f"The mean rule (unfloored) is within 1 of FIGARO's maxEE on each read, at FIGARO's cut, on at least 8 of the 10 datasets: {m13} of {n} run",
            m13 >= 8 if n == 10 else None),
           ("P14", f"On the 4 poor datasets the floored mean rule at AmpliTrim's cut is within 1 of FIGARO's pick on each read on at least 3: {m14} of {len(poor)} run",
            m14 >= 3 if len(poor) == 4 else None),
           ("P15", f"On the 6 clean datasets the floored mean rule at AmpliTrim's cut gives (2,2) on at least 5: {m15} of {len(clean)} run",
            m15 >= 5 if len(clean) == 6 else None),
           ("P16", f"The median rule matches FIGARO on no more datasets than the mean rule: median {m16}, mean {m13} (of {n} run)", m16 <= m13 if n == 10 else None)]
    if all("q" in results[k]["figaro"] for k in have):
        best = best_q(results)
        m17 = best["matches"]
        out.append(("P17", f"For at least one q of {', '.join(f'{q:g}' for q in QGRID)}, the unfloored percentile rule at FIGARO's cut is within 1 of FIGARO's maxEE on each read on at least 9 of the 10 "
                           f"datasets: the best is q = {best['q']:g} with {m17} of {n} run", m17 >= 9 if n == 10 else None))
        q = best["q"]
        clean22 = sum(ee_pair(results[k]["amplitrim"]["q"][q]["ee"]) == (2.0, 2.0) for k in clean)
        reach = sum(results[k]["amplitrim"]["q"][q]["pass"] >= 0.70 - 1e-9 for k in poor)
        ok = (clean22 >= 5 and reach >= 3) if len(clean) == 6 and len(poor) == 4 else None
        out.append(("P18", f"With that q, at AmpliTrim's cut, the floored rule gives (2,2) on at least 5 of the 6 clean datasets ({clean22} of {len(clean)} run) and reaches 70% of the pairs "
                           f"passing on at least 3 of the 4 poor ones ({reach} of {len(poor)} run)", ok))
    return out


def best_q(results: dict) -> dict:
    """The q of QGRID whose unfloored rule matches FIGARO on most datasets (within 1 on each read, at FIGARO's cut); ties go to the smaller q (the cleaner reads)."""
    rows = [{"q": q, "matches": sum(results[k]["figaro"]["q"][q]["match"] for k in results)} for q in QGRID]
    return max(rows, key=lambda r: (r["matches"], -r["q"]))


def q_window(results: dict) -> dict:
    """EXPLORATORY. A percentile rule gives FIGARO's pick on a read when q lies in (share of reads below pick - 1, share of reads below pick]; the intersection over every read run."""
    lo, hi, n = 0.0, 1.0, 0
    for r in results.values():
        for a, b in r["figaro"]["cdf"]:
            lo, hi, n = max(lo, a), min(hi, b), n + 1
    return {"lo": lo, "hi": hi, "n": n, "feasible": lo < hi}


def implied_windows(results: dict, stat: str) -> dict:
    """EXPLORATORY. FIGARO's maxEE is int(x) + 1, so its curve value x lies in [pick - 1, pick). If x were k times a statistic s of the reads at FIGARO's cut,
    k lies in [(pick - 1) / s, pick / s) for each read of each dataset; the k that would reproduce every pick is the intersection of those windows. `exact` counts
    the reads that the statistic alone (k = 1) reproduces."""
    lo, hi, n, exact = 0.0, float("inf"), 0, 0
    for r in results.values():
        for pick, s in zip(r["figaro"]["pick"], r["figaro"]["stats"][stat]):
            if s <= 0:
                continue
            lo, hi, n = max(lo, (pick - 1) / s), min(hi, pick / s), n + 1
            exact += figaro_pad(s) == int(pick)
    return {"lo": lo, "hi": hi, "n": n, "exact": exact, "feasible": lo < hi}


def exploratory_lines(results: dict) -> list:
    L = ["## Exploratory (not pre-registered; added after a first, partial run, and not part of the criteria)", "",
         "FIGARO's maxEE is `int(x) + 1`, so its curve value x at the cut lay in [pick - 1, pick). If x were k times a statistic of the reads at FIGARO's cut, which k would "
         "reproduce every pick (every read of every dataset run)? A k fitted here is fitted on the datasets that test it, so it would be post hoc and adds a parameter.", "",
         "| statistic of the reads at FIGARO's cut | k reproducing every pick | reads reproduced at k = 1 |", "|---|---|---|"]
    for name in STATS:
        w = implied_windows(results, name)
        k = f"{w['lo']:.2f} to {w['hi']:.2f}" if w["feasible"] else f"none (the windows do not overlap: {w['lo']:.2f} and above, below {w['hi']:.2f})"
        L.append(f"| {name} | {k} | {w['exact']} of {w['n']} |")
    return L + [""]


def variants_lines(results: dict) -> list:
    """EXPLORATORY: what each q does with and without a cap at AmpliTrim's cut, on every dataset where the rule leaves the default."""
    if not all("variants" in r["amplitrim"] for r in results.values()):
        return []
    L = ["## Variants with a cap (exploratory; added after the percentile run, at the researcher's request)", "",
         "At AmpliTrim's cut, the floored rule for q = 0.75 and 0.80, with and without a cap of 5 on each read (DADA2's tutorial gives c(2, 5) as an example for the reverse read). "
         "Only datasets where some variant leaves (2,2) are listed; the pairs passing and the expected errors of the kept reads are from the exact quality table.", "",
         "| dataset | default (2,2) | q 0.75 | q 0.75, cap 5 | q 0.80 | q 0.80, cap 5 | FIGARO's pick |", "|---|---|---|---|---|---|---|"]
    for k, r in results.items():
        a = r["amplitrim"]
        if all(ee_pair(v["ee"]) == (2.0, 2.0) for v in a["variants"].values()):
            continue
        cells = [f"{ee_text(v['ee'])}: {v['pass']:.1%}" for v in a["variants"].values()]
        d, f = a["points"]["default (2,2)"], a["points"]["FIGARO's pick"]
        L.append(f"| {k} | {d['pass']:.1%} | " + " | ".join(cells) + f" | {ee_text(f['ee'])}: {f['pass']:.1%} |")
    L += ["", "Expected errors per kept read (forward / reverse), same order:", "", "| dataset | q 0.75 | q 0.75, cap 5 | q 0.80 | q 0.80, cap 5 | FIGARO's pick |", "|---|---|---|---|---|---|"]
    for k, r in results.items():
        a = r["amplitrim"]
        if all(ee_pair(v["ee"]) == (2.0, 2.0) for v in a["variants"].values()):
            continue
        f = a["points"]["FIGARO's pick"]
        L.append(f"| {k} | " + " | ".join(f"{v['ee_f']:.2f} / {v['ee_r']:.2f}" for v in a["variants"].values()) + f" | {f['ee_f']:.2f} / {f['ee_r']:.2f} |")
    return L + [""]


def percentile_lines(results: dict) -> list:
    if not all("q" in r["figaro"] for r in results.values()):
        return []
    w = q_window(results)
    L = ["## The percentile rule (P17 and P18, section 17 of the register)", "",
         "Per read: `int(q-th percentile of the expected errors at the cut) + 1`, and the same never below 2. FIGARO's pick is reproduced on a read when q lies in "
         "(share of reads with fewer expected errors than pick - 1, share with fewer than pick]; the intersection over every read run is "
         + (f"**q from {w['lo']:.2f} to {w['hi']:.2f}**" if w["feasible"] else f"**empty** (above {w['lo']:.2f} and at most {w['hi']:.2f}): no single q reproduces every pick") + f" ({w['n']} reads).", "",
         "| q | datasets within 1 of FIGARO at FIGARO's cut (unfloored) | clean datasets at (2,2), AmpliTrim's cut | poor datasets reaching 70% passing, AmpliTrim's cut |", "|---|---|---|---|"]
    poor, clean = [k for k in POOR if k in results], [k for k in CLEAN if k in results]
    for q in QGRID:
        m = sum(results[k]["figaro"]["q"][q]["match"] for k in results)
        c = sum(ee_pair(results[k]["amplitrim"]["q"][q]["ee"]) == (2.0, 2.0) for k in clean)
        r = sum(results[k]["amplitrim"]["q"][q]["pass"] >= 0.70 - 1e-9 for k in poor)
        L.append(f"| {q:g} | {m} of {len(results)} | {c} of {len(clean)} | {r} of {len(poor)} |")
    return L + [""]


def ee_text(ee) -> str:
    return f"{ee[0]:g},{ee[1]:g}"


def dataset_lines(name: str, r: dict) -> list:
    f, a = r["figaro"], r["amplitrim"]
    L = [f"### {name}  ({'poor' if name in POOR else 'clean'}: FIGARO chose maxEE {ee_text(f['pick'])})", ""]
    gap = r["sanity_gap"]
    L += [f"Check on the reads: at FIGARO's cut {f['cut'][0]}/{f['cut'][1]} and maxEE {ee_text(f['pick'])}, AmpliTrim's filter passes {f['replica']['pass']:.1%} of the pairs and "
          f"DADA2 passed {f['dada2_pass']:.1%} in night 3" + (f" (**{gap * 100:.1f} points apart: are these the reads of the night 3 run?**)" if gap > SANITY_POINTS else ""), ""]
    L += ["| at | read-pair statistic | forward / reverse | mean rule | mean rule, floored | median rule | median rule, floored | FIGARO |", "|---|---|---|---|---|---|---|---|"]
    for label, block, cut in (("FIGARO's cut", f, f["cut"]), ("AmpliTrim's cut", a, a["cut"])):
        rl = block["rules"]
        L.append(f"| {label} {cut[0]}/{cut[1]} | mean expected errors | {block['mean_ee'][0]:.2f} / {block['mean_ee'][1]:.2f} | {ee_text(rl['mean'])} | {ee_text(rl['mean_floored'])} | "
                 f"{ee_text(rl['median'])} | {ee_text(rl['median_floored'])} | {ee_text(f['pick'])} |")
        L.append(f"| | median expected errors | {block['median_ee'][0]:.2f} / {block['median_ee'][1]:.2f} | | | | | |")
    st = f["stats"]
    L += ["", f"FIGARO's pick means its curve value was in [{f['pick'][0] - 1:g}, {f['pick'][0]:g}) forward and [{f['pick'][1] - 1:g}, {f['pick'][1]:g}) reverse. Expected errors of the reads at its cut "
          f"(forward / reverse): " + "; ".join(f"{name} {v[0]:.2f} / {v[1]:.2f}" for name, v in st.items()) + "."]
    L += ["", f"At AmpliTrim's cut ({a['cut'][0]}/{a['cut'][1]}), what each maxEE keeps (exact quality table):", "", "| maxEE | pairs passing | expected errors per kept read (F / R) |", "|---|---|---|"]
    for label, x in a["points"].items():
        L.append(f"| {label}: {ee_text(x['ee'])} | {x['pass']:.1%} | {x['ee_f']:.2f} / {x['ee_r']:.2f} |")
    for q, x in a.get("q", {}).items():
        L.append(f"| percentile rule q = {q:g}, floored: {ee_text(x['ee'])} (unfloored {ee_text(x['raw'])}; at FIGARO's cut {ee_text(f['q'][q]['raw'])}) | {x['pass']:.1%} | {x['ee_f']:.2f} / {x['ee_r']:.2f} |")
    return L + [""]


def load_config(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def reads_folder(cfg: dict, key: str):
    """The reads folder of a dataset, found the way night2.py finds it: mockA_high, v1v2 and its are top-level keys; every other dataset is under "extras"."""
    if cfg.get(key):
        return Path(os.path.expandvars(cfg[key]))
    extra = (cfg.get("extras") or {}).get(key) or {}
    return Path(os.path.expandvars(extra["ready"])) if extra.get("ready") else None


def write_tsv(path: Path, results: dict) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["dataset", "kind", "figaro_maxEE", "at", "cut_F", "cut_R", "mean_ee_F", "mean_ee_R", "median_ee_F", "median_ee_R", "mean", "mean_floored", "median", "median_floored",
                    "mean_matches_figaro", "mean_floored_matches_figaro", "sanity_gap_points"])
        for k, r in results.items():
            for at, block in (("figaro", r["figaro"]), ("amplitrim", r["amplitrim"])):
                rl = block["rules"]
                w.writerow([k, "poor" if k in POOR else "clean", ee_text(r["figaro"]["pick"]), at, *block["cut"], *(f"{v:.3f}" for v in (*block["mean_ee"], *block["median_ee"])),
                            ee_text(rl["mean"]), ee_text(rl["mean_floored"]), ee_text(rl["median"]), ee_text(rl["median_floored"]),
                            block["match"]["mean"], block["match"]["mean_floored"], f"{r['sanity_gap'] * 100:.2f}"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, help="amplitrim_paths.json (default: the one next to the package)")
    ap.add_argument("--only", nargs="+", help="run only these datasets (keys of amplitrim_paths.json)")
    ap.add_argument("--out", type=Path, default=Path("maxee_rules_out"), help="where to write the report and the table")
    ap.add_argument("--max-pairs", type=int, default=MAX_PAIRS, help="pairs pooled over all samples (default 60,000, as in the frontier check)")
    args = ap.parse_args(argv)
    cfg = load_config(args.config or ROOT / "amplitrim_paths.json")
    results, skipped = {}, []
    for key in ROWS:
        if args.only and key not in args.only:
            continue
        folder = reads_folder(cfg, key)
        if folder is None or not folder.is_dir():
            skipped.append(key)
            print(f"{key}: no reads folder in amplitrim_paths.json (top level or extras), skipped")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            info = build_subsample(folder, Path(tmp), max_pairs=args.max_pairs)
            pairs = ReadPairs.load(info.r1, info.r2)
            print(f"{key}: {pairs.n:,} pairs pooled", flush=True)
            try:
                results[key] = analyse(pairs, ROWS[key])
            except ValueError as err:
                skipped.append(key)
                print(f"{key}: SKIPPED, {err}")
    if not results:
        print("nothing was run: check amplitrim_paths.json")
        return 2
    args.out.mkdir(parents=True, exist_ok=True)
    L = ["# Night 4 (offline): a simple maxEE rule against FIGARO's", "", "The criteria are in TEST_REGISTER.md, section 16, and were fixed before this ran.", ""]
    L += ["## Verdicts", ""] + [f"- **{p}** {'MET' if ok else 'NOT MET' if ok is not None else 'not decided (not all datasets were run)'}: {txt}" for p, txt, ok in verdicts(results)] + [""]
    if skipped:
        L += [f"Skipped (no reads folder): {', '.join(skipped)}", ""]
    L += percentile_lines(results)
    L += variants_lines(results)
    L += exploratory_lines(results)
    L += ["## Datasets", ""]
    for k, r in results.items():
        L += dataset_lines(k, r)
    (args.out / "MAXEE_RULES_REPORT.md").write_text("\n".join(L), encoding="utf-8")
    write_tsv(args.out / "maxee_rules.tsv", results)
    print("\n".join(L[:12]))
    print(f"\nwritten: {args.out / 'MAXEE_RULES_REPORT.md'}\n         {args.out / 'maxee_rules.tsv'}\n\nNext: send MAXEE_RULES_REPORT.md.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
