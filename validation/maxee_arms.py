#!/usr/bin/env python
"""The arms and the measures of the maxEE head-to-head (night 2). Nothing here runs DADA2; night2.py calls compare_methods.py for that.

The question: when AmpliTrim chooses the maxEE as well as the truncation, does it end up where FIGARO ends up, with as many merged pairs and no noisier reads?

Arms (every arm is a truncation and a maxEE, run through DADA2 on the same reads with the same 12 nt overlap):
  A   AmpliTrim, maxEE 2 (its default)
  B   AmpliTrim at the maxEE of its own quality-first pick (--explore-max-ee), floor 50% of the pairs (the tool's default; the headline)
  B2  the same with a floor of 80%
  B4  AmpliTrim at the maxEE of 'the mean expected errors of the reads at the cut, plus one', rounded to whole numbers, at least 1 (the rule FIGARO
      is documented to use; the earlier check reproduced FIGARO's maxEE within 1.0 on 4 of 6 datasets)
  C   FIGARO's own pick (truncation and maxEE), given the longest amplicon AmpliTrim measured, and also the typical one when a tail makes them differ
  D   AmpliTrim's truncation at FIGARO's maxEE (a control: only the truncation differs)
  E   no truncation (maxEE 2, minimum length 50)
  F   a practitioner's cut: each read where its median quality falls below 30 (maxEE 2)
Night 3 adds (TEST_REGISTER.md, section 14):
  G   AmpliTrim's best truncation at maxEE (2,5), DADA2's own suggestion for the reverse read, on every dataset
  R   a rule, no extra run: G where arm A passes under 70% of the pairs, otherwise A (what the planned warning would offer)
  F2  a careful manual cut (post hoc): the median-Q30 cut, never longer than the reads (10th percentile of the lengths), never too short for the longest amplicon
      plus 12 nt (the missing nucleotides split between the reads), maxEE (2,2) relaxed to (2,5) when fewer than 60% of the pairs pass at (2,2)
Measures: the share of pairs DADA2 merges (paired per-sample difference), the mean expected errors per kept read (forward, reverse), and on mocks
the variants and strains recovered and the share of reads in exact ASVs. The margins are fixed here, before any run.
"""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

MARGIN_MERGED = 2.0                                          # points of merged share inside which two picks count as the same
MARGIN_NOISE = 0.02                                          # expected errors per read inside which two picks count as equally noisy (the margin of the frontier check)
MARGIN_ACCURACY = {"variants": 1, "strains": 1, "identical_points": 0.1}      # the margins already fixed for the overnight battery
FLOORS = (0.5, 0.8)
TAIL_RATIO = 0.95                                            # a tail exists when the typical amplicon is under this share of the longest
STABILITY_NT = 10                                            # the margin already fixed for the stability of a truncation
ALT_EE = "2,5"                                               # arm G: DADA2's suggestion for the reverse read when too few reads pass
POOR_BELOW = 0.70                                            # rule R: use ALT_EE where arm A passes under this share of the pairs (chosen after night 2)
F2_RELAX_BELOW = 0.60                                        # arm F2: relax to ALT_EE when fewer than this share of the pairs pass at (2,2)
F2_CAP_PERCENTILE = 10                                       # arm F2: a cut is never longer than the length that 90% of the reads reach


def ee_str(pair) -> str:
    f, r = pair
    return f"{f:g}" if f == r else f"{f:g},{r:g}"


def round_ee(x: float) -> int:
    """To the nearest whole number, halves up (Python's round() sends 2.5 to 2 and 3.5 to 4), and never below 1."""
    return max(1, int(math.floor(x + 0.5)))


# ------------------------------------------------------------------------------------------------------------------ picks
def explore(reads: Path, work: Path, marker: str = "16s", min_overlap: int = 12, min_len: int = 50, run_fn=None) -> dict:
    """AmpliTrim once with --explore-max-ee: its default truncation, the table of maxEE rows, the quality-first pick at each floor, and the amplicon lengths."""
    from amplitrim.config import ee_pair
    from amplitrim.pipeline import Settings, run
    from amplitrim.ranking import quality_first_suggestion
    res = (run_fn or run)(Settings(input_dir=reads, work_dir=work, max_ee=2.0, min_overlap=min_overlap, marker=marker, min_len=min_len, explore_max_ee=True))
    rec = res.recommendation
    rows = [{"ee": ee_pair(x["ee"]), "yield": x.get("yield_est"), "pass": x.get("pass_rate"), "trim_f": x.get("trim_f"), "trim_r": x.get("trim_r")} for x in (res.ee_search or [])]
    picks = {}
    for fl in FLOORS:
        row, reached = quality_first_suggestion(res.ee_search, fl) if res.ee_search else (None, False)
        picks[fl] = None if row is None else {"ee": ee_pair(row["ee"]), "reached": bool(reached), "yield": row.get("yield_est")}
    inserts = list(res.inserts or [])
    return {"table": rows, "picks": picks, "trim": (rec.trim_f, rec.trim_r) if rec is not None else None, "longest": res.longest_insert,
            "typical": int(np.median(inserts)) if inserts else None}


def load_pairs(reads: Path, tmp: Path, max_pairs: int = 60000):
    """The pooled quality table of a reads folder (what the frontier check uses)."""
    from amplitrim.subsample import build_subsample
    from amplitrim.trimmer import ReadPairs
    info = build_subsample(reads, tmp, max_pairs=max_pairs)
    return ReadPairs.load(info.r1, info.r2)


def rule_ee(pairs, trim_f: int, trim_r: int) -> tuple:
    """B4: the mean expected errors of the reads that reach the cut, plus one, rounded, at least 1 (forward, reverse)."""
    from round4_steps import mean_ee_at_cut
    ef, er = mean_ee_at_cut(pairs, trim_f, trim_r)
    return (round_ee(ef + 1.0), round_ee(er + 1.0))


def q30_cut(pairs, threshold: float = 30.0, minimum: int = 50) -> tuple:
    """F: where the median quality of the reads that reach each position first falls below the threshold, for the forward and the reverse reads."""
    def one(phred, lengths):
        longest = int(lengths.max())
        for j in range(longest):
            col = phred[lengths > j, j]
            if col.size and float(np.median(col)) < threshold:
                return max(j, minimum)
        return longest
    return one(pairs._phred1, pairs.len1), one(pairs._phred2, pairs.len2)


def f2_plan(pairs, longest, min_overlap: int = 12) -> dict:
    """F2: a careful manual cut, as a practitioner who reads the tutorial would make it (section 14 of the register).

    1. the median-Q30 cut of each read (arm F); 2. never longer than the read: capped at the length that 90% of the reads reach;
    3. never too short for the amplicon: if forward + reverse is under the longest amplicon + min_overlap, the missing nucleotides are added, half to each
       read (the forward takes the odd one), each up to its cap, and what one read cannot take goes to the other;
    4. maxEE (2,2), relaxed to (2,5) when fewer than F2_RELAX_BELOW of the pairs pass at (2,2) at that cut (AmpliTrim's exact pass rate)."""
    tf0, tr0 = q30_cut(pairs)
    cap_f = int(np.percentile(pairs.len1, F2_CAP_PERCENTILE, method="lower"))
    cap_r = int(np.percentile(pairs.len2, F2_CAP_PERCENTILE, method="lower"))
    tf, tr = min(tf0, cap_f), min(tr0, cap_r)
    need = (int(longest) + min_overlap) if longest else 0
    add = max(0, need - tf - tr)
    give_f = min(cap_f - tf, (add + 1) // 2)
    give_r = min(cap_r - tr, add - give_f)
    give_f = min(cap_f - tf, add - give_r)
    tf, tr = tf + give_f, tr + give_r
    pass22 = float(pairs.pass_table(2.0).rate(tf, tr))
    return {"cut": (tf, tr), "q30": (tf0, tr0), "caps": (cap_f, cap_r), "need": need, "reaches": tf + tr >= need, "pass22": pass22,
            "ee": ALT_EE if pass22 < F2_RELAX_BELOW else "2"}


def has_tail(longest, typical) -> bool:
    return bool(longest and typical and typical < TAIL_RATIO * longest)


def arm_plan(info: dict, rule: tuple = None) -> dict:
    """The maxEE of each AmpliTrim arm (as the strings compare_methods.py takes), in a fixed order, with the duplicates folded together."""
    arms = {"A": "2"}
    for label, fl in (("B", 0.5), ("B2", 0.8)):
        pick = info["picks"].get(fl)
        arms[label] = ee_str(pick["ee"]) if pick else None
    arms["B4"] = ee_str(rule) if rule else None
    return arms


def unique_ees(arms: dict) -> list:
    out = []
    for v in arms.values():
        if v and v not in out:
            out.append(v)
    return out


# ------------------------------------------------------------------------------------------------------------------ measures
def paired_diff(a: dict, b: dict, reps: int = 2000, seed: int = 1) -> dict:
    """Mean of (a - b) over the samples both have, in points; a bootstrap 95% interval over samples when there are at least 10."""
    names = sorted(set(a) & set(b))
    d = np.array([100.0 * (a[n] - b[n]) for n in names])
    if not len(d):
        return {"n": 0, "mean": None, "lo": None, "hi": None}
    out = {"n": len(d), "mean": float(d.mean()), "lo": None, "hi": None}
    if len(d) >= 10:
        rng = np.random.default_rng(seed)
        means = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(reps)])
        out["lo"], out["hi"] = float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))
    return out


def retention_class(diff, margin: float = MARGIN_MERGED) -> str:
    if diff is None:
        return "unknown"
    return "more" if diff > margin else ("fewer" if diff < -margin else "same")


def noise_class(a, b, margin: float = MARGIN_NOISE) -> str:
    """Mean expected errors per kept read (forward, reverse) of a against b: cleaner, about the same, noisier, or mixed (one read each way)."""
    if a is None or b is None or any(x != x for x in (*a, *b)):
        return "unknown"
    d = [x - y for x, y in zip(a, b)]
    lower, higher = any(x < -margin for x in d), any(x > margin for x in d)
    if lower and higher:
        return "mixed"
    return "cleaner" if lower else ("noisier" if higher else "about the same")


def verdict(retention: str, noise: str) -> str:
    """One word per pair of picks. 'better' and 'worse' need a gain or a loss on one axis and no loss or gain on the other."""
    if "unknown" in (retention, noise):
        return "no result"
    if noise == "mixed":
        return f"mixed noise ({retention} pairs)"
    gain = (retention == "more") + (noise == "cleaner")
    loss = (retention == "fewer") + (noise == "noisier")
    if gain and loss:
        return "trade-off: more pairs, noisier" if retention == "more" else "trade-off: fewer pairs, cleaner"
    if gain:
        return "better"
    if loss:
        return "worse"
    return "matches"


def accuracy_verdict(a: dict, b: dict) -> str:
    """On a mock: a against b by variants recovered, strains recovered and the share of reads in exact ASVs, with the margins fixed for the battery."""
    def num(row, key):
        try:
            return float(row[key])
        except (KeyError, TypeError, ValueError):
            return None

    def strains(row):
        try:
            return float(str(row["strains"]).split("/")[0])
        except (KeyError, TypeError, ValueError):
            return None
    checks = [(num(a, "found"), num(b, "found"), MARGIN_ACCURACY["variants"]), (strains(a), strains(b), MARGIN_ACCURACY["strains"])]
    ia, ib = num(a, "identical"), num(b, "identical")
    if ia is not None and ib is not None:
        checks.append((100 * ia, 100 * ib, MARGIN_ACCURACY["identical_points"]))
    better = sum(1 for x, y, m in checks if x is not None and y is not None and x - y > m)
    worse = sum(1 for x, y, m in checks if x is not None and y is not None and y - x > m)
    if not any(x is not None and y is not None for x, y, _ in checks):
        return "no truth"
    return "mixed" if better and worse else ("more accurate" if better else ("less accurate" if worse else "as accurate"))
