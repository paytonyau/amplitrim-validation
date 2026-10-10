"""
Candidate evaluation and the selection rule.

A candidate is a pair (trim_f, trim_r); its total length is trim_f + trim_r. Three things matter:

  coverage  Can the settings merge the amplicon lengths actually present? The longest insert they
            can merge is total length - min overlap, and at most `max_unmergeable_frac` of the
            observed inserts may be longer (otherwise long amplicons are systematically lost).
  yield     Fraction of read pairs that survive both the quality filter and merging.
  overlap   Mean overlap used. Less overlap means less reliance on the noisy ends of the reads.

Selection: among candidates that cover the amplicon range, take those whose yield is within
`yield_tolerance` of the best, and choose the one with the smallest overlap. Overlaps within 2 nt count as equal; ties go to the
highest yield, then the most balanced forward/reverse split, then the fewest observed mismatches per overlap base. No candidate covering the range with a usable yield is a result in
itself, reported as such.
"""

import math
from bisect import bisect_right
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import pandas as pd

from .config import ee_pair

NAN = float("nan")


@dataclass
class Candidate:
    trim_f: int
    trim_r: int
    ee_pass_rate: float            # share of pairs passing the quality filter
    merge_rate_strict: float       # raw reads, maxMismatch as in DADA2
    merge_rate_tolerant: float     # raw reads, vsearch default mismatch allowance
    mean_overlap: float
    mean_overlap_diffs: float      # observed mismatches per merged pair
    max_mergeable_insert: int
    covers_longest: bool
    yield_est: float               # ee_pass_rate * merge_rate_tolerant
    is_recommended: bool = False
    source: str = ""                # which method proposed the split
    balance: float = 0.0           # |trim_f / longest forward read - trim_r / longest reverse read|

    @property
    def total_length(self) -> int:
        return self.trim_f + self.trim_r

    @property
    def tag(self) -> str:
        return f"f{self.trim_f}_r{self.trim_r}"

    @property
    def diff_density(self) -> float:
        return self.mean_overlap_diffs / self.mean_overlap if self.mean_overlap else NAN


def frac_unmergeable(sorted_inserts: Sequence[int], max_mergeable: int) -> float:
    """Share of observed insert sizes longer than the longest mergeable insert."""
    if not sorted_inserts:
        return NAN
    return (len(sorted_inserts) - bisect_right(sorted_inserts, max_mergeable)) / len(sorted_inserts)


def coverage_floor(sorted_inserts: Sequence[int], min_overlap: int,
                   max_unmergeable_frac: float) -> Tuple[int, int]:
    """(shortest total length that covers the amplicon range, the longest amplicon it must merge)."""
    n = len(sorted_inserts)
    allowed_longer = int(max_unmergeable_frac * n + 1e-9)
    longest = sorted_inserts[n - allowed_longer - 1]
    return longest + min_overlap, longest


def _eligible(c: Candidate) -> bool:
    return (c.covers_longest and c.yield_est > 0
            and math.isfinite(c.mean_overlap) and math.isfinite(c.diff_density))


def select_recommendation(candidates: List[Candidate],
                          yield_tolerance: float = 0.01) -> Optional[Candidate]:
    for c in candidates:
        c.is_recommended = False
    eligible = [c for c in candidates if _eligible(c)]
    if not eligible:
        return None
    best = max(c.yield_est for c in eligible)
    near = [c for c in eligible if c.yield_est >= best - yield_tolerance]
    pick = min(near, key=lambda c: (int(c.mean_overlap // 2), -round(c.yield_est, 3), c.balance, c.diff_density))
    pick.is_recommended = True
    return pick


def to_dataframe(candidates: List[Candidate]) -> pd.DataFrame:
    """The table written to candidates.csv, one row per setting tested."""
    return pd.DataFrame([{
        "truncLenF": c.trim_f, "truncLenR": c.trim_r, "total_length": c.total_length,
        "estimated_yield": round(c.yield_est, 4), "pass_filter": round(c.ee_pass_rate, 4),
        "merge_rate": round(c.merge_rate_tolerant, 4),
        "merge_rate_strict": round(c.merge_rate_strict, 4),
        "mean_overlap": round(c.mean_overlap, 1) if math.isfinite(c.mean_overlap) else None,
        "longest_mergeable_insert": c.max_mergeable_insert,
        "covers_longest_amplicon": c.covers_longest, "suggested_by": c.source,
        "recommended": c.is_recommended,
    } for c in candidates])


MIN_RETENTION = 0.5     # the default floor for the quality-first maxEE suggestion; DADA2's author says even 40% of the reads is fine


def dominators(rows: list) -> list:
    """For each row, the rows that dominate it: they keep at least as many pairs AND at least as few expected errors per kept pair
    (forward plus reverse), and are strictly better on one of the two. Rows without data have none. There is no weight to choose: a
    row is dominated only when another is better on both counts at once."""
    def point(x):
        return (x["pass_rate"], x["ee_f"] + x["ee_r"]) if x.get("pass_rate") is not None and x.get("ee_f") is not None else None
    pts = [point(x) for x in rows]
    out = []
    for i, p in enumerate(pts):
        better = []
        if p is not None:
            for j, q in enumerate(pts):
                if j != i and q is not None and q[0] >= p[0] - 1e-12 and q[1] <= p[1] + 1e-12 and (q[0] > p[0] + 1e-12 or q[1] < p[1] - 1e-12):
                    better.append(rows[j])
        out.append(better)
    return out


def quality_first_suggestion(rows: list, min_retention: float = MIN_RETENTION):
    """Quality first. Among the rows that are not stricter than DADA2's default of 2 and are not beaten by another row, the one with the
    fewest expected errors per kept pair that still keeps at least `min_retention` of the pairs (estimated yield). If no row reaches
    it, the row that keeps the most. Returns (row, reached). It weighs retention only against read noise; it knows nothing about
    accuracy, and the default maxEE is never changed by it."""
    usable = [x for x, beaten in zip(rows, dominators(rows)) if not beaten and x.get("yield_est") is not None and min(ee_pair(x["ee"])) >= 2]
    if not usable:
        return None, False
    reaching = [x for x in usable if x["yield_est"] >= min_retention]
    if reaching:
        return min(reaching, key=lambda x: (x["ee_f"] + x["ee_r"], sum(ee_pair(x["ee"])))), True
    return max(usable, key=lambda x: x["yield_est"]), False
