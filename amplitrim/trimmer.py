"""
Pair-safe truncation and DADA2-style read filtering.

For a given (trim_f, trim_r), a read pair is kept only if, after truncation:
  - neither read is shorter than its truncation length,
  - no base has Q <= 2 (truncQ = 2),
  - neither read contains an N (maxN = 0),
  - both reads have at most `max_ee` expected errors.
Pairs are filtered jointly so R1 and R2 stay in sync for the merge test.

Every rule looks only at the first `trim` bases and can only get stricter as the truncation grows.
So each read has one maximum truncation length at which it still passes, and a pair passes at
(trim_f, trim_r) exactly when trim_f <= its forward maximum and trim_r <= its reverse maximum.
These maxima are computed once; after that the filter at any truncation, and the share of pairs
that pass at EVERY (trim_f, trim_r), are cheap.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from .subsample import read_fastq
from .config import ee_pair

_ERR = 10.0 ** (-np.arange(0, 94) / 10.0)      # Phred -> error probability
_CHUNK = 8192                                   # rows per block when accumulating expected errors


@dataclass
class TrimStats:
    n_in: int = 0
    n_too_short: int = 0
    n_low_q: int = 0
    n_has_n: int = 0
    n_ee_fail: int = 0
    n_out: int = 0

    @property
    def ee_pass_rate(self) -> float:
        """Fraction of input pairs that survive the filter."""
        return self.n_out / self.n_in if self.n_in else 0.0


def _phred(quals: List[str], width: int) -> np.ndarray:
    buf = "".join(q.ljust(width, "!") for q in quals).encode("ascii")
    return np.frombuffer(buf, dtype=np.uint8).reshape(len(quals), width) - 33


def _first_n(seqs: List[str], width: int) -> np.ndarray:
    """Index of the first 'N' in each read (= width when there is none)."""
    buf = "".join(s.ljust(width, "A") for s in seqs).encode("ascii")
    is_n = np.frombuffer(buf, dtype=np.uint8).reshape(len(seqs), width) == ord("N")
    return np.where(is_n.any(axis=1), is_n.argmax(axis=1), width)


def _first_low_quality(phred: np.ndarray) -> np.ndarray:
    low = phred <= 2
    return np.where(low.any(axis=1), low.argmax(axis=1), phred.shape[1])


def _ee_limit(phred: np.ndarray, max_ee: float) -> np.ndarray:
    """Longest truncation whose expected errors stay within max_ee."""
    out = np.empty(len(phred), dtype=np.int64)
    for i in range(0, len(phred), _CHUNK):
        out[i:i + _CHUNK] = (np.cumsum(_ERR[phred[i:i + _CHUNK]], axis=1) <= max_ee).sum(axis=1)
    return out


class PassTable:
    """Share of read pairs that pass the filter at every truncation pair (trim_f, trim_r)."""

    def __init__(self, max_f: np.ndarray, max_r: np.ndarray, width_f: int, width_r: int):
        counts = np.zeros((width_f + 2, width_r + 2))
        np.add.at(counts, (max_f, max_r), 1)
        # pairs with max_f >= trim_f and max_r >= trim_r
        self._rate = counts[::-1, ::-1].cumsum(axis=0).cumsum(axis=1)[::-1, ::-1] / max(len(max_f), 1)
        self.width_f, self.width_r = width_f, width_r

    def rate(self, trim_f: int, trim_r: int) -> float:
        if not (0 < trim_f <= self.width_f and 0 < trim_r <= self.width_r):
            return 0.0
        return float(self._rate[trim_f, trim_r])

    def best_split(self, total: int, prefer_f: Optional[int] = None) -> Optional[Tuple[int, int, float]]:
        """The (trim_f, trim_r) adding up to `total` with the highest pass rate. Ties go to the split
        closest to `prefer_f` (or to the longest forward read when not given)."""
        lo, hi = max(1, total - self.width_r), min(self.width_f, total - 1)
        if lo > hi:
            return None
        fs = np.arange(lo, hi + 1)
        rates = self._rate[fs, total - fs]
        near_best = fs[rates >= rates.max() - 1e-12]
        pick = int(near_best[np.argmin(np.abs(near_best - prefer_f))] if prefer_f is not None
                   else near_best.max())
        return pick, total - pick, float(self._rate[pick, total - pick])


    def near_best_splits(self, total: int, within: float = 0.005, k: int = 5) -> List[Tuple[int, int]]:
        """Up to k splits of `total`, spread evenly over all those whose pass rate is within `within`
        of the best. When pass rate cannot tell splits apart (all equal), this spans the whole range."""
        lo, hi = max(1, total - self.width_r), min(self.width_f, total - 1)
        if lo > hi:
            return []
        fs = np.arange(lo, hi + 1)
        rates = self._rate[fs, total - fs]
        ok = fs[rates >= rates.max() - within]
        if len(ok) > k:
            ok = ok[np.unique(np.linspace(0, len(ok) - 1, k).round().astype(int))]
        return [(int(f), total - int(f)) for f in ok]


def _cut_errors(phred: np.ndarray, lengths: np.ndarray, rows: np.ndarray, chunk: int = 4000) -> np.ndarray:
    """Expected errors of each selected read counted over its first `lengths` bases (in chunks, to bound the memory)."""
    idx = np.flatnonzero(rows)
    out = np.empty(len(idx))
    cols = np.arange(phred.shape[1])
    for a in range(0, len(idx), chunk):
        sel = idx[a:a + chunk]
        out[a:a + chunk] = np.where(cols[None, :] < lengths[sel][:, None], _ERR[phred[sel]], 0.0).sum(axis=1)
    return out


class ReadPairs:
    """A set of read pairs held in memory, filterable at any truncation lengths."""

    def __init__(self, recs1, recs2):
        self.recs1, self.recs2 = list(recs1), list(recs2)
        if len(self.recs1) != len(self.recs2):
            raise ValueError("R1 and R2 contain different numbers of reads")
        self.n = len(self.recs1)
        self.len1 = np.array([len(r[1]) for r in self.recs1], dtype=np.int64)
        self.len2 = np.array([len(r[1]) for r in self.recs2], dtype=np.int64)
        self.width_f, self.width_r = int(self.len1.max(initial=0)), int(self.len2.max(initial=0))
        self._phred1 = _phred([r[3] for r in self.recs1], self.width_f)
        self._phred2 = _phred([r[3] for r in self.recs2], self.width_r)
        self.low1, self.low2 = _first_low_quality(self._phred1), _first_low_quality(self._phred2)
        self.n1 = _first_n([r[1] for r in self.recs1], self.width_f)
        self.n2 = _first_n([r[1] for r in self.recs2], self.width_r)
        self._ee = {}

    @classmethod
    def load(cls, r1: Path, r2: Path) -> "ReadPairs":
        return cls(read_fastq(r1), read_fastq(r2))

    def _ee_limits(self, max_ee):
        """Per-read truncation limits for the forward and the reverse read; `max_ee` is one number or (forward, reverse)."""
        key = ee_pair(max_ee)
        if key not in self._ee:
            self._ee[key] = (_ee_limit(self._phred1, key[0]), _ee_limit(self._phred2, key[1]))
        return self._ee[key]

    def pass_table(self, max_ee: float = 2.0) -> PassTable:
        e1, e2 = self._ee_limits(max_ee)
        max_f = np.minimum.reduce([self.len1, self.low1, self.n1, e1])
        max_r = np.minimum.reduce([self.len2, self.low2, self.n2, e2])
        return PassTable(max_f, max_r, self.width_f, self.width_r)

    def untruncated(self, max_ee: float, min_len: int):
        """Outcome of DADA2's filterAndTrim with no truncLen, per pair. Each read is cut at its first base
        with Q <= 2 (truncQ = 2); the pair is dropped if either cut read is shorter than `min_len`, contains
        an N, or has more than `max_ee` expected errors. Returns the cut lengths and the failure masks."""
        e1, e2 = self._ee_limits(max_ee)
        l1, l2 = np.minimum(self.len1, self.low1), np.minimum(self.len2, self.low2)
        short = (l1 < min_len) | (l2 < min_len)
        has_n = ~short & ((self.n1 < l1) | (self.n2 < l2))
        ee_fail = ~short & ~has_n & ((e1 < l1) | (e2 < l2))
        return l1, l2, short, has_n, ee_fail

    def untruncated_load(self, max_ee, min_len: int, baseline=None) -> dict:
        """What DADA2's no-truncLen filter keeps at max_ee and how noisy and how long it is: the mean expected errors and the mean cut length
        (forward, reverse) of the kept pairs, and, with a baseline limit, the mean expected errors of the pairs kept at max_ee but not at the
        baseline (nan when there are none). Each read is cut at its first base with Q <= 2, and the errors are counted over the cut read."""
        l1, l2, short, has_n, ee_fail = self.untruncated(max_ee, min_len)
        kept = ~(short | has_n | ee_fail)
        added = None
        if baseline is not None:
            _, _, s0, n0, e0 = self.untruncated(baseline, min_len)
            added = kept & (s0 | n0 | e0)                             # kept now, but not at the baseline

        def mean(values):
            return float(values.mean()) if len(values) else float("nan")
        out = {"n_kept": int(kept.sum()), "ee_f": mean(_cut_errors(self._phred1, l1, kept)), "ee_r": mean(_cut_errors(self._phred2, l2, kept)),
               "len_f": mean(l1[kept]), "len_r": mean(l2[kept]), "add_f": float("nan"), "add_r": float("nan")}
        if added is not None:
            out["add_f"], out["add_r"] = mean(_cut_errors(self._phred1, l1, added)), mean(_cut_errors(self._phred2, l2, added))
        return out

    def filter_untruncated(self, out1: Path, out2: Path, max_ee: float = 2.0, min_len: int = 50) -> TrimStats:
        """Write the pairs that survive the no-truncLen filter, each read cut at its first Q <= 2 base."""
        l1, l2, short, has_n, ee_fail = self.untruncated(max_ee, min_len)
        idx = np.flatnonzero(~(short | has_n | ee_fail))
        st = TrimStats(n_in=self.n, n_too_short=int(short.sum()), n_has_n=int(has_n.sum()),
                       n_ee_fail=int(ee_fail.sum()), n_out=len(idx))
        out1.parent.mkdir(parents=True, exist_ok=True)
        with open(out1, "w", newline="\n") as o1, open(out2, "w", newline="\n") as o2:
            o1.write("".join(f"{h}\n{s[:l1[i]]}\n{p}\n{q[:l1[i]]}\n" for i in idx for h, s, p, q in (self.recs1[i],)))
            o2.write("".join(f"{h}\n{s[:l2[i]]}\n{p}\n{q[:l2[i]]}\n" for i in idx for h, s, p, q in (self.recs2[i],)))
        return st

    def _masks(self, trim_f: int, trim_r: int, max_ee):
        """Why each pair is dropped at this truncation (too short, a base at or below Q2, an N, too many expected errors) and which are kept."""
        e1, e2 = self._ee_limits(max_ee)
        short = (self.len1 < trim_f) | (self.len2 < trim_r)
        low_q = ~short & ((self.low1 < trim_f) | (self.low2 < trim_r))
        has_n = ~short & ~low_q & ((self.n1 < trim_f) | (self.n2 < trim_r))
        ee_fail = ~short & ~low_q & ~has_n & ((e1 < trim_f) | (e2 < trim_r))
        return short, low_q, has_n, ee_fail, ~(short | low_q | has_n | ee_fail)

    def error_profile(self, n: int = 4000) -> dict:
        """How expected errors (sum of 10^(-Q/10)) accumulate along each read: the 10th, 50th and 90th percentile over reads at every
        position, for the forward ('f') and reverse ('r') reads, from an even sample of at most n reads. Only reads long enough to
        reach a position count at it; a position reached by fewer than 30 reads is nan."""
        idx = np.linspace(0, self.n - 1, min(n, self.n)).astype(int)
        out = {}
        for key, phred, lens, width in (("f", self._phred1, self.len1, self.width_f), ("r", self._phred2, self.len2, self.width_r)):
            cum, ln = np.cumsum(_ERR[phred[idx]], axis=1), lens[idx]
            curves = np.full((3, width), np.nan)
            for j in range(width):
                v = cum[ln > j, j]
                if len(v) >= 30:
                    curves[:, j] = np.percentile(v, (10, 50, 90))
            out[key] = curves
        return out

    def added_error_load(self, trim_f: int, trim_r: int, max_ee, baseline=(2.0, 2.0)):
        """Mean expected errors, forward and reverse, of the pairs that a looser maxEE keeps beyond what the baseline limit keeps: how
        noisy the extra reads are. (nan, nan) when it adds none, for example when max_ee is not looser than the baseline."""
        keep, base = self._masks(trim_f, trim_r, max_ee)[4], self._masks(trim_f, trim_r, baseline)[4]
        added = keep & ~base
        if not added.any():
            return float("nan"), float("nan")
        return (float(_ERR[self._phred1[added, :trim_f]].sum(axis=1).mean()), float(_ERR[self._phred2[added, :trim_r]].sum(axis=1).mean()))

    def error_load(self, trim_f: int, trim_r: int, max_ee=2.0):
        """Mean expected errors (sum of 10^(-Q/10)) of the forward and the reverse read among the pairs that pass the filter. A looser
        maxEE keeps noisier reads, and this is how much noisier."""
        keep = self._masks(trim_f, trim_r, max_ee)[4]
        if not keep.any():
            return float("nan"), float("nan")
        return (float(_ERR[self._phred1[keep, :trim_f]].sum(axis=1).mean()), float(_ERR[self._phred2[keep, :trim_r]].sum(axis=1).mean()))

    def filter(self, trim_f: int, trim_r: int, out1: Path, out2: Path,
               max_ee: float = 2.0) -> TrimStats:
        """Truncate, filter and write the surviving pairs (plain FASTQ) to out1/out2."""
        short, low_q, has_n, ee_fail, keep = self._masks(trim_f, trim_r, max_ee)
        idx = np.flatnonzero(keep)
        st = TrimStats(n_in=self.n, n_too_short=int(short.sum()), n_low_q=int(low_q.sum()),
                       n_has_n=int(has_n.sum()), n_ee_fail=int(ee_fail.sum()), n_out=len(idx))
        out1.parent.mkdir(parents=True, exist_ok=True)
        with open(out1, "w", newline="\n") as o1, open(out2, "w", newline="\n") as o2:
            o1.write("".join(f"{h}\n{s[:trim_f]}\n{p}\n{q[:trim_f]}\n"
                             for h, s, p, q in (self.recs1[i] for i in idx)))
            o2.write("".join(f"{h}\n{s[:trim_r]}\n{p}\n{q[:trim_r]}\n"
                             for h, s, p, q in (self.recs2[i] for i in idx)))
        return st
