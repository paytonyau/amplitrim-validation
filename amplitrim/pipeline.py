"""
The calibration pipeline.

  1. Take one shared, deterministic subsample of the reads.
  2. Measure the amplicon lengths (pairs that overlap at full read length).
  3. Work out the shortest total length (truncLenF + truncLenR) that can merge the longest amplicons.
     That bound is exact, so shorter totals are never tested.
  4. Test totals upward from that bound. For each total, several forward/reverse splits with a high
     quality-filter pass rate (and, optionally, FIGARO's) get a DADA2-style filter and a vsearch merge
     test. The search stops once yield no longer improves.
  5. Recommend the shortest-overlap setting whose yield is within a tolerance of the best.
"""

import os
import shutil
import sys
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable, List, Optional

from .config import Config, ee_pair
from .figaro_wrapper import candidates_for_total
from .merge_wrapper import estimate_insert_sizes, run_vsearch
from .ranking import Candidate, coverage_floor, frac_unmergeable, select_recommendation
from .subsample import DEFAULT_MAX_PAIRS, SubsampleInfo, build_subsample
from .trimmer import ReadPairs

MIN_PAIRS_FOR_ESTIMATE = 100


class PipelineError(Exception):
    """A problem the user can fix (missing tool, unusable input)."""


@dataclass
class Settings:
    input_dir: Path
    work_dir: Path
    min_overlap: int = 12
    max_ee: object = 2.0              # one number, or (forward, reverse)
    max_pairs: int = DEFAULT_MAX_PAIRS
    subsample_every: int = 1
    longest_insert: Optional[int] = None
    step: int = 2
    yield_tolerance: Optional[float] = None
    exhaustive: bool = False
    # optional FIGARO cross-check
    figaro_script: Optional[Path] = None
    figaro_python: str = sys.executable
    top_n: int = 3
    splits: str = "quality"          # who proposes forward/reverse splits: quality, both or figaro
    max_unmergeable: float = 0.002  # share of the longest amplicons that may be left unmergeable
    explore_coverage: bool = False  # also show what covering fewer of the longest amplicons would allow
    min_retention: float = 0.5      # floor for the quality-first maxEE suggestion: the share of pairs it must keep
    explore_max_ee: bool = False    # also find the best truncation for several maxEE values and show where retention levels off
    marker: str = "16s"             # 16s: search truncation lengths; its: no truncation (DADA2's ITS workflow)
    min_len: int = 50               # minLen, used with --marker its


@dataclass
class Result:
    cfg: Config
    sub: SubsampleInfo
    inserts: List[int]                   # sorted observed amplicon (insert) lengths
    frac_overlapping: float              # share of pairs that overlap at full read length
    longest_insert: int
    floor: int                           # shortest total length that covers the amplicon range
    max_combined: int                    # longest total length the reads allow
    totals: List[int] = field(default_factory=list)
    candidates: List[Candidate] = field(default_factory=list)
    recommendation: Optional[Candidate] = None
    warnings: List[str] = field(default_factory=list)
    step: int = 2
    exhaustive: bool = False
    tolerance: float = 0.01
    splits: str = "quality"
    short_lost: float = 0.0              # share of pairs shorter than the recommended truncation
    marker: str = "16s"
    min_len: int = 50
    # per total: (total, best FIGARO pass rate, best pass rate over all splits, that split)
    split_checks: list = field(default_factory=list)
    coverage_search: list = field(default_factory=list)   # with --explore-coverage: the best setting for several coverage levels
    min_retention: float = 0.5
    error_profile: object = None                           # expected errors along the reads, for the report's figure
    ee_search: list = field(default_factory=list)     # with --explore-max-ee: the best truncation and yield for each maxEE
    ee_table: list = field(default_factory=list)     # what other maxEE values would keep, at the recommended truncation


# maxEE values shown next to the one in use. The tutorial value (2,2), the strict end, and the relaxed reverse-read
# settings the DADA2 tutorial gives as an example (2,5); nothing above 5 is offered, because values that high are rarely used.
EE_GRID = [(1, 1), (2, 2), (2, 3), (2, 4), (2, 5), (3, 3), (3, 5), (5, 5)]


def _max_ee_table(pairs: ReadPairs, rec: Candidate, cfg: Config, floor: int, work: Path, longest_reads) -> list:
    """Exact pass rate and merge-aware yield at the recommended truncation for each maxEE in EE_GRID (and the one in use)."""
    current = ee_pair(cfg.max_ee)
    grid = sorted({ee_pair(e) for e in EE_GRID} | {current}, key=lambda e: (e[0] + e[1], e[0], e[1]))
    rows = []
    for i, ee in enumerate(grid):
        if ee == current:
            c = rec
        else:
            c = _evaluate(pairs, rec.trim_f, rec.trim_r, replace(cfg, max_ee=ee), floor, work / f"maxee{i}", "maxEE table", longest_reads)
        ee_f, ee_r = pairs.error_load(rec.trim_f, rec.trim_r, ee)
        add_f, add_r = pairs.added_error_load(rec.trim_f, rec.trim_r, ee)
        rows.append({"add_f": add_f, "add_r": add_r, "ee": ee, "pass_rate": c.ee_pass_rate, "merge_rate": c.merge_rate_tolerant, "yield_est": c.yield_est,
                     "ee_f": ee_f, "ee_r": ee_r, "score": 100 * c.ee_pass_rate - ((ee[0] - 1) ** 2 + (ee[1] - 1) ** 2), "current": ee == current})
    return rows


def _evaluate(pairs: ReadPairs, trim_f: int, trim_r: int, cfg: Config, floor: int,
              work: Path, source: str = "", longest_reads=(0, 0)) -> Candidate:
    tag = f"f{trim_f}_r{trim_r}"
    o1, o2 = work / "trimmed" / f"{tag}_R1.fastq", work / "trimmed" / f"{tag}_R2.fastq"
    st = pairs.filter(trim_f, trim_r, o1, o2, cfg.max_ee)
    nan = float("nan")
    strict = tolerant = None
    if st.n_out > 0:
        mdir = work / "merge"
        strict = run_vsearch(o1, o2, mdir, tag, cfg.dada2_min_overlap, cfg.dada2_max_mismatch,
                             "strict", keep_merged=False)
        tolerant = run_vsearch(o1, o2, mdir, tag, cfg.dada2_min_overlap, None, "tolerant",
                               keep_merged=False)
    for p in (o1, o2):
        p.unlink(missing_ok=True)
    total = trim_f + trim_r
    merge = tolerant.merge_rate if tolerant else 0.0
    return Candidate(
        trim_f=trim_f, trim_r=trim_r, ee_pass_rate=st.ee_pass_rate,
        merge_rate_strict=strict.merge_rate if strict else 0.0, merge_rate_tolerant=merge,
        mean_overlap=(total - tolerant.mean_fragment_len) if tolerant else nan,
        mean_overlap_diffs=tolerant.mean_overlap_diffs if tolerant else nan,
        max_mergeable_insert=total - cfg.dada2_min_overlap, covers_longest=total >= floor,
        yield_est=st.ee_pass_rate * merge, source=source,
        balance=abs(trim_f / longest_reads[0] - trim_r / longest_reads[1]) if all(longest_reads) else 0.0)


def _warnings(res: Result, given_longest: bool) -> List[str]:
    w, rec, sub = [], res.recommendation, res.sub
    cfg = res.cfg
    qv = getattr(sub, "quality_values", 0)
    if 0 < qv <= 2:
        w.append(f"Quality scores take only {qv} distinct value{'s' if qv > 1 else ''}, so the per-base quality that "
                 f"DADA2 relies on is missing. In a test, DADA2 could not learn an error model from such data. This "
                 f"is consistent with a simplified-quality file such as the SRA 'lite' format; use the original "
                 f"FASTQ files instead.")
    if not given_longest and res.frac_overlapping < 0.9:
        w.append(f"Only {res.frac_overlapping:.0%} of read pairs could be merged at full read length. "
                 f"Some amplicons may be longer than your reads can span (about "
                 f"{res.max_combined - cfg.dada2_min_overlap} nt), some reads may be shorter than the rest, "
                 f"or the overlaps may be of poor quality or repetitive. The amplicon lengths above come "
                 f"from the pairs that did merge.")
    if res.floor <= res.max_combined and res.max_combined - res.floor < 6:
        w.append(f"The reads only just span your longest amplicons (spare length: "
                 f"{res.max_combined - res.floor} nt). Longer reads would be more reliable.")
    if sub.n_pairs < 2000:
        w.append(f"Only {sub.n_pairs:,} read pairs were analysed, so the estimates are rough.")
    if rec is None:
        if res.floor > res.max_combined:
            w.append(f"No setting can merge your longest amplicons ({res.longest_insert} nt): that needs "
                     f"{res.floor} nt of read in total, and the reads give {res.max_combined} nt.")
        elif res.marker == "its":
            w.append("No read pairs survived filtering and merging.")
        else:
            w.append("No setting that covers your amplicons kept any read pairs after filtering and merging.")
    else:
        if not sub.uniform_length and res.short_lost >= 0.01:
            w.append(f"Read lengths vary (forward {sub.min_len_f}-{sub.read_len_f} nt, reverse "
                     f"{sub.min_len_r}-{sub.read_len_r} nt). DADA2 discards reads shorter than truncLen, "
                     f"which loses {res.short_lost:.0%} of read pairs at these settings.")
        check = next((c for c in res.split_checks if c[0] == rec.total_length), None)
        if check and check[2] - check[1] >= 0.01:
            w.append(f"FIGARO's suggested splits of {rec.total_length} nt pass {check[1]:.1%} of pairs through the "
                     f"quality filter, against {check[2]:.1%} for {check[3][0]} + {check[3][1]}; that split "
                     f"was tested too.")
        if rec.yield_est < 0.7:
            w.append(f"Only {rec.yield_est:.0%} of read pairs are expected to pass filtering and merge.")
        elif rec.ee_pass_rate < 0.8:
            w.append(f"Quality filtering removes {1 - rec.ee_pass_rate:.0%} of read pairs at these settings.")
    return w


def _its_ee_rows(s: Settings, cfg: Config, pairs: ReadPairs, work: Path, current_pass: float, current_merge: float) -> list:
    """For each maxEE in EE_GRID (and the one in use): the share of pairs DADA2's no-truncLen filter keeps, the estimated yield (a vsearch merge of
    what is kept), how long and how noisy the kept reads are, and how noisy the reads are that this limit adds beyond maxEE 2. There is no
    truncation to re-optimise, so unlike the 16S exploration this is cheap and always done."""
    current = ee_pair(cfg.max_ee)
    grid = sorted({ee_pair(e) for e in EE_GRID} | {current}, key=lambda e: (e[0] + e[1], e[0], e[1]))
    rows = []
    for i, ee in enumerate(grid):
        if ee == current:
            pass_rate, merge = current_pass, current_merge
        else:
            o1, o2 = work / "trimmed" / f"its{i}_R1.fastq", work / "trimmed" / f"its{i}_R2.fastq"
            st = pairs.filter_untruncated(o1, o2, ee, s.min_len)
            merge = run_vsearch(o1, o2, work / "merge", f"its{i}", cfg.dada2_min_overlap, None, "tolerant", keep_merged=False).merge_rate if st.n_out > 0 else 0.0
            for p in (o1, o2):
                p.unlink(missing_ok=True)
            pass_rate = st.ee_pass_rate
        load = pairs.untruncated_load(ee, s.min_len, baseline=(2.0, 2.0))
        rows.append({"ee": ee, "current": ee == current, "trim_f": 0, "trim_r": 0, "pass_rate": pass_rate, "merge_rate": merge, "yield_est": pass_rate * merge,
                     "ee_f": load["ee_f"], "ee_r": load["ee_r"], "add_f": load["add_f"], "add_r": load["add_r"], "len_f": load["len_f"], "len_r": load["len_r"]})
    return rows


def _run_its(s: Settings, cfg: Config, sub: SubsampleInfo, inserts: List[int], frac: float, work: Path,
             say: Callable[[str], None]) -> Result:
    """ITS: amplicon lengths vary, so DADA2's workflow does not truncate. Evaluate exactly that: each read
    is cut at its first low-quality base, then kept if long enough, free of Ns and within maxEE."""
    pairs = ReadPairs.load(sub.r1, sub.r2)
    o1, o2 = work / "trimmed" / "its_R1.fastq", work / "trimmed" / "its_R2.fastq"
    st = pairs.filter_untruncated(o1, o2, cfg.max_ee, s.min_len)
    nan = float("nan")
    strict = tolerant = None
    if st.n_out > 0:
        strict = run_vsearch(o1, o2, work / "merge", "its", cfg.dada2_min_overlap, cfg.dada2_max_mismatch,
                             "strict", keep_merged=False)
        tolerant = run_vsearch(o1, o2, work / "merge", "its", cfg.dada2_min_overlap, None, "tolerant",
                               keep_merged=False)
    for p in (o1, o2):
        p.unlink(missing_ok=True)
    l1, l2, short, has_n, ee_fail = pairs.untruncated(cfg.max_ee, s.min_len)
    kept = ~(short | has_n | ee_fail)
    totals = (l1 + l2)[kept]
    merge = tolerant.merge_rate if tolerant else 0.0
    cand = Candidate(
        trim_f=0, trim_r=0, ee_pass_rate=st.ee_pass_rate,
        merge_rate_strict=strict.merge_rate if strict else 0.0, merge_rate_tolerant=merge,
        mean_overlap=(float(totals.mean()) - tolerant.mean_fragment_len) if tolerant else nan,
        mean_overlap_diffs=tolerant.mean_overlap_diffs if tolerant else nan,
        max_mergeable_insert=int(totals.max()) - cfg.dada2_min_overlap if len(totals) else 0,
        covers_longest=True, yield_est=st.ee_pass_rate * merge, source="no truncation")
    longest = s.longest_insert or (inserts[-1] if inserts else 0)
    tol = cfg.yield_tolerance if s.yield_tolerance is None else s.yield_tolerance
    res = Result(cfg, sub, inserts, frac, longest, 0, sub.read_len_f + sub.read_len_r, step=max(1, s.step),
                 exhaustive=s.exhaustive, tolerance=tol, splits=s.splits, marker="its", min_len=s.min_len)
    res.candidates = [cand]
    res.recommendation = cand if cand.yield_est > 0 else None
    say(f"  no truncation: {cand.yield_est:.1%} of pairs usable")
    say("Finding what other maxEE values would keep...")
    res.min_retention = s.min_retention
    res.ee_search = _its_ee_rows(s, cfg, pairs, work, st.ee_pass_rate, merge)
    res.warnings = _warnings(res, bool(s.longest_insert))
    return res


def check_vsearch_path(path: Path, encoding: Optional[str] = None) -> None:
    """The Windows build of vsearch opens files by their names in the system's ANSI code page, so a work folder
    with characters outside it cannot be used. Elsewhere there is nothing to check."""
    if encoding is None:
        if os.name != "nt":
            return
        encoding = "mbcs"
    try:
        str(path.resolve()).encode(encoding)
    except UnicodeEncodeError:
        raise PipelineError(f"The work folder {path} contains characters that the Windows version of vsearch cannot "
                            f"open (it reads file names in the system's ANSI code page). Use a work folder with "
                            f"plain characters, for example C:\\amplitrim_work.") from None


def _search(res: Result, pairs: ReadPairs, cfg: Config, s: Settings, work: Path, say) -> Optional[Candidate]:
    """Test total lengths from the coverage floor upwards under the filter in `cfg` and return the best setting. Fills
    res.candidates, res.totals and res.split_checks."""
    sub, floor, max_combined, tol = res.sub, res.floor, res.max_combined, res.tolerance
    totals = list(range(floor, max_combined + 1, max(1, s.step)))
    table = pairs.pass_table(cfg.max_ee)
    seen, best_so_far = {}, None
    for total in totals:
        group = []
        proposals = {}                    # (trim_f, trim_r) -> who proposed it
        if total == max_combined:         # no trimming at all: the full reads are the only pair
            proposals[(sub.read_len_f, sub.read_len_r)] = "full reads"
        else:
            figaro = [] if s.splits == "quality" else candidates_for_total(
                s.figaro_script, s.figaro_python, sub.dir, work / "figaro", total, cfg.dada2_min_overlap,
                s.top_n)
            for fc in figaro:
                proposals[(fc.trim_f, fc.trim_r)] = "FIGARO"
            best = table.best_split(total, prefer_f=figaro[0].trim_f if figaro else None)
            if s.splits != "figaro":
                for tf, tr in table.near_best_splits(total):
                    proposals.setdefault((tf, tr), "quality table")
                if best:
                    key = (best[0], best[1])
                    proposals[key] = "FIGARO + quality table" if proposals.get(key) == "FIGARO" else "quality table"
            if best and figaro:
                res.split_checks.append((total, max(table.rate(f.trim_f, f.trim_r) for f in figaro),
                                         best[2], (best[0], best[1])))
        for (tf, tr), source in proposals.items():
            if (tf, tr) not in seen:
                seen[(tf, tr)] = _evaluate(pairs, tf, tr, cfg, floor, work, source, (sub.read_len_f, sub.read_len_r))
                res.candidates.append(seen[(tf, tr)])
            group.append(seen[(tf, tr)])
        res.totals.append(total)
        if not group:
            continue
        top = max(group, key=lambda c: c.yield_est)
        say(f"  {total} nt: {top.yield_est:.1%} of pairs usable")
        gain = None if best_so_far is None else top.yield_est - best_so_far
        best_so_far = top.yield_est if best_so_far is None else max(best_so_far, top.yield_est)
        if not s.exhaustive and gain is not None and gain < tol:
            break

    return select_recommendation(res.candidates, tol)


COVERAGE_LEVELS = (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10)    # shares of the longest amplicons that may be left unmergeable


def _coverage_search(res: Result, pairs: ReadPairs, s: Settings, work: Path) -> list:
    """The best setting for several coverage levels: covering fewer of the longest amplicons lets the reads be shorter and cleaner.
    The level in use is the recommendation itself; the others get the same search with a lower coverage floor."""
    if s.longest_insert:                       # the longest amplicon was given, so there is no distribution to trade against
        return []
    cur, quick, rows = res.cfg.max_unmergeable_frac, replace(s, splits="quality"), []
    for i, lvl in enumerate(sorted(set(COVERAGE_LEVELS) | {cur})):
        floor, longest = coverage_floor(res.inserts, res.cfg.dada2_min_overlap, lvl)
        if lvl == cur:
            rec = res.recommendation
        elif floor > res.max_combined:
            rec = None
        else:
            res2 = replace(res, floor=floor, longest_insert=longest, candidates=[], totals=[], split_checks=[], recommendation=None)
            rec = _search(res2, pairs, res.cfg, quick, work / f"coverage{i}", lambda msg: None)
        ee_f, ee_r = pairs.error_load(rec.trim_f, rec.trim_r, res.cfg.max_ee) if rec else (None, None)
        rows.append({"level": lvl, "current": lvl == cur, "longest": rec.max_mergeable_insert if rec else None,
                     "lost": frac_unmergeable(res.inserts, rec.max_mergeable_insert) if rec else None,
                     "total": rec.total_length if rec else None, "trim_f": rec.trim_f if rec else None, "trim_r": rec.trim_r if rec else None,
                     "pass_rate": rec.ee_pass_rate if rec else None, "yield_est": rec.yield_est if rec else None, "ee_f": ee_f, "ee_r": ee_r})
    return rows


def _ee_search(res: Result, pairs: ReadPairs, s: Settings, work: Path) -> list:
    """For each maxEE in EE_GRID (and the one in use): the best truncation under that filter, how many pairs pass, and the yield.
    The quality-table search is used for every row but the one in use, which is the recommendation itself."""
    current = ee_pair(res.cfg.max_ee)
    grid = sorted({ee_pair(e) for e in EE_GRID} | {current}, key=lambda e: (e[0] + e[1], e[0], e[1]))
    quick = replace(s, splits="quality")
    rows = []
    for i, ee in enumerate(grid):
        if ee == current:
            rec = res.recommendation
        else:
            cfg2 = replace(res.cfg, max_ee=ee)
            rec = _search(replace(res, cfg=cfg2, candidates=[], totals=[], split_checks=[], recommendation=None),
                          pairs, cfg2, quick, work / f"eesearch{i}", lambda msg: None)
        ee_f, ee_r = pairs.error_load(rec.trim_f, rec.trim_r, ee) if rec else (None, None)
        add_f, add_r = pairs.added_error_load(rec.trim_f, rec.trim_r, ee) if rec else (None, None)
        rows.append({"add_f": add_f, "add_r": add_r, "ee": ee, "current": ee == current, "trim_f": rec.trim_f if rec else None, "trim_r": rec.trim_r if rec else None,
                     "ee_f": ee_f, "ee_r": ee_r,
                     "pass_rate": rec.ee_pass_rate if rec else None, "merge_rate": rec.merge_rate_tolerant if rec else None,
                     "yield_est": rec.yield_est if rec else None})
    return rows


def run(s: Settings, say: Callable[[str], None] = lambda msg: None) -> Result:
    cfg = Config(dada2_min_overlap=s.min_overlap, max_ee=s.max_ee, max_unmergeable_frac=s.max_unmergeable)
    if shutil.which("vsearch") is None:
        raise PipelineError("vsearch was not found on PATH. Install it (https://github.com/torognes/vsearch).")
    if not 0 < s.min_retention <= 1:
        raise PipelineError("--min-retention must be above 0 and at most 1 (a share of the pairs).")
    if not 0 <= s.max_unmergeable < 0.5:
        raise PipelineError("--allow-unmergeable must be at least 0 and below 0.5 (a share of the amplicons).")
    if s.min_overlap < 1 or min(ee_pair(s.max_ee)) <= 0:
        raise PipelineError("--min-overlap must be at least 1 and --max-ee greater than 0.")
    if s.marker not in ("16s", "its"):
        raise PipelineError(f"Unknown marker: {s.marker}")
    if s.splits not in ("both", "figaro", "quality"):
        raise PipelineError(f"Unknown splits mode: {s.splits}")
    if s.splits != "quality":
        if s.figaro_script is None:
            raise PipelineError("This mode needs --figaro-script.")
        if not s.figaro_script.exists():
            raise PipelineError(f"FIGARO script not found: {s.figaro_script}")
    work, tol = s.work_dir, (cfg.yield_tolerance if s.yield_tolerance is None else s.yield_tolerance)
    check_vsearch_path(work)

    if not s.input_dir.is_dir():
        raise PipelineError(f"Input folder not found: {s.input_dir}")
    try:
        sub = build_subsample(s.input_dir, work / "subsample", s.subsample_every, s.max_pairs)
    except ValueError as err:                      # no paired FASTQ files in the folder
        raise PipelineError(str(err)) from None
    if sub.n_pairs == 0:
        raise PipelineError("No usable read pairs were found in the input files.")
    if s.splits != "quality" and not sub.uniform_length:
        raise PipelineError("FIGARO needs reads of one length, and these differ. Run without --figaro-script "
                            "(reads of different lengths are fine then), or crop the reads to one length first.")
    stride = "" if sub.stride == 1 else f" (every {sub.stride}th)"
    lengths = (f"{sub.read_len_f} + {sub.read_len_r} nt" if sub.uniform_length else
               f"up to {sub.read_len_f} + {sub.read_len_r} nt")
    say(f"Reads: {sub.n_pairs:,} of {sub.n_input_pairs:,} pairs analysed{stride}; {lengths}")
    max_combined = sub.read_len_f + sub.read_len_r

    inserts = sorted(estimate_insert_sizes(sub.r1, sub.r2, work / "stage0", cfg.dada2_min_overlap))
    frac = len(inserts) / sub.n_pairs
    if s.marker == "its":
        return _run_its(s, cfg, sub, inserts, frac, work, say)
    if s.longest_insert:
        longest, floor = s.longest_insert, s.longest_insert + cfg.dada2_min_overlap
        say(f"Longest amplicon: {longest} nt (given); shortest total length that can merge it: {floor}")
    elif len(inserts) < MIN_PAIRS_FOR_ESTIMATE:
        raise PipelineError(
            f"Only {len(inserts)} of {sub.n_pairs:,} read pairs overlap at full read length, so the amplicon "
            f"lengths cannot be estimated. Give the longest amplicon you expect with --longest-insert N (nt).")
    else:
        floor, longest = coverage_floor(inserts, cfg.dada2_min_overlap, cfg.max_unmergeable_frac)
        say(f"Amplicons: {inserts[0]}-{inserts[-1]} nt, longest {longest}; "
            f"shortest total length that can merge them: {floor}")

    res = Result(cfg, sub, inserts, frac, longest, floor, max_combined, step=max(1, s.step),
                 exhaustive=s.exhaustive, tolerance=tol, splits=s.splits)
    totals = list(range(floor, max_combined + 1, max(1, s.step)))
    if not totals:
        res.warnings = _warnings(res, bool(s.longest_insert))
        return res

    pairs = ReadPairs.load(sub.r1, sub.r2)
    res.recommendation = _search(res, pairs, cfg, s, work, say)
    if res.recommendation is not None and not sub.uniform_length:
        r = res.recommendation
        res.short_lost = float(((pairs.len1 < r.trim_f) | (pairs.len2 < r.trim_r)).mean())
    res.error_profile = pairs.error_profile()
    if res.recommendation is not None:
        res.ee_table = _max_ee_table(pairs, res.recommendation, cfg, floor, work, (sub.read_len_f, sub.read_len_r))
        if s.explore_coverage:
            say("Finding the best setting for several coverage levels...")
            res.coverage_search = _coverage_search(res, pairs, s, work)
        if s.explore_max_ee:
            say("Finding the best truncation for several maxEE values...")
            res.min_retention = s.min_retention
            res.ee_search = _ee_search(res, pairs, s, work)
    res.warnings = _warnings(res, bool(s.longest_insert))
    return res
