"""
Stage 2: sequence-based merge testing with vsearch.

vsearch exposes the two settings that mirror DADA2's mergePairs (minimum overlap, maximum
mismatches) and reports the mismatches observed in each overlap.

Two vsearch passes per candidate:
  strict    - maxdiffs = DADA2's maxMismatch (default 0) applied to RAW reads.
              A LOWER bound on DADA2's merge success: DADA2 merges *denoised*
              sequences, from which sequencing errors have been removed.
  tolerant  - vsearch's default diff allowance. Closer proxy for "do these two
              reads geometrically overlap", i.e. what DADA2 can merge after
              denoising. Used for ranking.

Both tools are invoked as subprocesses only; no code is vendored.
"""

import logging
import os
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)
NAN = float("nan")


@dataclass
class MergeOutcome:
    engine: str
    setting: str
    total_pairs: int
    merged_pairs: int
    mean_fragment_len: float = NAN
    mean_overlap_diffs: float = NAN       # mean observed mismatches in merged region
    fail_reasons: Dict[str, int] = field(default_factory=dict)
    merged_path: Optional[Path] = None

    @property
    def merge_rate(self) -> float:
        return self.merged_pairs / self.total_pairs if self.total_pairs else 0.0


# ---------------------------------------------------------------- vsearch ---

def parse_vsearch_log(text: str, setting: str = "") -> MergeOutcome:
    def num(pattern, cast=float, default=NAN):
        m = re.search(pattern, text, re.MULTILINE)
        return cast(m.group(1)) if m else default

    total = num(r"^\s*(\d+)\s+Pairs", int, 0)
    merged = num(r"^\s*(\d+)\s+Merged", int, 0)
    reasons: Dict[str, int] = {}
    block = re.search(r"failed merging due to various reasons:\n(.*?)\n\s*\n", text, re.DOTALL)
    if block:
        for line in block.group(1).splitlines():
            m = re.match(r"^\s*(\d+)\s+(.+?)\s*$", line)
            if m:
                reasons[m.group(2)] = int(m.group(1))
    return MergeOutcome(
        engine="vsearch", setting=setting, total_pairs=total, merged_pairs=merged,
        mean_fragment_len=num(r"^\s*([\d.]+)\s+Mean fragment length"),
        mean_overlap_diffs=num(r"^\s*([\d.]+)\s+Mean observed errors in merged region\s*$"),
        fail_reasons=reasons,
    )


def run_vsearch(r1: Path, r2: Path, out_dir: Path, tag: str, min_overlap: int,
                max_diffs: Optional[int], setting: str, keep_merged: bool = True,
                threads: int = 1) -> MergeOutcome:
    out_dir.mkdir(parents=True, exist_ok=True)
    merged = out_dir / f"{tag}.{setting}.merged.fasta"
    log = out_dir / f"{tag}.{setting}.vsearch.log"
    cmd = ["vsearch", "--fastq_mergepairs", str(r1), "--reverse", str(r2),
           "--fastaout", str(merged if keep_merged else os.devnull),
           "--fastq_minovlen", str(min_overlap), "--fastq_allowmergestagger",
           "--threads", str(threads), "--log", str(log), "--quiet"]
    if max_diffs is not None:
        cmd += ["--fastq_maxdiffs", str(max_diffs)]
    res = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if res.returncode != 0:
        raise RuntimeError(f"vsearch failed ({' '.join(cmd)}):\n{res.stderr}")
    out = parse_vsearch_log(log.read_bytes().decode("utf-8", errors="replace"), setting)   # paths in it may be in the ANSI code page
    out.merged_path = merged if keep_merged else None
    return out


def merged_lengths(fasta: Path) -> List[int]:
    """Lengths of merged sequences (= observed insert/fragment sizes)."""
    lengths, cur = [], 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    lengths.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
    if cur:
        lengths.append(cur)
    return lengths


# ------------------------------------------------------------------ Stage 0 ---

def estimate_insert_sizes(r1: Path, r2: Path, out_dir: Path, min_overlap: int) -> List[int]:
    """Observed insert (fragment) sizes from a tolerant merge of FULL-length
    reads. The upper tail of this distribution is the data-driven answer to
    'what is the longest amplicon?' - the quantity FIGARO's -a cannot discover
    for itself. Caveat: only pairs that overlap at full length are seen, so
    inserts approaching (len_f + len_r - min_overlap) are under-represented."""
    out = run_vsearch(r1, r2, out_dir, "stage0", min_overlap, None, "tolerant")
    return merged_lengths(out.merged_path)
