"""
One shared, deterministic subsample of the reads.

Reads are expected to have adapters and primers already removed, and not to be quality-trimmed.
Their lengths may differ (for example after adapter read-through trimming).
"""

import gzip
import logging
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

DEFAULT_MAX_PAIRS = 50_000

logger = logging.getLogger(__name__)

# Accepts sample_R1.fastq.gz, sample_1.fq, AND standard Illumina names such as
# sample_S1_L001_R1_001.fastq.gz (the optional _001 chunk suffix).
PAIR_RE = re.compile(
    r"^(?P<sample>.+?)_R?(?P<mate>[12])(?:_\d{3})?[._](?:fastq|fq)(?:\.gz)?$")


def find_pairs(input_dir: Path) -> List[Tuple[str, Path, Path]]:
    found = {}
    for p in sorted(Path(input_dir).iterdir()):
        m = PAIR_RE.match(p.name)
        if m:
            found.setdefault(m["sample"], {})[m["mate"]] = p
    pairs = [(s, d["1"], d["2"]) for s, d in sorted(found.items()) if "1" in d and "2" in d]
    if not pairs:
        seen = sorted(p.name for p in Path(input_dir).iterdir() if p.is_file())[:5]
        raise ValueError(
            f"No paired FASTQ files found in {input_dir}. Expected names like "
            f"sample_R1.fastq.gz / sample_R2.fastq.gz or sample_S1_L001_R1_001.fastq.gz "
            f"(both mates present). Files seen: {seen if seen else 'none'}")
    return pairs


def read_fastq(path: Path) -> Iterator[Tuple[str, str, str, str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().rstrip("\n")
            p = fh.readline().rstrip("\n")
            q = fh.readline().rstrip("\n")
            yield h.rstrip("\n"), s, p, q


@dataclass
class SubsampleInfo:
    n_pairs: int             # pairs kept
    read_len_f: int          # longest forward / reverse read
    read_len_r: int
    uniform_length: bool
    dir: Path                # folder holding the subsample (also FIGARO's input)
    r1: Path
    r2: Path
    n_input_pairs: int = 0   # pairs in the input files
    stride: int = 1          # every Nth pair was kept
    min_len_f: int = 0       # shortest forward / reverse read
    min_len_r: int = 0
    quality_values: int = 0  # distinct quality scores seen in a sample of the reads (0 = not measured)


def count_pairs(path: Path) -> int:
    """Number of reads in a FASTQ file (counts newlines in large blocks)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            n += chunk.count(b"\n")
    return n // 4


def build_subsample(input_dir: Path, out_dir: Path, every_nth: int = 1,
                    max_pairs: Optional[int] = DEFAULT_MAX_PAIRS) -> SubsampleInfo:
    """Write a deterministic subsample: every Nth pair across ALL samples, where N is the larger of
    `every_nth` and the stride needed to stay within `max_pairs` (None or 0 = no cap)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    r1_path, r2_path = out_dir / "subsample_R1.fastq.gz", out_dir / "subsample_R2.fastq.gz"
    pairs = find_pairs(input_dir)
    n_input = sum(count_pairs(p1) for _s, p1, _p2 in pairs)
    stride = max(every_nth, math.ceil(n_input / max_pairs)) if max_pairs else every_nth
    n_seen = n_kept = 0
    lens_f, lens_r, qvals = set(), set(), set()
    with gzip.open(r1_path, "wt", newline="\n") as o1, \
            gzip.open(r2_path, "wt", newline="\n") as o2:
        for _sample, p1, p2 in pairs:
            for (h1, s1, pl1, q1), (h2, s2, pl2, q2) in zip(read_fastq(p1), read_fastq(p2)):
                n_seen += 1
                if (n_seen - 1) % stride or not s1 or not s2:
                    continue
                o1.write(f"{h1}\n{s1}\n{pl1}\n{q1}\n")
                o2.write(f"{h2}\n{s2}\n{pl2}\n{q2}\n")
                lens_f.add(len(s1))
                lens_r.add(len(s2))
                n_kept += 1
                if n_kept <= 200 or n_kept % 25 == 0:
                    qvals.update(q1)
                    qvals.update(q2)
    uniform = len(lens_f) == 1 and len(lens_r) == 1
    return SubsampleInfo(n_kept, max(lens_f, default=0), max(lens_r, default=0), uniform, out_dir,
                         r1_path, r2_path, n_input, stride, min(lens_f, default=0), min(lens_r, default=0), len(qvals))
