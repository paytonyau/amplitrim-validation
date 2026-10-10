"""
Optional cross-check: FIGARO's quality-based trim positions.

FIGARO is called as a subprocess; its algorithm is not reimplemented or vendored. Its candidates
depend only on the sum of its -a and -m options, which is the total retained length
(truncLenF + truncLenR), so a run is requested by that total alone. Reads are expected to carry no
primers, so FIGARO is told the primer lengths are 0.

Notes on stock FIGARO:
  * It needs reads of one length.
  * It calls `collections.Iterable`, which Python removed in 3.10. FIGARO is therefore launched with
    a small shim (amplitrim/_compat/sitecustomize.py) on PYTHONPATH that restores the removed
    aliases; the FIGARO script itself is not modified.
  * Its default filename parser rejects Illumina-style names (sample_R1_001.fastq.gz); the
    subsample written by AmpliTrim uses names it accepts.
"""

import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List

COMPAT_DIR = Path(__file__).parent / "_compat"


@dataclass
class FigaroCandidate:
    trim_f: int                # DADA2 truncLenF
    trim_r: int                # DADA2 truncLenR


@dataclass
class FigaroRecommendation:
    """FIGARO's own top entry, as a user following FIGARO would use it: trim positions together with the maxEE it pairs them with."""
    trim_f: int
    trim_r: int
    max_ee: tuple              # (forward, reverse) maxExpectedError
    retention: float           # FIGARO's predicted share of reads kept by the filter (0-1); it does not model merging
    score: float


def parse_figaro_recommendation(path: Path) -> FigaroRecommendation:
    """The first entry of FIGARO's output, which it lists best first."""
    entries = json.loads(path.read_text(encoding="utf-8"))
    if not entries:
        raise RuntimeError("FIGARO listed no candidates for these reads.")
    e = entries[0]
    return FigaroRecommendation(trim_f=e["trimPosition"][0], trim_r=e["trimPosition"][1],
                                max_ee=(e["maxExpectedError"][0], e["maxExpectedError"][1]),
                                retention=e["readRetentionPercent"] / 100, score=e["score"])


def figaro_env(compat: bool = True) -> dict:
    """Environment for the FIGARO subprocess: the current one, plus the compatibility shim first
    on PYTHONPATH."""
    env = os.environ.copy()
    if compat:
        parts = [str(COMPAT_DIR)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])
        env["PYTHONPATH"] = os.pathsep.join(parts)
    return env


def run_figaro(figaro_script: Path, python_exe: str, input_dir: Path, out_dir: Path,
               a: int, m: int, compat: bool = True) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [python_exe, str(figaro_script), "-i", str(input_dir), "-o", str(out_dir),
           "-f", "0", "-r", "0", "-a", str(a), "-m", str(m)]
    res = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace", env=figaro_env(compat))
    path = out_dir / "trimParameters.json"
    if res.returncode != 0 or not path.exists():
        hint = ""
        if "No module named 'scipy'" in res.stderr:
            hint = '\nFIGARO needs scipy. Install it with: pip install "amplitrim[figaro]"'
        raise RuntimeError(f"FIGARO failed (exit {res.returncode}).\n{res.stderr[-600:]}{hint}")
    return path


def parse_figaro_json(path: Path, top_n: int) -> List[FigaroCandidate]:
    """The first `top_n` candidates; FIGARO lists them best first."""
    return [FigaroCandidate(trim_f=e["trimPosition"][0], trim_r=e["trimPosition"][1])
            for e in json.loads(path.read_text(encoding="utf-8"))[:top_n]]


def candidates_for_total(figaro_script: Path, python_exe: str, input_dir: Path, work_dir: Path,
                         total: int, min_overlap: int, top_n: int = 3,
                         compat: bool = True) -> List[FigaroCandidate]:
    """FIGARO's best trim pairs whose lengths add up to `total` (truncLenF + truncLenR)."""
    path = run_figaro(figaro_script, python_exe, input_dir, work_dir / f"total{total}",
                      total - min_overlap, min_overlap, compat)
    return parse_figaro_json(path, top_n)
