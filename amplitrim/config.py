"""
Defaults for the settings that should mirror the DADA2 call you will actually run (for example
`qiime dada2 denoise-paired --p-max-ee-f 2`), because merging is tested under those rules.
`min_overlap` and `max_ee` can be set on the command line.
"""

import argparse
from dataclasses import dataclass


def ee_pair(max_ee):
    """(forward, reverse) maxEE from one number (both reads) or a pair, as DADA2's maxEE = c(F, R)."""
    if isinstance(max_ee, (tuple, list)):
        return float(max_ee[0]), float(max_ee[1])
    return float(max_ee), float(max_ee)


def parse_ee(text):
    """Command-line value for maxEE: '2' (both reads) or '2,5' (forward, reverse). Equal values collapse to one number."""
    try:
        parts = [float(x) for x in str(text).replace(";", ",").split(",")]
    except ValueError:
        parts = []
    if len(parts) not in (1, 2) or any(p <= 0 for p in parts):
        raise argparse.ArgumentTypeError(f"'{text}' is not a positive number or a pair such as 2,5 (forward, reverse)")
    return parts[0] if len(parts) == 1 or parts[0] == parts[1] else (parts[0], parts[1])


@dataclass
class Config:
    dada2_min_overlap: int = 12        # mergePairs minOverlap (DADA2's default)
    dada2_max_mismatch: int = 0        # mergePairs maxMismatch (DADA2's default)
    max_ee: object = 2.0              # one number, or (forward, reverse) as in maxEE = c(F, R)

    # At most this share of the observed amplicons may be longer than the longest the settings
    # can merge; otherwise long amplicons are systematically lost.
    max_unmergeable_frac: float = 0.002

    # Settings whose usable-pair yield is within this (absolute) of the best are treated as tied;
    # the one with the smallest overlap then wins. The search also stops once a longer total
    # length improves the yield by less than this.
    yield_tolerance: float = 0.01
