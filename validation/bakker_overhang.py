#!/usr/bin/env python
"""Bakker's ASVs: exact reference sequences, reference sequences carrying primer overhangs, or other (register sections 38 and 39).

  python validation\\bakker_overhang.py --run results\\28_its_bakker --references validation\\references_bakker

Where an amplicon is shorter than the reads, each read runs past its far end; when the pair is merged, the overhangs put the forward primer (or its last
bases) back at the start of the sequence and the reverse primer's reverse complement at its end. Such a sequence is a true member's sequence that the
2-edit rule counts as spurious. Every ASV of the arm (std by default) of each Bakker dataset is classified, weighted by its reads. Writes OVERHANG_CHECK.md to --run.
"""
import argparse
import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bakker_primers as bp  # noqa: E402
from bakker_prep import AMPLICONS, rc  # noqa: E402


def read_refs(path: Path) -> list:
    out = []
    for l in Path(path).read_text(encoding="utf-8").splitlines():
        if l and not l.startswith(">"):
            out.append(l.strip())
    return out


def overhanging(seq: str, refs: list, fwd: str, rev_rc: str, max_mm: int = 2) -> bool:
    """A reference inside seq, with at most the forward primer's last bases before it and the reverse primer's first bases (reverse-complemented) after it."""
    for r in refs:
        i = seq.find(r)
        if i < 0:
            continue
        head, tail = seq[:i], seq[i + len(r):]
        if (head or tail) and len(head) <= len(fwd) and len(tail) <= len(rev_rc) \
                and (not head or bp.mismatches(head, fwd[-len(head):]) <= max_mm) and (not tail or bp.mismatches(tail, rev_rc[:len(tail)]) <= max_mm):
            return True
    return False


def classify(seqs: dict, refs: list, fwd: str, rev: str) -> dict:
    """seqs: {sequence: reads}. Returns the shares of the reads that are exact, overhang and other."""
    exact, rr, c = set(refs), rc(rev), {"exact": 0, "overhang": 0, "other": 0}
    for s, n in seqs.items():
        c["exact" if s in exact else "overhang" if overhanging(s, refs, fwd, rr) else "other"] += n
    tot = sum(c.values()) or 1
    return {k: v / tot for k, v in c.items()} | {"reads": sum(c.values())}


def run_check(run: Path, references: Path, arm: str = "std") -> dict:
    out = {}
    for label, (_, f, r) in AMPLICONS.items():
        p = Path(run) / label / "dada2" / arm / "seqs.tsv"
        if not p.is_file():
            continue
        seqs = {}
        with open(p, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                seqs[row["sequence"]] = seqs.get(row["sequence"], 0) + int(float(row["abundance"]))
        out[label] = classify(seqs, read_refs(Path(references) / f"references_{label}.fasta"), bp.PRIMERS[f], bp.PRIMERS[r])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--arm", default="std")
    a = ap.parse_args(argv)
    res = run_check(a.run, a.references, a.arm)
    L = ["# Bakker's ASVs: exact, with primer overhangs, other", "", f"Arm `{a.arm}`; shares of the reads.", "",
         "| Dataset | Reads | Exact reference | Reference with primer overhangs | Other |", "|---|---|---|---|---|"]
    L += [f"| {lab} | {x['reads']:,} | {x['exact']:.1%} | {x['overhang']:.1%} | {x['other']:.1%} |" for lab, x in res.items()]
    (a.run / "OVERHANG_CHECK.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0 if res else 2


if __name__ == "__main__":
    sys.exit(main())
