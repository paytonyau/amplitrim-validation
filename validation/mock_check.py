#!/usr/bin/env python3
"""Compare the ASVs DADA2 found in a mock community with the sequences that should be there.

    python validation/mock_check.py asvs_120x118.fasta HMP_MOCK.v35.fasta

The ASV file is the one real_data.py --dada2 writes (asvs_<F>x<R>.fasta; ">ASV1;size=N" headers). The reference is a FASTA of the
expected sequences; an ASV may be shorter than its reference (a sub-region), because the terminal gaps are not counted.
Every ASV is matched to its best reference with vsearch and is counted as
  exact  - identical to part of a reference
  near   - within a few bases of a reference but not identical (--min-id 0.99, about 4-5 bases over 470 nt, as in the rule of the paper that
           defined the mock A reference): residual errors, or a real variant the reference does not list
  none   - no reference that close (chimeras, contaminants, or reads that are mostly errors)
All shares are of the reads in the ASVs. Needs vsearch on PATH.

Variants: one strain often has several 16S sequences (operon copies that differ by a few bases), so a reference may list several sequences per
strain. Name them >ID|STRAIN|label and the tool also reports how many strains have at least one variant recovered, next to the variants.
"""
import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def read_fasta(path):
    records, name, seq = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name, seq = line[1:].strip(), []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def size_of(header: str) -> int:
    m = re.search(r"size=(\d+)", header)
    return int(m.group(1)) if m else 1


def read_asvs(asv_fasta) -> list:
    """The ASVs in a FASTA file; none when DADA2 wrote no file at all, which is what happens when no pair merged."""
    return read_fasta(asv_fasta) if Path(asv_fasta).exists() else []


def strain_of(name: str) -> str:
    """The strain of a reference named ID|STRAIN|label; any other name is its own strain."""
    parts = name.split("|")
    return parts[1] if len(parts) >= 3 else name


def strain_counts(ref_names: list, recovered: list):
    """(strains with at least one variant recovered, strains), or (None, None) when no strain has more than one variant."""
    strains = {}
    for r in ref_names:
        strains.setdefault(strain_of(r), [0, 0])[1] += 1
    if len(strains) == len(ref_names):
        return None, None
    for r in recovered:
        strains[strain_of(r)][0] += 1
    return sum(1 for found, _ in strains.values() if found), len(strains)


def summarise(asv_fasta, reference_fasta, min_id: float = 0.99, vsearch: str = "vsearch") -> dict:
    asvs, refs = read_asvs(asv_fasta), read_fasta(reference_fasta)
    ref_names = [n.split()[0] for n, _ in refs]
    if not asvs:                     # nothing to score: every expected sequence is missing
        found, total = strain_counts(ref_names, [])
        return {"asvs": 0, "reads": 1, "categories": {"exact": (0, 0), "near": (0, 0), "none": (0, 0)}, "per_reference": {}, "recovered": [],
                "only_near": [], "missing": ref_names, "references": len(ref_names), "split": {}, "strains_found": found, "strains_total": total}
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "hits.tsv"
        res = subprocess.run([vsearch, "--usearch_global", str(asv_fasta), "--db", str(reference_fasta), "--id", str(min_id), "--strand", "both",
                              "--maxaccepts", "0", "--maxrejects", "0", "--userout", str(out), "--userfields", "query+target+id", "--quiet"],
                             capture_output=True, text=True)
        if res.returncode != 0:
            raise SystemExit(f"vsearch failed: {res.stderr.strip()[-300:]}")
        best = {}
        for line in out.read_text().splitlines():
            q, t, ident = line.split("\t")
            if q not in best or float(ident) > best[q][1]:
                best[q] = (t.split()[0], float(ident))
    cats = {"exact": [0, 0], "near": [0, 0], "none": [0, 0]}
    per_ref, unmatched = {}, []
    for name, seq in asvs:
        size, hit = size_of(name), best.get(name.split()[0] if name.split() else name)
        if hit is None:
            cat = "none"
            unmatched.append((size, len(seq), name.split(";")[0]))
        else:
            cat = "exact" if hit[1] >= 99.999 else "near"
            per_ref.setdefault(hit[0], {"exact": 0, "near": 0})[cat] += 1
        cats[cat][0] += 1
        cats[cat][1] += size
    total_reads = sum(c[1] for c in cats.values()) or 1
    recovered = [r for r in ref_names if per_ref.get(r, {}).get("exact")]
    found, total = strain_counts(ref_names, recovered)
    return {"asvs": len(asvs), "reads": total_reads, "categories": cats, "per_reference": per_ref, "strains_found": found, "strains_total": total,
            "recovered": recovered,
            "only_near": [r for r in ref_names if per_ref.get(r) and not per_ref[r]["exact"]],
            "missing": [r for r in ref_names if r not in per_ref], "references": len(ref_names),
            "split": {r: v["exact"] + v["near"] for r, v in per_ref.items() if v["exact"] + v["near"] > 1},
            "unmatched": sorted(unmatched, reverse=True)[:5]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("asvs")
    ap.add_argument("reference")
    ap.add_argument("--min-id", type=float, default=0.99, help="identity below which an ASV counts as having no reference; 0.99 is about 4-5 bases over 470 nt (default 0.99)")
    args = ap.parse_args()
    vsearch = shutil.which("vsearch")
    if not vsearch:
        raise SystemExit("vsearch was not found on PATH.")
    s = summarise(args.asvs, args.reference, args.min_id, vsearch)
    print(f"{s['asvs']} ASVs carrying {s['reads']:,} reads; {s['references']} reference sequences\n")
    for label, key in (("identical to a reference", "exact"), (f"within {1 - args.min_id:.1%} of a reference, not identical", "near"),
                       ("no reference within that distance", "none")):
        n, reads = s["categories"][key]
        print(f"  {label + ':':44s} {n:4d} ASVs, {reads / s['reads']:6.1%} of the reads")
    print(f"\nreference sequences recovered by an identical ASV: {len(s['recovered'])} of {s['references']}")
    if s["strains_total"]:
        print(f"strains with at least one variant recovered: {s['strains_found']} of {s['strains_total']}")
    if s["only_near"]:
        print(f"  only a near match: {', '.join(s['only_near'])}")
    if s["missing"]:
        print(f"  not found at all:  {', '.join(s['missing'])}")
    if s["split"]:
        print("references represented by more than one ASV (residual errors or real variants): "
              + ", ".join(f"{r} ({n})" for r, n in sorted(s["split"].items(), key=lambda kv: -kv[1])[:8]))
    if s["unmatched"]:
        print("largest ASVs without a reference: " + ", ".join(f"{n} (size {size}, {length} nt)" for size, length, n in s["unmatched"]))


if __name__ == "__main__":
    sys.exit(main())
