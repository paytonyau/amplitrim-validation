#!/usr/bin/env python
"""Rebuild the expected sequences of the staggered mock community A (Huelpuesch, Rauer et al. 2023, BMC Biology 21:269) from the paper's repository.

The mock has 15 cultured strains. The authors Sanger-sequenced each strain and defined the mock ASVs of their DADA2 run as the ASVs that contain the
middle of a strain's Sanger read (positions 68-430; ambiguous calls match either base), plus ASVs within 4 edits of those. This script repeats that
from the raw chromatograms and the authors' ASV table, and writes the matching ASVs as a reference FASTA, named >ASV|STRAIN|species so that the
variants (operon copies) of one strain are grouped by mock_check.py. It reproduces the paper's 52 mock ASVs to within one ASV (53 at tolerance 3).

  python validation/build_mockA_reference.py --out mockA_reference.fasta [--tolerance 3] [--cache folder] [--samples REGEX [--min-reads N]]

The full reference holds every variant found in any of the 23 mock samples, so the variants that exist only in other dilutions (a third of them) cannot be
recovered from a subset of samples and would be counted as missed. With --samples (a regular expression on the authors' sample names, for example
"^(103|104|111|112|119|120)-KitA" for the 10^8 and 10^9 cell samples) the reference keeps only the variants with at least --min-reads reads in those samples.

Needs internet (GitHub), pandas, biopython and rapidfuzz. The source repository, https://github.com/LuiseRauer/Decontamination_benchmarking, is GPL-3.0:
this script downloads from it and does not redistribute its data. The reference is the authors' DADA2 ASVs checked against Sanger sequences, not independent
truth, and it favours the authors' own DADA2 settings (truncLen 299/280, trimLeft 20/17, maxEE 2,6).
"""

import argparse
import csv
import html
import re
import subprocess
import sys
from pathlib import Path

REPO = "Input/Mock_staggered-A_Huelpuesch-Rauer-et-al"
RAW = "https://raw.githubusercontent.com/LuiseRauer/Decontamination_benchmarking/main/" + REPO
PAGE = "https://github.com/LuiseRauer/Decontamination_benchmarking/tree/main/" + REPO
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
BASES = {c: i for i, c in enumerate("ACGT")}


def region_mask(region: str):
    """A (length, 4) boolean array: which bases each position of a Sanger region allows (an ambiguity code allows several)."""
    import numpy as np
    mask = np.zeros((len(region), 4), bool)
    for i, c in enumerate(region):
        for b in IUPAC.get(c, "ACGT"):
            mask[i, BASES[b]] = True
    return mask


def best_mismatches(asv: str, region: str):
    """Fewest positions at which the ASV disagrees with the region over every place the region fits in it; None if the ASV is shorter."""
    import numpy as np
    n = len(region)
    if len(asv) < n:
        return None
    enc = np.array([BASES.get(c, 0) for c in asv])
    windows = np.lib.stride_tricks.sliding_window_view(enc, n)
    allowed = region_mask(region)[np.arange(n)[None, :], windows]
    return int((~allowed).sum(axis=1).min())


def match_asvs(asvs: dict, regions: dict, tolerance: int) -> dict:
    """{strain: {asv id: mismatches}} for the ASVs that contain a strain's region with at most `tolerance` mismatches."""
    out = {}
    for strain, region in regions.items():
        hits = {}
        for asv_id, seq in asvs.items():
            m = best_mismatches(seq, region)
            if m is not None and m <= tolerance:
                hits[asv_id] = m
        out[strain] = hits
    return out


def add_close(asvs: dict, direct: set, max_edits: int = 4) -> set:
    """The ASVs not in `direct` that are within `max_edits` edits of one that is."""
    from rapidfuzz.distance import Levenshtein
    return {a for a, s in asvs.items() if a not in direct and direct and min(Levenshtein.distance(s, asvs[d]) for d in direct) <= max_edits}


def build(asvs: dict, regions: dict, species: dict, tolerance: int = 3, max_edits: int = 4):
    """(FASTA text, rows of the mapping table, strains without a match)."""
    matches = match_asvs(asvs, regions, tolerance)
    direct = {a for hits in matches.values() for a in hits}
    final = sorted(direct | add_close(asvs, direct, max_edits))
    fasta, rows = [], []
    for a in final:
        candidates = [(strain, hits[a]) for strain, hits in matches.items() if a in hits]
        strain, mismatches = min(candidates, key=lambda x: x[1]) if candidates else ("near", None)       # found only by its edit distance
        label = species.get(strain, strain).replace(" ", "_")
        fasta.append(f">{a}|{strain}|{label}\n{asvs[a]}\n")
        rows.append((a, strain, species.get(strain, ""), mismatches, len(asvs[a])))
    return "".join(fasta), rows, [s for s, hits in matches.items() if not hits]


def present_in(counts, ids: list, sample_regex: str, min_reads: int = 1) -> list:
    """The ids with at least `min_reads` reads in total over the samples (columns of the table `counts`) whose names match the regular expression."""
    columns = [c for c in counts.columns if re.search(sample_regex, str(c))]
    if not columns:
        raise SystemExit(f"--samples {sample_regex!r} matches no sample in the authors' ASV table")
    totals = counts.loc[ids, columns].apply(lambda col: col.astype(float)).sum(axis=1)
    return [i for i in ids if totals[i] >= min_reads]


def fetch(url: str, dest: Path) -> None:
    res = subprocess.run(["curl", "-s", "-m", "90", "-L", "-f", "-o", str(dest), url], capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit(f"could not download {url} (curl exit {res.returncode}); check the internet connection")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=Path("mockA_reference.fasta"))
    ap.add_argument("--tolerance", type=int, default=3, help="mismatches allowed between a Sanger region and an ASV (default 3)")
    ap.add_argument("--samples", help="keep only the variants present in the authors' samples whose names match this regular expression")
    ap.add_argument("--min-reads", type=int, default=1, help="with --samples: the least reads (summed over those samples) a variant needs (default 1)")
    ap.add_argument("--cache", type=Path, default=Path("mockA_reference_downloads"), help="folder for the downloaded files")
    args = ap.parse_args()
    import importlib.util
    missing = [m for m in ("pandas", "Bio", "rapidfuzz") if importlib.util.find_spec(m) is None]
    if missing:
        raise SystemExit(f"Missing: {', '.join(missing)}. Install what is missing: pip install pandas biopython rapidfuzz")
    import pandas as pd
    from Bio import SeqIO
    sanger = args.cache / "sanger"
    sanger.mkdir(parents=True, exist_ok=True)
    table, overview = args.cache / "asvs.csv", args.cache / "overview.csv"
    fetch(f"{RAW}/Mock_staggered-A_ASVs.csv", table)
    fetch(f"{RAW}/Mock_staggered_A_community-overview_short.csv", overview)
    page = subprocess.run(["curl", "-s", "-m", "60", "-L", "-A", "Mozilla/5.0", f"{PAGE}/SangerSequences"], capture_output=True, text=True).stdout
    names = list(dict.fromkeys(html.unescape(n) for n in re.findall(r"SangerSequences/([^\"]+\.ab1)\"", page)))
    if not names:
        raise SystemExit("could not list the Sanger files on GitHub")
    regions = {}
    for n in names:
        m = re.match(r"\w\d_(S\d+)_", n)
        if not m:
            continue                                                            # the 27F primer file
        fetch(f"{RAW}/SangerSequences/{n}", sanger / n)
        regions[m.group(1)] = str(SeqIO.read(sanger / n, "abi").seq)[67:430]
    species = {}
    with open(overview, encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh, delimiter=";"):
            species[r["ID Sanger sequencing"]] = r["Species"].strip()
    df = pd.read_csv(table, sep=";", dtype=str)
    fasta, rows, unmatched = build(dict(zip(df["OTU_ID"], df["Sequence"])), regions, species, args.tolerance)
    if args.samples:
        counts = df.set_index("OTU_ID")
        keep = set(present_in(counts, [r[0] for r in rows], args.samples, args.min_reads))
        dropped = len(rows) - len(keep)
        rows = [r for r in rows if r[0] in keep]
        fasta = "".join(f">{rec}" for rec in fasta.split(">")[1:] if rec.split("|")[0] in keep)
        print(f"kept the {len(keep)} variants with at least {args.min_reads} read(s) in the samples matching {args.samples!r}; left out {dropped} that are not in them")
        unmatched = [s for s in regions if s not in {r[1] for r in rows}]
    args.out.write_text(fasta, encoding="utf-8")
    with open(args.out.with_suffix(".tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["asv", "strain_id", "species", "sanger_mismatches", "length"])
        w.writerows(rows)
    print(f"{len(rows)} reference sequences for {len(regions) - len(unmatched)} of {len(regions)} strains written to {args.out} (map: {args.out.with_suffix('.tsv')})")
    if unmatched:
        print("strains with no variant in the reference (the paper reports Streptococcus pneumoniae as not detected; with --samples, strains absent from those samples are listed too): " + ", ".join(f"{s} {species.get(s, '')}" for s in unmatched))
    return 0


if __name__ == "__main__":
    sys.exit(main())
