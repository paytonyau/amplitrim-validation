#!/usr/bin/env python
"""Recovery against input for the staggered mock community A: one DADA2 run and one scoring per dilution.

The mock (15 strains, 18% to 0.18%) was diluted tenfold from 10^9 to 10^2 cells in triplicate, and contaminants take over as the input falls (under 2% of
the reads at 10^8 and 10^9, over 96% below 10^4 cells). For every dilution this runs DADA2 on the samples of that dilution only, scores the ASVs against the
variants that are present at that dilution (the reference variants with reads in those samples in the authors' own ASV table), and prints how many of those
variants and strains are recovered, the share of reads in identical ASVs, and the ASVs that are not in the reference (contaminants, chimeras, errors).

  python validation/mockA_dilutions.py READS --runinfo runs_info.tsv --asv-table Mock_staggered-A_ASVs.csv --reference mockA_reference.fasta --trunc 257,260 --max-ee 2,6

READS is the folder of primer-free reads (the `ready` folder of real_data.py, samples named by run accession); --runinfo is the ENA table (it gives each run's
sample name, such as 103-KitA-1-1-8, whose last digit is log10 of the input cells); --asv-table is the authors' table from
https://github.com/LuiseRauer/Decontamination_benchmarking (Input/Mock_staggered-A_Huelpuesch-Rauer-et-al/Mock_staggered-A_ASVs.csv, GPL-3.0);
--reference is the file build_mockA_reference.py writes (all 53 variants: the presence at each dilution is taken from the table).
"""

import argparse
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

SAMPLE_PATTERN = r"^\d+-KitA-1-\d-(\d)_"          # the group is the dilution: log10 of the cells put in


def group_samples(design: dict, available: set, pattern: str = SAMPLE_PATTERN) -> dict:
    """{dilution digit: sorted run accessions} for the runs in `available` whose sample name matches the pattern (design: run -> (name, pairs))."""
    out = {}
    for run, (name, _pairs) in design.items():
        m = re.search(pattern, name)
        if m and run in available:
            out.setdefault(m.group(1), []).append(run)
    return {d: sorted(runs) for d, runs in out.items()}


def restricted_reference(reference: list, counts, dilution: str, min_reads: int = 1) -> list:
    """The reference records (name, sequence) whose variant has at least `min_reads` reads in the authors' samples of this dilution."""
    from build_mockA_reference import present_in
    ids = [name.split()[0].split("|")[0] for name, _ in reference]
    keep = set(present_in(counts, ids, rf"-KitA-\d-\d-{dilution}$", min_reads))
    return [(name, seq) for name, seq in reference if name.split()[0].split("|")[0] in keep]


def strain_of(name: str) -> str:
    parts = name.split()[0].split("|")
    return parts[1] if len(parts) >= 3 else name


def format_table(rows: list) -> list:
    head = (f"{'cells':>6s} {'samples':>7s} {'pass':>6s} {'merged':>7s} {'ASVs':>5s} | {'variants':>9s} {'strains':>8s} | {'identical':>9s} {'near':>6s} {'none':>6s} | "
            f"{'ASVs outside the reference':>26s}")
    out = [head]
    for r in rows:
        if r.get("error"):
            out.append(f"{'1e' + r['dilution']:>6s} {r['samples']:>7d} {r['error']}")
            continue
        out.append(f"{'1e' + r['dilution']:>6s} {r['samples']:>7d} {r['pass']:6.1%} {r['merged']:7.1%} {r['asvs']:5d} | {r['found']:>4d}/{r['present']:<4d} "
                   f"{r['strains_found']:>3d}/{r['strains_present']:<4d} | {r['identical']:9.2%} {r['near']:6.2%} {r['none']:6.2%} | {r['outside']:>26d}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reads", type=Path)
    ap.add_argument("--runinfo", type=Path, required=True, help="the ENA table of the project (run_accession, sample_alias, read_count)")
    ap.add_argument("--asv-table", type=Path, required=True, help="the authors' ASV table, Mock_staggered-A_ASVs.csv")
    ap.add_argument("--reference", type=Path, required=True, help="mockA_reference.fasta from build_mockA_reference.py")
    ap.add_argument("--trunc", required=True, metavar="F,R", help="truncLenF,truncLenR")
    ap.add_argument("--max-ee", default="2,6", help="maxEE, one value or two (default 2,6)")
    ap.add_argument("--dilutions", nargs="+", metavar="D", help="log10 of the cells to include (default: all in the reads folder), for example 9 8 7")
    ap.add_argument("--near-id", type=float, default=0.99, help="identity at or above which a non-identical ASV counts as near a reference (default 0.99)")
    ap.add_argument("--min-overlap", type=int, default=12)
    ap.add_argument("--out", type=Path, default=Path("mockA_dilutions"))
    args = ap.parse_args()
    for path, what in ((args.runinfo, "--runinfo"), (args.asv_table, "--asv-table"), (args.reference, "--reference")):
        if not path.is_file():
            raise SystemExit(f"{what} {path}: no such file")
    if not args.reads.is_dir():
        raise SystemExit(f"{args.reads} is not a folder of reads")
    if not shutil.which("vsearch"):
        raise SystemExit("vsearch was not found on PATH (needed to score the ASVs against the reference).")
    try:
        tf, tr = (int(x) for x in args.trunc.split(","))
    except ValueError:
        raise SystemExit("--trunc must be two whole numbers, for example 257,260")

    import mock_check
    import pandas as pd
    from amplitrim.subsample import find_pairs
    from compare_methods import subset_reads
    from per_sample import load_design
    from real_data import Dada2Error, dada2_run, ee_text
    available = {s for s, _a, _b in find_pairs(args.reads)}
    groups = group_samples(load_design(args.runinfo), available)
    if args.dilutions:
        groups = {d: runs for d, runs in groups.items() if d in set(args.dilutions)}
    if not groups:
        raise SystemExit("no sample of the reads folder matches a mock dilution (names like 103-KitA-1-1-8 are expected in --runinfo)")
    counts = pd.read_csv(args.asv_table, sep=";", dtype=str).set_index("OTU_ID")
    reference = mock_check.read_fasta(args.reference)
    ee = ee_text(args.max_ee)
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    print(f"truncation {tf}/{tr}, maxEE {ee}; one DADA2 run per dilution\n")
    print(format_table([])[0])
    for d in sorted(groups, reverse=True):
        runs = groups[d]
        ref = restricted_reference(reference, counts, d)
        ref_path = args.out / f"reference_1e{d}.fasta"
        ref_path.write_text("".join(f">{n}\n{s}\n" for n, s in ref), encoding="utf-8")
        sub = subset_reads(args.reads, args.out / f"reads_1e{d}", ",".join(runs))
        asv = args.out / f"asvs_1e{d}.fasta"
        row = {"dilution": d, "samples": len(runs)}
        try:
            res = dada2_run(sub, tf, tr, ee, args.min_overlap, 0, None, asv)
        except Dada2Error as err:
            row["error"] = f"DADA2 could not run: {err}"
            rows.append(row)
            print(format_table([row])[1])
            continue
        s = mock_check.summarise(asv, ref_path, args.near_id, "vsearch")
        c, total = s["categories"], max(s["reads"], 1)
        strains = {strain_of(n) for n, _ in ref}
        found_strains = {strain_of(r) for r in s["recovered"]}
        row.update({"pass": res["filtered"], "merged": res["merged"], "asvs": res["asvs"], "present": len(ref), "found": len(s["recovered"]),
                    "strains_present": len(strains), "strains_found": len(found_strains), "identical": c["exact"][1] / total,
                    "near": c["near"][1] / total, "none": c["none"][1] / total, "outside": c["near"][0] + c["none"][0], "missing": s["missing"]})
        rows.append(row)
        print(format_table([row])[1])
        shutil.rmtree(args.out / f"reads_1e{d}", ignore_errors=True)
    lines = format_table(rows)
    lines += ["", "variants/strains: recovered exactly / present at that dilution in the authors' table. ASVs outside the reference are contaminants, chimeras or errors; at low input they are mostly contaminants."]
    for r in rows:
        if not r.get("error") and r["missing"]:
            lines.append(f"1e{r['dilution']}: present but not recovered: " + ", ".join(m.split("|")[0] + "|" + m.split("|")[-1] for m in r["missing"][:12]) + (" ..." if len(r["missing"]) > 12 else ""))
    (args.out / "mockA_dilutions.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n" + "\n".join(lines[len(rows) + 1:]))
    print(f"\nSaved: {args.out / 'mockA_dilutions.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
