#!/usr/bin/env python3
"""Does maxEE change what DADA2 finds? Runs real DADA2 at ONE truncation with several maxEE values and, if you give the expected
sequences of a mock community, scores every result against them.

    python validation/max_ee_sweep.py C:\\data\\results\\mockA\\ready --trunc 275,242 --pairs 6000 --reference expected.fasta ^
           --max-ee 1,1 2,2 2,3 2,5 2,8 --out C:\\data\\sweep

Use the `ready` folder real_data.py made (primers already removed). Each maxEE is one value ("2") or a forward,reverse pair ("2,5").
Without --reference the table shows retention and the ASV count; with it, also how many expected sequences were found, the share of
reads in ASVs identical to them, and how many ASVs are not. The truncation is fixed so that only maxEE differs between rows.
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))


def abundance_summary(asv_fasta) -> dict:
    """What the ASVs' read counts look like, without any reference: the share of reads in the 15 most abundant ASVs, the number of
    ASVs with at least 0.1% of the reads, and the share of reads in ASVs below 0.1% (errors, contaminants, rare variants). A mock
    community of N strains should show about N dominant ASVs and a small tail."""
    from mock_check import read_asvs, size_of
    sizes = sorted((size_of(n) for n, _ in read_asvs(asv_fasta)), reverse=True)
    total = max(sum(sizes), 1)
    return {"top15": sum(sizes[:15]) / total, "n_01": sum(1 for x in sizes if x / total >= 0.001),
            "tail": sum(x for x in sizes if x / total < 0.001) / total}


def format_row(ee: str, d: dict, score: dict = None, spectrum: dict = None) -> str:
    row = f"{ee:>8s} {d['filtered']:6.1%} {d['merged']:7.1%} {d['nonchim']:8.1%} {d['asvs']:5d}"
    if spectrum:
        row += f" | {spectrum['top15']:7.2%} {spectrum['n_01']:6d} {spectrum['tail']:6.2%}"
    if score:
        c, tot = score["categories"], max(score["reads"], 1)
        row += (f" | {c['exact'][1] / tot:9.2%} {c['near'][1] / tot:6.2%} {c['none'][1] / tot:6.2%} | "
                f"{len(score['recovered']):>4d} of {score['references']:<3d} {c['near'][0] + c['none'][0]:>5d}")
        if score.get("strains_total"):
            row += f" {score['strains_found']:>3d}/{score['strains_total']:<3d}"
    return row


def header(with_reference: bool, with_spectrum: bool = True, with_strains: bool = False) -> str:
    h = f"{'maxEE':>8s} {'pass':>6s} {'merged':>7s} {'no chim':>8s} {'ASVs':>5s}"
    if with_spectrum:
        h += f" | {'top 15':>7s} {'>=0.1%':>6s} {'tail':>6s}"
    if with_reference:
        h += f" | {'identical':>9s} {'near':>6s} {'none':>6s} | {'expected found':>14s} {'extra':>5s}" + (f" {'strains':>7s}" if with_strains else "")
    return h


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reads", type=Path, help="folder of primer-free paired FASTQs (the `ready` folder of real_data.py)")
    ap.add_argument("--trunc", required=True, metavar="F,R", help="truncLenF,truncLenR, the same for every row")
    ap.add_argument("--max-ee", nargs="+", default=["1,1", "2,2", "2,3", "2,5", "2,8"], metavar="EE")
    ap.add_argument("--pairs", type=int, default=0, metavar="N", help="give DADA2 at most N evenly spaced pairs per sample (default all)")
    ap.add_argument("--min-overlap", type=int, default=12)
    ap.add_argument("--reference", type=Path, help="FASTA of the expected sequences (a mock community); name variants of one strain >ID|STRAIN|label")
    ap.add_argument("--near-id", type=float, default=0.99, help="identity at or above which a non-identical ASV counts as near a reference (default 0.99, about 4-5 bases)")
    ap.add_argument("--out", type=Path, default=Path("max_ee_sweep"))
    args = ap.parse_args()
    from real_data import Dada2Error, dada2_run, ee_text, thin_reads
    tf, tr = (int(x) for x in args.trunc.split(","))
    values = [ee_text(v) for v in args.max_ee]
    args.out.mkdir(parents=True, exist_ok=True)
    reads = thin_reads(args.reads, args.out / "reads", args.pairs) if args.pairs else args.reads
    if args.reference:
        import shutil
        import mock_check
        if not shutil.which("vsearch"):
            raise SystemExit("vsearch was not found on PATH (needed to score the ASVs against the reference).")
    has_strains = bool(args.reference) and mock_check.strain_counts([n.split()[0] for n, _ in mock_check.read_fasta(args.reference)], [])[1] is not None
    lines = [f"truncation {tf}/{tr}, one DADA2 run per maxEE", "", header(bool(args.reference), True, has_strains)]
    print("\n".join(lines))
    for ee in values:
        asv = args.out / f"asvs_maxEE_{ee.replace(',', '_')}.fasta"
        try:
            d = dada2_run(reads, tf, tr, ee, args.min_overlap, 0, None, asv)
        except Dada2Error as e:
            line = f"{ee:>8s} DADA2 could not run: {e}"
            print(line)
            lines.append(line)
            continue
        score = mock_check.summarise(asv, args.reference, args.near_id, vsearch="vsearch") if args.reference else None
        line = format_row(ee, d, score, abundance_summary(asv))
        print(line)
        lines.append(line)
    (args.out / "max_ee_sweep.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nTable saved to {args.out / 'max_ee_sweep.txt'}; the ASVs of each row are in the same folder.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
