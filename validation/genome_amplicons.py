"""The amplicons that a mock community should give, found by searching the strains' genomes with the two primers (an in silico PCR).

A folder holds one genome (FASTA, optionally gzipped) per strain. Each genome is searched on both strands for the forward primer and for the reverse primer
(ambiguity codes allowed, a few mismatches allowed) and every pair of sites that lies a plausible distance apart is one copy of the amplicon. The copies of a
strain that are identical between the primers are one expected variant. The output is the form that make_reference.py writes and mock_check.py reads
(>STRAIN_vN|STRAIN|N_copies), and the table adds what a table of shares of genomic DNA cannot give: the expected share of the READS, which depends on how many
copies of the gene each strain carries (ATCC's own FAQ says the 16S shares differ from the genomic ones for that reason).

The table (a tab-separated file with the columns atcc_number, name_on_atcc_page, share_percent, accession and, optionally, tier) says which genome file belongs to
which strain: the file name must contain the accession, as the files of download_msa1003_genomes.bat do.

    python validation\\genome_amplicons.py C:\\data\\longread\\pacbio_kinnex16s\\ref\\genomes --table msa1003_genomes.tsv
        --primer-f AGRGTTYGATYMTGGCTCAG --primer-r RGYTACCTTGTTACGACTT --out C:\\data\\runs\\genome_amplicons_msa1003

Writes OUT/expected_amplicons.fasta, OUT/genome_amplicons.tsv (one line per strain), OUT/expected_variants.tsv (one line per variant) and OUT/EXPECTED_AMPLICONS.md.
Needs Python 3.9+ and numpy.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_reference as mr  # noqa: E402
from longread_describe import Primer, open_text, revcomp  # noqa: E402

LUT = np.full(256, 4, dtype=np.uint8)                      # A, C, G, T -> 0..3; anything else (N, other ambiguity codes) -> 4, which never matches
for _i, _c in enumerate("ACGT"):
    LUT[ord(_c)] = _i
GENOME_SUFFIXES = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
ACCESSION = re.compile(r"GC[AF]_\d{9}\.\d+")


def read_fasta(path: Path):
    """(name, sequence in capitals) of every record of a FASTA file."""
    name, parts = None, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                fields = line[1:].split()
                name, parts = (fields[0] if fields else "unnamed"), []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts).upper()


def sites(codes: np.ndarray, primer: Primer, max_mismatch: int):
    """(positions, mismatches) of every window of the sequence that matches the primer with at most max_mismatch mismatches."""
    width = len(primer)
    n = len(codes) - width + 1
    if n <= 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    mm = np.zeros(n, dtype=np.int16)
    for j in range(width):
        mm += ~primer.allowed[j][codes[j:j + n]]
    hit = np.flatnonzero(mm <= max_mismatch)
    return hit, mm[hit].astype(np.int64)


def pairs(seq: str, forward: Primer, reverse_site: Primer, max_mismatch: int, min_length: int, max_length: int) -> list:
    """(start, end, forward mismatches, reverse mismatches) of every amplicon on this strand, the primer sites included. The reverse site is the first one that
    leaves at least min_length between the primers; the pair counts when that is at most max_length. Where two forward sites reach the same reverse site, the better one is kept."""
    codes = LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
    fi, fm = sites(codes, forward, max_mismatch)
    ri, rm = sites(codes, reverse_site, max_mismatch)
    best = {}
    for i, m in zip(fi.tolist(), fm.tolist()):
        a = i + len(forward)
        idx = int(np.searchsorted(ri, a + min_length))
        if idx >= len(ri):
            continue
        j = int(ri[idx])
        if j - a > max_length:
            continue
        if j not in best or m < best[j][1]:
            best[j] = (i, m, int(rm[idx]))
    return [(i, j + len(reverse_site), mf, mr_) for j, (i, mf, mr_) in sorted(best.items())]


def scan_genome(path: Path, primer_f: str, primer_r: str, max_mismatch: int, min_length: int, max_length: int) -> dict:
    """Every amplicon of one genome, in the orientation of the gene, with the primer sites included."""
    forward, reverse_site = Primer(primer_f), Primer(revcomp(primer_r.upper()))
    hits, contigs, length = [], 0, 0
    for name, seq in read_fasta(path):
        contigs += 1
        length += len(seq)
        for strand, s in (("+", seq), ("-", revcomp(seq))):
            for start, end, mf, mr_ in pairs(s, forward, reverse_site, max_mismatch, min_length, max_length):
                hits.append({"contig": name, "strand": strand, "start": start, "mm_f": mf, "mm_r": mr_, "segment": s[start:end]})
    return {"hits": hits, "contigs": contigs, "length": length}


def find_genomes(folder: Path) -> list:
    return [p for p in sorted(Path(folder).iterdir()) if p.is_file() and p.name.lower().endswith(GENOME_SUFFIXES)]


def read_table(path: Path) -> dict:
    """{accession: row} of the strain table."""
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["accession"].strip(): r for r in csv.DictReader(fh, delimiter="\t") if r.get("accession", "").strip()}


def label_of(row: dict) -> str:
    return re.sub(r"[^A-Za-z0-9_\-]", "_", f"{row['name_on_atcc_page'].strip().replace(' ', '_')}_{row['atcc_number'].replace(' ', '')}")


def parse_variants(fasta_text: str) -> dict:
    """{strain: [(variant name, copies, length), ...]} from the FASTA that make_reference.build writes."""
    out, name, seq = {}, None, []
    def flush():
        if name is not None:
            ident, strain, copies = name.split("|")
            out.setdefault(strain, []).append((ident, int(copies.split("_")[0]), len("".join(seq))))
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            flush()
            name, seq = line[1:], []
        else:
            seq.append(line.strip())
    flush()
    return out


def expected_shares(strains: dict) -> None:
    """Adds 'read_share' to each strain: its share of genomic DNA times its copies, as a share of the total. Strains without a share or without copies stay at None."""
    weights = {k: (v["share"] * v["copies"]) for k, v in strains.items() if v["share"] is not None and v["copies"]}
    total = sum(weights.values())
    for k, v in strains.items():
        v["read_share"] = weights[k] / total if k in weights and total else None


def run(folder: Path, table_path, primer_f: str, primer_r: str, out: Path, max_mismatch: int = 3, min_length: int = 1000, max_length: int = 2000,
        depth: int = 30000, log=print) -> dict:
    files = find_genomes(folder)
    if not files:
        raise SystemExit(f"No genome files ({', '.join(GENOME_SUFFIXES)}) in {folder}")
    table = read_table(Path(table_path)) if table_path else {}
    records, strains, warnings, used = [], {}, [], set()
    for path in files:
        m = ACCESSION.search(path.name)
        row = table.get(m.group(0)) if m else None
        if table and row is None:
            warnings.append(f"{path.name}: no row of the table has this accession, so the file is named by itself and has no ATCC share")
        label = label_of(row) if row else re.sub(r"[^A-Za-z0-9_\-]", "_", re.sub(r"\.(fna|fa|fasta)(\.gz)?$", "", path.name))
        if row:
            used.add(row["accession"].strip())
        scan = scan_genome(path, primer_f, primer_r, max_mismatch, min_length, max_length)
        for h in scan["hits"]:
            records.append((f"{label}.fna", f"{h['contig']}_{h['strand']}_{h['start']}", h["segment"]))
        tier = (row or {}).get("tier", "")
        share = float(row["share_percent"]) if row and row.get("share_percent") not in (None, "") else None
        strains[label] = {"file": path.name, "tier": tier, "share": share, "contigs": scan["contigs"], "length": scan["length"], "copies": len(scan["hits"]),
                          "mm_f": max((h["mm_f"] for h in scan["hits"]), default=None), "mm_r": max((h["mm_r"] for h in scan["hits"]), default=None), "variants": []}
        log(f"{label}: {scan['contigs']} sequence(s), {scan['length']:,} nt, {len(scan['hits'])} copy(ies) of the amplicon")
        if not scan["hits"]:
            warnings.append(f"{label}: NO amplicon found; the primers do not match this genome within {max_mismatch} mismatches at a distance of {min_length} to {max_length} nt")
        if "draft" in tier.lower():
            warnings.append(f"{label}: a draft genome; copies of the rRNA operon may be collapsed, so its copy number and its expected read share are uncertain")
        if strains[label]["copies"] and max(strains[label]["mm_f"], strains[label]["mm_r"]) >= 2:
            warnings.append(f"{label}: a primer site has up to {max(strains[label]['mm_f'], strains[label]['mm_r'])} mismatches, so this strain may amplify less than the others")
    for acc, row in table.items():
        if acc not in used:
            warnings.append(f"{label_of(row)} ({acc}): no genome file with this accession in the folder")
    fasta_text, report = mr.build(records, primer_f, primer_r, max_mismatch, "file", None, min_length=50)
    variants = parse_variants(fasta_text)
    for label, v in strains.items():
        v["variants"] = variants.get(label, [])
        if v["variants"]:
            v["copies"] = sum(c for _, c, _ in v["variants"])
    expected_shares(strains)
    out.mkdir(parents=True, exist_ok=True)
    (out / "expected_amplicons.fasta").write_text(fasta_text, encoding="utf-8")
    write_tables(strains, out, depth)
    write_report(strains, report, warnings, out, primer_f, primer_r, max_mismatch, min_length, max_length, depth)
    log(f"Saved: {out / 'EXPECTED_AMPLICONS.md'}")
    return {"strains": strains, "warnings": warnings, "fasta": fasta_text, "report": report}


def write_tables(strains: dict, out: Path, depth: int) -> None:
    with open(out / "genome_amplicons.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["strain", "file", "tier", "sequences", "genome_length", "copies", "variants", "length_min", "length_max", "max_mismatches_forward", "max_mismatches_reverse",
                    "share_genomic_percent", "expected_read_share_percent", f"expected_reads_per_{depth}"])
        for label, v in strains.items():
            lengths = [n for _, _, n in v["variants"]]
            rs = v["read_share"]
            w.writerow([label, v["file"], v["tier"], v["contigs"], v["length"], v["copies"], len(v["variants"]), min(lengths, default=""), max(lengths, default=""), v["mm_f"] if v["mm_f"] is not None else "",
                        v["mm_r"] if v["mm_r"] is not None else "", v["share"] if v["share"] is not None else "", f"{100 * rs:.4f}" if rs is not None else "", f"{rs * depth:.1f}" if rs is not None else ""])
    with open(out / "expected_variants.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["variant", "strain", "copies", "length", "expected_read_share_percent"])
        for label, v in strains.items():
            for name, copies, n in v["variants"]:
                rs = v["read_share"]
                w.writerow([name, label, copies, n, f"{100 * rs * copies / v['copies']:.4f}" if rs is not None and v["copies"] else ""])


def write_report(strains: dict, report: list, warnings: list, out: Path, primer_f: str, primer_r: str, max_mismatch: int, min_length: int, max_length: int, depth: int) -> None:
    found = [k for k, v in strains.items() if v["variants"]]
    lines = ["# Expected amplicons of the mock community, from the genomes", "",
             f"Primers {primer_f} and {primer_r}, at most {max_mismatch} mismatches each, between {min_length:,} and {max_length:,} nt between the primers; "
             f"{len(strains)} genome(s), {len(found)} with at least one amplicon, {sum(len(strains[k]['variants']) for k in found)} expected variants "
             f"from {sum(strains[k]['copies'] for k in found)} copies of the gene.", "",
             "The expected read share is the strain's share of genomic DNA times its number of copies of the gene, as a share of the total over all strains; "
             f"the last column is that share of {depth:,} reads.", "",
             "| strain | tier | genome (sequences) | copies | variants | length (nt) | primer mismatches (5', 3') | genomic share % | expected read share % | reads per " + f"{depth:,} |", "|---|---|---|---|---|---|---|---|---|---|"]
    for label, v in sorted(strains.items(), key=lambda kv: -(kv[1]["read_share"] or 0)):
        lengths = [n for _, _, n in v["variants"]]
        span = "-" if not lengths else (str(lengths[0]) if min(lengths) == max(lengths) else f"{min(lengths)}-{max(lengths)}")
        mm = "-" if v["mm_f"] is None else f"{v['mm_f']}, {v['mm_r']}"
        rs = v["read_share"]
        lines.append(f"| {label} | {v['tier'] or '-'} | {v['length']:,} ({v['contigs']}) | {v['copies']} | {len(v['variants'])} | {span} | {mm} | "
                     f"{'-' if v['share'] is None else format(v['share'], 'g')} | {'-' if rs is None else f'{100 * rs:.3f}'} | {'-' if rs is None else f'{rs * depth:.1f}'} |")
    if warnings:
        lines += ["", "## To check", ""] + [f"- {w}" for w in warnings]
    lines += ["", "## What make_reference.py reported", "", "```"] + report + ["```"]
    (out / "EXPECTED_AMPLICONS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("folder", type=Path, help="folder of genome files, one per strain")
    ap.add_argument("--table", type=Path, help="tab-separated strain table (atcc_number, name_on_atcc_page, share_percent, accession, tier)")
    ap.add_argument("--primer-f", required=True)
    ap.add_argument("--primer-r", required=True, help="5'->3' as ordered")
    ap.add_argument("--out", type=Path, default=Path("genome_amplicons"))
    ap.add_argument("--max-mismatch", type=int, default=3)
    ap.add_argument("--min-length", type=int, default=1000, help="fewest nt between the primers (default 1000)")
    ap.add_argument("--max-length", type=int, default=2000, help="most nt between the primers (default 2000)")
    ap.add_argument("--depth", type=int, default=30000, help="reads per sample for the column of expected reads (default 30000)")
    a = ap.parse_args(argv)
    run(a.folder, a.table, a.primer_f, a.primer_r, a.out, a.max_mismatch, a.min_length, a.max_length, a.depth)
    return 0


if __name__ == "__main__":
    sys.exit(main())
