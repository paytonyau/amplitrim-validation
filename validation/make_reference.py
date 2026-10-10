#!/usr/bin/env python
"""Expected sequences of a mock community, in the form the accuracy tests use, from full-length reference rRNA sequences.

A reference such as ZymoBIOMICS.STD.refseq (a folder or zip with the 16S and 18S sequences of the community) holds whole genes. The reads, after primer removal, cover only
the amplicon between the primers, so each reference sequence is cut between the primers of the dataset (an in silico PCR, primers with ambiguity codes and a few mismatches
allowed, either orientation), and copies that are identical inside the amplicon are collapsed into one expected variant. Sequences without both primer sites (the 18S
sequences of the yeasts, for a 16S primer pair) are left out and reported.

The output names every variant >STRAINNAME_vN|STRAINNAME|N copies, the form mock_check.py reads to report variants and strains separately.

    python validation\\make_reference.py ZymoBIOMICS.STD.refseq.v2.zip --primer-f AGAGTTTGATYMTGGCTCAG --primer-r ATTACCGCGGCTGCTGG --out mockB_reference.fasta

or, for a dataset of the catalog, python validation\\check_build.py --make-reference mockB=ZymoBIOMICS.STD.refseq.v2.zip (it takes the primers of the dataset and registers the result).
"""
import argparse
import re
import sys
import zipfile
from pathlib import Path

CODES = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
FASTA_SUFFIXES = (".fa", ".fas", ".fasta", ".fna", ".ffn", ".txt")
GENE_PATH = re.compile(r"ssrrna|16s|rrna", re.I)          # a reference folder often holds the genomes as well as the rRNA genes: when some files are named for the genes, only those are read
GENERIC_NAMES = {"", "consensus", "reference", "representative", "sequence", "seq", "unnamed"}
MAX_GENE_LENGTH = 10_000                                   # a record longer than this is a genome or a contig, not an rRNA gene


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def mismatches(window: str, primer: str) -> int:
    return sum(1 for base, code in zip(window, primer) if base not in CODES.get(code, ""))


def find_site(seq: str, primer: str, max_mismatch: int, start: int = 0):
    """(mismatches, position) of the best match of the primer (ambiguity codes allowed) at or after `start`, the leftmost of equally good ones; None when none has at most max_mismatch."""
    best = None
    for i in range(start, len(seq) - len(primer) + 1):
        mm = mismatches(seq[i:i + len(primer)], primer)
        if mm <= max_mismatch and (best is None or mm < best[0]):
            best = (mm, i)
            if mm == 0:
                break
    return best


def read_fasta_text(text: str) -> list:
    records, header, parts = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header, parts = line[1:].strip(), []
        elif line.strip():
            parts.append(re.sub(r"[^A-Za-z]", "", line).upper().replace("U", "T"))
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def prefer_gene_files(names: list) -> list:
    """The files named for the genes (a path with ssrRNA, 16S or rRNA) when there are any, otherwise all of them."""
    genes = [n for n in names if GENE_PATH.search(str(n))]
    return genes or names


def read_inputs(paths: list) -> list:
    """[(source file name, header, sequence)] from FASTA files, folders (any depth) and zip archives."""
    out = []
    for p in map(Path, paths):
        if not p.exists():
            raise SystemExit(f"{p} does not exist: give the real location of the file or folder")
        if p.is_dir():
            files = prefer_gene_files(sorted(f for f in p.rglob("*") if f.is_file() and f.suffix.lower() in FASTA_SUFFIXES))
            for f in files:
                out += [(f.name, h, s) for h, s in read_fasta_text(f.read_text(encoding="utf-8", errors="replace"))]
        elif p.suffix.lower() == ".zip":
            with zipfile.ZipFile(p) as z:
                members = prefer_gene_files([n for n in sorted(z.namelist()) if n.lower().endswith(FASTA_SUFFIXES) and not n.endswith("/")])
                for name in members:
                    out += [(Path(name).name, h, s) for h, s in read_fasta_text(z.read(name).decode("utf-8", errors="replace"))]
        elif p.is_file():
            out += [(p.name, h, s) for h, s in read_fasta_text(p.read_text(encoding="utf-8", errors="replace"))]
        else:
            raise SystemExit(f"{p} does not exist")
    return out


def strain_name(header: str, source: str, mode: str = "header") -> str:
    """The strain of a sequence: the first word of its header without the gene and copy number (Escherichia_coli_16S_3 -> Escherichia_coli), or the name of its file."""
    stem = re.sub(r"\.(fa|fas|fasta|fna|ffn|txt)$", "", source, flags=re.I)
    if mode == "file":
        name = stem
    else:
        name = re.split(r"[|\s]", header)[0]
        name = re.sub(r"^(consensus|representative|reference)[_.\-]+(?=[A-Za-z])", "", name, flags=re.I)         # Consensus_Imtechella_halotolerans_16S -> Imtechella_halotolerans
        name = re.sub(r"[_.\-]?(16S|18S|28S|SSU|rRNA|ssrRNA|rrn).*$", "", name, flags=re.I)
        name = re.sub(r"[_.\-]?\d+$", "", name)
        if name.lower() in GENERIC_NAMES:                                                                         # a header that names no organism: the file does, cleaned the same way
            name = re.sub(r"[_.\-]?\d+$", "", re.sub(r"[_.\-]?(16S|18S|28S|SSU|rRNA|ssrRNA|rrn).*$", "", stem, flags=re.I))
    name = re.sub(r"[^A-Za-z0-9_\-]", "_", (name or stem).replace(".", "_")).strip("_")                   # Allobacillus.halotolerans and Allobacillus_halotolerans are one strain
    name = re.sub(r"_+", "_", name)
    return name or "unnamed"


def best_mismatches(seq: str, primer: str) -> int:
    """The fewest mismatches the primer has anywhere in the sequence (a primer longer than the sequence counts as all mismatches)."""
    if len(seq) < len(primer):
        return len(primer)
    return min(mismatches(seq[i:i + len(primer)], primer) for i in range(len(seq) - len(primer) + 1))


def why_not(seq: str, primer_f: str, primer_r: str, max_mismatch: int) -> str:
    """For a sequence that could not be cut: which primer site is missing and how close the best site came, in the orientation that fits best."""
    reverse_site = revcomp(primer_r.upper())
    best = None
    for s in (seq, revcomp(seq)):
        f, r = best_mismatches(s, primer_f.upper()), best_mismatches(s, reverse_site)
        if best is None or f + r < best[0] + best[1]:
            best = (f, r)
    f, r = best
    parts = [f"forward primer {'found' if f <= max_mismatch else f'not found (best site {f} mismatches)'}", f"reverse primer {'found' if r <= max_mismatch else f'not found (best site {r} mismatches)'}"]
    return "; ".join(parts) + (" (the sequence may start or end inside the amplicon)" if (f <= max_mismatch) != (r <= max_mismatch) else "")


def amplicon(seq: str, primer_f: str, primer_r: str, max_mismatch: int = 3):
    """(sequence between the primers, orientation, forward mismatches, reverse mismatches) or None."""
    reverse_site = revcomp(primer_r.upper())
    for orientation, s in (("+", seq), ("-", revcomp(seq))):
        f = find_site(s, primer_f.upper(), max_mismatch)
        if f is None:
            continue
        r = find_site(s, reverse_site, max_mismatch, start=f[1] + len(primer_f))
        if r is None:
            continue
        return s[f[1] + len(primer_f):r[1]], orientation, f[0], r[0]
    return None


def build(records: list, primer_f: str, primer_r: str, max_mismatch: int = 3, strain_mode: str = "header", exclude: str = None, min_length: int = 50):
    """(fasta text, report lines). Identical amplicons of one strain are one variant; the report says what was kept, collapsed and left out."""
    by_strain, consensus, skipped = {}, {}, []
    for source, header, seq in records:
        if len(seq) > MAX_GENE_LENGTH:
            skipped.append((header, f"{len(seq):,} nt: a genome or contig, not an rRNA gene"))
            continue
        if exclude and re.search(exclude, f"{source} {header}", re.I):
            skipped.append((header, "excluded by --exclude"))
            continue
        hit = amplicon(seq, primer_f, primer_r, max_mismatch)
        if hit is None:
            skipped.append((header, "no site for both primers: " + why_not(seq, primer_f, primer_r, max_mismatch)))
            continue
        inner = hit[0]
        if len(inner) < min_length:
            skipped.append((header, f"only {len(inner)} nt between the primers"))
            continue
        target = consensus if re.match(r"(?i)consensus", header.split()[0]) else by_strain
        target.setdefault(strain_name(header, source, strain_mode), {}).setdefault(inner, []).append(header)
    for strain, variants in consensus.items():                                                                   # a consensus of a strain's copies was never seen in a genome: it is left out when real copies exist
        if strain in by_strain:
            for inner, heads in variants.items():
                skipped += [(h, f"a consensus of {strain}'s other sequences, not an observed variant") for h in heads]
        else:
            by_strain[strain] = variants
    fasta, report, owners = [], [], {}
    total_seq = total_var = 0
    for strain in sorted(by_strain):
        variants = sorted(by_strain[strain].items(), key=lambda kv: (-len(kv[1]), kv[0]))
        copies = sum(len(h) for _, h in variants)
        total_seq, total_var = total_seq + copies, total_var + len(variants)
        heads = sorted({h.split()[0][:48] for _, hs in variants for h in hs})
        report.append(f"  {strain}: {copies} sequence(s) -> {len(variants)} distinct amplicon(s) of {', '.join(str(len(s)) for s, _ in variants)} nt   [from {', '.join(heads[:3])}{', ...' if len(heads) > 3 else ''}]")
        for i, (inner, heads) in enumerate(variants, 1):
            fasta.append(f">{strain}_v{i}|{strain}|{len(heads)}_copies\n{inner}")
            owners.setdefault(inner, []).append(strain)
    shared = sorted({s for strains in owners.values() if len(set(strains)) > 1 for s in set(strains)})
    lines = [f"{len(records)} sequences read; {total_seq} cut between the primers ({primer_f} / {primer_r}, at most {max_mismatch} mismatches), {len(skipped)} left out;",
             f"{total_var} expected variants of {len(by_strain)} strains:"] + report
    if shared:
        lines.append(f"  note: strains with an amplicon identical to another strain's (scored on the same reads, so one of them cannot be told apart): {', '.join(shared)}")
    if skipped:
        lines.append("left out:\n" + "\n".join(f"  {h.split()[0][:48]}: {why}" for h, why in skipped[:12]) + (f"\n  and {len(skipped) - 12} more" if len(skipped) > 12 else ""))
    return "\n".join(fasta) + ("\n" if fasta else ""), lines


def strain_of_name(name: str) -> str:
    """The strain of an expected sequence named ID|STRAIN|label; any other name is its own strain."""
    parts = name.split("|")
    return parts[1] if len(parts) >= 3 else name


def compare_references(candidate_fasta: Path, reference_fasta: Path, min_id: float = 0.99, vsearch: str = "vsearch", labels: dict = None) -> list:
    """Report lines: which sequences of one reference have a counterpart (at least min_id identity, either strand) in another, and which strains on each side have none.
    Used to ask whether a published reference (the Zymo file) describes the same community as the expected sequences already in use."""
    import shutil
    import subprocess
    import tempfile
    if not shutil.which(vsearch):
        raise SystemExit("vsearch is needed for the comparison and was not found on the PATH")
    labels = labels or {}
    cand = [(h.split()[0], q) for h, q in read_fasta_text(Path(candidate_fasta).read_text(encoding="utf-8"))]
    ref = [(h.split()[0], q) for h, q in read_fasta_text(Path(reference_fasta).read_text(encoding="utf-8"))]
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "hits.tsv"
        res = subprocess.run([vsearch, "--usearch_global", str(candidate_fasta), "--db", str(reference_fasta), "--id", str(min_id), "--strand", "both", "--maxaccepts", "0", "--maxrejects", "0",
                              "--userout", str(out), "--userfields", "query+target+id", "--quiet"], capture_output=True, text=True)
        if res.returncode != 0:
            raise SystemExit(f"vsearch failed: {res.stderr.strip()[-300:]}")
        matched, reached = set(), set()
        for line in out.read_text().splitlines():
            q, t, _ident = line.split("\t")
            matched.add(q.split()[0])
            reached.add(t.split()[0])
    def name_of(strain: str) -> str:
        return f"{strain} ({labels[strain]})" if labels.get(strain) else strain
    cand_strains, ref_strains = {}, {}
    for n, _ in cand:
        cand_strains.setdefault(strain_of_name(n), []).append(n in matched)
    for n, _ in ref:
        ref_strains.setdefault(strain_of_name(n), []).append(n in reached)
    lost_cand = sorted(s for s, flags in cand_strains.items() if not any(flags))
    lost_ref = sorted(s for s, flags in ref_strains.items() if not any(flags))
    lines = [f"{len(cand)} sequences of {len(cand_strains)} strains in the candidate reference; {len(matched)} of them match a variant of {Path(reference_fasta).name} at at least {min_id:.0%}, {len(cand) - len(matched)} do not.",
             f"{len(ref)} expected variants of {len(ref_strains)} strains in {Path(reference_fasta).name}; {len(reached)} are matched by a candidate sequence, {len(ref) - len(reached)} are not.",
             "strains of the candidate with no counterpart: " + (", ".join(name_of(s) for s in lost_cand) if lost_cand else "none"),
             "strains of the expected sequences with no counterpart in the candidate: " + (", ".join(name_of(s) for s in lost_ref) if lost_ref else "none")]
    if lost_ref or lost_cand:
        lines.append("The two describe different communities, or the same one only in part: the candidate can check the strains they share and cannot replace the expected sequences.")
    else:
        lines.append("Every strain on each side has a counterpart: the candidate describes the same community.")
    return lines


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="FASTA files, folders or zip archives of full-length reference sequences")
    ap.add_argument("--primer-f")
    ap.add_argument("--primer-r", help="the reverse primer as written 5' to 3' (its reverse complement is searched)")
    ap.add_argument("--out", help="the FASTA of expected sequences to write")
    ap.add_argument("--list", action="store_true", help="only list the files, headers and lengths found in the input (no primers needed)")
    ap.add_argument("--max-mismatch", type=int, default=3)
    ap.add_argument("--strain-from", choices=("header", "file"), default="header", help="name the strains from the sequence headers (default) or from the file names")
    ap.add_argument("--exclude", metavar="REGEX", help="leave out sequences whose file name or header matches (for example 18S|ITS)")
    args = ap.parse_args(argv)
    if args.list:
        for source, header, seq in read_inputs(args.inputs):
            print(f"{source}\t{header.split()[0][:60]}\t{len(seq)} nt\t-> strain {strain_name(header, source, args.strain_from)}")
        return 0
    if not (args.primer_f and args.primer_r and args.out):
        raise SystemExit("--primer-f, --primer-r and --out are needed (or --list to see what the input holds)")
    fasta, lines = build(read_inputs(args.inputs), args.primer_f, args.primer_r, args.max_mismatch, args.strain_from, args.exclude)
    print("\n".join(lines))
    if not fasta:
        print("no sequence has sites for both primers: check the primers and the orientation of the reference", file=sys.stderr)
        return 1
    Path(args.out).write_text(fasta, encoding="utf-8")
    print(f"written: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
