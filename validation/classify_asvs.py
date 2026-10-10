#!/usr/bin/env python
"""What are the ASVs of a mock community run that are not the expected sequences? Chimeras, contaminants, error variants, or something else.

  python validation/classify_asvs.py ASVS.fasta --reference mockA_reference.fasta [--present mockA_high_input_reference.fasta]
         [--asv-table Mock_staggered-A_ASVs.csv --controls-regex "NK|Wasser"] [--out classes.tsv]

ASVS.fasta is the output of DADA2 as the validation scripts write it (>ASV1;size=900). Every ASV falls in the first class that fits:

  expected          identical to an expected sequence (to part of one, as in mock_check.py) present in the samples
  listed elsewhere  identical to an expected sequence that is not present in these samples (--present narrows what counts as present)
  contaminant       identical to (or within 2 edits of) an ASV that the authors' table finds in the negative controls (--asv-table, columns
                    matching --controls-regex), with the taxonomy of the table when it has one
  chimera           flagged by vsearch --uchime_denovo (parents among the ASVs) or --uchime_ref (parents among the expected sequences), with
                    parents that differ from each other (identity below 99.5%), so that two variants of one strain do not make a "chimera"
  error variant     within 3 edits of an ASV at least 5 times as abundant, or of an expected sequence
  unlisted variant  97% to 99.9% identical to an expected sequence but none of the above: a real variant the reference does not list, or an error
  unexplained       nothing close to it

The classes are evidence, not proof: an error variant can be a real variant, a chimera call can be wrong, and a contaminant call needs the
negative controls. Needs vsearch on PATH and rapidfuzz (conda install -c conda-forge rapidfuzz).
"""

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

CLASSES = ["expected", "listed elsewhere", "contaminant", "chimera", "error variant", "unlisted variant", "unexplained"]
SIZE = re.compile(r";?size=(\d+);?")


def read_sized(path) -> list:
    """[{'id', 'size', 'seq'}] from a FASTA with ;size=N in the headers (size 1 when there is none)."""
    out, name, seq = [], None, []
    def flush():
        if name is not None:
            m = SIZE.search(name)
            out.append({"id": SIZE.sub("", name).strip().split()[0], "size": int(m[1]) if m else 1, "seq": "".join(seq).upper()})
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            flush()
            name, seq = line[1:].strip(), []
        elif line.strip():
            seq.append(line.strip())
    flush()
    return out


def write_sized(records: list, path: Path) -> None:
    ordered = sorted(records, key=lambda r: -r["size"])
    Path(path).write_text("".join(f">{r['id']};size={r['size']}\n{r['seq']}\n" for r in ordered), encoding="utf-8")


def best_hits(query: Path, db: Path, min_id: float, vsearch: str = "vsearch") -> dict:
    """{query id: (target id, identity in percent)} for the best hit at or above min_id (identity as vsearch reports it, terminal gaps not counted)."""
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        hits = Path(tmp) / "hits.tsv"
        res = subprocess.run([vsearch, "--usearch_global", str(query), "--db", str(db), "--id", str(min_id), "--strand", "both", "--maxaccepts", "0", "--maxrejects", "0",
                              "--userout", str(hits), "--userfields", "query+target+id", "--quiet"], capture_output=True, text=True)
        if res.returncode != 0:
            raise SystemExit(f"vsearch failed: {res.stderr.strip()[-300:]}")
        for line in hits.read_text().splitlines():
            q, t, i = line.split("\t")
            q, t = SIZE.sub("", q).strip(), SIZE.sub("", t).strip()
            if q not in out or float(i) > out[q][1]:
                out[q] = (t, float(i))
    return out


def uchime_flags(asvs: list, reference: Path = None, vsearch: str = "vsearch", max_parent_identity: float = 99.5) -> dict:
    """{asv id: 'denovo' | 'reference'} for the ASVs vsearch calls chimeric with parents that differ from each other."""
    flags = {}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fa = tmp / "asvs.fasta"
        write_sized(asvs, fa)
        runs = [("denovo", ["--uchime_denovo", str(fa), "--sizein"])]
        if reference is not None:
            runs.append(("reference", ["--uchime_ref", str(fa), "--db", str(reference), "--sizein"]))
        for kind, args in runs:
            out = tmp / f"{kind}.tsv"
            res = subprocess.run([vsearch] + args + ["--uchimeout", str(out), "--quiet"], capture_output=True, text=True)
            if res.returncode != 0:
                raise SystemExit(f"vsearch --{args[0].lstrip('-')} failed: {res.stderr.strip()[-300:]}")
            for line in out.read_text().splitlines():
                f = line.split("\t")
                if len(f) >= 18 and f[17].strip() == "Y" and float(f[8]) < max_parent_identity:
                    flags.setdefault(SIZE.sub("", f[1]).strip(), kind)
    return flags


def load_controls(asv_table: Path, controls_regex: str) -> list:
    """[{'seq', 'id', 'reads', 'taxon'}] for the ASVs of the authors' table with reads in the negative-control columns."""
    import pandas as pd
    t = pd.read_csv(asv_table, sep=";", dtype=str)
    cols = [c for c in t.columns if re.search(controls_regex, c, re.I) and c not in ("OTU_ID", "Sequence")]
    if not cols:
        raise SystemExit(f"no column of {asv_table} matches --controls-regex {controls_regex!r}")
    tax = [c for c in ("Genus", "Family", "Order", "Class", "Phylum") if c in t.columns]
    out = []
    for _, r in t.iterrows():
        reads = sum(float(r[c]) for c in cols if str(r[c]) not in ("nan", ""))
        if reads > 0:
            taxon = next((str(r[c]) for c in tax if str(r[c]) not in ("nan", "")), "")
            out.append({"seq": str(r["Sequence"]).upper(), "id": str(r["OTU_ID"]), "reads": int(reads), "taxon": taxon})
    return out


def classify(asv_fasta: Path, reference: Path, present: Path = None, asv_table: Path = None, controls_regex: str = "NK|Wasser",
             error_edits: int = 3, error_ratio: float = 5.0, vsearch: str = "vsearch") -> list:
    """One dict per ASV: id, size, length, class, detail."""
    try:
        from rapidfuzz.distance import Levenshtein
    except ImportError:
        raise SystemExit("rapidfuzz is needed: conda install -c conda-forge rapidfuzz (or pip install rapidfuzz)")
    import mock_check
    asvs = read_sized(asv_fasta)
    refs = [(n.split()[0], s.upper()) for n, s in mock_check.read_fasta(reference)]
    present_ids = {n.split()[0] for n, _ in mock_check.read_fasta(present)} if present else {n for n, _ in refs}
    hits = best_hits(Path(asv_fasta), Path(reference), 0.97, vsearch)
    chim = uchime_flags(asvs, reference, vsearch)
    controls = load_controls(asv_table, controls_regex) if asv_table else []
    rows = []
    for a in asvs:
        row = {"id": a["id"], "size": a["size"], "length": len(a["seq"]), "class": "", "detail": ""}
        hit = hits.get(a["id"])
        if hit and hit[1] >= 99.999:
            row["class"] = "expected" if hit[0] in present_ids else "listed elsewhere"
            row["detail"] = hit[0]
            rows.append(row)
            continue
        near_control = next(((c, Levenshtein.distance(a["seq"], c["seq"], score_cutoff=2)) for c in controls if Levenshtein.distance(a["seq"], c["seq"], score_cutoff=2) <= 2), None)
        if near_control:
            c, d = near_control
            row["class"], row["detail"] = "contaminant", f"{c['id']} in the negative controls ({c['reads']} reads, {d} edits)" + (f"; {c['taxon']}" if c["taxon"] else "")
        elif a["id"] in chim:
            row["class"], row["detail"] = "chimera", f"vsearch uchime {chim[a['id']]}"
        else:
            parent = None
            for other in asvs:
                if other["id"] != a["id"] and other["size"] >= error_ratio * a["size"] and Levenshtein.distance(a["seq"], other["seq"], score_cutoff=error_edits) <= error_edits:
                    parent = (other["id"], Levenshtein.distance(a["seq"], other["seq"]), f"ASV of {other['size']} reads")
                    break
            if parent is None:
                for rid, rseq in refs:
                    if Levenshtein.distance(a["seq"], rseq, score_cutoff=error_edits) <= error_edits:
                        parent = (rid, Levenshtein.distance(a["seq"], rseq), "expected sequence")
                        break
            if parent:
                row["class"], row["detail"] = "error variant", f"{parent[1]} edit(s) from {parent[0]} ({parent[2]})"
            elif hit:
                row["class"], row["detail"] = "unlisted variant", f"{hit[1]:.1f}% identical to {hit[0]}"
            else:
                row["class"], row["detail"] = "unexplained", "no expected sequence within 3%"
        rows.append(row)
    return rows


def summarise(rows: list) -> dict:
    """{class: (ASVs, reads)} and the total reads."""
    total = sum(r["size"] for r in rows) or 1
    out = {c: [0, 0] for c in CLASSES}
    for r in rows:
        out[r["class"]][0] += 1
        out[r["class"]][1] += r["size"]
    return {"classes": {c: tuple(v) for c, v in out.items()}, "reads": total, "asvs": len(rows)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("asvs", type=Path)
    ap.add_argument("--reference", type=Path, required=True, help="all expected sequences (FASTA)")
    ap.add_argument("--present", type=Path, help="the expected sequences present in these samples (default: all of --reference)")
    ap.add_argument("--asv-table", type=Path, help="the authors' ASV table (Mock_staggered-A_ASVs.csv) for the negative controls")
    ap.add_argument("--controls-regex", default="NK|Wasser", help="columns of that table that are negative controls (default NK|Wasser)")
    ap.add_argument("--out", type=Path, help="write one line per ASV to this TSV")
    args = ap.parse_args()
    if not shutil.which("vsearch"):
        raise SystemExit("vsearch was not found on PATH.")
    rows = classify(args.asvs, args.reference, args.present, args.asv_table, args.controls_regex)
    s = summarise(rows)
    print(f"{s['asvs']} ASVs carrying {s['reads']:,} reads\n")
    for c in CLASSES:
        n, reads = s["classes"][c]
        print(f"  {c:18s} {n:4d} ASVs  {reads / s['reads']:7.2%} of the reads")
    outside = [r for r in rows if r["class"] not in ("expected", "listed elsewhere")]
    print("\nthe largest ASVs outside the expected sequences:")
    for r in sorted(outside, key=lambda r: -r["size"])[:12]:
        print(f"  {r['id']:8s} {r['size']:6d} reads {r['length']} nt  {r['class']}: {r['detail']}")
    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["asv", "reads", "length", "class", "detail"])
            for r in rows:
                w.writerow([r["id"], r["size"], r["length"], r["class"], r["detail"]])
        print(f"\nwritten to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
