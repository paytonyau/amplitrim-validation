#!/usr/bin/env python
"""Reference ITS sequences for the mock communities of Arfken et al. 2023, built from NCBI genomes by in-silico PCR with the paper's primers.

Arfken AM, Frey JF, Carrillo NI, Dike NI, Onyeachonamm O, Rivera DN, Davies CP, Summers KL (2023). Porcine fungal mock community analyses: Implications for
mycobiome investigations. Frontiers in Cellular and Infection Microbiology 13:928353. doi:10.3389/fcimb.2023.928353 (BioProject PRJNA693350).

  python validation\\ncbi_references.py --out C:\\data\\references_its --list
        shows, for every member of the two mocks, the genome assemblies NCBI has and the one that would be chosen; nothing is downloaded
  python validation\\ncbi_references.py --out C:\\data\\references_its --compare-isolates C:\\data\\runs\\25_its_followup_arfken
        downloads the chosen genomes, cuts out the ITS1 and ITS2 amplicons, and compares them with the sequences of the isolates seen in run 25
  --use Af=GCF_000002655.1   uses that assembly for that member (after the review)          --only Af,Cn   only those members

The mocks: the ATCC mycobiome standard MSA-1010 (ten strains with ATCC numbers, Table 3 of the paper) and the five isolates from piglet faeces (C, K, M, P, T).
Primers (Table 1 of the paper): ITS1 amplicon ITS1F CTTGGTCATTTAGAGGAAGTAA and ITS2 GCTGCGTTCTTCATCGATGC; ITS2 amplicon ITS3 GCATCGATGAAGAACGCAGC and ITS4
TCCTCCGCTTATTGATATGC. The products are written without the primers, like the reads, which had their primers removed. A product is accepted with at most --max-mismatch
mismatches per primer, anywhere in it (default 3). When a genome gives no product for a marker, NCBI's nucleotide records of ITS of the species are tried, and the report says so.
Writes references_ITS1.fasta, references_ITS2.fasta and REFERENCES_REPORT.md to --out (genomes are kept in --out\\genomes, so a second run downloads nothing).
How good a reference is, is said for each member: 'strain' (an assembly of the very strain), or 'species' (another strain of the species; real strain differences can then
look like errors). The piglet isolates were never sequenced as genomes in public, except perhaps K. slooffiae (Davies et al. 2021, the same laboratory), so they are 'species'.
"""
import argparse
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PRIMERS = {"ITS1": ("CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"), "ITS2": ("GCATCGATGAAGAACGCAGC", "TCCTCCGCTTATTGATATGC")}
# member id, mock, names to search at NCBI (in order), strings that identify the strain in an assembly's strain, isolate or name
TARGETS = [
    ("Af", "ATCC", ["Aspergillus fumigatus"], ["MYA-4609", "Af293", "CBS 101355"]),
    ("Cn", "ATCC", ["Cryptococcus neoformans var. grubii", "Cryptococcus neoformans"], ["208821", "H99"]),
    ("Ti", "ATCC", ["Trichophyton interdigitale"], ["9533"]),
    ("Pc", "ATCC", ["Penicillium chrysogenum", "Penicillium rubens"], ["10106"]),
    ("Fk", "ATCC", ["Fusarium keratoplasticum"], ["36031"]),
    ("Ca", "ATCC", ["Candida albicans"], ["10231"]),
    ("Cg", "ATCC", ["Nakaseomyces glabratus", "Candida glabrata"], ["2001", "CBS138", "CBS 138"]),
    ("Mg", "ATCC", ["Malassezia globosa"], ["MYA-4612", "CBS 7966"]),
    ("Sc", "ATCC", ["Saccharomyces cerevisiae"], ["201390", "BY4743", "S288C"]),
    ("Cd", "ATCC", ["Cutaneotrichosporon dermatis", "Trichosporon dermatis"], ["204094"]),
    ("C", "piglet", ["Candida albicans"], []),
    ("K", "piglet", ["Kazachstania slooffiae"], []),
    ("M", "piglet", ["Lichtheimia corymbifera"], []),
    ("P", "piglet", ["Pichia fermentans"], []),
    ("T", "piglet", ["Trichosporon asahii"], []),
]
LEVEL = {"complete genome": 4, "chromosome": 3, "scaffold": 2, "contig": 1}
CATEGORY = {"reference genome": 2, "representative genome": 1}


# ---------------------------------------------------------------- the network (one function, replaceable in tests)
def fetch(url: str, tries: int = 3, sleep=time.sleep) -> bytes:
    last = None
    for attempt in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "amplitrim-validation (ncbi_references.py)"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:                                           # NCBI is sometimes busy; try again
            last = e
            if attempt < tries:
                sleep(3 * attempt)
    raise RuntimeError(f"NCBI did not answer {url}: {last}")


# ---------------------------------------------------------------- sequences
def rc(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))


def mismatches(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y and x != "N" and y != "N")


def sites(seq: str, primer: str, max_mm: int) -> list:
    """Every start position where `primer` matches `seq` with at most max_mm mismatches, anywhere in the primer. The primer is cut into max_mm + 1 pieces; a site with at
    most max_mm mismatches matches at least one piece exactly, so searching the pieces finds every such site (a mismatch at the 3' end included: ITS1F has one against
    Saccharomycetes, which Arfken et al. 2023 still amplified), and the search stays fast on a whole genome."""
    n, k = len(primer), max_mm + 1
    bounds = [round(i * n / k) for i in range(k + 1)]
    found = set()
    for lo, hi in zip(bounds, bounds[1:]):
        piece = primer[lo:hi]
        i = seq.find(piece)
        while i != -1:
            start = i - lo
            if start >= 0 and start + n <= len(seq) and start not in found and mismatches(seq[start:start + n], primer) <= max_mm:
                found.add(start)
            i = seq.find(piece, i + 1)
    return sorted(found)


def pcr(seq: str, fwd: str, rev: str, max_mm: int = 3, min_len: int = 50, max_len: int = 1500) -> list:
    """Products between the two primers, on both strands, without the primers: each forward site with the nearest reverse site downstream."""
    out = []
    for strand, s in (("+", seq.upper()), ("-", rc(seq.upper()))):
        rev_rc = rc(rev)
        fsites = sites(s, fwd, max_mm)
        rsites = sites(s, rev_rc, max_mm)                                  # the reverse primer as it reads on this strand
        for f in fsites:
            end = f + len(fwd)
            nxt = [r for r in rsites if min_len <= r - end <= max_len]
            if nxt:
                out.append((strand, f, nxt[0], s[end:nxt[0]]))
    return out


def fasta_records(handle):
    name, parts = None, []
    for line in handle:
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(parts)
            name, parts = line[1:].split()[0] if len(line) > 1 else "", []
        elif line:
            parts.append(line)
    if name is not None:
        yield name, "".join(parts)


def zip_records(path: Path):
    with zipfile.ZipFile(path) as z:
        for member in z.namelist():
            if member.endswith((".fna", ".fa", ".fasta")):
                with z.open(member) as fh:
                    yield from fasta_records(io.TextIOWrapper(fh, encoding="ascii", errors="replace"))


def levenshtein(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


# ---------------------------------------------------------------- NCBI
def norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s or "").lower())


def assemblies(name: str, get=fetch) -> list:
    """The genome assemblies NCBI Datasets lists for a taxon name (the reference ones first, then up to 1000 others)."""
    found, seen = [], set()
    for extra in ("&filters.reference_only=true", ""):
        url = f"{DATASETS}/genome/taxon/{urllib.parse.quote(name)}/dataset_report?page_size=1000{extra}"
        try:
            reports = json.loads(get(url).decode("utf-8")).get("reports", [])
        except (RuntimeError, ValueError):
            reports = []
        for r in reports:
            acc = r.get("accession") or r.get("current_accession")
            if not acc or acc in seen:
                continue
            seen.add(acc)
            org, info = r.get("organism", {}), r.get("assembly_info", {})
            names = org.get("infraspecific_names", {}) or {}
            found.append({"accession": acc, "organism": org.get("organism_name", ""), "strain": names.get("strain", "") or names.get("isolate", ""),
                          "assembly_name": info.get("assembly_name", ""), "level": info.get("assembly_level", ""), "category": info.get("refseq_category", ""),
                          "n50": (r.get("assembly_stats", {}) or {}).get("contig_n50") or 0})
    return found


def choose(reports: list, hints: list):
    """(assembly, how): an assembly of the very strain if one exists ('strain'), else the reference or representative genome, else the most complete ('species')."""
    if not reports:
        return None, "none"
    hn = [norm(h) for h in hints if h]

    def is_strain(r):
        text = norm(r["strain"]) + "|" + norm(r["assembly_name"])
        return any(h and h in text for h in hn)
    rank = lambda r: (CATEGORY.get(str(r["category"]).lower(), 0), LEVEL.get(str(r["level"]).lower(), 0), r["accession"].startswith("GCF_"), float(r["n50"] or 0))
    strain = [r for r in reports if is_strain(r)]
    if strain:
        return max(strain, key=rank), "strain"
    return max(reports, key=rank), "species"


def download(acc: str, folder: Path, get=fetch) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    dest = folder / f"{acc}.zip"
    if not (dest.is_file() and dest.stat().st_size > 0 and zipfile.is_zipfile(dest)):
        data = get(f"{DATASETS}/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA")
        part = dest.with_suffix(".part")
        part.write_bytes(data)
        if not zipfile.is_zipfile(part):
            part.unlink()
            raise RuntimeError(f"the download of {acc} is not a zip file")
        part.replace(dest)
    return dest


def its_records(name: str, get=fetch, retmax: int = 20, strain: str = "") -> list:
    """(accession, sequence) of NCBI nucleotide records of ITS of the species (of the strain, when `strain` is given), for a member whose genome gives no product."""
    term = f'"{name}"[Organism] AND (internal transcribed spacer[Title] OR ITS1[Title] OR ITS2[Title])'
    if strain:
        term += f' AND "{strain}"[All Fields]'
    ids = json.loads(get(f"{EUTILS}/esearch.fcgi?db=nucleotide&retmode=json&retmax={retmax}&term={urllib.parse.quote(term)}").decode("utf-8")).get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    text = get(f"{EUTILS}/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id={','.join(ids)}").decode("utf-8", errors="replace")
    return list(fasta_records(io.StringIO(text)))


def products(records, max_mm, min_len, max_len) -> dict:
    """marker -> {product sequence: copies}"""
    out = {m: {} for m in PRIMERS}
    for _, seq in records:
        for marker, (f, r) in PRIMERS.items():
            for _, _, _, p in pcr(seq, f, r, max_mm, min_len, max_len):
                out[marker][p] = out[marker].get(p, 0) + 1
    return out


# ---------------------------------------------------------------- the isolates of run 25
def isolate_sequences(folder: Path) -> dict:
    """marker -> isolate letter -> the most abundant merged sequence of that isolate at the default (run 25: arfken_its1_seqs.tsv, arfken_its2_seqs.tsv)."""
    out = {}
    for marker, f in (("ITS1", "arfken_its1_seqs.tsv"), ("ITS2", "arfken_its2_seqs.tsv")):
        p = Path(folder) / f
        if not p.is_file():
            continue
        tot = {}
        for line in p.read_text(encoding="utf-8").splitlines()[1:]:
            s, k, q, a = line.split("\t")
            if int(float(k)) == 0:
                iso = s.split("-")[0]
                tot.setdefault(iso, {})
                tot[iso][q] = tot[iso].get(q, 0) + float(a)
        out[marker] = {iso: max(d, key=d.get) for iso, d in tot.items()}
    return out


# ---------------------------------------------------------------- the run
def main(argv=None, get=fetch) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--list", action="store_true", help="show the assemblies NCBI has and the choice; download nothing")
    ap.add_argument("--use", action="append", default=[], metavar="ID=ACCESSION", help="use this assembly for this member")
    ap.add_argument("--only", default="", help="members to do, e.g. Af,Cn (default: all)")
    ap.add_argument("--compare-isolates", type=Path, help="run 25's folder with arfken_its1_seqs.tsv and arfken_its2_seqs.tsv")
    ap.add_argument("--max-mismatch", type=int, default=3, help="mismatches allowed per primer, anywhere in it (default 3)")
    ap.add_argument("--min-length", type=int, default=50)
    ap.add_argument("--max-length", type=int, default=1500)
    a = ap.parse_args(argv)
    forced = {}
    for item in a.use:
        k, sep, v = item.partition("=")
        if not sep or not v:
            raise SystemExit(f"--use wants ID=ACCESSION, not {item!r}")
        forced[k] = v
    only = {x for x in a.only.split(",") if x}
    targets = [t for t in TARGETS if not only or t[0] in only]
    unknown = (only | set(forced)) - {t[0] for t in TARGETS}
    if unknown:
        raise SystemExit(f"unknown member(s): {sorted(unknown)}; the members are {[t[0] for t in TARGETS]}")
    a.out.mkdir(parents=True, exist_ok=True)
    isolates = isolate_sequences(a.compare_isolates) if a.compare_isolates else {}
    rows, fasta = [], {m: [] for m in PRIMERS}
    for tid, mock, names, hints in targets:
        reports, used_name = [], names[0]
        for name in names:
            reports = assemblies(name, get)
            if reports:
                used_name = name
                break
        if tid in forced:
            pick = next((r for r in reports if r["accession"] == forced[tid]), {"accession": forced[tid], "organism": used_name, "strain": "(given with --use)", "assembly_name": "",
                                                                                 "level": "", "category": "", "n50": 0})
            how = "strain" if hints and any(norm(h) in norm(pick["strain"]) + norm(pick["assembly_name"]) for h in hints) else "given"
        else:
            pick, how = choose(reports, hints)
        if mock == "piglet" and how == "strain":
            how = "species"                                              # the piglet isolates were never sequenced as genomes in public
        print(f"{tid:3s} {mock:6s} {used_name:38s} {len(reports):4d} assemblies; chosen: "
              + (f"{pick['accession']} ({pick['strain'] or 'no strain'}, {pick['level'] or '?'}, {pick['category'] or '-'}) [{how}]" if pick else "none"))
        if a.list:
            for r in sorted(reports, key=lambda r: (-CATEGORY.get(str(r["category"]).lower(), 0), -LEVEL.get(str(r["level"]).lower(), 0)))[:5]:
                print(f"        {r['accession']:18s} {str(r['strain'])[:30]:30s} {r['level']:16s} {r['category']}")
            continue
        found, source = {m: {} for m in PRIMERS}, ""
        kind = {m: "" for m in PRIMERS}
        if pick:
            try:
                found = products(zip_records(download(pick["accession"], a.out / "genomes", get)), a.max_mismatch, a.min_length, a.max_length)
                source = f"genome {pick['accession']}"
                for m in PRIMERS:
                    if found[m]:
                        kind[m] = f"{how} (genome)"
            except RuntimeError as e:
                print(f"    {e}")
        for label, strain in ([(h, h) for h in hints] if mock == "ATCC" else []) + [("species", "")]:
            missing = [m for m in PRIMERS if not found[m]]
            if not missing:
                break
            extra = products(its_records(used_name, get, strain=strain), a.max_mismatch, a.min_length, a.max_length)
            for m in missing:
                if extra[m]:
                    found[m] = extra[m]
                    kind[m] = "strain (record)" if strain else "species (record)"
                    source += f"; {m} from NCBI ITS records" + (f" of {strain}" if strain else " of the species")
        row = {"id": tid, "mock": mock, "species": used_name, "assembly": pick["accession"] if pick else "", "strain": pick["strain"] if pick else "", "how": how, "source": source.strip("; ")}
        for m in PRIMERS:
            row[f"{m}_kind"] = kind[m] or "none"
        for m in PRIMERS:
            seqs = sorted(found[m].items(), key=lambda kv: -kv[1])
            row[m] = seqs
            for i, (q, n) in enumerate(seqs, 1):
                fasta[m].append((f"{tid}|{used_name.replace(' ', '_')}|{row['assembly']}|{re.sub(r'[^a-z]+', '_', kind[m]).strip('_')}|copy{i}x{n}|{len(q)}nt", q))
            if mock == "piglet" and m in isolates and tid in isolates[m] and seqs:
                row[f"{m}_edits"] = min(levenshtein(q, isolates[m][tid]) for q, _ in seqs)
                row[f"{m}_isolate_len"] = len(isolates[m][tid])
        rows.append(row)
    if a.list:
        print("\nNothing was downloaded. Without --list the chosen genomes are downloaded; --use ID=ACCESSION changes a choice.")
        return 0
    for m, recs in fasta.items():
        (a.out / f"references_{m}.fasta").write_text("".join(f">{h}\n{q}\n" for h, q in recs), encoding="utf-8")
    L = ["# Reference ITS sequences from NCBI", "", "In-silico PCR with the primers of Arfken et al. 2023 (Table 1); products without the primers. 'strain': an assembly of the very strain; "
         "'species': another strain of the species; 'given': chosen with --use.", "",
         "| Member | Mock | Species | Assembly | Strain in the assembly | ITS1 reference is | ITS2 reference is | ITS1 products (nt x copies) | ITS2 products (nt x copies) | Edits to the isolate (ITS1 / ITS2) | Source |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        prod = lambda m: ", ".join(f"{len(q)}x{n}" for q, n in r[m]) or "none"
        ed = " / ".join(f"{r[f'{m}_edits']} (isolate {r[f'{m}_isolate_len']} nt)" if f"{m}_edits" in r else "-" for m in PRIMERS) if r["mock"] == "piglet" else ""
        L.append(f"| {r['id']} | {r['mock']} | {r['species']} | {r['assembly'] or '-'} | {r['strain'] or '-'} | {r['ITS1_kind']} | {r['ITS2_kind']} | {prod('ITS1')} | {prod('ITS2')} | {ed} | {r['source'] or '-'} |")
    none = [r["id"] for r in rows if not r["ITS1"] and not r["ITS2"]]
    atcc = [r for r in rows if r["mock"] == "ATCC"]
    exact = {m: [r["id"] for r in atcc if r[f"{m}_kind"].startswith("strain")] for m in PRIMERS}
    L += ["", f"Members without any product: {', '.join(none) if none else 'none'}.", ""]
    if atcc:
        L += ["ATCC members with a reference of the very strain (from its genome or its own ITS records): "
              + "; ".join(f"{m} {len(exact[m])} of {len(atcc)} ({', '.join(exact[m]) or 'none'})" for m in PRIMERS) + ".", ""]
    L += [
          "Read before the plan uses these: a member with no product, a product far from the length the paper or run 25 shows, or a piglet member many edits from its isolate "
          "needs another assembly (--use) or stays out of the truth. Next: send this file."]
    (a.out / "REFERENCES_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"\nwritten: {a.out / 'REFERENCES_REPORT.md'}, references_ITS1.fasta, references_ITS2.fasta")
    return 2 if none else 0


if __name__ == "__main__":
    sys.exit(main())
