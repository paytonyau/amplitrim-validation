#!/usr/bin/env python
"""Bakker (2018): ready reads and reference sequences for the ITS tests.

  python validation\\bakker_prep.py --base C:\\data --appendix C:\\data\\bakker\\men12760-sup-0004-appendixs2.txt --out-references C:\\data\\references_bakker

1. Reads. Each read pair of the raw folders is assigned to an amplicon by its primers: the forward primer at the start of one read and the reverse primer at the
   start of the other (after a spacer of up to --max-offset bases, with up to --max-mismatch mismatches; IUPAC codes honoured). A pair whose forward primer is on
   R2 is swapped, so that R1 always starts at the forward primer. The spacer and the primer are cut off both reads, like the other datasets' reads. Pairs with
   neither primer pair are dropped and counted.
     22_mock_ITS1_bakker_pcr_PRJNA377530\\raw         -> ...\\prepared\\its1\\ready        (label bakker_pcr_its1)
     23_mock_ITSmix_bakker_polymerase_PRJNA377530\\raw -> ...\\prepared\\its1\\ready and ...\\prepared\\its2\\ready   (bakker_pol_its1, bakker_pol_its2)
2. References. The ITS1 and ITS2 sequences of the 19 isolates (Appendix S2 of the paper) are cut with the same primers, without the primers (in-silico PCR),
   and written in the form ncbi_references.py writes (member|species|isolate|strain_record|copy|length): they are of the very strains, so every member is judged.
The primers are those of the paper's methods (AMPLICONS): BITS and B58S3 on the first run (an ITS1 fragment), ITS1F and ITS2 (ITS1) and ITS3-KYO2 and
ITS4-KYO3 (ITS2) on the second; the primer census of the reads confirmed them. Each dataset has its own reference file, references_<label>.fasta.
Writes PREP_REPORT.md to --base and REFERENCES_REPORT.md to --out-references.
"""
import argparse
import gzip
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bakker_primers as bp  # noqa: E402

# The primers of the paper's methods (sections 2.3 and 2.4), confirmed by the primer census of the reads: label -> (locus, forward, reverse, raw run, ready folder)
AMPLICONS = {"bakker_pcr_its1": ("ITS1", "BITS", "B58S3"), "bakker_pol_its1": ("ITS1", "ITS1F", "ITS2"), "bakker_pol_its2": ("ITS2", "ITS3_KYO2", "ITS4_KYO3")}
# Members the marker cannot tell apart under the 2-edit rule count as one: on both ITS1 fragments F. oxysporum and F. verticillioides are 1 edit apart
# (on ITS2, 6 edits). Recall is then over 18 members on the ITS1 datasets and 19 on ITS2.
ONE_MEMBER = {"bakker_pcr_its1": {"Fusarium.oxysporum": "Fusarium.oxysporum+verticillioides", "Fusarium.verticillioides": "Fusarium.oxysporum+verticillioides"},
              "bakker_pol_its1": {"Fusarium.oxysporum": "Fusarium.oxysporum+verticillioides", "Fusarium.verticillioides": "Fusarium.oxysporum+verticillioides"}}
RAW = {"pcr": "22_mock_ITS1_bakker_pcr_PRJNA377530", "polymerase": "23_mock_ITSmix_bakker_polymerase_PRJNA377530"}
TARGETS = [("bakker_pcr_its1", "pcr", "its1"), ("bakker_pol_its1", "polymerase", "its1"), ("bakker_pol_its2", "polymerase", "its2")]
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def find(seq: str, primer: str, max_mm: int = 2, max_offset: int = 10):
    """(offset, mismatches) of `primer` at the start of seq (after up to max_offset bases), the fewest mismatches first; None if not there."""
    best = None
    for off in range(0, max_offset + 1):
        part = seq[off:off + len(primer)]
        if len(part) < len(primer):
            break
        mm = bp.mismatches(part, primer)
        if mm <= max_mm and (best is None or mm < best[1]):
            best = (off, mm)
    return best


def fastq(path: Path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s, p, q = fh.readline().rstrip("\n"), fh.readline().rstrip("\n"), fh.readline().rstrip("\n")
            yield h.rstrip("\n"), s, p, q


def trim3(seq: str, primer_rc: str, rate: float = 0.1, min_overlap: int = 3) -> int:
    """Where a read runs past the end of a short amplicon into the opposite primer: the length to keep. The opposite primer appears reverse-complemented,
    in full (up to `rate` mismatches) or, at the very end of the read, its first bases only (at least `min_overlap`), as cutadapt's -a/-A do in DADA2's
    ITS workflow. Returns len(seq) when the read does not run into it. (Bakker's reads did; run 28 left those bases on, see register section 38.)"""
    s, p = seq.upper(), primer_rc
    for i in range(0, len(s) - len(p) + 1):
        if bp.mismatches(s[i:i + len(p)], p) <= int(rate * len(p)):
            return i
    for L in range(min(len(p) - 1, len(s)), min_overlap - 1, -1):
        if bp.mismatches(s[-L:], p[:L]) <= int(rate * L):
            return len(s) - L
    return len(s)


def assign(r1: tuple, r2: tuple, fwd: str, rev: str, max_mm: int, max_offset: int):
    """The pair trimmed and in the forward orientation, or None when the pair does not carry this primer pair. r1, r2: (header, seq, plus, qual)."""
    for a, b, swapped in ((r1, r2, False), (r2, r1, True)):
        f, r = find(a[1].upper(), fwd, max_mm, max_offset), find(b[1].upper(), rev, max_mm, max_offset)
        if f and r:
            ka, kb = f[0] + len(fwd), r[0] + len(rev)
            sa, qa, sb, qb = a[1][ka:], a[3][ka:], b[1][kb:], b[3][kb:]
            ea, eb = trim3(sa, rc(rev)), trim3(sb, rc(fwd))                       # read-through into the opposite primer at short amplicons
            return (a[0], sa[:ea], a[2], qa[:ea]), (b[0], sb[:eb], b[2], qb[:eb]), swapped, (ea < len(sa)) or (eb < len(sb))
    return None


def prepare(raw: Path, out: Path, fwd: str, rev: str, max_mm: int = 2, max_offset: int = 10) -> list:
    """Writes <sample>_R1/_R2.fastq.gz of the pairs that carry fwd/rev to out; returns (sample, pairs, kept, swapped, read through the far primer)."""
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for f1 in sorted(Path(raw).glob("*_R1*.fastq.gz")):
        f2 = f1.with_name(f1.name.replace("_R1", "_R2", 1))
        if not f2.is_file():
            continue
        sample = f1.name.split("_R1")[0]
        n = kept = swapped = through = 0
        p1, p2 = out / f"{sample}_R1.fastq.gz.part", out / f"{sample}_R2.fastq.gz.part"
        with gzip.open(p1, "wt", newline="\n") as o1, gzip.open(p2, "wt", newline="\n") as o2:
            for a, b in zip(fastq(f1), fastq(f2)):
                n += 1
                got = assign(a, b, fwd, rev, max_mm, max_offset)
                if got:
                    x, y, sw, th = got
                    kept += 1
                    swapped += sw
                    through += th
                    o1.write(f"{x[0]}\n{x[1]}\n{x[2]}\n{x[3]}\n")
                    o2.write(f"{y[0]}\n{y[1]}\n{y[2]}\n{y[3]}\n")
        p1.replace(out / f"{sample}_R1.fastq.gz")
        p2.replace(out / f"{sample}_R2.fastq.gz")
        rows.append((sample, n, kept, swapped, through))
    return rows


def appendix(path: Path) -> dict:
    """{'ITS1': [(name, seq)], 'ITS2': [...]} from Appendix S2 (sections '## ITS1:' and '## ITS2:')."""
    out, sec, cur = {"ITS1": [], "ITS2": []}, None, None
    for line in Path(path).read_text(encoding="utf-8", errors="replace").replace("\r", "").splitlines():
        m = re.match(r"^##\s*(ITS[12])", line)
        if m:
            sec = m.group(1)
        elif line.startswith(">") and sec:
            cur = [line[1:].strip(), ""]
            out[sec].append(cur)
        elif cur and line.strip() and not line.startswith("#"):
            cur[1] += line.strip().upper()
    return {k: [tuple(x) for x in v] for k, v in out.items()}


def amplicon(seq: str, fwd: str, rev: str, max_mm: int = 3) -> str:
    """The sequence between the forward primer and the reverse primer's site, without the primers (the first forward site, the nearest reverse site after it)."""
    rrc = rc(rev)
    fsites = [i for i in range(0, len(seq) - len(fwd) + 1) if bp.mismatches(seq[i:i + len(fwd)], fwd) <= max_mm]
    for f in fsites:
        start = f + len(fwd)
        for j in range(start, len(seq) - len(rrc) + 1):
            if bp.mismatches(seq[j:j + len(rrc)], rrc) <= max_mm:
                return seq[start:j]
    return ""


def expand(seq: str, limit: int = 64) -> list:
    """Every concrete sequence an IUPAC-coded sequence stands for (R -> A, G; Y -> C, T; ...): a Sanger reference with ambiguity codes is judged by the variants
    it allows, so that a perfect read is not counted as edits away from it. More than `limit` variants is refused."""
    out = [""]
    for b in seq:
        out = [x + c for x in out for c in bp.IUPAC.get(b, b)]
        if len(out) > limit:
            raise ValueError(f"more than {limit} variants")
    return out


def member_of(name: str) -> tuple:
    """('Alternaria.alternata', 'NRRL.A-15700') from 'Alternaria.alternata.NRRL.A-15700'."""
    parts = name.split(".")
    return ".".join(parts[:2]), ".".join(parts[2:]) or "-"


def references(appx: dict, amplicons: dict, out: Path) -> dict:
    """One reference file per dataset (references_<label>.fasta), cut from the dataset's locus with the dataset's primers."""
    out.mkdir(parents=True, exist_ok=True)
    report = {}
    for label, (locus, f, r) in amplicons.items():
        recs, missing = [], []
        for name, seq in appx.get(locus, []):
            prod = amplicon(seq, bp.PRIMERS[f], bp.PRIMERS[r])
            member, isolate = member_of(name)
            member = ONE_MEMBER.get(label, {}).get(member, member)
            if prod:
                recs.append((member, isolate, prod))
            else:
                missing.append(name)
        copies, lines, coded = {}, [], []
        for member, isolate, prod in recs:
            variants = expand(prod)
            if len(variants) > 1:
                coded.append(f"{member} {isolate} ({len(variants)} variants)")
            for v in variants:
                copies[member] = copies.get(member, 0) + 1
                lines.append(f">{member}|{member.replace('.', '_')}|{isolate}|strain_record|copy{copies[member]}x1|{len(v)}nt\n{v}\n")
        (out / f"references_{label}.fasta").write_text("".join(lines), encoding="utf-8")
        report[label] = {"members": sorted(copies), "sequences": len(recs), "written": len(lines), "coded": coded, "missing": missing,
                         "lengths": sorted(len(p) for _, _, p in recs)}
    L = ["# Reference sequences of Bakker (2018), cut from Appendix S2", ""]
    for label, r in report.items():
        locus, f, rv = amplicons[label]
        L += [f"## {label}: {locus} ({f} / {rv})", "", f"{r['sequences']} sequences of {len(r['members'])} members; amplicons {r['lengths'][0] if r['lengths'] else '-'} to "
              f"{r['lengths'][-1] if r['lengths'] else '-'} nt without the primers.", f"Without a product: {', '.join(r['missing']) or 'none'}.",
              f"With ambiguity codes, written as every variant they allow ({r['written']} sequences in the file): {', '.join(r['coded']) or 'none'}.",
              f"Counted as one member (1 edit apart): {', '.join(sorted(set(ONE_MEMBER.get(label, {}).values()))) or 'none'}.", ""]
    (out / "REFERENCES_REPORT.md").write_text("\n".join(L) + "Next: send this file.\n", encoding="utf-8")
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", type=Path, required=True)
    ap.add_argument("--appendix", type=Path, required=True, help="men12760-sup-0004-appendixs2.txt")
    ap.add_argument("--out-references", type=Path, required=True)
    ap.add_argument("--max-mismatch", type=int, default=2)
    ap.add_argument("--max-offset", type=int, default=10)
    ap.add_argument("--references-only", action="store_true")
    a = ap.parse_args(argv)
    amps = dict(AMPLICONS)
    rep = references(appendix(a.appendix), amps, a.out_references)
    for label, r in rep.items():
        print(f"references {label}: {r['sequences']} amplicons of {len(r['members'])} members ({r['written']} sequences written); without a product: {len(r['missing'])}")
    if a.references_only:
        return 0 if all(len(r["members"]) == (18 if lab in ONE_MEMBER else 19) for lab, r in rep.items()) else 2
    L = ["# Bakker (2018): the reads assigned by primers", "", "The far primer is cut where a read runs into it (short amplicons).", "",
         "| Dataset | Primers | Sample | Read pairs | Kept | Of them swapped | Of them cut at the far primer |", "|---|---|---|---|---|---|---|"]
    for label, src, sub in TARGETS:
        raw = a.base / RAW[src] / "raw"
        if not raw.is_dir():
            print(f"not found: {raw}")
            return 2
        _, f, r = amps[label]
        rows = prepare(raw, a.base / RAW[src] / "prepared" / sub / "ready", bp.PRIMERS[f], bp.PRIMERS[r], a.max_mismatch, a.max_offset)
        L += [f"| {label} | {f} / {r} | {s} | {n:,} | {k:,} ({k / n:.0%}) | {w:,} | {t:,} ({t / k:.0%}) |" if k else f"| {label} | {f} / {r} | {s} | {n:,} | 0 | 0 | 0 |"
              for s, n, k, w, t in rows]
        print(f"{label}: {sum(x[2] for x in rows):,} of {sum(x[1] for x in rows):,} pairs kept")
    (a.base / "PREP_REPORT.md").write_text("\n".join(L) + "\n\nNext: send this file and REFERENCES_REPORT.md.\n", encoding="utf-8")
    print(f"written: {a.base / 'PREP_REPORT.md'}, {a.out_references / 'REFERENCES_REPORT.md'}")
    return 0 if all(len(r["members"]) == (18 if lab in ONE_MEMBER else 19) for lab, r in rep.items()) else 2


if __name__ == "__main__":
    sys.exit(main())
