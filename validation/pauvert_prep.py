#!/usr/bin/env python
"""The fourth mock community (Pauvert et al. 2019, Fungal Ecology 41:23-33): its reads made ready for the v0.4 night, and its truth from the authors' Table S1.

  python validation\\pauvert_prep.py --raw C:\\data\\24_mock_ITS1_pauvert_dataverse\\raw --table C:\\data\\24_mock_ITS1_pauvert_dataverse\\mmc2.csv
      --out C:\\data\\24_mock_ITS1_pauvert_dataverse

Reads (`Forward_R1_ITSEnz.DikA.fastq.gz` with `Reverse_R2_ITSEnz.DikA.fastq.gz`, and DikB, DikC): a pair is kept when ITS1F starts one read and ITS2 the other
(bakker_prep.assign: spacer up to 10 bases, up to 2 mismatches, swapped mates turned), both primers are cut at the start and, where a read runs past the amplicon,
at its end (bakker_prep.trim3); written as prepared/its1/ready/DikA_R1.fastq.gz and so on.

Truth (TEST_REGISTER.md, the section of the v0.4 night): each strain's Sanger sequence starts right after ITS1F. A complete one, holding the
reverse-complemented ITS2 site (up to 3 mismatches, at least 50 bases in), is cut before it. A partial one (the Sanger read stopped before the site, often
within ITS1) is kept whole and judged as a prefix: a sequence matches it when the sequence's first bases, as many as the reference has, are within 2 edits
of it (partial references shorter than MIN_PARTIAL are left out). Ambiguity codes are written as every variant they allow, up to 64, and a sequence with more
is left out. Strains whose sequences are identical or within 2 edits (a partial one against the other's first bases) count as one member. Expected: the members holding a strain
the authors found in the raw reads (Present_In_Raw_MiSeq_Data = Yes). Writes references_pauvert_its1.fasta (every strain kept, as a reference for judging),
expected_pauvert_its1.tsv and PAUVERT_REFERENCES.md.
"""
import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bakker_prep as B  # noqa: E402

ITS1F, ITS2 = "CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"
LABEL = "pauvert_its1"
MIN_START, MAX_MM, LIMIT, SAME, MIN_PARTIAL = 50, 3, 64, 2, 150


def within(a: str, b: str, k: int = SAME) -> bool:
    """Whether the edit distance of a and b is at most k (a banded Levenshtein: exact at that bound)."""
    if abs(len(a) - len(b)) > k:
        return False
    prev = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        lo, hi = max(1, i - k), min(len(b), i + k)
        cur = [k + 1] * (len(b) + 1)
        cur[0] = i if i <= k else k + 1
        for j in range(lo, hi + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
        if min(cur[max(0, lo - 1):hi + 1]) > k:
            return False
        prev = cur
    return prev[len(b)] <= k


def prefix_within(part: str, seq: str, k: int = SAME) -> bool:
    """Whether the partial reference `part` is within k edits of the first bases of `seq` (as many as `part` has, give or take k, for indels)."""
    return any(within(part, seq[:len(part) + d], k) for d in range(-k, k + 1) if len(part) + d <= len(seq))


def same(a: dict, b: dict) -> bool:
    """Two strains are one member: identical or within SAME edits, a partial reference against the first bases of the other."""
    for x in a["variants"]:
        for y in b["variants"]:
            if a["partial"] and b["partial"]:
                n = min(len(x), len(y))
                if within(x[:n], y[:n]):
                    return True
            elif a["partial"]:
                if prefix_within(x, y):
                    return True
            elif b["partial"]:
                if prefix_within(y, x):
                    return True
            elif within(x, y):
                return True
    return False


def cut_its1(seq: str) -> str:
    """The sequence before the reverse-complemented ITS2 site (the first site at least MIN_START bases in, up to MAX_MM mismatches); '' when there is none."""
    site = B.rc(ITS2)
    for j in range(MIN_START, len(seq) - len(site) + 1):
        if B.bp.mismatches(seq[j:j + len(site)], site) <= MAX_MM:
            return seq[:j]
    return ""


def table(path: Path) -> list:
    """Table S1 as dicts (semicolon-separated)."""
    with open(path, encoding="utf-8-sig", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=";"))


def strain_name(r: dict) -> str:
    return re.sub(r"[^A-Za-z0-9.-]+", "_", r["Mock_Strain"].strip()).strip("_")


def truth(rows: list, out: Path) -> dict:
    """Writes the references, the expected members and the report; returns the counts."""
    out.mkdir(parents=True, exist_ok=True)
    kept, no_site, too_coded = [], [], []
    for r in rows:
        name, seq = strain_name(r), re.sub(r"[^A-Za-z]", "", r["Sanger_Sequence"]).upper()
        cut = cut_its1(seq)
        partial = not cut
        if partial:
            if len(seq) < MIN_PARTIAL:
                no_site.append(f"{name} ({len(seq)} nt)")
                continue
            cut = seq
        try:
            variants = B.expand(cut, LIMIT)
        except ValueError:
            too_coded.append(f"{name} ({len(re.findall('[^ACGT]', cut))} codes)")
            continue
        kept.append({"name": name, "present": r["Present_In_Raw_MiSeq_Data"].strip().lower() == "yes", "variants": variants, "len": len(cut), "partial": partial})
    parent = list(range(len(kept)))                                   # strains within SAME edits (any variant) are one member
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(len(kept)):
        for j in range(i + 1, len(kept)):
            if find(i) != find(j) and same(kept[i], kept[j]):
                parent[find(j)] = find(i)
    groups = {}
    for i, k in enumerate(kept):
        groups.setdefault(find(i), []).append(k)
    members, lines, expected = {}, [], []
    for g in groups.values():
        names = sorted(k["name"] for k in g)
        member = names[0] if len(names) == 1 else f"{names[0]}+{len(names) - 1}"
        members[member] = g
        if any(k["present"] for k in g):
            expected.append(member)
        n = 0
        for k in g:
            for v in k["variants"]:
                n += 1
                lines.append(f">{member}|{member.replace('.', '_')}|{k['name']}|{'partial' if k['partial'] else 'complete'}|copy{n}x1|{len(v)}nt\n{v}\n")
    (out / f"references_{LABEL}.fasta").write_text("".join(lines), encoding="utf-8")
    (out / f"expected_{LABEL}.tsv").write_text("member\tstrains\n" + "".join(f"{m}\t{', '.join(k['name'] for k in members[m])}\n" for m in sorted(expected)), encoding="utf-8")
    merged = {m: [k["name"] for k in g] for m, g in members.items() if len(g) > 1}
    lens = sorted(k["len"] for k in kept if not k["partial"])
    plens = sorted(k["len"] for k in kept if k["partial"])
    absent_only = sorted(m for m in members if m not in expected)
    res = {"strains": len(rows), "kept": len(kept), "partial": len(plens), "no_site": no_site, "too_coded": too_coded, "members": len(members), "expected": len(expected),
           "merged": merged, "absent_members": absent_only, "lengths": (lens[0], lens[len(lens) // 2], lens[-1]) if lens else None, "written": len(lines)}
    L = ["# The truth of Pauvert et al. (2019), from Table S1", "",
         f"{res['strains']} strains. Complete (cut before the ITS2 site): {len(lens)}, {lens[0]} to {lens[-1]} nt, median {lens[len(lens) // 2]}, without the primers.",
         f"Partial (the Sanger read stopped before the ITS2 site; judged as a prefix): {len(plens)}, {plens[0] if plens else '-'} to {plens[-1] if plens else '-'} nt.",
         f"Left out, partial and shorter than {MIN_PARTIAL} nt: {len(no_site)} ({', '.join(no_site) or 'none'}).",
         f"Left out, more than {LIMIT} variants: {len(too_coded)} ({', '.join(too_coded) or 'none'}).",
         f"Members (strains identical or within {SAME} edits counted as one): {res['members']}; {len(merged)} of them hold more than one strain:",
         ""] + [f"- {m}: {', '.join(v)}" for m, v in sorted(merged.items())] + [
         "", f"**Expected members** (holding a strain the authors found in the raw reads): **{res['expected']}** (the gate: at least 120).",
         f"Members whose strains the authors did not find in the raw reads (references only, not expected): {len(absent_only)}.",
         f"References written: {res['written']} sequences (ambiguity codes expanded).", "", "Next: send this file."]
    (out / "PAUVERT_REFERENCES.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    return res


def reads(raw: Path, out: Path) -> list:
    """Pairs Forward_R1_*.DikX with Reverse_R2_*.DikX, keeps the pairs that carry the primers and cuts them (bakker_prep.assign); one row per sample."""
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for f1 in sorted(Path(raw).glob("Forward_R1_*.fastq.gz")):
        f2 = f1.with_name(f1.name.replace("Forward_R1_", "Reverse_R2_", 1))
        if not f2.is_file():
            continue
        m = re.search(r"(Dik[A-Z])", f1.name)
        sample = m.group(1) if m else f1.name.split(".")[0]
        n = kept = swapped = through = 0
        p1, p2 = out / f"{sample}_R1.fastq.gz.part", out / f"{sample}_R2.fastq.gz.part"
        with gzip.open(p1, "wt", newline="\n") as o1, gzip.open(p2, "wt", newline="\n") as o2:
            for a, b in zip(B.fastq(f1), B.fastq(f2)):
                n += 1
                got = B.assign(a, b, ITS1F, ITS2, 2, 10)
                if got:
                    x, y, sw, th = got
                    kept, swapped, through = kept + 1, swapped + sw, through + th
                    o1.write(f"{x[0]}\n{x[1]}\n{x[2]}\n{x[3]}\n")
                    o2.write(f"{y[0]}\n{y[1]}\n{y[2]}\n{y[3]}\n")
        p1.replace(out / f"{sample}_R1.fastq.gz")
        p2.replace(out / f"{sample}_R2.fastq.gz")
        rows.append((sample, n, kept, swapped, through))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", type=Path, help="the folder of the six FASTQ files")
    ap.add_argument("--table", type=Path, required=True, help="Table S1 (mmc2.csv)")
    ap.add_argument("--out", type=Path, required=True, help="the dataset folder (24_mock_ITS1_pauvert_dataverse)")
    ap.add_argument("--references-only", action="store_true")
    a = ap.parse_args(argv)
    res = truth(table(a.table), a.out / "references")
    print((a.out / "references" / "PAUVERT_REFERENCES.md").read_text(encoding="utf-8"))
    if res["expected"] < 120:
        print(f"The truth gate fails: {res['expected']} expected members (at least 120).")
        return 2
    if a.references_only:
        return 0
    if not a.raw or not a.raw.is_dir():
        ap.error("--raw must be the folder of the FASTQ files")
    rows = reads(a.raw, a.out / "prepared" / "its1" / "ready")
    if not rows:
        print(f"No Forward_R1_*/Reverse_R2_* pairs in {a.raw}")
        return 2
    L = ["# Pauvert's reads, prepared", "", "| Sample | Pairs | Kept (ITS1F/ITS2) | Mates swapped | Read past the amplicon |", "|---|---|---|---|---|"]
    L += [f"| {s} | {n:,} | {k:,} ({k / max(1, n):.1%}) | {sw:,} | {th:,} |" for s, n, k, sw, th in rows]
    (a.out / "PREP_REPORT.md").write_text("\n".join(L) + "\n\nNext: send this file.\n", encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
