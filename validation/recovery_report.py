#!/usr/bin/env python
"""How many of the known members each mock community gave back (TEST_REGISTER.md, the recovery report; described, no criteria).

  python validation\\recovery_report.py --its-run results\\31_v04_night --its-references validation\\references_its
      --bakker-references validation\\references_bakker --pauvert-references C:\\data\\24_mock_ITS1_pauvert_dataverse\\references
      --sixteen-run results\\26_its_night --base C:\\data --out C:\\data\\runs\\recovery

For each mock community: the members expected in each sample, and the members found, counted two ways: holding at least 0.1% of the sample's reads (the
rule of sections 33 to 44) and holding at least one read (as most papers count); then the members found in any sample and in all samples.
ITS: run 31's arms std (DADA2's defaults) and sw15 (the trim), judged as in sections 33 to 44 (within 2 edits; Pauvert with its partial references).
16S: run 26's 16S controls (DADA2 at AmpliTrim's cut, maxEE 2,2), judged as mock_check.py judges (an ASV identical to part of a reference, vsearch);
a strain is found when one of its variants is. Writes RECOVERY.md and recovery.tsv.
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import its_night_report as N  # noqa: E402
import v04_report as V  # noqa: E402
import mock_check as MC  # noqa: E402

SEEN = N.SEEN_SHARE                                   # 0.1% of a sample's reads
ITS_SETS = ["atcc_its1", "atcc_its2", "arfken_mix_its1", "arfken_mix_its2", "bakker_pcr_its1", "bakker_pol_its1", "bakker_pol_its2", V.PAUVERT]
SIXTEEN = {"mockA_high": "01_mock_V1-V3_mockA_PRJEB67797/reference/mockA_reference.fasta",
           "mockB": "02_mock_V1-V3_mockB_PRJEB67827/reference/reference.fasta",
           "d6300_V3-V4": "05_mock_V3-V4_zymoD6300_PRJNA1261676/reference/d6300_V3-V4_expected_sequences.fasta",
           "d6300_V4-V5": "06_mock_V4-V5_zymoD6300_miseq/reference/d6300_V4-V5_expected_sequences.fasta",
           "miseqsop_mock": "14_mock_V4_miseqsop/HMP_MOCK.v35.strains.fasta"}
PUBLISHED = {V.PAUVERT: "75.8% of the 189 strains per replicate (Pauvert et al. 2019, their recommended DADA2 approach, Si5)"}


def count(found_by_sample: dict, expected_by_sample: dict) -> dict:
    """found_by_sample: {sample: {member: share of the sample's reads}}; expected_by_sample: {sample: [members]}."""
    rows = []
    for s, exp in expected_by_sample.items():
        f = {m: x for m, x in found_by_sample.get(s, {}).items() if m in exp}
        rows.append((len(exp), sum(1 for x in f.values() if x >= SEEN), sum(1 for x in f.values() if x > 0)))
    exp_any = set().union(*[set(e) for e in expected_by_sample.values()]) if expected_by_sample else set()
    any01, any1 = set(), set()
    all01 = all1 = None
    for s, exp in expected_by_sample.items():
        f = found_by_sample.get(s, {})
        a01, a1 = {m for m in exp if f.get(m, 0) >= SEEN}, {m for m in exp if f.get(m, 0) > 0}
        any01, any1 = any01 | a01, any1 | a1
    common = set.intersection(*[set(e) for e in expected_by_sample.values()]) if expected_by_sample else set()
    all01 = {m for m in common if all(found_by_sample.get(s, {}).get(m, 0) >= SEEN for s in expected_by_sample)}
    all1 = {m for m in common if all(found_by_sample.get(s, {}).get(m, 0) > 0 for s in expected_by_sample)}
    n = len(rows) or 1
    return {"samples": len(rows), "expected": sum(r[0] for r in rows) / n, "found01": sum(r[1] for r in rows) / n, "found1": sum(r[2] for r in rows) / n,
            "rate01": sum(r[1] / r[0] for r in rows if r[0]) / n, "rate1": sum(r[2] / r[0] for r in rows if r[0]) / n,
            "min01": min((r[1] for r in rows), default=0), "max01": max((r[1] for r in rows), default=0),
            "any_expected": len(exp_any), "any01": len(any01), "any1": len(any1), "all_expected": len(common), "all01": len(all01), "all1": len(all1),
            "found_any1": any1}


def its_found(seqs: dict, truth, mix: bool) -> tuple:
    found, expected = {}, {}
    for s, d in sorted(seqs.items()):
        tot = sum(d.values()) or 1
        expected[s] = N.expected_members(s, truth, mix)
        f = {}
        for q, a in d.items():
            kind, m, _ = truth.classify(q)
            if kind == "true" and m:
                f[m] = f.get(m, 0) + a / tot
        found[s] = f
    return found, expected


def its_rows(run: Path, its_refs: Path, bakker_refs: Path, pauvert_refs: Path) -> list:
    out = []
    refs = {m: N.read_references(its_refs / f"references_{m}.fasta") for m in ("ITS1", "ITS2") if (its_refs / f"references_{m}.fasta").is_file()}
    truths = {}
    for marker, lab in N.PRIMARY.items():
        if marker in refs:
            truths[lab] = (N.fasta_truth(refs[marker]), False)
    for marker, (iso, mix) in N.SECOND.items():
        p = run / iso / "dada2" / "std" / "seqs.tsv"
        if p.is_file():
            truths[mix] = (N.isolate_truth(N.read_seqs(p)), True)
    for lab in ("bakker_pcr_its1", "bakker_pol_its1", "bakker_pol_its2"):
        f = bakker_refs / f"references_{lab}.fasta"
        if f.is_file():
            truths[lab] = (N.fasta_truth(N.read_references(f)), False)
    if (pauvert_refs / f"references_{V.PAUVERT}.fasta").is_file():
        truths[V.PAUVERT] = (V.pauvert_truth(pauvert_refs), False)
    strains = {}
    exp = pauvert_refs / f"expected_{V.PAUVERT}.tsv"
    if exp.is_file():
        with open(exp, encoding="utf-8") as fh:
            strains = {r["member"]: len(r["strains"].split(", ")) for r in csv.DictReader(fh, delimiter="\t")}
    for lab in ITS_SETS:
        if lab not in truths:
            continue
        truth, mix = truths[lab]
        for arm in ("std", "sw15"):
            p = run / lab / "dada2" / arm / "seqs.tsv"
            if not p.is_file():
                continue
            seqs = N.read_seqs(p)
            samples = [r["sample"] for r in N.read_samples(run / lab / "dada2" / arm / "samples.tsv")]
            for s in samples:
                seqs.setdefault(s, {})
            c = count(*its_found(seqs, truth, mix))
            c["strains_any1"] = sum(strains.get(m, 1) for m in c["found_any1"]) if lab == V.PAUVERT else None
            out.append({"marker": "ITS", "dataset": lab, "arm": arm, **c})
    return out


def best_hits(seqs: list, reference: Path, vsearch: str = "vsearch") -> dict:
    """{sequence: (reference name, identity)} as mock_check.summarise matches (vsearch, terminal gaps not counted)."""
    with tempfile.TemporaryDirectory() as tmp:
        q, out = Path(tmp) / "q.fasta", Path(tmp) / "hits.tsv"
        q.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)), encoding="utf-8")
        res = subprocess.run([vsearch, "--usearch_global", str(q), "--db", str(reference), "--id", "0.99", "--strand", "both", "--maxaccepts", "0",
                              "--maxrejects", "0", "--userout", str(out), "--userfields", "query+target+id", "--quiet"], capture_output=True, text=True)
        if res.returncode != 0:
            raise SystemExit(f"vsearch failed: {res.stderr.strip()[-300:]}")
        best = {}
        for line in out.read_text().splitlines():
            qq, t, ident = line.split("\t")
            i = int(qq[1:])
            if i not in best or float(ident) > best[i][1]:
                best[i] = (t.split()[0], float(ident))
    return {seqs[i]: v for i, v in best.items()}


def sixteen_rows(run: Path, base: Path, hits=best_hits) -> list:
    out = []
    for name, rel in SIXTEEN.items():
        ref, p = Path(base) / rel, Path(run) / f"16s_{name}" / "dada2" / "std" / "seqs.tsv"
        if not ref.is_file() or not p.is_file():
            out.append({"marker": "16S", "dataset": name, "arm": "std", "missing": str(ref if not ref.is_file() else p)})
            continue
        ref_names = [n.split()[0] for n, _ in MC.read_fasta(ref)]
        strains = sorted({MC.strain_of(r) for r in ref_names})
        seqs = N.read_seqs(p)
        for r in N.read_samples(Path(run) / f"16s_{name}" / "dada2" / "std" / "samples.tsv"):
            seqs.setdefault(r["sample"], {})
        best = hits(sorted({q for d in seqs.values() for q in d}), ref)
        found, expected = {}, {}
        for s, d in seqs.items():
            tot, f = sum(d.values()) or 1, {}
            for q, a in d.items():
                h = best.get(q)
                if h and h[1] >= 99.999:                                          # identical to part of a reference (mock_check's "exact")
                    st = MC.strain_of(h[0])
                    f[st] = f.get(st, 0) + a / tot
            found[s], expected[s] = f, strains
        exact_variants = {best[q][0] for d in seqs.values() for q in d if best.get(q) and best[q][1] >= 99.999}
        out.append({"marker": "16S", "dataset": name, "arm": "std", "variants": f"{len(exact_variants)} of {len(ref_names)}", **count(found, expected)})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--its-run", type=Path, required=True)
    ap.add_argument("--its-references", type=Path, required=True)
    ap.add_argument("--bakker-references", type=Path, required=True)
    ap.add_argument("--pauvert-references", type=Path, required=True)
    ap.add_argument("--sixteen-run", type=Path, required=True)
    ap.add_argument("--base", type=Path, required=True, help="the data folder holding the 16S mock communities and their references")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    rows = its_rows(a.its_run, a.its_references, a.bakker_references, a.pauvert_references) + sixteen_rows(a.sixteen_run, a.base)
    a.out.mkdir(parents=True, exist_ok=True)
    L = ["# Recovery of the known members of the mock communities", "",
         "Members found per sample: holding at least 0.1% of the sample's reads (the rule of sections 33 to 44), and holding at least one read. "
         "Then the members found in any sample and in all samples (of those expected in every sample). ITS: run 31 (std: DADA2's defaults; sw15: the trim). "
         "16S: run 26 (DADA2 at AmpliTrim's cut, maxEE 2,2), a strain found when one of its variants is (exact, as mock_check.py).", "",
         "| Marker | Dataset | Arm | Samples | Expected per sample | Found per sample, at least 0.1% | Found per sample, any read | In any sample (0.1% / any read) | In all samples (0.1% / any read) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if "missing" in r:
            L.append(f"| {r['marker']} | {r['dataset']} | {r['arm']} | not found: {r['missing']} | | | | | |")
            continue
        extra = f" ({r['strains_any1']} strains)" if r.get("strains_any1") else ""
        L.append(f"| {r['marker']} | {r['dataset']} | {r['arm']} | {r['samples']} | {r['expected']:.1f} | {r['found01']:.1f} ({r['rate01']:.1%}; {r['min01']} to {r['max01']}) "
                 f"| {r['found1']:.1f} ({r['rate1']:.1%}) | {r['any01']} / {r['any1']} of {r['any_expected']}{extra} | {r['all01']} / {r['all1']} of {r['all_expected']} |")
    L += ["", "Published, for comparison: " + "; ".join(f"{k}: {v}" for k, v in PUBLISHED.items()) + ".",
          "16S exact variants recovered (pooled): " + "; ".join(f"{r['dataset']} {r['variants']}" for r in rows if r.get("variants")) + ".", "", "Next: send this file."]
    (a.out / "RECOVERY.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    cols = ["marker", "dataset", "arm", "samples", "expected", "found01", "found1", "rate01", "rate1", "any_expected", "any01", "any1", "all_expected", "all01", "all1"]
    (a.out / "recovery.tsv").write_text("\t".join(cols) + "\n" + "".join("\t".join(str(r.get(c, "")) for c in cols) + "\n" for r in rows if "missing" not in r), encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
