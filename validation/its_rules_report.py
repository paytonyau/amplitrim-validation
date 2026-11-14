#!/usr/bin/env python
"""ITS rules night: the analysis (TEST_REGISTER.md, the section of the ITS rules night; criteria P94 to P104).

  python validation\\its_rules_report.py --out C:\\data\\runs\\its_rules --references validation\\references_its --previous results\\26_its_night
  (--truth mock3=C:\\data\\references_mock3   adds a further mock community with its own references_ITS1/ITS2.fasta; it must then be safe there too)

Reads what validation\\its_night.py --plan rules left in --out (per dataset: dada2\\<arm>\\samples.tsv and seqs.tsv, eer.json). The truths, the accuracy and the
flows are those of validation\\its_night_report.py (section 33), unchanged. Added here:
- Safe (P94, P98, P102): against std, the spurious share rises by at most 1 point and recall falls by at most 2 points, on the ATCC standard AND on the piglet mixes,
  for both markers (and on every mock given with --truth).
- Useful (P95, P99, P103): on the ATCC standard's ITS2 the yield (after chimera removal, of the input pairs) rises by at least 5 points.
- The mechanism of question A (P96): on the ATCC standard's ITS2 the share of denoised pairs that merge rises by at least 2 points.
- The flags (P97 on the first ITS study; P100, P104 on both): the median ASVs per sample rise by at most 25%.
- The consistency check: std must reproduce run 26's counts on every dataset of both runs (--previous); if not, no criterion is decided.
  P101: sw15 must reproduce run 26's sw15 likewise.
Writes ITS_RULES_SUMMARY.md, ITS_RULES_REPORT.md and its_rules_metrics.tsv to --out.
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import its_night_report as N  # noqa: E402  (the truths, accuracy and flows of section 33)

TQ_ARMS = ["tq3", "tq4", "tq5", "tq6", "tq8"]
SW_ARMS = ["sw10", "sw20", "sw25"]
PRIMARY = {"ITS1": "atcc_its1", "ITS2": "atcc_its2"}
MIXES = {"ITS1": "arfken_mix_its1", "ITS2": "arfken_mix_its2"}
STUDIES = ("its2", "its")
COUNTS = ("input", "filtered", "denoised", "merged", "nonchim")
SAFE_SPURIOUS, SAFE_RECALL, USEFUL, MECHANISM, FLAG = 1.0, 2.0, 5.0, 2.0, 0.25


def counts(path: Path) -> dict:
    with open(path, newline="", encoding="utf-8") as fh:
        return {r["sample"]: tuple(round(float(r[c])) for c in COUNTS) for r in csv.DictReader(fh, delimiter="\t")}


def merge_share(path: Path) -> float:
    rows = N.read_samples(path)
    den = sum(r["denoised"] for r in rows)
    return sum(r["merged"] for r in rows) / den if den else float("nan")


def analyse(out: Path, refs: dict, extra: dict = None) -> dict:
    """Section 33's analysis, plus the merge share per arm, the truths of further mocks and the computed maxEE."""
    R = N.analyse(out, refs)
    R["merge"], R["eer"] = {}, {}
    for lab_dir in sorted(p for p in Path(out).iterdir() if (p / "dada2").is_dir()):
        for arm_dir in sorted((lab_dir / "dada2").iterdir()):
            if (arm_dir / "samples.tsv").is_file():
                R["merge"][(lab_dir.name, arm_dir.name)] = merge_share(arm_dir / "samples.tsv")
        if (lab_dir / "eer.json").is_file():
            R["eer"][lab_dir.name] = json.loads((lab_dir / "eer.json").read_text(encoding="utf-8"))
    R["truth_labels"] = [PRIMARY["ITS1"], PRIMARY["ITS2"], MIXES["ITS1"], MIXES["ITS2"]]
    for label, folder in (extra or {}).items():                        # a further mock: its FASTA references judge it, as the ATCC standard's do
        files = [Path(folder) / f"references_{m}.fasta" for m in ("ITS1", "ITS2")]
        t = N.fasta_truth([x for f in files if f.is_file() for x in N.read_references(f)])
        for arm_dir in sorted((Path(out) / label / "dada2").iterdir()) if (Path(out) / label / "dada2").is_dir() else []:
            if (arm_dir / "seqs.tsv").is_file():
                R["acc"][(label, arm_dir.name)] = N.accuracy(N.read_seqs(arm_dir / "seqs.tsv"), t, False)
        R["truth_labels"].append(label)
    return R


def consistency(out: Path, previous: Path, arm: str) -> tuple:
    """(datasets compared, differences) for one arm: every sample's counts must equal the earlier run's."""
    compared, diffs = 0, []
    for lab_dir in sorted(p for p in Path(out).iterdir() if (p / "dada2" / arm / "samples.tsv").is_file()):
        old = Path(previous) / lab_dir.name / "dada2" / arm / "samples.tsv"
        if not old.is_file():
            continue
        compared += 1
        a, b = counts(lab_dir / "dada2" / arm / "samples.tsv"), counts(old)
        if a != b:
            bad = sorted(s for s in set(a) | set(b) if a.get(s) != b.get(s))
            diffs.append(f"{lab_dir.name}: {len(bad)} sample(s) differ, for example {bad[0]}: {a.get(bad[0])} against {b.get(bad[0])}")
    return compared, diffs


def verdicts(R: dict, gate_ok, sw15_check) -> list:
    F, A, M = R["flow"], R["acc"], R["merge"]

    def safe(arm):
        parts, ok = [], True
        for lab in R["truth_labels"]:
            if (lab, arm) not in A or (lab, "std") not in A:
                return None, f"not decided: {arm} or std is missing on {lab}"
            ds = (A[(lab, arm)]["spurious"] - A[(lab, "std")]["spurious"]) * 100
            dr = (A[(lab, "std")]["recall"] - A[(lab, arm)]["recall"]) * 100
            ok = ok and ds <= SAFE_SPURIOUS and dr <= SAFE_RECALL
            parts.append(f"{lab} spurious {ds:+.2f}, recall {-dr:+.1f}")
        return ok, "safe on every truth: " + "; ".join(parts)

    def useful(arm):
        lab = PRIMARY["ITS2"]
        if (lab, arm) not in F or (lab, "std") not in F:
            return None, "not decided: an arm is missing"
        g = (F[(lab, arm)]["yield"] - F[(lab, "std")]["yield"]) * 100
        return g >= USEFUL, f"the ATCC standard's ITS2 yield rises by at least 5 points: {g:+.1f}"

    def mechanism(arm):
        lab = PRIMARY["ITS2"]
        if (lab, arm) not in M or (lab, "std") not in M:
            return None, "not decided: an arm is missing"
        g = (M[(lab, arm)] - M[(lab, "std")]) * 100
        return g >= MECHANISM, f"the share of denoised pairs that merge rises by at least 2 points on the ATCC standard's ITS2: {g:+.2f}"

    def flag(arm, studies):
        rises = [(s, F[(s, arm)]["asv_median"] / F[(s, "std")]["asv_median"] - 1) for s in studies
                 if (s, arm) in F and (s, "std") in F and F[(s, "std")]["asv_median"]]
        return (all(r <= FLAG for _, r in rises) if rises else None), "a flag: the median ASVs per sample rise by at most 25%: " + (", ".join(f"{s} {r:+.0%}" for s, r in rises) or "not run")

    V = []
    for arm in TQ_ARMS:
        V += [(f"P94 {arm}", *safe(arm)), (f"P95 {arm}", *useful(arm)), (f"P96 {arm}", *mechanism(arm)), (f"P97 {arm}", *flag(arm, ("its",)))]
    for arm in SW_ARMS:
        V += [(f"P98 {arm}", *safe(arm)), (f"P99 {arm}", *useful(arm)), (f"P100 {arm}", *flag(arm, STUDIES))]
    V.append(("P101", sw15_check[0], sw15_check[1]))
    V += [("P102", *safe("eer")), ("P103", *useful("eer")), ("P104", *flag("eer", STUDIES))]
    if not gate_ok:
        V = [(c, None, "not decided (the consistency check failed): " + t) for c, _, t in V]
    # the advice of section 35
    ok = lambda arm: safe(arm)[0] and useful(arm)[0]
    tq_ok = [a for a in TQ_ARMS if gate_ok and ok(a)]
    sw_ok = [a for a in SW_ARMS + ["sw15"] if gate_ok and ok(a)]
    lab2 = PRIMARY["ITS2"]
    R["advice_tq"] = max(tq_ok, key=lambda a: int(a[2:])) if tq_ok else None
    R["advice_sw"] = max(sw_ok, key=lambda a: F[(lab2, a)]["yield"]) if sw_ok else None
    both = [a for a in (R["advice_tq"], R["advice_sw"]) if a]
    R["advice"] = max(both, key=lambda a: F[(lab2, a)]["yield"]) if both else None
    R["advice_eer"] = bool(gate_ok and safe("eer")[0] and useful("eer")[0])
    return V


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--previous", type=Path, required=True, help="run 26's folder (results/26_its_night or C:\\data\\runs\\26_its_night)")
    ap.add_argument("--truth", action="append", default=[], metavar="LABEL=FOLDER", help="a further mock community and the folder of its references")
    a = ap.parse_args(argv)
    refs = {m: N.read_references(a.references / f"references_{m}.fasta") for m in ("ITS1", "ITS2") if (a.references / f"references_{m}.fasta").is_file()}
    extra = dict(t.split("=", 1) for t in a.truth)
    if not a.out.is_dir() or len(refs) < 2 or not a.previous.is_dir():
        print("Not found: --out, --references (references_ITS1/ITS2.fasta) or --previous")
        return 2
    R = analyse(a.out, refs, extra)
    n_std, d_std = consistency(a.out, a.previous, "std")
    n_sw, d_sw = consistency(a.out, a.previous, "sw15")
    gate_ok = n_std > 0 and not d_std
    sw15_check = (None, "not decided: no dataset to compare") if not n_sw else (not d_sw, f"sw15 reproduces run 26's sw15 on {n_sw} datasets" + ("" if not d_sw else ": " + "; ".join(d_sw)))
    V = verdicts(R, gate_ok, sw15_check)
    S = ["# ITS rules night: summary", "", "The criteria (P94 to P104) are in TEST_REGISTER.md and were fixed before this ran.", "",
         f"**The consistency check:** std reproduces run 26 on {n_std} datasets" + (" exactly." if gate_ok else ": FAILED, so no criterion is decided: " + "; ".join(d_std or ["no dataset to compare"])), ""]
    S += [N.line(c, ok, t) for c, ok, t in V]
    S += ["", "**What follows (section 35):** "
          + (f"`truncQ` option: `{R['advice_tq']}`. " if R["advice_tq"] else "no truncQ value from 3 to 8 is safe and useful: truncQ stays 2. ")
          + (f"The trim strength: `{R['advice_sw']}`. " if R["advice_sw"] else "No trim strength is safe and useful. ")
          + (f"Advised for poor read ends: `{R['advice']}`. " if R["advice"] else "")
          + ("The computed maxEE (`eer`) may be used for ITS." if R["advice_eer"] else "ITS keeps maxEE 2,2."), "", "Next: send this file."]
    F, A, Mg = R["flow"], R["acc"], R["merge"]
    L = ["# ITS rules night: report", "", "## The computed maxEE per dataset (the arm eer)", "", "| Dataset | P75 of EE, forward / reverse | maxEE | Reads used |", "|---|---|---|---|"]
    L += [f"| {lab} | {e['percentile'][0]:.2f} / {e['percentile'][1]:.2f} | {e['max_ee'][0]},{e['max_ee'][1]} | {e['reads'][0]:,} / {e['reads'][1]:,} |" for lab, e in sorted(R["eer"].items())]
    L += ["", "## Flows and merging (shares of the input pairs; merge share: of the denoised pairs)", "",
          "| Dataset | Arm | Input | Pass | Merged | After chimera removal | Merge share | ASVs |", "|---|---|---|---|---|---|---|---|"]
    L += [f"| {lab} | {arm} | {f['input']:,.0f} | {N.pct(f['pass'])} | {N.pct(f['merged'])} | {N.pct(f['yield'])} | {N.pct(Mg.get((lab, arm), float('nan')))} | {f['asv_median']:.0f} |"
          for (lab, arm), f in sorted(F.items())]
    L += ["", "## Accuracy (means over samples)", "", "| Dataset | Arm | Spurious | Recall | Exact | Members found | Unjudged |", "|---|---|---|---|---|---|---|"]
    L += [f"| {lab} | {arm} | {N.pct(x['spurious'])} | {N.pct(x['recall'])} | {N.pct(x['exact'])} | {x['found']:.1f} | {N.pct(x['unjudged'])} |" for (lab, arm), x in sorted(A.items())]
    L += ["", "## Verdicts", ""] + [N.line(c, ok, t) for c, ok, t in V]
    Mt = ["dataset\tarm\tinput\tpass\tmerged\tyield\tmerge_share\tasv_median\tspurious\trecall\texact\tfound\tunjudged"]
    for (lab, arm), f in sorted(F.items()):
        x = A.get((lab, arm), {})
        Mt.append("\t".join([lab, arm] + [f"{f[k]:.6g}" for k in ("input", "pass", "merged", "yield")] + [f"{Mg.get((lab, arm), float('nan')):.6g}", f"{f['asv_median']:.6g}"]
                            + [f"{x[k]:.6g}" if k in x else "" for k in ("spurious", "recall", "exact", "found", "unjudged")]))
    (a.out / "ITS_RULES_SUMMARY.md").write_text("\n".join(S) + "\n", encoding="utf-8")
    (a.out / "ITS_RULES_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    (a.out / "its_rules_metrics.tsv").write_text("\n".join(Mt) + "\n", encoding="utf-8")
    print("\n".join(S))
    print(f"written: {a.out / 'ITS_RULES_SUMMARY.md'}, ITS_RULES_REPORT.md, its_rules_metrics.tsv")
    return 0 if gate_ok and all(ok is not None for _, ok, _ in V) else 2


if __name__ == "__main__":
    sys.exit(main())
