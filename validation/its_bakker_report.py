#!/usr/bin/env python
"""The third mock community (Bakker 2018): the analysis (TEST_REGISTER.md, the section of the Bakker night; criteria P105 to P114).

  python validation\\its_bakker_report.py --out C:\\data\\runs\\its_bakker --references validation\\references_its --bakker-references validation\\references_bakker
        --previous results\\26_its_night --previous-rules results\\27_its_rules

Reads what validation\\its_night.py --plan bakker left in --out. The truths, accuracy and flows are those of its_night_report.py (section 33); the consistency check
is that of its_rules_report.py (section 35). Truths: the ATCC standard (NCBI references), the piglet mixes (the isolates' own sequences), and Bakker's three
datasets (bakker_pcr_its1, bakker_pol_its1, bakker_pol_its2) against the sequences of the 19 isolates (Appendix S2 of the paper, cut with each dataset's primers
by bakker_prep.py into references_<label>.fasta); every one of the 19 members is expected in every Bakker sample. Bakker's samples are also reported by condition: the standard workflow, a high
PCR cycle number (HC) and a long extension (LE) on the first run; GoTaq, Phusion and PreAmp on the second.
"Safe" and "useful" are those of sections 33 and 35 (spurious up at most 1 point, recall down at most 2 points; yield up at least 5 points).
Writes ITS_BAKKER_SUMMARY.md, ITS_BAKKER_REPORT.md and its_bakker_metrics.tsv to --out.
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import its_night_report as N  # noqa: E402
import its_rules_report as RR  # noqa: E402

ATCC = {"ITS1": "atcc_its1", "ITS2": "atcc_its2"}
MIXES = ["arfken_mix_its1", "arfken_mix_its2"]
BAKKER = {"bakker_pcr_its1": "ITS1", "bakker_pol_its1": "ITS1", "bakker_pol_its2": "ITS2"}
STUDIES = ("its2", "its")
OLD = ["atcc_its1", "atcc_its2", "arfken_its1", "arfken_its2", "arfken_mix_its1", "arfken_mix_its2", "its2", "its"]
SAFE_SPURIOUS, SAFE_RECALL, USEFUL, LOW_PASS, TRIGGER_LOSS, EER_LOSS, COMBO_GAIN, COMBO_LOSS, FLAG = 1.0, 2.0, 5.0, 0.5, 2.0, 1.0, 2.0, 1.0, 0.25


def condition(sample: str) -> str:
    s = sample.lower()
    for token, name in (("_hc", "high cycles"), ("_le", "long extension"), ("gotaq", "GoTaq"), ("phusion", "Phusion"), ("preamp", "PreAmp")):
        if token in s:
            return name
    return "standard"


def analyse(out: Path, refs: dict, bakker_refs: dict) -> dict:
    R = N.analyse(out, refs)
    R["by_condition"], R["merge"] = {}, {}
    for lab_dir in sorted(p for p in Path(out).iterdir() if (p / "dada2").is_dir()):
        for arm_dir in sorted((lab_dir / "dada2").iterdir()):
            if (arm_dir / "samples.tsv").is_file():
                R["merge"][(lab_dir.name, arm_dir.name)] = RR.merge_share(arm_dir / "samples.tsv")   # the share of denoised pairs that merge
    for lab in BAKKER:
        if lab not in bakker_refs:
            continue
        truth = N.fasta_truth(bakker_refs[lab])
        d = Path(out) / lab / "dada2"
        for arm_dir in sorted(d.iterdir()) if d.is_dir() else []:
            if not (arm_dir / "seqs.tsv").is_file():
                continue
            seqs = N.read_seqs(arm_dir / "seqs.tsv")
            samples = [r["sample"] for r in N.read_samples(arm_dir / "samples.tsv")]
            for s in samples:
                seqs.setdefault(s, {})                                         # a sample that kept nothing counts, with recall 0
            R["acc"][(lab, arm_dir.name)] = N.accuracy(seqs, truth, False)
            for cond in sorted({condition(s) for s in seqs}):
                R["by_condition"][(lab, arm_dir.name, cond)] = N.accuracy({s: v for s, v in seqs.items() if condition(s) == cond}, truth, False)
    R["truth_labels"] = [ATCC["ITS1"], ATCC["ITS2"]] + MIXES + list(BAKKER)
    return R


def verdicts(R: dict, gate_ok: bool, reproduce: tuple, previous_ok=True) -> list:
    F, A = R["flow"], R["acc"]
    labels = R["truth_labels"]

    def safe(arm, only=None):
        parts, ok = [], True
        for lab in (only or labels):
            if (lab, arm) not in A or (lab, "std") not in A:
                return None, f"not decided: {arm} or std is missing on {lab}"
            ds = (A[(lab, arm)]["spurious"] - A[(lab, "std")]["spurious"]) * 100
            dr = (A[(lab, "std")]["recall"] - A[(lab, arm)]["recall"]) * 100
            ok = ok and ds <= SAFE_SPURIOUS and dr <= SAFE_RECALL
            parts.append(f"{lab} spurious {ds:+.2f}, recall {-dr:+.1f}")
        return ok, "; ".join(parts)

    def gain(lab, arm, base="std"):
        return (F[(lab, arm)]["yield"] - F[(lab, base)]["yield"]) * 100 if (lab, arm) in F and (lab, base) in F else None

    V = []
    ok, txt = safe("sw15")
    V.append(("P105", ok, "the trim (sw15) is safe on every truth: " + txt))
    low = [lab for lab in labels if (lab, "std") in F and F[(lab, "std")]["pass"] < LOW_PASS]
    gains = [(lab, gain(lab, "sw15")) for lab in low]
    V.append(("P106", (all(g is not None and g >= USEFUL for _, g in gains) if gains else None),
              "on every truth dataset where fewer than 50% of the pairs pass at maxEE 2,2, the trim raises the yield by at least 5 points: "
              + (", ".join(f"{lab} {g:+.1f}" for lab, g in gains) if gains else "not decided: no such dataset")))
    high = [lab for lab in labels if (lab, "std") in F and F[(lab, "std")]["pass"] >= LOW_PASS]
    parts, ok = [], True if high else None
    for lab in high:
        s_ok, _ = safe("sw15", [lab])
        g = gain(lab, "sw15")
        if s_ok is None or g is None:
            ok = None
            break
        ok = ok and s_ok and g >= -TRIGGER_LOSS
        parts.append(f"{lab} {'safe' if s_ok else 'NOT safe'}, yield {g:+.1f}")
    V.append(("P107", ok, "the 50% trigger: where at least 50% of the pairs pass at maxEE 2,2, the trim is safe and loses at most 2 points of yield: " + ("; ".join(parts) or "no such dataset")))
    ok, txt = safe("eer")
    V.append(("P108", ok, "the computed maxEE (eer) is safe on every truth: " + txt))
    losses = [(lab, gain(lab, "eer")) for lab in labels if gain(lab, "eer") is not None]
    V.append(("P109", (all(g >= -EER_LOSS for _, g in losses) if losses else None),
              "the computed maxEE loses at most 1 point of yield on every truth dataset: " + ", ".join(f"{lab} {g:+.1f}" for lab, g in losses)))
    ok, txt = safe("sweer")
    V.append(("P110", ok, "the trim with maxEE computed on the trimmed reads (sweer) is safe on every truth: " + txt))
    g_atcc = gain(ATCC["ITS2"], "sweer", "sw15")
    worse = [(lab, gain(lab, "sweer", "sw15")) for lab in labels if gain(lab, "sweer", "sw15") is not None]
    V.append(("P111", (g_atcc is not None and g_atcc >= COMBO_GAIN and all(g >= -COMBO_LOSS for _, g in worse)) if g_atcc is not None else None,
              f"sweer against the trim alone: the ATCC standard's ITS2 yield {g_atcc:+.1f} points (at least +2), and nowhere lower than -1: "
              + ", ".join(f"{lab} {g:+.1f}" for lab, g in worse) if g_atcc is not None else "not decided: an arm is missing"))
    for arm in ("sw15", "eer", "sweer"):
        rises = [(s, F[(s, arm)]["asv_median"] / F[(s, "std")]["asv_median"] - 1) for s in STUDIES if (s, arm) in F and (s, "std") in F and F[(s, "std")]["asv_median"]]
        V.append((f"P112 {arm}", (all(r <= FLAG for _, r in rises) if rises else None), "a flag: the median ASVs per sample rise by at most 25%: "
                  + (", ".join(f"{s} {r:+.0%}" for s, r in rises) or "not run")))
    V.append(("P113", reproduce[0], reproduce[1]))
    unsafe = safe("tq10")
    V.append(("P114", (unsafe[0] is False) if unsafe[0] is not None else None, "the known-unsafe reference (truncQ 10) is still unsafe somewhere: " + unsafe[1]))
    if not gate_ok:
        V = [(c, None, "not decided (the consistency check failed): " + t) for c, _, t in V]
    R["verdict"] = {c: ok for c, ok, _ in V}
    return V


def follows(R: dict) -> list:
    v = R["verdict"]
    L = []
    if v.get("P105") and v.get("P106"):
        L.append("The trim advice (window 4, mean quality 15) holds on a third mock community, with other read lengths and error-prone PCR conditions.")
    elif v.get("P105") is False:
        L.append("The trim was not safe on every truth: the advice is restricted to what it was safe on, and the report names where it was not.")
    if v.get("P107"):
        L.append("The trim did not harm data that pass well at maxEE 2,2: the 50% trigger is not needed for safety; it may stay as a guide to when the trim helps.")
    elif v.get("P107") is False:
        L.append("The trim harmed some data that pass well at maxEE 2,2: the 50% trigger stays.")
    if v.get("P108") and v.get("P109"):
        L.append("The computed maxEE for ITS holds on a third mock community.")
    elif v.get("P108") is False or v.get("P109") is False:
        L.append("The computed maxEE did not hold on every truth: it is not used in 0.3.")
    if v.get("P110") and v.get("P111"):
        L.append("0.3 may combine the trim with maxEE computed on the trimmed reads.")
    elif v.get("P110") and v.get("P111") is False:
        L.append("The combination is safe but not better than the trim alone: the trim and the computed maxEE stay separate.")
    elif v.get("P110") is False:
        L.append("The combination is not safe: the trim and the computed maxEE are never combined.")
    return L or ["No criterion was decided."]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--bakker-references", type=Path, required=True)
    ap.add_argument("--previous", type=Path, required=True, help="run 26 (std and sw15)")
    ap.add_argument("--previous-rules", type=Path, required=True, help="run 27 (eer)")
    ap.add_argument("--overhang-max", type=float, default=None, help="the re-run's gate (section 39): no criterion is decided if more than this share of the reads of "
                    "std on a Bakker dataset carries primer overhangs (bakker_overhang.py)")
    a = ap.parse_args(argv)
    refs = {m: N.read_references(a.references / f"references_{m}.fasta") for m in ("ITS1", "ITS2") if (a.references / f"references_{m}.fasta").is_file()}
    brefs = {lab: N.read_references(a.bakker_references / f"references_{lab}.fasta") for lab in BAKKER if (a.bakker_references / f"references_{lab}.fasta").is_file()}
    if not a.out.is_dir() or len(refs) < 2 or len(brefs) < len(BAKKER) or not a.previous.is_dir() or not a.previous_rules.is_dir():
        print("Not found: --out, --references (references_ITS1/ITS2.fasta), --bakker-references (references_<label>.fasta of the three Bakker datasets), --previous or --previous-rules")
        return 2
    R = analyse(a.out, refs, brefs)
    n_std, d_std = RR.consistency(a.out, a.previous, "std")
    n_sw, d_sw = RR.consistency(a.out, a.previous, "sw15")
    n_ee, d_ee = RR.consistency(a.out, a.previous_rules, "eer")
    gate_ok = n_std > 0 and not d_std
    over_line = ""
    if a.overhang_max is not None:
        import bakker_overhang as BO
        oh = BO.run_check(a.out, a.bakker_references)
        bad = {k: v["overhang"] for k, v in oh.items() if v["overhang"] > a.overhang_max}
        over_line = (f"**The overhang check:** reads of std carrying primer overhangs: " + ", ".join(f"{k} {v['overhang']:.2%}" for k, v in oh.items())
                     + (f" (at most {a.overhang_max:.0%}: passed)." if not bad and len(oh) == len(BAKKER) else f": FAILED (more than {a.overhang_max:.0%}, or a dataset missing), so no criterion is decided."))
        gate_ok = gate_ok and not bad and len(oh) == len(BAKKER)
    reproduce = ((not d_sw and not d_ee) if (n_sw and n_ee) else None,
                 f"sw15 reproduces run 26 on {n_sw} datasets and eer reproduces run 27 on {n_ee}" + ("" if not (d_sw or d_ee) else ": " + "; ".join(d_sw + d_ee)))
    V = verdicts(R, gate_ok, reproduce)
    S = ["# The Bakker night: summary", "", "The criteria (P105 to P114) are in TEST_REGISTER.md and were fixed before this ran.", "",
         f"**The consistency check:** std reproduces run 26 on {n_std} datasets" + (" exactly." if (n_std > 0 and not d_std) else ": FAILED, so no criterion is decided: " + "; ".join(d_std or ["nothing to compare"])), ""]
    if over_line:
        S += [over_line, ""]
    S += [N.line(c, ok, t) for c, ok, t in V]
    S += ["", "**What follows (section 37):** " + " ".join(follows(R)), "", "Next: send this file."]
    F, A = R["flow"], R["acc"]
    L = ["# The Bakker night: report", "", "## Bakker's datasets by condition (means over the samples of the condition)", "",
         "| Dataset | Arm | Condition | Spurious | Recall | Members found (of 19) | Samples |", "|---|---|---|---|---|---|---|"]
    L += [f"| {lab} | {arm} | {cond} | {N.pct(x['spurious'])} | {N.pct(x['recall'])} | {x['found']:.1f} | {x['samples']} |" for (lab, arm, cond), x in sorted(R["by_condition"].items())]
    L += ["", "## Flows (shares of the input pairs; merge share: of the denoised pairs, which falls when a trim takes the overlap away)", "",
          "| Dataset | Arm | Input | Pass | Merged | After chimera removal | Merge share | ASVs |", "|---|---|---|---|---|---|---|---|"]
    L += [f"| {lab} | {arm} | {f['input']:,.0f} | {N.pct(f['pass'])} | {N.pct(f['merged'])} | {N.pct(f['yield'])} | {N.pct(R['merge'].get((lab, arm), float('nan')))} | {f['asv_median']:.0f} |"
          for (lab, arm), f in sorted(F.items())]
    L += ["", "## Accuracy (means over samples)", "", "| Dataset | Arm | Spurious | Recall | Exact | Members found | Unjudged |", "|---|---|---|---|---|---|---|"]
    L += [f"| {lab} | {arm} | {N.pct(x['spurious'])} | {N.pct(x['recall'])} | {N.pct(x['exact'])} | {x['found']:.1f} | {N.pct(x['unjudged'])} |" for (lab, arm), x in sorted(A.items())]
    L += ["", "## Verdicts", ""] + [N.line(c, ok, t) for c, ok, t in V]
    M = ["dataset\tarm\tinput\tpass\tmerged\tyield\tasv_median\tspurious\trecall\texact\tfound\tunjudged"]
    for (lab, arm), f in sorted(F.items()):
        x = A.get((lab, arm), {})
        M.append("\t".join([lab, arm] + [f"{f[k]:.6g}" for k in ("input", "pass", "merged", "yield", "asv_median")] + [f"{x[k]:.6g}" if k in x else "" for k in ("spurious", "recall", "exact", "found", "unjudged")]))
    (a.out / "ITS_BAKKER_SUMMARY.md").write_text("\n".join(S) + "\n", encoding="utf-8")
    (a.out / "ITS_BAKKER_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    (a.out / "its_bakker_metrics.tsv").write_text("\n".join(M) + "\n", encoding="utf-8")
    print("\n".join(S))
    print(f"written: {a.out / 'ITS_BAKKER_SUMMARY.md'}, ITS_BAKKER_REPORT.md, its_bakker_metrics.tsv")
    return 0 if gate_ok and all(ok is not None for _, ok, _ in V) else 2


if __name__ == "__main__":
    sys.exit(main())
