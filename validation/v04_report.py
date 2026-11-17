#!/usr/bin/env python
"""The v0.4 night (TEST_REGISTER.md, the section of the v0.4 night): gates, criteria P117 to P130, and what follows.

  python validation\\v04_report.py --out C:\\data\\runs\\v04_night --references validation\\references_its --bakker-references validation\\references_bakker
      --pauvert-references C:\\data\\24_mock_ITS1_pauvert_dataverse\\references --previous results\\29_bakker_rerun --cli-16s results\\24_night10_closing_2026-11-06\\cli

Writes V04_SUMMARY.md, V04_REPORT.md and v04_metrics.tsv to --out.
"""
import argparse
import csv
import gzip
import json
import math
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import its_night_report as N  # noqa: E402
import its_rules_report as RR  # noqa: E402
import its_bakker_report as BR  # noqa: E402
import bakker_overhang as BO  # noqa: E402
import bakker_prep as BP  # noqa: E402

PAUVERT = "pauvert_its1"
ITS1F, ITS2 = "CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"
TRIMS = ["sw15"] + [f"w{w}" for w in (2, 3, 5, 6, 8, 10, 12, 16, 20)] + ["q12", "q18"]
MINLENS, TQS, FWD = ["ml20", "ml100", "ml150"], ["tq3", "tq5", "tq7"], ["fwd", "fwd15"]
TRUTH = ["atcc_its1", "atcc_its2", "arfken_mix_its1", "arfken_mix_its2", "bakker_pcr_its1", "bakker_pol_its1", "bakker_pol_its2", PAUVERT]
ITS_ALL = list(BR.OLD) + list(BR.BAKKER) + [PAUVERT]
EARLIER = ITS_ALL[:-1]                                   # the eleven datasets of run 29
STUDIES = ("its2", "its")
SP_UP, REC_DOWN, TIE, FLAG = 0.01, 0.02, 0.005, 0.25
EPS = 1e-9                                               # a result exactly at a limit is within it (binary fractions are not exact)
PRED_ARM = {("truncQ", "2"): "std", ("window", "4"): "sw15", ("quality", "15"): "sw15"}


# ------------------------------------------------------------------------------------------------ Pauvert's truth: complete and partial references
class PartialTruth(N.Truth):
    """Every kept strain judges (true within 2 edits); a partial reference is compared with the sequence's first bases, as many as it has (give or take 2);
    recall counts only the expected members."""

    def __init__(self, refs: list, expected: list):
        super().__init__([(m, s) for m, k, s in refs])
        self.partial = {(m, s) for m, k, s in refs if k == "partial"}
        self.members = sorted(set(expected))

    def dist(self, q: str, m: str, r: str, prefix: bool) -> int:
        E = N.EDITS
        partial = (m, r) in self.partial
        if not partial:
            return N.within(q, r[:len(q)] if prefix else r, E)
        if prefix and len(q) < len(r):                                  # a forward read shorter than the partial reference
            return N.within(q, r[:len(q)], E)
        return min(N.within(q[:len(r) + d], r, E) for d in range(-E, E + 1) if 0 < len(r) + d <= len(q))  if len(q) >= len(r) - E else E + 1

    def classify(self, q: str, prefix: bool = False):
        key = (q, prefix)
        if key not in self._cache:
            best = None
            for m, r in self.expected:
                d = self.dist(q, m, r, prefix)
                if d <= N.EDITS and (best is None or d < best[1]):
                    best = (m, d)
            self._cache[key] = ("true", best[0], best[1] == 0) if best else ("spurious", None, False)
        return self._cache[key]


def pauvert_truth(folder: Path) -> PartialTruth:
    refs = N.read_references(Path(folder) / f"references_{PAUVERT}.fasta")
    with open(Path(folder) / f"expected_{PAUVERT}.tsv", encoding="utf-8") as fh:
        expected = [r["member"] for r in csv.DictReader(fh, delimiter="\t")]
    return PartialTruth(refs, expected)


# ------------------------------------------------------------------------------------------------ the analysis
def fastq_bases(folder: Path) -> tuple:
    """(reads, bases) over every FASTQ of a folder (both mates)."""
    n = b = 0
    for f in sorted(Path(folder).glob("*.fastq*")):
        op = gzip.open if f.name.endswith(".gz") else open
        with op(f, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    n += 1
                    b += len(line.rstrip("\n"))
    return n, b


def bases_cut(d: Path) -> dict:
    """Mean bases cut per read by each trim of the night (both mates): from cut.json, which the night writes before deleting the trimmed copies, or
    counted from the copies when they are still there."""
    if (d / "cut.json").is_file():
        c = json.loads((d / "cut.json").read_text(encoding="utf-8"))["mean_bases_cut"]
        return {arm: (float("nan") if v is None else v) for arm, v in c.items() if arm in TRIMS}
    n0, b0 = fastq_bases(d / "reads")
    out = {}
    for arm in TRIMS:
        src = d / ("reads_sw15" if arm == "sw15" else f"reads_{arm}")
        if src.is_dir() and n0:
            n1, b1 = fastq_bases(src)
            out[arm] = (b0 - b1) / n0 if n1 == n0 else float("nan")
    return out


def judge(R: dict, out: Path, lab: str, truth, mix: bool, arms=None):
    d = Path(out) / lab / "dada2"
    for arm_dir in sorted(d.iterdir()) if d.is_dir() else []:
        if arms is not None and arm_dir.name not in arms:
            continue
        if (arm_dir / "seqs.tsv").is_file() and (arm_dir / "samples.tsv").is_file():
            seqs = N.read_seqs(arm_dir / "seqs.tsv")
            for r in N.read_samples(arm_dir / "samples.tsv"):
                seqs.setdefault(r["sample"], {})                               # a sample that kept nothing counts, with recall 0
            R["acc"][(lab, arm_dir.name)] = N.accuracy(seqs, truth, mix, prefix=arm_dir.name in FWD)


def analyse(out: Path, refs: dict, bakker_refs: dict, pauvert_dir: Path) -> dict:
    R = BR.analyse(out, refs, bakker_refs)                              # flows, merge shares, ATCC, mixes, Bakker (run 29's)
    out = Path(out)
    for marker, lab in N.PRIMARY.items():                               # forward-only arms as prefixes, on every truth dataset
        if marker in refs:
            judge(R, out, lab, N.fasta_truth(refs[marker]), False, FWD)
    for marker, (iso, mix) in N.SECOND.items():
        p = out / iso / "dada2" / "std" / "seqs.tsv"
        if p.is_file():
            judge(R, out, mix, N.isolate_truth(N.read_seqs(p)), True, FWD)
    for lab in BR.BAKKER:
        if lab in bakker_refs:
            judge(R, out, lab, N.fasta_truth(bakker_refs[lab]), False, FWD)
    if (Path(pauvert_dir) / f"references_{PAUVERT}.fasta").is_file():
        R["pauvert_expected"] = len(pauvert_truth(pauvert_dir).members)
        judge(R, out, PAUVERT, pauvert_truth(pauvert_dir), False)
    R["truth_labels"] = TRUTH
    R["cut"] = {lab: bases_cut(out / lab) for lab in ITS_ALL if (out / lab / "reads").is_dir() or (out / lab / "cut.json").is_file()}
    R["est"] = {}
    for lab in ITS_ALL:
        for arm, name in (("std", "estimates.tsv"), ("sw15", "estimates_sw15.tsv")):
            e, fl = out / lab / name, R["flow"].get((lab, arm))
            if e.is_file() and fl and fl["input"]:
                rows = N.read_samples(e)
                R["est"][(lab, arm)] = {k: sum(r[f"merged_d{k}"] for r in rows) / fl["input"] for k in N.DIFFS}
    R["pred"] = {lab: json.loads((out / lab / "predict.json").read_text(encoding="utf-8")) for lab in ITS_ALL if (out / lab / "predict.json").is_file()}
    return R


# ------------------------------------------------------------------------------------------------ the gates
def overhang_gate(out: Path, bakker_dir: Path, pauvert_dir: Path, limit: float = 0.01, min_bases: int = 3) -> tuple:
    shares = {lab: v["overhang"] for lab, v in BO.run_check(out, bakker_dir, "std", min_bases).items()}
    p = Path(out) / PAUVERT / "dada2" / "std" / "seqs.tsv"
    if p.is_file():
        seqs = {}
        for s in N.read_seqs(p).values():
            for q, a in s.items():
                seqs[q] = seqs.get(q, 0) + a
        shares[PAUVERT] = BO.classify(seqs, BO.read_refs(Path(pauvert_dir) / f"references_{PAUVERT}.fasta"), ITS1F, ITS2, min_bases)["overhang"]
    need = list(BR.BAKKER) + [PAUVERT]
    ok = all(lab in shares and shares[lab] <= limit for lab in need)
    return ok, shares


def consistency_gate(out: Path, previous: Path) -> tuple:
    texts, ok = [], True
    for arm in ("std", "sw15", "eer"):
        n, diffs = RR.consistency(out, previous, arm)
        ok = ok and n == len(EARLIER) and not diffs
        texts.append(f"{arm}: {n} datasets compared" + (f", differences: {'; '.join(diffs)}" if diffs else ", identical"))
    return ok, texts


# ------------------------------------------------------------------------------------------------ the criteria
def acc(R, lab, arm):
    return R["acc"].get((lab, arm))


def yld(R, lab, arm):
    f = R["flow"].get((lab, arm))
    return f["yield"] if f else float("nan")


def safe(R, lab, arm, ref) -> tuple:
    """(safe, text) of arm against ref on one truth dataset: spurious up at most 1 point, recall down at most 2."""
    a, b = acc(R, lab, arm), acc(R, lab, ref)
    if not a or not b or any(math.isnan(x) for x in (a["spurious"], a["recall"], b["spurious"], b["recall"])):
        return None, f"{lab} missing"
    ds, dr = (a["spurious"] - b["spurious"]) * 100, (a["recall"] - b["recall"]) * 100          # in points, as runs 27 to 29 computed them
    return ds <= SP_UP * 100 + EPS and dr >= -REC_DOWN * 100 - EPS, f"{lab} spurious {ds:+.2f}, recall {dr:+.1f}"


def safe_all(R, arm, ref, labs=TRUTH) -> tuple:
    res = [safe(R, lab, arm, ref) for lab in labs]
    if any(ok is None for ok, _ in res):
        return None, "; ".join(t for _, t in res)
    return all(ok for ok, _ in res), "; ".join(t for _, t in res)


def differs(R, lab, arm, ref) -> bool:
    a, b = R["flow"].get((lab, arm)), R["flow"].get((lab, ref))
    return bool(a and b and (abs(a["pass"] - b["pass"]) > 1e-9 or abs(a["yield"] - b["yield"]) > 1e-9))


def mean_cut(R, arm) -> float:
    v = [R["cut"][lab][arm] for lab in TRUTH if lab in R["cut"] and arm in R["cut"][lab] and not math.isnan(R["cut"][lab][arm])]
    return statistics.fmean(v) if v else float("nan")


def best_trim(R) -> dict:
    """P130: among the trims safe against std on every truth dataset, the highest mean yield over the truth datasets; within TIE of it, the fewest bases cut."""
    rows = []
    for arm in TRIMS:
        ok, _ = safe_all(R, arm, "std")
        ys = [yld(R, lab, arm) for lab in TRUTH]
        rows.append({"arm": arm, "safe": ok, "mean_yield": statistics.fmean(ys) if not any(math.isnan(y) for y in ys) else float("nan"),
                     "atcc_its2": yld(R, "atcc_its2", arm), "cut": mean_cut(R, arm),
                     "merge": statistics.fmean([R["merge"].get((lab, arm), float("nan")) for lab in TRUTH])})
    cands = [r for r in rows if r["safe"] and not math.isnan(r["mean_yield"])]
    if not cands:
        return {"rows": rows, "best": None}
    top = max(r["mean_yield"] for r in cands)
    tied = [r for r in cands if r["mean_yield"] >= top - TIE]
    best = min(tied, key=lambda r: (r["cut"] if not math.isnan(r["cut"]) else float("inf"), -r["mean_yield"]))
    return {"rows": rows, "best": best["arm"], "tied": [r["arm"] for r in tied]}


def verdicts(R: dict, gates_ok: bool) -> list:
    """[(id, ok, text)]; ok is None when not decided."""
    V = []
    def add(cid, ok, text):
        V.append((cid, ok if gates_ok else None, text))
    add("P117", *safe_all(R, "sw15", "std", [PAUVERT]))
    add("P118", *safe_all(R, "eer", "std", [PAUVERT]))
    for arm in TRIMS[1:]:
        add(f"P119 {arm}", *safe_all(R, arm, "std"))
        a2 = yld(R, "atcc_its2", arm) - yld(R, "atcc_its2", "sw15")
        worst = min(yld(R, lab, arm) - yld(R, lab, "sw15") for lab in TRUTH)
        add(f"P120 {arm}", None if math.isnan(a2) or math.isnan(worst) else (a2 >= 0.02 and worst >= -0.01),
            f"ATCC ITS2 yield {a2 * 100:+.1f} points against sw15 (at least +2), lowest {worst * 100:+.1f} (at least -1)")
    for arm in MINLENS:
        add(f"P121 {arm}", *safe_all(R, arm, "sw15"))
    for arm in TQS:
        labs = [lab for lab in TRUTH if differs(R, lab, arm, "std")]
        if not labs:
            V.append((f"P122 {arm}", None, "changes nothing on any truth dataset (not judged)"))
            V.append((f"P123 {arm}", None, "changes nothing on any truth dataset (not judged)"))
            continue
        add(f"P122 {arm}", *safe_all(R, arm, "std", labs))
        gain = max(yld(R, lab, arm) - yld(R, lab, "std") for lab in labs)
        add(f"P123 {arm}", gain >= 0.05, f"the largest yield gain where it differs: {gain * 100:+.1f} points ({', '.join(labs)})")
    ok, text = safe_all(R, "fwd", "std")
    rec = [(lab, acc(R, lab, "fwd"), acc(R, lab, "std")) for lab in TRUTH]
    rec_ok = all(a and b and a["recall"] >= b["recall"] for _, a, b in rec) if all(a and b for _, a, b in rec) else None
    add("P124", None if ok is None or rec_ok is None else (ok and rec_ok), f"{text}; recall not below std's everywhere: {rec_ok}")
    ok, text = safe_all(R, "fwd15", "fwd")
    ratio = []
    for lab in TRUTH:
        a, b, fa, fb = acc(R, lab, "fwd15"), acc(R, lab, "fwd"), R["flow"].get((lab, "fwd15")), R["flow"].get((lab, "fwd"))
        if a and b and fa and fb and a["found"] and b["found"]:
            ratio.append((lab, (fa["asv_median"] / a["found"]) / (fb["asv_median"] / b["found"]) - 1))
    flag = all(r <= FLAG for _, r in ratio) if len(ratio) == len(TRUTH) else None
    add("P129", None if ok is None or flag is None else (ok and flag),
        f"{text}; ASVs per member found against fwd: {', '.join(f'{lab} {r * 100:+.0f}%' for lab, r in ratio)} (at most +25%)")
    for arm in ("std", "sw15"):
        diffs = {k: [R["est"][(lab, arm)][k] - R["flow"][(lab, arm)]["merged"] for lab in ITS_ALL if (lab, arm) in R["est"] and (lab, arm) in R["flow"]] for k in (0, 10)}
        if len(diffs[0]) < len(ITS_ALL):
            V.append((f"P125 {arm}", None, f"estimates for {len(diffs[0])} of {len(ITS_ALL)} datasets"))
            continue
        mae = {k: statistics.fmean(abs(x) for x in v) for k, v in diffs.items()}
        add(f"P125 {arm}", mae[0] <= mae[10] - 0.05 and max(abs(x) for x in diffs[0]) <= 0.10,
            f"mean absolute difference from DADA2's merged share: 0 differences {mae[0] * 100:.1f} points, 10 differences {mae[10] * 100:.1f}; largest at 0: "
            f"{max(abs(x) for x in diffs[0]) * 100:.1f} (at most 10)")
    worst, missing = (0.0, ""), []
    for lab in ITS_ALL:
        p = R["pred"].get(lab)
        if not p:
            missing.append(lab)
            continue
        for kind, vals in p.items():
            if kind == "pairs":
                continue
            for v, pred in vals.items():
                arm = PRED_ARM.get((kind, v)) or {"truncQ": f"tq{v}", "window": f"w{v}", "quality": f"q{v}", "minLen": f"ml{v}"}[kind]
                f = R["flow"].get((lab, arm))
                if pred is None or not f:
                    missing.append(f"{lab} {kind} {v}")
                    continue
                if abs(pred - f["pass"]) >= abs(worst[0]):
                    worst = (pred - f["pass"], f"{lab}, {kind} {v}")
    add("P126", None if missing else abs(worst[0]) <= 0.01, f"the largest difference {worst[0] * 100:+.2f} points ({worst[1]})" + (f"; missing: {', '.join(missing[:5])}" if missing else ""))
    add("P127", *R.get("h", (None, "the 16S comparison was not run")))
    for arm in ["sw15"] + TRIMS[1:] + ["eer"] + TQS + ["fwd"]:
        rises = []
        for lab in STUDIES:
            a, b = R["flow"].get((lab, arm)), R["flow"].get((lab, "std"))
            if a and b and b["asv_median"]:
                rises.append((lab, a["asv_median"] / b["asv_median"] - 1))
        add(f"P128 {arm}", None if len(rises) < 2 else all(r <= FLAG for _, r in rises), "a flag: the median ASVs per sample rise by at most 25%: " + ", ".join(f"{lab} {r * 100:+.0f}%" for lab, r in rises))
    b = best_trim(R)
    add("P130", None if b["best"] is None else True, f"the best trim: {b['best']} (tied within 0.5 points: {', '.join(b.get('tied', []))})")
    return V


def compare_16s(out: Path, cli: Path) -> tuple:
    """H: AmpliTrim 0.4 (defaults, and every new setting at its default) against run 24's 0.2.0 on each 16S dataset."""
    rows, ok, n = [], True, 0
    for d in sorted(p for p in Path(out).iterdir() if p.name.startswith("16s_")):
        name = d.name[len("16s_"):]
        old_p, a, b = Path(cli) / name / "report" / "summary.json", d / "a04" / "summary.json", d / "a04_explicit" / "summary.json"
        if not (old_p.is_file() and a.is_file() and b.is_file()):
            rows.append(f"{name}: missing"); ok = False
            continue
        n += 1
        o, x, y = (json.loads(p.read_text(encoding="utf-8"))["recommended"] for p in (old_p, a, b))
        keys = ("truncLenF", "truncLenR", "maxEE_f", "maxEE_r")
        same = all(o[k] == x[k] for k in keys) and abs(o["estimated_yield"] - x["estimated_yield"]) <= 1e-4 and x == y
        ok = ok and same
        rows.append(f"{name}: {'same' if same else 'DIFFERENT'} ({x['truncLenF']},{x['truncLenR']}, maxEE {x['maxEE_f']},{x['maxEE_r']})")
    return (ok and n > 0) if n else None, f"{n} 16S datasets: " + "; ".join(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True)
    ap.add_argument("--bakker-references", type=Path, required=True)
    ap.add_argument("--pauvert-references", type=Path, required=True)
    ap.add_argument("--previous", type=Path, required=True, help="run 29's results, for the consistency gate")
    ap.add_argument("--cli-16s", type=Path, required=True, help="run 24's command-line reports, for H")
    a = ap.parse_args(argv)
    refs = {m: N.read_references(a.references / f"references_{m}.fasta") for m in ("ITS1", "ITS2") if (a.references / f"references_{m}.fasta").is_file()}
    bakker = {lab: N.read_references(a.bakker_references / f"references_{lab}.fasta") for lab in BR.BAKKER if (a.bakker_references / f"references_{lab}.fasta").is_file()}
    R = analyse(a.out, refs, bakker, a.pauvert_references)
    R["h"] = compare_16s(a.out, a.cli_16s)
    c_ok, c_text = consistency_gate(a.out, a.previous)
    o_ok, shares = overhang_gate(a.out, a.bakker_references, a.pauvert_references)
    t_ok = R.get("pauvert_expected", 0) >= 120
    gates = c_ok and o_ok and t_ok
    V = verdicts(R, gates)
    word = lambda ok: "MET" if ok else "NOT MET" if ok is False else "not decided"
    L = ["# The v0.4 night: summary", "", "The criteria (P117 to P130) are in TEST_REGISTER.md and were fixed before this ran.", "",
         f"**Gate, consistency** with run 29: {'passed' if c_ok else 'FAILED'} ({'; '.join(c_text)}).",
         f"**Gate, overhangs** of at least 3 bases in std: {', '.join(f'{k} {v:.2%}' for k, v in shares.items())} (at most 1%: {'passed' if o_ok else 'FAILED'}).",
         f"**Gate, Pauvert's truth:** {R.get('pauvert_expected', 0)} expected members (at least 120: {'passed' if t_ok else 'FAILED'}).", ""]
    if not gates:
        L += ["A gate failed, so no criterion is decided until the failure is explained.", ""]
    L += [f"- **{c}** {word(ok)}: {t}" for c, ok, t in V]
    b = best_trim(R)
    L += ["", "## The trims (B)", "", "| Trim | Safe on every truth | Mean yield (truth) | ATCC ITS2 yield | Mean bases cut per read | Mean merge share |", "|---|---|---|---|---|---|"]
    L += [f"| {r['arm']}{' (best)' if r['arm'] == b['best'] else ''} | {r['safe']} | {r['mean_yield']:.2%} | {r['atcc_its2']:.2%} | {r['cut']:.1f} | {r['merge']:.2%} |" for r in b["rows"]]
    L += ["", "Next: send this file."]
    (a.out / "V04_SUMMARY.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    cols = ["dataset", "arm", "input", "pass", "merged", "yield", "asv_median", "spurious", "recall", "exact", "found"]
    lines = ["\t".join(cols)]
    for (lab, arm), f in sorted(R["flow"].items()):
        x = R["acc"].get((lab, arm), {})
        lines.append("\t".join([lab, arm] + [f"{f[k]:.6g}" for k in ("input", "pass", "merged", "yield", "asv_median")] + [f"{x[k]:.6g}" if k in x else "" for k in ("spurious", "recall", "exact", "found")]))
    (a.out / "v04_metrics.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0 if gates else 2


if __name__ == "__main__":
    sys.exit(main())
