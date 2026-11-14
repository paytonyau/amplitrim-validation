#!/usr/bin/env python
"""ITS overnight run: the analysis (TEST_REGISTER.md, the section of the ITS overnight run; criteria P83 to P93).

  python validation\\its_night_report.py --out C:\\data\\runs\\its_night --references C:\\data\\references_its

Reads what validation\\its_night.py left in --out (per dataset: dada2\\<arm>\\samples.tsv and seqs.tsv, estimates.tsv, amplitrim\\summary.json) and the reference
sequences of validation\\ncbi_references.py (references_ITS1.fasta, references_ITS2.fasta). Writes ITS_NIGHT_REPORT.md (every table), ITS_NIGHT_SUMMARY.md (one verdict
line per criterion, '- **P83** MET: ...') and its_night_metrics.tsv (one row per dataset and arm).

Truth, as section 33 fixed it before the run:
- The ATCC standard (atcc_its1, atcc_its2) is the primary truth of both markers. Expected: the reference sequences of the members whose reference is of the very strain
  ('strain_genome' or 'strain_record' in the FASTA headers); every copy counts. A sequence is true within 2 edits of an expected sequence, unjudged within 3% of its
  length of a reference of another member (species-level references), spurious otherwise. Unjudged reads are left out of the spurious share.
- The piglet mixes (arfken_mix_its1, arfken_mix_its2) are the second view: the isolate reference is the merged sequences of the isolates in the std arm, after
  chimera removal, that hold at least 1% of an isolate's reads (replicates pooled); true within 2 edits; isolate M is flagged in the report.
- For the forward-read arm (fwd) a sequence is compared with the start of each reference sequence, of the same length.
- Recall: the share of the expected members (the strain-exact members; the isolates of a mix) found in a sample, a member being found when a true sequence of it holds
  at least 0.1% of the sample's reads; the mean over the samples. Reads in exact sequences: the share of judged reads in sequences equal to an expected sequence.
- Yield: the share of the input pairs left after chimera removal; the merged share (before chimera removal) is what the estimates of Q1 are compared with, as in
  Supplementary Figure S2.
"""
import argparse
import csv
import json
import math
import statistics
import sys
from pathlib import Path

ITS_GROUPS = ["atcc_its1", "atcc_its2", "arfken_its1", "arfken_its2", "arfken_mix_its1", "arfken_mix_its2", "its2", "its"]
PRIMARY = {"ITS1": "atcc_its1", "ITS2": "atcc_its2"}
SECOND = {"ITS1": ("arfken_its1", "arfken_mix_its1"), "ITS2": ("arfken_its2", "arfken_mix_its2")}
MIXES = {"MixA": "KCTPM", "MixB": "KC", "MixC": "KT"}
Q2_ARMS = ["ee23", "ee24", "ee25", "ee33", "ee55", "tq10", "tq20", "sw15"]
MAXEE_ARMS = {"std": (2, 2), "ee23": (2, 3), "ee24": (2, 4), "ee25": (2, 5), "ee33": (3, 3), "ee55": (5, 5)}
DIFFS = (10, 5, 2, 0)
EDITS, UNJUDGED_SHARE, REF_SHARE, SEEN_SHARE = 2, 0.03, 0.01, 0.001
P83_POINTS, P84_POINTS = 10.0, 2.0
P85_SPURIOUS, P85_RECALL, P86_GAIN, P87_RISE = 1.0, 2.0, 5.0, 0.25
P90_GAIN, P91_RHO, P92_POINTS, P93_POINTS = 10.0, 0.8, 1.0, 0.5


# ------------------------------------------------------------------------------------------------ reading
def read_samples(path: Path) -> list:
    with open(path, newline="", encoding="utf-8") as fh:
        return [{k: (v if k == "sample" else float(v)) for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t")]


def read_seqs(path: Path) -> dict:
    out = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.setdefault(r["sample"], {})
            out[r["sample"]][r["sequence"]] = out[r["sample"]].get(r["sequence"], 0) + float(r["abundance"])
    return out


def read_references(path: Path) -> list:
    """(member, kind, sequence) from a FASTA of ncbi_references.py; the header is member|species|assembly|kind|copy|length."""
    out, head = [], None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            head = line[1:].split("|")
        elif line.strip() and head:
            out.append((head[0], head[3], line.strip().upper()))
    return out


# ------------------------------------------------------------------------------------------------ edit distance (banded: only small distances matter)
def within(a: str, b: str, k: int) -> int:
    """The edit distance when it is at most k, else k + 1."""
    if abs(len(a) - len(b)) > k:
        return k + 1
    prev = {j: j for j in range(0, min(len(b), k) + 1)}
    for i in range(1, len(a) + 1):
        cur, lo, hi = {}, max(0, i - k), min(len(b), i + k)
        for j in range(lo, hi + 1):
            if j == 0:
                cur[j] = i
                continue
            best = prev.get(j - 1, k + 1) + (a[i - 1] != b[j - 1])
            best = min(best, prev.get(j, k + 1) + 1, cur.get(j - 1, k + 1) + 1)
            cur[j] = min(best, k + 1)
        if min(cur.values()) > k:
            return k + 1
        prev = cur
    return min(prev.get(len(b), k + 1), k + 1)


class Truth:
    """expected: [(member, sequence)]; other: [(member, sequence)] that make a sequence unjudged (within 3% of its length)."""

    def __init__(self, expected: list, other: list = ()):
        self.expected, self.other = list(expected), list(other)
        self.members = sorted({m for m, _ in self.expected})
        self._cache = {}

    def classify(self, q: str, prefix: bool = False):
        """('true', member, exact), ('unjudged', None, False) or ('spurious', None, False)."""
        key = (q, prefix)
        if key not in self._cache:
            cut = (lambda r: r[:len(q)]) if prefix else (lambda r: r)
            best = None
            for m, r in self.expected:
                d = within(q, cut(r), EDITS)
                if d <= EDITS and (best is None or d < best[1]):
                    best = (m, d)
            if best:
                self._cache[key] = ("true", best[0], best[1] == 0)
            else:
                k = max(EDITS, math.ceil(UNJUDGED_SHARE * len(q)))
                near = any(within(q, cut(r), k) <= k for _, r in self.other)
                self._cache[key] = ("unjudged", None, False) if near else ("spurious", None, False)
        return self._cache[key]


def fasta_truth(refs: list) -> Truth:
    exact = [(m, s) for m, k, s in refs if k.startswith("strain")]
    return Truth(exact, [(m, s) for m, k, s in refs if not k.startswith("strain")])


def isolate_truth(iso_seqs: dict) -> Truth:
    """The isolate reference: per isolate letter (the sample name's start, e.g. 'K' of 'K-ITS1a'), the sequences holding at least 1% of its reads, replicates pooled."""
    pooled = {}
    for sample, seqs in iso_seqs.items():
        letter = sample.split("-")[0]
        for q, a in seqs.items():
            pooled.setdefault(letter, {})
            pooled[letter][q] = pooled[letter].get(q, 0) + a
    expected = []
    for letter, seqs in sorted(pooled.items()):
        tot = sum(seqs.values())
        expected += [(letter, q) for q, a in seqs.items() if tot > 0 and a >= REF_SHARE * tot]
    return Truth(expected)


def expected_members(sample: str, truth: Truth, mix: bool) -> list:
    if not mix:
        return truth.members
    name = sample.split("-")[0]
    return [m for m in MIXES.get(name, "") if m in truth.members]


def accuracy(seqs: dict, truth: Truth, mix: bool, prefix: bool = False) -> dict:
    """The means over the samples: spurious share, recall, reads in exact sequences, members found; the share of unjudged reads."""
    rows = []
    for sample, s in sorted(seqs.items()):
        tot = sum(s.values())
        if tot <= 0:
            continue
        sp = un = ex = 0.0
        found = set()
        for q, a in s.items():
            cls, member, exact = truth.classify(q, prefix)
            if cls == "spurious":
                sp += a
            elif cls == "unjudged":
                un += a
            else:
                ex += a if exact else 0
                if a >= SEEN_SHARE * tot:
                    found.add(member)
        judged = tot - un
        exp = expected_members(sample, truth, mix)
        rows.append({"spurious": sp / judged if judged else 0.0, "exact": ex / judged if judged else 0.0, "unjudged": un / tot,
                     "recall": len(found & set(exp)) / len(exp) if exp else float("nan"), "found": len(found & set(exp))})
    if not rows:
        return {"spurious": float("nan"), "recall": float("nan"), "exact": float("nan"), "found": float("nan"), "unjudged": float("nan"), "samples": 0}
    mean = lambda k: statistics.fmean([r[k] for r in rows if not math.isnan(r[k])]) if any(not math.isnan(r[k]) for r in rows) else float("nan")
    return {"spurious": mean("spurious"), "recall": mean("recall"), "exact": mean("exact"), "found": mean("found"), "unjudged": mean("unjudged"), "samples": len(rows)}


def flows(samples: list) -> dict:
    s = lambda k: sum(r[k] for r in samples)
    inp = s("input")
    f = lambda k: s(k) / inp if inp else float("nan")
    return {"input": inp, "pass": f("filtered"), "merged": f("merged"), "yield": f("nonchim"), "asv_median": statistics.median([r["asvs"] for r in samples]) if samples else float("nan")}


def spearman(x: list, y: list) -> float:
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for t in range(i, j + 1):
                r[order[t]] = (i + j) / 2 + 1
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.fmean(rx), statistics.fmean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


# ------------------------------------------------------------------------------------------------ the analysis
def pts(x):
    return f"{x * 100:+.2f}" if not math.isnan(x) else "n/a"


def pct(x):
    return f"{x:.2%}" if not (isinstance(x, float) and math.isnan(x)) else "n/a"


def analyse(out: Path, refs: dict) -> dict:
    out = Path(out)
    R = {"flow": {}, "acc": {}, "est": {}, "amp": {}, "missing": []}
    labels = [p.name for p in sorted(out.iterdir()) if p.is_dir() and (p / "dada2").is_dir()]
    seqs_cache = {}
    for lab in labels:
        for arm_dir in sorted((out / lab / "dada2").iterdir()):
            if (arm_dir / "DONE").is_file():
                R["flow"][(lab, arm_dir.name)] = flows(read_samples(arm_dir / "samples.tsv"))
                seqs_cache[(lab, arm_dir.name)] = arm_dir / "seqs.tsv"
        e = out / lab / "estimates.tsv"
        if e.is_file() and (lab, "std") in R["flow"]:
            rows = read_samples(e)
            inp = R["flow"][(lab, "std")]["input"]
            R["est"][lab] = {k: sum(r[f"merged_d{k}"] for r in rows) / inp if inp else float("nan") for k in DIFFS}
        s = out / lab / "amplitrim" / "summary.json"
        if s.is_file():
            table = json.loads(s.read_text(encoding="utf-8")).get("maxEE_table", [])
            R["amp"][lab] = {tuple(t["maxEE"]): (t["pass_rate"], t["estimated_yield"]) for t in table}
    truths = {}
    for marker, lab in PRIMARY.items():
        if marker in refs:
            truths[lab] = (fasta_truth(refs[marker]), False)
    for marker, (iso, mix) in SECOND.items():
        if (iso, "std") in seqs_cache:
            truths[mix] = (isolate_truth(read_seqs(seqs_cache[(iso, "std")])), True)
    R["isolate_reference"] = {lab: {m: sum(1 for x, _ in t.expected if x == m) for m in t.members} for lab, (t, mixflag) in truths.items() if mixflag}
    for (lab, arm), path in seqs_cache.items():
        if lab in truths:
            t, mixflag = truths[lab]
            R["acc"][(lab, arm)] = accuracy(read_seqs(path), t, mixflag, prefix=(arm == "fwd"))
    return R


def verdicts(R: dict) -> list:
    """[(criterion, met: True/False/None, text)]"""
    V, F, A, E = [], R["flow"], R["acc"], R["est"]
    # Q1
    if "its2" in E and ("its2", "std") in F:
        dada = F[("its2", "std")]["merged"]
        gaps = {k: E["its2"][k] - dada for k in DIFFS}
        best = min((5, 2, 0), key=lambda k: abs(gaps[k]))
        V.append(("P83", abs(gaps[best]) * 100 <= P83_POINTS, "On the second ITS study, at least one stricter estimate is within 10 points of DADA2's merged share: "
                  + ", ".join(f"{k} differences {pts(gaps[k])}" for k in DIFFS) + f" (closest: {best})"))
        s16 = [lab for lab in E if lab.startswith("16s_") and (lab, "std") in F]
        if s16:
            mad = lambda k: statistics.fmean([abs(E[lab][k] - F[(lab, "std")]["merged"]) for lab in s16]) * 100
            V.append(("P84", mad(best) - mad(10) <= P84_POINTS, f"On the {len(s16)} 16S controls, the estimate with {best} differences is not farther from DADA2 than "
                      f"today's (10) by more than 2 points: mean absolute difference {mad(best):.2f} against {mad(10):.2f} points"))
        else:
            V.append(("P84", None, "not decided: no 16S control was run"))
        R["q1_best"] = best
    else:
        V += [("P83", None, "not decided: the second ITS study's std arm or estimates are missing"), ("P84", None, "not decided: P83 could not be computed")]
    # Q2
    safe, useful = {}, {}
    for arm in Q2_ARMS:
        ok, parts = True, []
        for marker, lab in PRIMARY.items():
            if (lab, arm) not in A or (lab, "std") not in A:
                ok = None
                break
            ds = (A[(lab, arm)]["spurious"] - A[(lab, "std")]["spurious"]) * 100
            dr = (A[(lab, "std")]["recall"] - A[(lab, arm)]["recall"]) * 100
            ok = ok and ds <= P85_SPURIOUS and dr <= P85_RECALL
            parts.append(f"{marker} spurious {ds:+.2f} points, recall {-dr:+.1f} points")
        safe[arm] = ok
        V.append((f"P85 {arm}", ok, "safe against std on the ATCC standard: " + ("; ".join(parts) if parts else "not decided: an arm is missing")))
        lab = PRIMARY["ITS2"]
        if (lab, arm) in F and (lab, "std") in F:
            g = (F[(lab, arm)]["yield"] - F[(lab, "std")]["yield"]) * 100
            useful[arm] = g >= P86_GAIN
            V.append((f"P86 {arm}", useful[arm], f"on the ATCC standard's ITS2 the yield rises by at least 5 points: {g:+.1f} points"))
        else:
            useful[arm] = None
            V.append((f"P86 {arm}", None, "not decided: an arm is missing"))
        rises = []
        for real in ("its2", "its"):
            if (real, arm) in F and (real, "std") in F and F[(real, "std")]["asv_median"]:
                rises.append((real, F[(real, arm)]["asv_median"] / F[(real, "std")]["asv_median"] - 1))
        V.append((f"P87 {arm}", all(r <= P87_RISE for _, r in rises) if rises else None,
                  "a flag: the median ASVs per sample rise by at most 25% on the two ITS studies: " + (", ".join(f"{n} {r:+.0%}" for n, r in rises) or "not run")))
    good = [a for a in Q2_ARMS if safe.get(a) and useful.get(a)]
    lab2 = PRIMARY["ITS2"]
    R["advice"] = max(good, key=lambda a: F[(lab2, a)]["yield"]) if good else None
    # Q3 and Q4
    for crit, arm, rule in (("P88", "fwd", "forward reads alone"), ("P89", "mfirst", "merging first")):
        parts, ok = [], True
        for marker, lab in PRIMARY.items():
            if (lab, arm) not in A or (lab, "std") not in A:
                ok = None
                break
            ds = (A[(lab, arm)]["spurious"] - A[(lab, "std")]["spurious"]) * 100
            dr = (A[(lab, arm)]["recall"] - A[(lab, "std")]["recall"]) * 100
            ok = ok and ds <= 1.0 and dr >= -1e-9
            parts.append(f"{marker} spurious {ds:+.2f} points, recall {dr:+.1f} points")
        V.append((crit, ok, f"{rule}: spurious at most 1 point above std and recall not below it, on the ATCC standard: " + ("; ".join(parts) or "not decided")))
    if ("its2", "mfirst") in F and ("its2", "std") in F:
        g = (F[("its2", "mfirst")]["yield"] - F[("its2", "std")]["yield"]) * 100
        V.append(("P90", g >= P90_GAIN, f"merging first raises the yield of the second ITS study by at least 10 points: {g:+.1f} points"))
    else:
        V.append(("P90", None, "not decided: an arm is missing"))
    # P91, P92
    rhos, gaps = [], []
    for lab in ITS_GROUPS:
        amp = R["amp"].get(lab, {})
        arms = [a for a, m in MAXEE_ARMS.items() if m in amp and (lab, a) in F]
        if len(arms) >= 3:
            rhos.append((lab, spearman([amp[MAXEE_ARMS[a]][1] for a in arms], [F[(lab, a)]["merged"] for a in arms]), len(arms)))
        gaps += [(lab, a, (amp[MAXEE_ARMS[a]][0] - F[(lab, a)]["pass"]) * 100) for a in arms]
    V.append(("P91", all(r >= P91_RHO for _, r, _ in rhos) if rhos else None,
              "AmpliTrim's predicted yields order the maxEE arms as DADA2's merged shares do, rank correlation at least 0.8 on each ITS group: "
              + (", ".join(f"{lab} {r:.2f} ({n} arms)" for lab, r, n in rhos) or "not run")))
    worst = max(gaps, key=lambda g: abs(g[2])) if gaps else None
    V.append(("P92", all(abs(g) <= P92_POINTS for _, _, g in gaps) if gaps else None,
              f"AmpliTrim's predicted pass rate is within 1 point of DADA2's at every maxEE arm on every ITS group: {len(gaps)} comparisons"
              + (f", the largest {worst[2]:+.2f} points ({worst[0]}, {worst[1]})" if worst else "")))
    # P93
    parts, ok = [], True
    for marker, lab in PRIMARY.items():
        if not all((lab, a) in A for a in ("std", "mo20", "mo30")):
            ok = None
            break
        d20 = (A[(lab, "std")]["spurious"] - A[(lab, "mo20")]["spurious"]) * 100
        d30 = (A[(lab, "std")]["spurious"] - A[(lab, "mo30")]["spurious"]) * 100
        ok = ok and d20 <= P93_POINTS and d30 <= P93_POINTS
        parts.append(f"{marker} 12 nt against 20 nt {d20:+.2f}, against 30 nt {d30:+.2f} points")
    V.append(("P93", ok, "the spurious share at 12 nt is at most 0.5 point above 20 and 30 nt on the ATCC standard: " + ("; ".join(parts) or "not decided")))
    return V


def line(cid, ok, text):
    return f"- **{cid}** " + ("MET: " if ok is True else "NOT MET: " if ok is False else "") + (text if ok is not None else (text if text.startswith("not decided") else "not decided: " + text))


def report(R: dict, V: list) -> tuple:
    F, A = R["flow"], R["acc"]
    S = ["# ITS overnight run: summary", "", "The criteria (P83 to P93) are in TEST_REGISTER.md and were fixed before this ran.", ""]
    S += [line(*v) for v in V]
    S += ["", f"**The advice for ITS (section 33):** " + (f"`{R['advice']}`, the arm that is safe on both markers and useful on ITS2 with the largest ITS2 yield." if R.get("advice")
                                                       else "keep maxEE 2,2: no arm was both safe on both markers and useful on ITS2."), "", "Next: send this file."]
    L = ["# ITS overnight run: report", "", "## Flows (shares of the input pairs; ASVs: median per sample)", "",
         "| Dataset | Arm | Input pairs | Pass | Merged | After chimera removal | ASVs |", "|---|---|---|---|---|---|---|"]
    L += [f"| {lab} | {arm} | {f['input']:,.0f} | {pct(f['pass'])} | {pct(f['merged'])} | {pct(f['yield'])} | {f['asv_median']:.0f} |" for (lab, arm), f in sorted(F.items())]
    L += ["", "## Accuracy (means over samples; the ATCC standard against NCBI references, the piglet mixes against the isolates' own sequences)", "",
          "| Dataset | Arm | Spurious | Recall | Reads in exact sequences | Members found | Unjudged reads | Samples |", "|---|---|---|---|---|---|---|---|"]
    L += [f"| {lab} | {arm} | {pct(a['spurious'])} | {pct(a['recall'])} | {pct(a['exact'])} | {a['found']:.1f} | {pct(a['unjudged'])} | {a['samples']} |" for (lab, arm), a in sorted(A.items())]
    L += ["", "The isolate reference (sequences per isolate): " + "; ".join(f"{lab}: " + ", ".join(f"{m} {n}" for m, n in v.items()) for lab, v in R["isolate_reference"].items())
          + ". Isolate M is flagged: its sequences were 56 edits from the NCBI record of *L. corymbifera* (section 33).", "",
          "## Q1: the estimates against DADA2's merged share (points)", "", "| Dataset | " + " | ".join(f"{k} differences" for k in DIFFS) + " | DADA2 merged |", "|---|" + "---|" * (len(DIFFS) + 1)]
    for lab, e in sorted(R["est"].items()):
        if (lab, "std") in F:
            d = F[(lab, "std")]["merged"]
            L.append(f"| {lab} | " + " | ".join(pts(e[k] - d) for k in DIFFS) + f" | {pct(d)} |")
    L += ["", "## Verdicts", ""] + [line(*v) for v in V] + ["", "Next: send ITS_NIGHT_SUMMARY.md."]
    M = ["dataset\tarm\tinput\tpass\tmerged\tyield\tasv_median\tspurious\trecall\texact\tfound\tunjudged"]
    for (lab, arm), f in sorted(F.items()):
        a = A.get((lab, arm), {})
        M.append("\t".join([lab, arm] + [f"{f[k]:.6g}" for k in ("input", "pass", "merged", "yield", "asv_median")]
                           + [f"{a[k]:.6g}" if k in a else "" for k in ("spurious", "recall", "exact", "found", "unjudged")]))
    return S, L, M


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--references", type=Path, required=True, help="the folder of references_ITS1.fasta and references_ITS2.fasta (ncbi_references.py)")
    a = ap.parse_args(argv)
    refs = {m: read_references(a.references / f"references_{m}.fasta") for m in ("ITS1", "ITS2") if (a.references / f"references_{m}.fasta").is_file()}
    if not a.out.is_dir() or not refs:
        print(f"Not found: {a.out if not a.out.is_dir() else a.references / 'references_ITS1.fasta'}")
        return 2
    R = analyse(a.out, refs)
    V = verdicts(R)
    S, L, M = report(R, V)
    (a.out / "ITS_NIGHT_SUMMARY.md").write_text("\n".join(S) + "\n", encoding="utf-8")
    (a.out / "ITS_NIGHT_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    (a.out / "its_night_metrics.tsv").write_text("\n".join(M) + "\n", encoding="utf-8")
    print("\n".join(S))
    print(f"written: {a.out / 'ITS_NIGHT_SUMMARY.md'}, ITS_NIGHT_REPORT.md, its_night_metrics.tsv")
    return 0 if all(ok is not None for _, ok, _ in V) else 2


if __name__ == "__main__":
    sys.exit(main())
