#!/usr/bin/env python
"""Over-merging on real reads: how many pairs does vsearch join at the wrong place when the reads are cut to a candidate truncation and merged at the minimum overlap?

  python validation\\overmerge_real.py                                   soil V4 (key v4), the cuts of the study, AmpliTrim's recommendation now
  python validation\\overmerge_real.py --dataset v4 --recommended 124/115 --out C:\\data\\runs\\overmerge_v4
  python validation\\overmerge_real.py --ref-overlap 20                  for reads that overlap less at full length

Real reads have no truth, so the reference is the pair merged from its FULL-length reads, where the overlap is long and the join is reliable. Every pair that overlaps by at
least --ref-overlap (40 nt) at full length is a reference pair, with its reference merged length L. Each cut is then applied (the DADA2-style filter at maxEE 2) and the
surviving reference pairs are merged by vsearch at the minimum overlap of 12 nt. Truncated reads that still span the amplicon must give the same length L; a pair that
comes back at another length was joined at the wrong offset (a false join). Also measured, for each cut: the expected errors that the two reads carry in their overlap, which
is where DADA2 refuses a pair (mismatches). The cuts are the 3 x 3 grid of the study around 120/118 (AmpliTrim's cut in round B6) and around AmpliTrim's recommendation now.
Needs vsearch only (no R). The criteria (P31 to P33) are in TEST_REGISTER.md, section 22, and were fixed before the run. Writes OVERMERGE_REPORT.md and overmerge.tsv.
"""
import argparse
import csv
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT)]

KEY = "v4"
EARLIER_RECOMMENDED, EARLIER_BEST = (120, 118), (130, 108)       # round B6 on soil V4: AmpliTrim's cut of that build, and the best cell of its grid
MIN_OVERLAP, REF_OVERLAP, MAX_EE, STEP = 12, 40, 2.0, 10
WRONG_LIMIT = 0.01                                                # P31: at most this share of the merges is at the wrong length
SAME_POINTS = 1.0                                                 # P32: the wrong-length merges of the two cells differ by at most this many points of the reference pairs


def parse_cell(text: str) -> tuple:
    a, _, b = str(text).partition("/")
    return int(a), int(b)


def grid(center: tuple, step: int = STEP) -> list:
    return [(center[0] + df, center[1] + dr) for df in (-step, 0, step) for dr in (-step, 0, step)]


def cell_text(cell) -> str:
    return f"{cell[0]}/{cell[1]}"


def read_lengths(path) -> dict:
    """{pair number: merged length} from a FASTA whose labels are p<number>."""
    out, label, n = {}, None, 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if label is not None:
                    out[label] = n
                label, n = int(line[1:].split()[0][1:]), 0
            else:
                n += len(line)
    if label is not None:
        out[label] = n
    return out


def write_pairs(pairs, idx, tf, tr, out1: Path, out2: Path) -> None:
    """The pairs `idx`, cut to tf and tr bases (None: the full reads), with the labels p<number>, as plain FASTQ."""
    out1.parent.mkdir(parents=True, exist_ok=True)
    with open(out1, "w", newline="\n") as o1, open(out2, "w", newline="\n") as o2:
        for i in idx:
            i = int(i)
            _, s1, _, q1 = pairs.recs1[i]
            _, s2, _, q2 = pairs.recs2[i]
            a, b = (len(s1) if tf is None else tf), (len(s2) if tr is None else tr)
            o1.write(f"@p{i}\n{s1[:a]}\n+\n{q1[:a]}\n")
            o2.write(f"@p{i}\n{s2[:b]}\n+\n{q2[:b]}\n")


def build_reference(pairs, merge, work: Path, ref_overlap: int = REF_OVERLAP, min_overlap: int = MIN_OVERLAP) -> dict:
    """{pair number: merged length} for the pairs whose full-length reads merge with an overlap of at least `ref_overlap`."""
    write_pairs(pairs, range(pairs.n), None, None, work / "ref_R1.fastq", work / "ref_R2.fastq")
    m = merge(work / "ref_R1.fastq", work / "ref_R2.fastq", work, "ref", min_overlap, None, "tolerant", keep_merged=True)
    lengths = read_lengths(m.merged_path) if getattr(m, "merged_path", None) else {}
    return {i: L for i, L in lengths.items() if int(pairs.len1[i]) + int(pairs.len2[i]) - L >= ref_overlap}


def analyse_cell(pairs, ref: dict, cell: tuple, merge, work: Path, min_overlap: int = MIN_OVERLAP, max_ee: float = MAX_EE) -> dict:
    """What the cut does to the reference pairs: passing the filter, merged at the reference length, merged at another length, and the errors in the overlap."""
    tf, tr = cell
    keep = pairs._masks(tf, tr, max_ee)[4]
    idx = np.flatnonzero(keep)
    n_ref = len(ref)
    passing = [int(i) for i in idx if int(i) in ref]
    row = {"cell": cell, "total": tf + tr, "pass": float(keep.mean()) if pairs.n else 0.0, "reference_pairs": n_ref, "pass_ref": len(passing) / n_ref if n_ref else 0.0}
    d = work / f"cell_{tf}x{tr}"
    for name, diffs in (("tol", None), ("strict", 0)):
        lengths = {}
        if len(idx):
            write_pairs(pairs, idx, tf, tr, d / "R1.fastq", d / "R2.fastq")
            m = merge(d / "R1.fastq", d / "R2.fastq", d, f"{tf}x{tr}", min_overlap, diffs, "tolerant" if diffs is None else "strict", keep_merged=True)
            lengths = read_lengths(m.merged_path) if getattr(m, "merged_path", None) else {}
        correct = sum(1 for i in passing if lengths.get(i) == ref[i])
        wrong = sum(1 for i in passing if i in lengths and lengths[i] != ref[i])
        row[f"correct_{name}"] = correct / n_ref if n_ref else 0.0
        row[f"wrong_{name}"] = wrong / n_ref if n_ref else 0.0
        row[f"wrong_of_merged_{name}"] = wrong / (correct + wrong) if correct + wrong else 0.0
    ees = []
    for i in passing:                                             # the last `overlap` bases of each truncated read lie in the overlap, whatever the split
        ov = tf + tr - ref[i]
        if 1 <= ov <= min(tf, tr):
            ees.append(float((10.0 ** (-pairs._phred1[i, tf - ov:tf].astype(np.float64) / 10.0)).sum() + (10.0 ** (-pairs._phred2[i, tr - ov:tr].astype(np.float64) / 10.0)).sum()))      # the table is uint8: negate a float
    row["overlap_pairs"] = len(ees)
    row["overlap_ee_mean"] = float(np.mean(ees)) if ees else float("nan")
    row["overlap_ee_ge1"] = float(np.mean(np.array(ees) >= 1.0)) if ees else float("nan")
    return row


def verdicts(rows: dict, now) -> list:
    """P31 to P33 from the rows by cell; `now` is AmpliTrim's recommendation now (or None)."""
    a, b = rows.get(EARLIER_RECOMMENDED), rows.get(EARLIER_BEST)
    named = [c for c in (now, EARLIER_RECOMMENDED, EARLIER_BEST) if c is not None and c in rows]
    out = []
    if named and all(rows[c]["reference_pairs"] for c in named):
        worst = max(named, key=lambda c: rows[c]["wrong_of_merged_tol"])
        out.append(("P31", f"At each cut that matters ({', '.join(cell_text(c) for c in named)}) at most {WRONG_LIMIT:.0%} of the merges are at the wrong length: the highest is "
                           f"{rows[worst]['wrong_of_merged_tol']:.2%} at {cell_text(worst)}", all(rows[c]["wrong_of_merged_tol"] <= WRONG_LIMIT for c in named)))
    else:
        out.append(("P31", "not decided: no reference pairs, or a cut was not run", None))
    if a and b and a["reference_pairs"]:
        gap = abs(a["wrong_tol"] - b["wrong_tol"]) * 100
        out.append(("P32", f"Wrong-length merges do not tell {cell_text(EARLIER_RECOMMENDED)} from {cell_text(EARLIER_BEST)} (they differ by at most {SAME_POINTS:g} point of the reference pairs): "
                           f"{a['wrong_tol']:.2%} against {b['wrong_tol']:.2%}, {gap:.2f} points apart", gap <= SAME_POINTS))
        ok = a["overlap_pairs"] and b["overlap_pairs"]
        out.append(("P33", f"(exploratory) The overlap of {cell_text(EARLIER_BEST)} carries fewer expected errors than that of {cell_text(EARLIER_RECOMMENDED)}: "
                           f"{b['overlap_ee_mean']:.3f} against {a['overlap_ee_mean']:.3f} per pair", (b["overlap_ee_mean"] < a["overlap_ee_mean"]) if ok else None))
    else:
        out += [("P32", "not decided: a cut was not run", None), ("P33", "not decided: a cut was not run", None)]
    return out


def report_lines(rows: dict, now, n_pairs: int, n_ref: int, ref_overlap: int) -> list:
    L = ["# Over-merging on real reads", "", "The criteria (P31 to P33) are in TEST_REGISTER.md, section 22, and were fixed before this ran.", "", "## Verdicts", ""]
    for p, txt, ok in verdicts(rows, now):
        L.append(f"- **{p}** {'MET' if ok else 'NOT MET' if ok is not None else 'not decided'}: {txt}")
    L += ["", f"## The cuts ({n_pairs:,} pairs analysed; {n_ref:,} reference pairs, whose full-length reads overlap by at least {ref_overlap} nt)", "",
          "Shares of the reference pairs. 'Correct' is merged at the reference length; 'wrong' is merged at another length (a false join). 'No mismatch' is vsearch with no difference allowed in the "
          "raw overlap. Overlap errors are the expected errors in the overlap of the two reads, per pair.", "",
          "| cut | total | pass the filter | correct | wrong | wrong, of the merged | correct, no mismatch | wrong, no mismatch | overlap errors, mean | overlap errors >= 1 |", "|---|---|---|---|---|---|---|---|---|---|"]
    tags = {EARLIER_RECOMMENDED: " (round B6 recommendation)", EARLIER_BEST: " (best cell of the B6 grid)"}
    for cell, r in sorted(rows.items()):
        label = cell_text(cell) + tags.get(cell, "") + (" (recommended now)" if cell == now else "")
        L.append(f"| {label} | {r['total']} | {r['pass_ref']:.1%} | {r['correct_tol']:.1%} | {r['wrong_tol']:.2%} | {r['wrong_of_merged_tol']:.2%} | {r['correct_strict']:.1%} | {r['wrong_strict']:.2%} | "
                 f"{r['overlap_ee_mean']:.3f} | {r['overlap_ee_ge1']:.1%} |")
    a, b = rows.get(EARLIER_RECOMMENDED), rows.get(EARLIER_BEST)
    if a and b:
        L += ["", f"- What vsearch sees of the gap of round B6: {cell_text(EARLIER_BEST)} merged {b['correct_tol'] * 100:.1f}% of the reference pairs correctly against {a['correct_tol'] * 100:.1f}% for "
                 f"{cell_text(EARLIER_RECOMMENDED)} ({(b['correct_tol'] - a['correct_tol']) * 100:+.1f} points); DADA2 measured 66.4% against 59.6% (+6.8 points)."]
    return L + ["", "Next: send this file."]


def main(argv=None, merge=None, build=None, recommend=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, help="amplitrim_paths.json (default: the one next to the package)")
    ap.add_argument("--dataset", default=KEY, help="the dataset key (default: v4, soil V4)")
    ap.add_argument("--recommended", default="auto", help="AmpliTrim's cut now, as 124/115 (default: run AmpliTrim to find it)")
    ap.add_argument("--ref-overlap", type=int, default=REF_OVERLAP, help="the least overlap of the full-length reads for a reference pair (default 40)")
    ap.add_argument("--out", type=Path, default=Path("overmerge_real_out"))
    a = ap.parse_args(argv)
    from maxee_rules import reads_folder
    cfg_path = a.config or ROOT / "amplitrim_paths.json"
    cfg = json.loads(cfg_path.read_text(encoding="utf-8")) if cfg_path.exists() else {}
    folder = reads_folder(cfg, a.dataset)
    if folder is None or not folder.is_dir():
        print(f"No reads folder for '{a.dataset}' in {cfg_path}.")
        return 2
    if merge is None:
        import shutil
        if shutil.which("vsearch") is None:
            print("vsearch was not found on the PATH.")
            return 2
        from amplitrim.merge_wrapper import run_vsearch as merge
    if build is None:
        from amplitrim.subsample import build_subsample as build
    from amplitrim.trimmer import ReadPairs
    a.out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="overmerge_") as tmp:
        work = Path(tmp)
        now = None
        if a.recommended != "auto":
            now = parse_cell(a.recommended)
        else:
            if recommend is None:
                from amplitrim import pipeline

                def recommend(path):
                    res = pipeline.run(pipeline.Settings(input_dir=path, work_dir=work / "amplitrim", marker="16s", max_ee=2.0))
                    return None if res.recommendation is None else (int(res.recommendation.trim_f), int(res.recommendation.trim_r))
            now = recommend(folder)
        info = build(folder, work / "sub")
        pairs = ReadPairs.load(info.r1, info.r2)
        ref = build_reference(pairs, merge, work, a.ref_overlap)
        if not ref:
            print(f"No pair overlaps by {a.ref_overlap} nt at full length: give --ref-overlap a smaller number (the reads here overlap less).")
            return 2
        cells = sorted(set(grid(EARLIER_RECOMMENDED)) | (set(grid(now)) if now else set()))
        rows = {}
        for cell in cells:
            if cell[0] < 50 or cell[1] < 50 or cell[0] > pairs.width_f or cell[1] > pairs.width_r:
                continue
            rows[cell] = analyse_cell(pairs, ref, cell, merge, work)
            print(f"  {cell_text(cell)}: {rows[cell]['correct_tol']:.1%} correct, {rows[cell]['wrong_tol']:.2%} wrong", flush=True)
        L = report_lines(rows, now, pairs.n, len(ref), a.ref_overlap)
    (a.out / "OVERMERGE_REPORT.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    with open(a.out / "overmerge.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["cell", "total", "pass_ref", "correct_tol", "wrong_tol", "wrong_of_merged_tol", "correct_strict", "wrong_strict", "overlap_ee_mean", "overlap_ee_ge1"])
        for cell, r in sorted(rows.items()):
            w.writerow([cell_text(cell), r["total"], f"{r['pass_ref']:.5f}", f"{r['correct_tol']:.5f}", f"{r['wrong_tol']:.5f}", f"{r['wrong_of_merged_tol']:.5f}", f"{r['correct_strict']:.5f}",
                        f"{r['wrong_strict']:.5f}", f"{r['overlap_ee_mean']:.5f}", f"{r['overlap_ee_ge1']:.5f}"])
    print("\n".join(L[:9]))
    print(f"\nwritten: {a.out / 'OVERMERGE_REPORT.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
