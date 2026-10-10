#!/usr/bin/env python
"""Overnight validation: every test the sequence data on hand can support, in one command, with one summary at the end.

  python validation/overnight.py --out C:\\data\\overnight --mockA-high READY --mockA-all READY --runinfo runs_info.tsv --asv-table asvs.csv
         --reference mockA_reference.fasta [--v1v2 READY] [--its READY] [--figaro-script figaro.py] [--quick]

Run it first with --quick (about 30 minutes): the same steps on small slices, to see that every path and tool works on this machine. Then run it
without --quick (about 8 to 10 hours). A step whose inputs are missing is skipped with a message, a step that fails is recorded and the run
goes on, and a step that finished is not repeated when you start the command again (use --fresh to repeat everything).

Steps, by tier (see --list for the inputs each needs):
  Tier 1  mock_grid         mock A, six high-input samples: truncations (AmpliTrim's, the authors', two longer ones) x maxEE (2,6 and 3,5), scored against
                            the expected sequences: accuracy of AmpliTrim's truncation against the authors', and the coverage rule against truth
          mock_dilutions    the dilution series at 257/260, 265/270 and 279/263
          mock_clip250      mock A reads clipped to 2x250: a tight-overlap case with known truth (AmpliTrim, FIGARO, fixed settings)
          figaro_mockA      AmpliTrim against FIGARO on the six mock A samples, with a paired per-sample interval
          figaro_v1v2       the same on the V1-V2 study (PRJEB42528)
  Tier 2  stability         does the recommendation change with the subsample size?
          per_sample_rank   per-sample rank agreement between the estimate and DADA2 (V1-V2 and the mock A series)
          mock_pooled       DADA2 with pooled inference on mock A
          v1v2_fulldepth    calibration of the yield estimate at full depth on V1-V2 (8 samples)
          asv_classify      what the ASVs outside the expected sequences are (chimeras, contaminants, error variants, unlisted variants), per setting of mock_grid
          coverage_levels   the cost and the benefit of leaving out 1%, 0.5%, 0.2%, 0.1% and 0% of the longest amplicons, on every dataset (the evidence for the default)
          landscape         is AmpliTrim's recommendation near the best of a DADA2 grid of truncations around it?
          extra_datasets    new data given with --extra NAME=READY|REFERENCE (an independent mock, V3-V4, V4-V5): recommendation, coverage levels, comparison with FIGARO
          battery           the standard battery on mock A, V1-V2 and ITS (and, as extra_datasets, on every new dataset): filter replica, yield estimate, against no truncation and FIGARO
  Tier 3  its_chimeras      chimera removal in the ITS study under q2-dada2-like settings (the unexplained 0.1% against 13.9%)
          mock_miss_trace   where the two unexplained mock A misses are lost (reads present? filtered? denoised?)
  Tier 4  round 4, run by `check_build.py --round4` in phases (criteria in TEST_REGISTER.md, section 9; see round4_steps.py):
          split_tiebreak, missing_variants, figaro_usable_pairs, asv_length_check, figaro_rule_check   from results already on disk (minutes)
          richness_sim      does the excess of the yield estimate grow with the richness of a simulated community? (hours)
          overmerge_sim     merges at the wrong length, and the 12 nt minimum overlap, with repeats, chimeras and indels (about an hour)
Not testable with the data on hand (listed in the summary): an independent mock community, V3-V4 2x250 and V4-V5 reads, other platforms.

The margins that decide "met" and "not met" are in MARGINS below; they are printed in the summary. Edit them BEFORE the run, not after seeing results.
"""

import argparse
import contextlib
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

MARGINS = {
    "recovered_variants": 1,         # AmpliTrim may recover this many variants fewer than the authors' truncation
    "recovered_strains": 1,
    "identical_points": 0.10,        # and this many points fewer reads in identical ASVs
    "extra_asvs": 2,                 # and this many more ASVs outside the reference
    "figaro_equivalence_points": 1.0,   # AmpliTrim equals FIGARO at equal maxEE when the 95% interval of the per-sample difference lies within +/- this
    "stability_truncation_nt": 10,   # recommended truncation of either read varies by at most this across subsample sizes
    "stability_yield_points": 1.0,
    "per_sample_rho": 0.80,          # rank correlation of the estimated and the observed per-sample merged share
    "landscape_regret_points": 1.0,  # AmpliTrim's recommendation is within this many points of the best DADA2 setting that meets the coverage floor
}
QUICK_PAIRS = 1500
ITS_DEEP_SAMPLES = "ERR13140519,ERR13140520,ERR13140522,ERR13140523,ERR13140521,ERR13140525,ERR13140527,ERR13140528"
HIGH_INPUT_PATTERN = r"-KitA-\d-\d-(8|9)$"


# ---------------------------------------------------------------------------------------------------------------- plumbing
class Ctx:
    def __init__(self, args):
        self.args = args
        self.out = Path(args.out)
        self.quick = args.quick
        self.cur = self.out
        self.paths = {k: (Path(v) if v else None) for k, v in (
            ("mockA_high", args.mockA_high), ("mockA_all", args.mockA_all), ("runinfo", args.runinfo), ("asv_table", args.asv_table),
            ("reference", args.reference), ("v1v2", args.v1v2), ("its", args.its), ("figaro", args.figaro_script))}
        self.out.mkdir(parents=True, exist_ok=True)
        self._log = open(self.out / "overnight.log", "a", encoding="utf-8")

    def log(self, msg: str) -> None:
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        self._log.write(line + "\n")
        self._log.flush()

    def p(self, name: str):
        return self.paths.get(name)

    def sh(self, cmd: list, name: str) -> str:
        """Run a command, keep its output in the step folder, and raise with the tail of the error if it fails."""
        cmd = [str(c) for c in cmd]
        self.log("  $ " + " ".join(cmd))
        res = subprocess.run(cmd, capture_output=True, text=True)
        (self.cur / f"{name}.log").write_text(res.stdout + "\n--- stderr ---\n" + res.stderr, encoding="utf-8")
        if res.returncode != 0:
            raise RuntimeError(f"{Path(cmd[1]).name} failed (exit {res.returncode}): {(res.stderr or res.stdout).strip()[-500:]}")
        return res.stdout


def verdict(criterion: str, ok, detail: str) -> dict:
    return {"criterion": criterion, "outcome": "info" if ok is None else ("met" if ok else "NOT MET"), "detail": detail}


@contextlib.contextmanager
def dada2_env(pool=None, chimera=None, minfold=None):
    keys = {"AMPLITRIM_DADA2_POOL": pool, "AMPLITRIM_DADA2_CHIMERA": chimera, "AMPLITRIM_DADA2_MINFOLD": minfold}
    old = {k: os.environ.get(k) for k in keys}
    try:
        for k, v in keys.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = str(v)
        yield
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def pick_even(names: list, k: int) -> list:
    """k names spread evenly over the sorted list."""
    names = sorted(names)
    if k >= len(names):
        return names
    return [names[int(i * len(names) / k)] for i in range(k)]


def clip_reads(src: Path, dst: Path, n: int) -> int:
    """Copy every pair keeping only the first n bases of each read (an emulation of shorter reads); returns the number of pairs."""
    from amplitrim.subsample import find_pairs, read_fastq
    dst.mkdir(parents=True, exist_ok=True)
    total = 0
    for sample, r1, r2 in find_pairs(src):
        with gzip.open(dst / f"{sample}_R1.fastq.gz", "wt", newline="\n") as o1, gzip.open(dst / f"{sample}_R2.fastq.gz", "wt", newline="\n") as o2:
            for (h1, s1, p1, q1), (h2, s2, p2, q2) in zip(read_fastq(r1), read_fastq(r2)):
                o1.write(f"{h1}\n{s1[:n]}\n{p1}\n{q1[:n]}\n")
                o2.write(f"{h2}\n{s2[:n]}\n{p2}\n{q2[:n]}\n")
                total += 1
    return total


def longest_read(folder: Path) -> int:
    from amplitrim.subsample import find_pairs, read_fastq
    best = 0
    for _s, r1, r2 in find_pairs(folder):
        for path in (r1, r2):
            for i, (_h, seq, _p, _q) in enumerate(read_fastq(path)):
                best = max(best, len(seq))
                if i >= 500:
                    break
    return best


def table_lines(rows: list, cols: list) -> list:
    """cols = [(header, key or callable, width)]"""
    def cell(r, key):
        return key(r) if callable(key) else r[key]
    out = ["  ".join(f"{h:>{w}s}" for h, _k, w in cols)]
    for r in rows:
        out.append("  ".join(f"{str(cell(r, k)):>{w}s}" for _h, k, w in cols))
    return out


# ---------------------------------------------------------------------------------------------------------------- shared test pieces
def high_input_reference(ctx: Ctx) -> Path:
    """The expected sequences present in the 10^8 and 10^9 cell samples (from the authors' table), written once."""
    out = ctx.out / "_reference" / "mockA_high_input_reference.fasta"
    if out.exists():
        return out
    import pandas as pd
    import mock_check
    from build_mockA_reference import present_in
    counts = pd.read_csv(ctx.p("asv_table"), sep=";", dtype=str).set_index("OTU_ID")
    ref = mock_check.read_fasta(ctx.p("reference"))
    ids = [n.split()[0].split("|")[0] for n, _ in ref]
    keep = set(present_in(counts, ids, HIGH_INPUT_PATTERN))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(f">{n}\n{s}\n" for n, s in ref if n.split()[0].split("|")[0] in keep), encoding="utf-8")
    return out


def run_and_score(ctx: Ctx, reads: Path, tf: int, tr: int, ee: str, ref: Path, tag: str, min_len: int = 0) -> dict:
    """DADA2 at one setting, then the ASVs scored against the expected sequences."""
    import mock_check
    import real_data
    asv = ctx.cur / f"asvs_{tag}.fasta"
    track = ctx.cur / f"track_{tag}.tsv"
    t0 = time.time()
    d = real_data.dada2_run(reads, tf, tr, real_data.ee_text(ee), 12, min_len, track, asv)
    s = mock_check.summarise(asv, ref, ctx.args.near_id, "vsearch")
    c, total = s["categories"], max(s["reads"], 1)
    ref_records = mock_check.read_fasta(ref)
    lengths = {n.split()[0]: len(q) for n, q in ref_records}
    limit = tf + tr - 12
    strains = {mock_check.strain_of(n) for n, _ in ref_records}
    ctx.log(f"  {tag}: pass {d['filtered']:.1%}, merged {d['merged']:.1%}, recovered {len(s['recovered'])}/{len(ref_records)} ({time.time() - t0:.0f} s)")
    return {"tag": tag, "tf": tf, "tr": tr, "ee": ee, "limit": limit, "pass": d["filtered"], "merged": d["merged"], "nonchim": d["nonchim"], "asvs": d["asvs"],
            "variants": len(s["recovered"]), "present": len(ref_records), "strains_found": len({mock_check.strain_of(r) for r in s["recovered"]}),
            "strains_present": len(strains), "identical": c["exact"][1] / total, "near": c["near"][1] / total, "none": c["none"][1] / total,
            "extra": c["near"][0] + c["none"][0], "recovered": s["recovered"], "missing": s["missing"],
            "longer_than_limit": [n for n in lengths if lengths[n] > limit], "lengths": lengths}


def noninferior(a: dict, b: dict):
    d_var, d_str = a["variants"] - b["variants"], a["strains_found"] - b["strains_found"]
    d_id, d_ex = 100 * (a["identical"] - b["identical"]), a["extra"] - b["extra"]
    ok = (d_var >= -MARGINS["recovered_variants"] and d_str >= -MARGINS["recovered_strains"] and d_id >= -MARGINS["identical_points"]
          and d_ex <= MARGINS["extra_asvs"])
    return ok, f"variants {d_var:+d}, strains {d_str:+d}, reads in identical ASVs {d_id:+.2f} points, ASVs outside the reference {d_ex:+d}"


def short(name: str) -> str:
    return name.split("|")[-1].replace("_", " ")


PAIR_LINE = re.compile(r"^\s*(?P<name>.+?): (?P<mean>[+-]?\d+\.\d) points \(95% interval (?P<lo>[+-]?\d+\.\d) to (?P<hi>[+-]?\d+\.\d)\)")


def parse_pairwise(text: str) -> list:
    out = []
    for line in text.splitlines():
        m = PAIR_LINE.match(line)
        if m:
            out.append({"name": m["name"].strip(), "mean": float(m["mean"]), "lo": float(m["lo"]), "hi": float(m["hi"])})
    return out


def read_tsv(path: Path) -> list:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def compare_with_figaro(ctx: Ctx, reads: Path, out: Path, extra: list, label: str, with_numbers: bool = False) -> tuple:
    """compare_methods.py with FIGARO, then the paired per-sample comparison with FIGARO as the baseline. Returns (lines, verdicts), and with_numbers also
    {'like_for_like': {...} or None, 'end_to_end': {...} or None}: AmpliTrim minus FIGARO in points with the 95% interval of the per-sample differences, at FIGARO's
    maxEE (like for like) and with AmpliTrim's own default setting against FIGARO as it ships (end to end)."""
    numbers = {"like_for_like": None, "end_to_end": None}
    def done(lines, verdicts):
        return (lines, verdicts, numbers) if with_numbers else (lines, verdicts)
    cmd = [sys.executable, HERE / "compare_methods.py", reads, "--out", out, "--figaro-script", ctx.p("figaro")] + extra
    ctx.sh(cmd, f"compare_{label}")
    rows = read_tsv(out / "compare_methods.tsv")
    lines = [ln for ln in (out / "compare_methods.txt").read_text(encoding="utf-8").splitlines()]
    figaro = next((r for r in rows if r["name"].lower().startswith("figaro")), None)
    verdicts = []
    if figaro is None:
        return done(lines, [verdict(f"{label}: FIGARO result", None, "no result: FIGARO failed or refused these reads (it is not a difference in outcome); the reason is in the notes of the comparison")])
    pair_text = ctx.sh([sys.executable, HERE / "compare_methods.py", "--pairwise", out, "--baseline", figaro["name"]], f"pairwise_{label}")
    lines += ["", f"Paired per-sample differences against {figaro['name']}:"] + [ln for ln in pair_text.splitlines() if "points" in ln]
    if not parse_pairwise(pair_text):
        return done(lines, [verdict(f"{label}: paired comparison with FIGARO", None, "no paired per-sample comparison was printed (it needs at least 6 samples)")])
    first_amplitrim = next((r["name"] for r in rows if r["name"].startswith("AmpliTrim")), None)            # AmpliTrim's default, with its own maxEE
    for pr in parse_pairwise(pair_text):
        if first_amplitrim is not None and pr["name"] == first_amplitrim:
            numbers["end_to_end"] = {"mean": pr["mean"], "lo": pr["lo"], "hi": pr["hi"]}
        if pr["name"].lower().startswith("amplitrim") and f"maxEE {figaro['maxEE']}" in pr["name"]:
            m = MARGINS["figaro_equivalence_points"]
            numbers["like_for_like"] = {"mean": pr["mean"], "lo": pr["lo"], "hi": pr["hi"]}
            verdicts.append(verdict(f"{label}: AmpliTrim equals FIGARO at FIGARO's maxEE ({figaro['maxEE']})", pr["lo"] >= -m and pr["hi"] <= m,
                                    f"{pr['mean']:+.1f} points (95% interval {pr['lo']:+.1f} to {pr['hi']:+.1f}); margin +/-{m}"))
    if not verdicts:
        verdicts.append(verdict(f"{label}: like-for-like row found", False, f"no AmpliTrim row at FIGARO's maxEE {figaro['maxEE']} in the pairwise output"))
    return done(lines, verdicts)


# ---------------------------------------------------------------------------------------------------------------- Tier 1
def step_mock_grid(ctx: Ctx) -> dict:
    import real_data
    ref = high_input_reference(ctx)
    reads = ctx.p("mockA_high")
    if ctx.quick:
        reads = real_data.thin_reads(reads, ctx.cur / "reads_quick", QUICK_PAIRS)
    truncs = [("AmpliTrim", 257, 260), ("AmpliTrim alt", 268, 249), ("authors", 279, 263), ("longer 524", 257, 267), ("longer 535", 265, 270)]
    ees = ["2,6", "3,5"]
    if ctx.quick:
        truncs, ees = truncs[:3], ees[:1]
    rows = []
    for name, tf, tr in truncs:
        for ee in ees:
            r = run_and_score(ctx, reads, tf, tr, ee, ref, f"{tf}_{tr}_ee{ee.replace(',', '_')}")
            r["name"] = name
            rows.append(r)
    cols = [("setting", lambda r: f"{r['name']} {r['tf']}/{r['tr']}", 22), ("maxEE", "ee", 5), ("limit", "limit", 5), ("pass", lambda r: f"{r['pass']:.1%}", 6),
            ("merged", lambda r: f"{r['merged']:.1%}", 6), ("ASVs", "asvs", 4), ("variants", lambda r: f"{r['variants']}/{r['present']}", 8),
            ("strains", lambda r: f"{r['strains_found']}/{r['strains_present']}", 7), ("identical", lambda r: f"{r['identical']:.2%}", 9), ("outside", "extra", 7)]
    lines = table_lines(rows, cols)
    with open(ctx.cur / "mock_grid.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["setting", "truncF", "truncR", "maxEE", "pass", "merged", "nonchim", "asvs", "variants", "present", "strains", "strains_present", "identical", "near", "none", "outside"])
        for r in rows:
            w.writerow([r["name"], r["tf"], r["tr"], r["ee"], round(r["pass"], 4), round(r["merged"], 4), round(r["nonchim"], 4), r["asvs"], r["variants"], r["present"],
                        r["strains_found"], r["strains_present"], round(r["identical"], 5), round(r["near"], 5), round(r["none"], 5), r["extra"]])
    by = {(r["name"], r["ee"]): r for r in rows}
    verdicts = []
    for ee in ees:
        for mine in ("AmpliTrim", "AmpliTrim alt"):
            if (mine, ee) in by and ("authors", ee) in by:
                ok, detail = noninferior(by[(mine, ee)], by[("authors", ee)])
                verdicts.append(verdict(f"accuracy: {mine} {by[(mine, ee)]['tf']}/{by[(mine, ee)]['tr']} not inferior to the authors' 279/263 at maxEE {ee}", ok, detail))
    broken = [(r["tag"], n) for r in rows for n in r["recovered"] if n in r["longer_than_limit"]]
    verdicts.append(verdict("coverage rule: no variant longer than (total - 12) nt is ever recovered", not broken,
                            "none recovered" if not broken else f"recovered although too long: {broken}"))
    for key, label in (("Enterococcus_faecalis", "E. faecalis (508 nt)"), ("Lactobacillus_reuteri", "L. reuteri (523 nt)")):
        seen = {f"{r['tf']}/{r['tr']}": any(key in n for n in r["recovered"]) for r in rows if r["ee"] == ees[0]}
        verdicts.append(verdict(f"coverage: {label} recovered once the total allows it", None, ", ".join(f"{k}: {'yes' if v else 'no'}" for k, v in seen.items())))
    return {"title": "Mock A, six high-input samples: truncation x maxEE against the expected sequences", "lines": lines, "verdicts": verdicts,
            "files": ["mock_grid.tsv"]}


LONG_ROW = re.compile(r"^\s*1e(?P<d>\d)\s+(?P<n>\d+)\s+(?P<pass>[\d.]+)%\s+(?P<merged>[\d.]+)%\s+(?P<asvs>\d+)\s*\|\s*(?P<found>\d+)/(?P<present>\d+)\s+(?P<sf>\d+)/(?P<sp>\d+)")
MISS_ROW = re.compile(r"^1e(?P<d>\d): present but not recovered: (?P<names>.+)$")


def parse_dilution_text(text: str) -> tuple:
    rows = [m.groupdict() for m in (LONG_ROW.match(ln) for ln in text.splitlines()) if m]
    misses = {m["d"]: [x.strip().rstrip(" ...") for x in m["names"].split(",")] for m in (MISS_ROW.match(ln) for ln in text.splitlines()) if m}
    return rows, misses


def step_mock_dilutions(ctx: Ctx) -> dict:
    import pandas as pd
    truncs = [(257, 260), (265, 270), (279, 263)]
    extra = []
    if ctx.quick:
        truncs, extra = truncs[:1], ["--dilutions", "9", "5"]
    counts = pd.read_csv(ctx.p("asv_table"), sep=";", dtype=str).set_index("OTU_ID")
    lens = {k: len(v) for k, v in counts["Sequence"].items()}
    lines, verdicts = [], []
    for tf, tr in truncs:
        out = ctx.cur / f"t{tf}_{tr}"
        text = ctx.sh([sys.executable, HERE / "mockA_dilutions.py", ctx.p("mockA_all"), "--runinfo", ctx.p("runinfo"), "--asv-table", ctx.p("asv_table"),
                       "--reference", ctx.p("reference"), "--trunc", f"{tf},{tr}", "--max-ee", "2,6", "--out", out] + extra, f"dilutions_{tf}_{tr}")
        rows, misses = parse_dilution_text(text)
        found = sum(int(r["found"]) for r in rows)
        present = sum(int(r["present"]) for r in rows)
        missed = [(d, a.split("|")[0]) for d, names in misses.items() for a in names]
        limit = tf + tr - 12
        long_missed = [x for x in missed if lens.get(x[1], 0) > limit]
        lines += [f"--- {tf}/{tr} (merge limit {limit} nt): {found} of {present} present variants recovered ({found / max(present, 1):.0%}); {len(missed)} misses, "
                  f"{len(long_missed)} of them longer than the limit"]
        lines += [f"    1e{r['d']}: {r['found']}/{r['present']} variants, {r['sf']}/{r['sp']} strains, merged {r['merged']}%" for r in rows]
        verdicts.append(verdict(f"dilutions {tf}/{tr}: misses longer than the merge limit", None,
                                f"{len(long_missed)} of {len(missed)} misses; recovered {found}/{present} ({found / max(present, 1):.0%})"))
    return {"title": "Mock A dilution series at several truncations", "lines": lines, "verdicts": verdicts, "files": []}


def step_mock_clip250(ctx: Ctx) -> dict:
    ref = high_input_reference(ctx)
    clipped = ctx.cur / "reads_2x250"
    n_read = longest_read(ctx.p("mockA_high"))
    if n_read < 251:
        raise RuntimeError(f"the reads are at most {n_read} nt, so there is nothing to clip to 250 nt")
    pairs = clip_reads(ctx.p("mockA_high"), clipped, 250)
    ctx.log(f"  clipped {pairs:,} pairs to 250 nt each")
    extra = ["--amplitrim-max-ee", "2,6", "--setting", "full250=250,250,2,6", "--setting", "t240=240,240,2,6", "--reference", ref, "--pairs", str(QUICK_PAIRS if ctx.quick else 0)]
    out = ctx.cur / "compare"
    verdicts = []
    if ctx.p("figaro"):
        lines, verdicts = compare_with_figaro(ctx, clipped, out, extra, "clip250")
    else:
        ctx.sh([sys.executable, HERE / "compare_methods.py", clipped, "--out", out] + extra, "compare_clip250")
        lines = (out / "compare_methods.txt").read_text(encoding="utf-8").splitlines()
        verdicts.append(verdict("clip250: FIGARO comparison", None, "skipped: no --figaro-script"))
    rows = read_tsv(out / "compare_methods.tsv")
    mine = next((r for r in rows if r["name"].startswith("AmpliTrim")), None)
    full = next((r for r in rows if r["name"] == "full250"), None)
    if mine and full:
        verdicts.append(verdict("clip250: AmpliTrim's recommendation recovers at least as many variants as no truncation (250/250)",
                                int(mine["found"]) >= int(full["found"]) - MARGINS["recovered_variants"], f"AmpliTrim {mine['found']}, 250/250 {full['found']} variants"))
    return {"title": "Mock A clipped to 2x250 (tight overlap, known truth)", "lines": lines, "verdicts": verdicts, "files": ["compare/compare_methods.txt"]}


def step_figaro_mockA(ctx: Ctx) -> dict:
    ref = high_input_reference(ctx)
    extra = ["--amplitrim-max-ee", "2", "2,6", "3,5", "--setting", "authors=279,263,2,6", "--reference", ref, "--pairs", str(QUICK_PAIRS if ctx.quick else 5000)]
    lines, verdicts = compare_with_figaro(ctx, ctx.p("mockA_high"), ctx.cur / "compare", extra, "mockA")
    return {"title": "AmpliTrim against FIGARO on mock A (six high-input samples)", "lines": lines, "verdicts": verdicts, "files": ["compare/compare_methods.txt"]}


def step_figaro_v1v2(ctx: Ctx) -> dict:
    extra = ["--amplitrim-max-ee", "2", "--setting", "more_overlap=230,200,2", "--pairs", str(QUICK_PAIRS if ctx.quick else 5000)]
    lines, verdicts = compare_with_figaro(ctx, ctx.p("v1v2"), ctx.cur / "compare", extra, "v1v2")
    return {"title": "AmpliTrim against FIGARO on the V1-V2 study", "lines": lines, "verdicts": verdicts, "files": ["compare/compare_methods.txt"]}


# ---------------------------------------------------------------------------------------------------------------- Tier 2
def recommended(ctx: Ctx, reads: Path, tag: str, **kw):
    from amplitrim.pipeline import Settings, run
    res = run(Settings(input_dir=reads, work_dir=ctx.cur / f"work_{tag}", **kw))
    return res


def step_stability(ctx: Ctx) -> dict:
    sets = [(k, ctx.p(k)) for k in ("mockA_high", "mockA_all", "v1v2") if ctx.p(k)] + [(n, ready) for n, (ready, _ref) in parse_extra(ctx.args.extra).items()]
    sizes, everys = ([5000, 10000], [1]) if ctx.quick else ([20000, 50000, 100000], [1, 2])
    lines, verdicts, data = [], [], []
    for name, reads in sets:
        rows = []
        for n in sizes:
            for k in everys:
                t0 = time.time()
                res = recommended(ctx, reads, f"{name}_{n}_{k}", max_pairs=n, subsample_every=k)
                rec = res.recommendation
                rows.append({"pairs": n, "every": k, "tf": rec.trim_f if rec else None, "tr": rec.trim_r if rec else None,
                             "yield": rec.yield_est if rec else None, "floor": res.floor, "sec": time.time() - t0})
                ctx.log(f"  {name}: max {n:,} pairs, every {k}: {rows[-1]['tf']}/{rows[-1]['tr']} ({rows[-1]['sec']:.0f} s)")
        ok_rows = [r for r in rows if r["tf"] is not None]
        lines += [f"--- {name}"] + table_lines(rows, [("max pairs", "pairs", 9), ("every", "every", 5), ("truncF", lambda r: r["tf"], 6), ("truncR", lambda r: r["tr"], 6),
                                                       ("yield", lambda r: "-" if r["yield"] is None else f"{r['yield']:.1%}", 7), ("floor", "floor", 6), ("seconds", lambda r: f"{r['sec']:.0f}", 7)])
        if ok_rows:
            sf, sr = max(r["tf"] for r in ok_rows) - min(r["tf"] for r in ok_rows), max(r["tr"] for r in ok_rows) - min(r["tr"] for r in ok_rows)
            sy = 100 * (max(r["yield"] for r in ok_rows) - min(r["yield"] for r in ok_rows))
            data.append({"name": name, "truncation_spread": max(sf, sr), "yield_spread": sy})
            verdicts.append(verdict(f"stability ({name}): recommendation does not depend on the subsample", max(sf, sr) <= MARGINS["stability_truncation_nt"] and sy <= MARGINS["stability_yield_points"],
                                    f"truncF varies by {sf} nt, truncR by {sr} nt, estimated yield by {sy:.2f} points over {len(ok_rows)} runs"))
        else:
            verdicts.append(verdict(f"stability ({name})", False, "AmpliTrim made no recommendation"))
    return {"title": "Stability of the recommendation across subsample sizes", "lines": lines, "verdicts": verdicts, "files": [], "data": data}


def per_sample_rank_one(ctx: Ctx, name: str, reads: Path, tf: int, tr: int, ee: str):
    """Per-sample rank agreement of AmpliTrim's estimated yield with DADA2's merged share, for one dataset at one setting: (lines, verdicts, data or None)."""
    import real_data
    from amplitrim.config import parse_ee
    from compare_methods import sample_merged_share
    from per_sample import per_sample, spearman
    thin = real_data.thin_reads(reads, ctx.cur / f"reads_rank_{name}", QUICK_PAIRS if ctx.quick else 5000)
    track = ctx.cur / f"track_{name}.tsv"
    real_data.dada2_run(thin, tf, tr, real_data.ee_text(ee), 12, 0, track, None)
    observed = sample_merged_share(track)
    est = {r["sample"]: r["estimated_yield"] for r in per_sample(reads, tf, tr, 1000 if ctx.quick else 4000, 12, parse_ee(ee))}
    common = sorted(set(observed) & set(est))
    if len(common) < 10:
        return [], [verdict(f"per-sample rank ({name})", None, f"only {len(common)} samples in common; at least 10 are needed")], None
    a, b = [est[x] for x in common], [observed[x] for x in common]
    rho = spearman(a, b)
    offs = [100 * (x - y) for x, y in zip(a, b)]
    lines = [f"--- {name} at {tf}/{tr}, maxEE {ee}: {len(common)} samples; offset (estimate - DADA2) mean {sum(offs) / len(offs):+.1f} points, range {min(offs):+.1f} to {max(offs):+.1f}; rank correlation {rho:+.2f}"]
    return (lines, [verdict(f"per-sample rank agreement ({name})", rho >= MARGINS["per_sample_rho"], f"Spearman {rho:+.2f} over {len(common)} samples (margin {MARGINS['per_sample_rho']})")],
            {"name": name, "rho": rho, "samples": len(common), "offset_mean": sum(offs) / len(offs)})


def step_per_sample_rank(ctx: Ctx) -> dict:
    from amplitrim.config import parse_ee
    jobs = []                                          # each dataset at the truncation AmpliTrim itself recommends for its maxEE
    lines, verdicts, data = [], [], []
    for name, ee in (("mockA_all", "2,6"), ("v1v2", "2")):
        if ctx.p(name):
            res = guarded(name, "per-sample rank", lambda: recommended(ctx, ctx.p(name), name, max_pairs=20000 if ctx.quick else 100000, max_ee=parse_ee(ee)), lines, verdicts)
            if res is not None and res.recommendation:
                jobs.append((name, ctx.p(name), (res.recommendation.trim_f, res.recommendation.trim_r), ee))
    for name, (ready, _ref) in parse_extra(ctx.args.extra).items():                    # new datasets are tested like the built-in ones
        res = guarded(name, "per-sample rank", lambda: recommended(ctx, ready, name, max_pairs=20000 if ctx.quick else 100000, max_ee=parse_ee("2")), lines, verdicts)
        if res is not None and res.recommendation:
            jobs.append((name, ready, (res.recommendation.trim_f, res.recommendation.trim_r), "2"))
    for name, reads, (tf, tr), ee in jobs:
        out = guarded(name, "per-sample rank", lambda: per_sample_rank_one(ctx, name, reads, tf, tr, ee), lines, verdicts)
        if out is None:
            continue
        l, v, d = out
        lines += l
        verdicts += v
        if d:
            data.append(d)
    return {"title": "Per-sample rank agreement of the yield estimate with DADA2", "lines": lines, "verdicts": verdicts, "files": [], "data": data}


def step_mock_pooled(ctx: Ctx) -> dict:
    import real_data
    ref = high_input_reference(ctx)
    reads = ctx.p("mockA_high")
    if ctx.quick:
        reads = real_data.thin_reads(reads, ctx.cur / "reads_quick", QUICK_PAIRS)
    modes = [("FALSE", "per sample (default)"), ("pseudo", "pseudo-pooled")] + ([] if ctx.quick else [("TRUE", "fully pooled")])
    rows = []
    for mode, label in modes:
        with dada2_env(pool=mode):
            r = run_and_score(ctx, reads, 257, 260, "2,6", ref, f"pool_{mode}")
        r["name"] = label
        rows.append(r)
    lines = table_lines(rows, [("inference", "name", 22), ("merged", lambda r: f"{r['merged']:.1%}", 6), ("ASVs", "asvs", 5), ("variants", lambda r: f"{r['variants']}/{r['present']}", 8),
                               ("strains", lambda r: f"{r['strains_found']}/{r['strains_present']}", 7), ("identical", lambda r: f"{r['identical']:.2%}", 9), ("outside", "extra", 7)])
    base = rows[0]
    verdicts = [verdict(f"pooled inference ({r['name']}) does not lose variants against per-sample inference", r["variants"] >= base["variants"] - MARGINS["recovered_variants"],
                        f"variants {r['variants'] - base['variants']:+d}, strains {r['strains_found'] - base['strains_found']:+d}, outside the reference {r['extra'] - base['extra']:+d}")
                for r in rows[1:]]
    return {"title": "Mock A with pooled DADA2 inference", "lines": lines, "verdicts": verdicts, "files": []}


def step_v1v2_fulldepth(ctx: Ctx) -> dict:
    from amplitrim.subsample import find_pairs
    from compare_methods import subset_reads
    names = pick_even([s for s, _a, _b in find_pairs(ctx.p("v1v2"))], 3 if ctx.quick else 8)
    sub = subset_reads(ctx.p("v1v2"), ctx.cur / "reads", ",".join(names))
    extra = ["--amplitrim-max-ee", "2", "3", "--setting", "more_overlap=230,200,2", "--pairs", str(QUICK_PAIRS if ctx.quick else 0)]
    out = ctx.cur / "compare"
    if ctx.p("figaro"):
        lines, verdicts = compare_with_figaro(ctx, sub, out, extra, "v1v2_full")
    else:
        ctx.sh([sys.executable, HERE / "compare_methods.py", sub, "--out", out] + extra, "compare_v1v2_full")
        lines, verdicts = (out / "compare_methods.txt").read_text(encoding="utf-8").splitlines(), []
    rows = read_tsv(out / "compare_methods.tsv")
    ratios = []
    for r in rows:
        try:
            ratios.append(float(r["merged"]) / float(r["est_yield"]))
        except (KeyError, ValueError, ZeroDivisionError):
            pass
    if ratios:
        verdicts.append(verdict("calibration at full depth: merged / estimated yield (V1-V2, 8 samples)", None,
                                f"mean {sum(ratios) / len(ratios):.2f}, range {min(ratios):.2f} to {max(ratios):.2f} over {len(ratios)} settings (V1-V2 at 5,000 pairs: 0.92; ITS: 0.65 thinned, 0.88 full)"))
    return {"title": "V1-V2 at full depth (8 samples): calibration of the yield estimate", "lines": lines, "verdicts": verdicts, "files": ["compare/compare_methods.txt"]}


# ---------------------------------------------------------------------------------------------------------------- Tier 3
def step_its_chimeras(ctx: Ctx) -> dict:
    import real_data
    from amplitrim.subsample import find_pairs
    from compare_methods import subset_reads
    have = [s for s, _a, _b in find_pairs(ctx.p("its"))]
    wanted = [n for n in ctx.args.its_deep_samples.split(",") if n]
    names = wanted if all(n in have for n in wanted) else pick_even(have, 8)
    if names is not wanted:
        ctx.log(f"  the named ITS samples are not all in the folder; using {len(names)} evenly spread samples instead")
    if ctx.quick:
        names = names[:3]
    sub = subset_reads(ctx.p("its"), ctx.cur / "reads", ",".join(names))
    if ctx.quick:
        sub = real_data.thin_reads(sub, ctx.cur / "reads_quick", QUICK_PAIRS)
    settings = [("DADA2 R defaults (as used so far)", dict()), ("consensus, minFold 1.0 (q2-dada2 default)", dict(minfold=1.0)), ("pooled, minFold 1.0", dict(chimera="pooled", minfold=1.0))]
    rows = []
    for label, kw in settings:
        with dada2_env(**kw):
            d = real_data.dada2_run(sub, 0, 0, "2", 12, 50, None, ctx.cur / f"asvs_{len(rows)}.fasta")
        removed = (d["merged"] - d["nonchim"]) / d["merged"] if d["merged"] else float("nan")
        rows.append({"label": label, "merged": d["merged"], "nonchim": d["nonchim"], "asvs": d["asvs"], "removed": removed})
        ctx.log(f"  {label}: chimeras removed {removed:.1%}")
    lines = table_lines(rows, [("chimera removal", "label", 42), ("merged", lambda r: f"{r['merged']:.1%}", 7), ("non-chimeric", lambda r: f"{r['nonchim']:.1%}", 12),
                               ("removed", lambda r: f"{r['removed']:.1%}", 8), ("ASVs", "asvs", 6)])
    lines.append("The original authors' pipeline removed 13.9% of merged reads on these samples.")
    best = max(rows, key=lambda r: r["removed"])
    return {"title": "ITS: chimera removal under different settings", "lines": lines,
            "verdicts": [verdict("ITS chimera discrepancy explained by the chimera settings", None,
                                 f"largest removal {best['removed']:.1%} with '{best['label']}' against 13.9% in the authors' pipeline")], "files": []}


def step_mock_miss_trace(ctx: Ctx) -> dict:
    import mock_check
    import pandas as pd
    import real_data
    from compare_methods import subset_reads
    from per_sample import load_design
    from amplitrim.subsample import find_pairs
    from mockA_dilutions import group_samples, restricted_reference
    from rapidfuzz.distance import Levenshtein
    counts = pd.read_csv(ctx.p("asv_table"), sep=";", dtype=str).set_index("OTU_ID")
    reference = mock_check.read_fasta(ctx.p("reference"))
    available = {s for s, _a, _b in find_pairs(ctx.p("mockA_all"))}
    groups = group_samples(load_design(ctx.p("runinfo")), available)
    lines = []
    n_traced = 0
    for d in (["9"] if ctx.quick else ["9", "6"]):
        if d not in groups:
            continue
        runs = groups[d]
        ref = restricted_reference(reference, counts, d)
        ref_path = ctx.cur / f"reference_1e{d}.fasta"
        ref_path.write_text("".join(f">{n}\n{s}\n" for n, s in ref), encoding="utf-8")
        sub = subset_reads(ctx.p("mockA_all"), ctx.cur / f"reads_1e{d}", ",".join(runs))
        asv = ctx.cur / f"asvs_1e{d}.fasta"
        real_data.dada2_run(sub, 257, 260, "2,6", 12, 0, None, asv)
        s = mock_check.summarise(asv, ref_path, ctx.args.near_id, "vsearch")
        asvs = mock_check.read_asvs(asv)
        merged_seqs = {}
        for run in runs:                                          # the merged full-length reads, before any filtering or denoising
            r1, r2 = next((a, b) for s_, a, b in find_pairs(ctx.p("mockA_all")) if s_ == run)
            fa = ctx.cur / f"merged_{run}.fasta"
            subprocess.run(["vsearch", "--fastq_mergepairs", str(r1), "--reverse", str(r2), "--fastq_minovlen", "12", "--fastq_allowmergestagger", "--fastaout", str(fa), "--quiet"],
                           check=True, capture_output=True)
            seq = []
            for line in fa.read_text().splitlines():
                if line.startswith(">"):
                    if seq:
                        merged_seqs["".join(seq)] = merged_seqs.get("".join(seq), 0) + 1
                    seq = []
                else:
                    seq.append(line.strip())
            if seq:
                merged_seqs["".join(seq)] = merged_seqs.get("".join(seq), 0) + 1
            fa.unlink()
        refseq = dict((n.split()[0], q) for n, q in ref)
        lines.append(f"--- 10^{d} cells ({len(runs)} samples): {len(s['missing'])} present variant(s) not recovered exactly")
        for name in s["missing"]:
            target = refseq[name]
            exact = merged_seqs.get(target, 0)
            near1 = sum(c for q, c in merged_seqs.items() if abs(len(q) - len(target)) <= 1 and Levenshtein.distance(q, target, score_cutoff=1) <= 1)
            near2 = sum(c for q, c in merged_seqs.items() if abs(len(q) - len(target)) <= 2 and Levenshtein.distance(q, target, score_cutoff=2) <= 2)
            nearest = min(((Levenshtein.distance(seq, target), name_, int((re.search(r"size=(\d+)", name_) or [0, 0])[1])) for name_, seq in asvs), default=None)
            lines.append(f"    {short(name)} ({len(target)} nt): reads identical to it before filtering: {exact}; within 1 edit: {near1}; within 2 edits: {near2}; "
                         + (f"nearest ASV: {nearest[0]} edits, {nearest[2]} reads" if nearest else "no ASVs"))
            n_traced += 1
    return {"title": "Where the unexplained mock A misses are lost", "lines": lines, "verdicts": [verdict("misses traced", None, f"{n_traced} missed variants traced; read counts before filtering are listed above")], "files": []}



# ---------------------------------------------------------------------------------------------------------------- the default coverage setting, the landscape, ASV classes
LEVELS = [0.01, 0.005, 0.002, 0.001, 0.0]


def dada_on(ctx: Ctx, reads: Path, tf: int, tr: int, ee: str, ref, tag: str) -> dict:
    """DADA2 at one setting; scored against the expected sequences when there is a reference."""
    import real_data
    if ref is not None:
        return run_and_score(ctx, reads, tf, tr, ee, ref, tag)
    asv = ctx.cur / f"asvs_{tag}.fasta"
    d = real_data.dada2_run(reads, tf, tr, real_data.ee_text(ee), 12, 0, ctx.cur / f"track_{tag}.tsv", asv)
    ctx.log(f"  {tag}: pass {d['filtered']:.1%}, merged {d['merged']:.1%}, {d['asvs']} ASVs")
    return {"tag": tag, "tf": tf, "tr": tr, "ee": ee, "limit": tf + tr - 12, "pass": d["filtered"], "merged": d["merged"], "nonchim": d["nonchim"], "asvs": d["asvs"],
            "variants": None, "present": None, "strains_found": None, "strains_present": None, "identical": None, "extra": None}


def _thin(ctx: Ctx, reads: Path, name: str, pairs):
    import real_data
    return real_data.thin_reads(reads, ctx.cur / f"reads_{name}", pairs) if pairs else reads


def coverage_dataset(ctx: Ctx, name: str, reads: Path, ee: str, ref=None, thin=5000, levels=None) -> dict:
    """AmpliTrim's recommendation at each coverage level (the share of the longest amplicons that may be left out), then DADA2 at each distinct truncation."""
    from amplitrim.config import parse_ee
    from amplitrim.ranking import frac_unmergeable
    from classify_asvs import read_sized
    levels = levels or ([0.01, 0.001] if ctx.quick else LEVELS)
    dreads = _thin(ctx, reads, name, thin)
    cache, rows = {}, []
    for level in levels:
        res = recommended(ctx, reads, f"{name}_cov{level}", max_unmergeable=level, max_ee=parse_ee(ee), **({"max_pairs": 20000} if ctx.quick else {}))
        rec = res.recommendation
        if rec is None:
            rows.append({"level": level, "tf": None})
            continue
        key = (rec.trim_f, rec.trim_r)
        if key not in cache:
            cache[key] = dada_on(ctx, dreads, rec.trim_f, rec.trim_r, ee, ref, f"{name}_{rec.trim_f}_{rec.trim_r}")
        limit = rec.trim_f + rec.trim_r - 12
        rows.append({"level": level, "tf": rec.trim_f, "tr": rec.trim_r, "limit": limit, "left": frac_unmergeable(res.inserts, limit), "est": rec.yield_est, "d": cache[key]})
    done = [r for r in rows if r["tf"] is not None]
    census = None
    if done:                                                        # what the most inclusive truncation finds beyond each limit
        widest = max(done, key=lambda r: r["limit"])
        asv = ctx.cur / f"asvs_{widest['d']['tag']}.fasta"
        census = read_sized(asv) if asv.exists() else None
    for r in done:
        if census:
            beyond = [a for a in census if len(a["seq"]) > r["limit"]]
            tot = sum(a["size"] for a in census) or 1
            r["beyond"] = (len(beyond), sum(a["size"] for a in beyond) / tot)
        else:
            r["beyond"] = None
    cols = [("level", lambda r: f"{r['level']:.1%}", 6), ("truncation", lambda r: "-" if r["tf"] is None else f"{r['tf']}/{r['tr']}", 10), ("limit", lambda r: r.get("limit", "-"), 5),
            ("left out", lambda r: "-" if r["tf"] is None else f"{r['left']:.2%}", 8), ("estimate", lambda r: "-" if r["tf"] is None else f"{r['est']:.1%}", 8),
            ("merged", lambda r: "-" if r["tf"] is None else f"{r['d']['merged']:.1%}", 7)]
    if ref is not None:
        cols += [("variants", lambda r: "-" if r["tf"] is None else f"{r['d']['variants']}/{r['d']['present']}", 8), ("strains", lambda r: "-" if r["tf"] is None else f"{r['d']['strains_found']}/{r['d']['strains_present']}", 7),
                 ("outside", lambda r: "-" if r["tf"] is None else r["d"]["extra"], 7)]
    cols.append(("ASVs beyond the limit", lambda r: "-" if r.get("beyond") is None else f"{r['beyond'][0]} ({r['beyond'][1]:.2%} of reads)", 24))
    return {"name": name, "rows": rows, "lines": [f"--- {name} (maxEE {ee})"] + table_lines(rows, cols), "data": coverage_data(name, ref, rows)}


def coverage_data(name: str, ref, rows: list) -> dict:
    """The numbers of a coverage-level run in a form that can be saved and read again (check_build.py --decide reads them)."""
    out = []
    for r in rows:
        if r["tf"] is None:
            continue
        d = r["d"]
        out.append({"level": r["level"], "tf": r["tf"], "tr": r["tr"], "limit": r["limit"], "left": r["left"], "est": r["est"], "merged": d["merged"],
                    "variants": d.get("variants"), "present": d.get("present"), "strains_found": d.get("strains_found"), "strains_present": d.get("strains_present"),
                    "outside": d.get("extra")})
    return {"name": name, "has_truth": ref is not None, "rows": out}


def cost_benefit(res: dict) -> str:
    rows = {r["level"]: r for r in res["rows"] if r["tf"] is not None}
    a, b = rows.get(0.01), rows.get(0.001)
    if not a or not b:
        return f"{res['name']}: not enough levels to compare"
    text = f"{res['name']}: covering the longest 0.1% instead of 1% changes the truncation from {a['tf']}/{a['tr']} to {b['tf']}/{b['tr']} and the merged share by {100 * (b['d']['merged'] - a['d']['merged']):+.1f} points"
    if a["d"]["variants"] is not None:
        text += f"; variants {b['d']['variants'] - a['d']['variants']:+d}, strains {b['d']['strains_found'] - a['d']['strains_found']:+d}, ASVs outside the reference {b['d']['extra'] - a['d']['extra']:+d}"
    if a.get("beyond"):
        text += f"; the 1% setting leaves {a['beyond'][0]} ASVs ({a['beyond'][1]:.2%} of reads) beyond its limit"
    return text


def _coverage_jobs(ctx: Ctx) -> list:
    jobs = []
    if ctx.p("mockA_high") and ctx.p("reference") and ctx.p("asv_table"):
        jobs.append(("mockA_high", ctx.p("mockA_high"), "2,6", high_input_reference(ctx), 1500 if ctx.quick else None))
    if ctx.p("v1v2"):
        jobs.append(("v1v2", ctx.p("v1v2"), "2", None, 1500 if ctx.quick else 5000))
    for name, (ready, ref) in parse_extra(ctx.args.extra).items():
        jobs.append((name, ready, "2", ref, 1500 if ctx.quick else 5000))
    return jobs


def parse_extra(items) -> dict:
    """{name: (reads folder, reference or None)} from --extra NAME=READY|REFERENCE."""
    out = {}
    for item in items or []:
        name, _, rest = item.partition("=")
        ready, _, ref = rest.partition("|")
        if not name or not ready:
            raise SystemExit(f"--extra {item!r}: use NAME=READY or NAME=READY|REFERENCE")
        out[name] = (Path(ready), Path(ref) if ref else None)
    return out


def guarded(name: str, what: str, fn, lines: list, verdicts: list):
    """fn() for one dataset. A failure is written into the report as a criterion that was not met and the other datasets go on (a failing dataset used to lose the whole step)."""
    try:
        return fn()
    except KeyboardInterrupt:
        raise
    except Exception as err:
        detail = f"{type(err).__name__}: {str(err)[:200]}"
        lines += [f"--- {name}: {what} FAILED: {detail}", ""]
        verdicts.append(verdict(f"{what} ({name}) ran", False, detail))
        return None


def step_coverage_levels(ctx: Ctx) -> dict:
    jobs = _coverage_jobs(ctx)
    lines, verdicts, data = [], [], []
    for name, reads, ee, ref, thin in jobs:
        res = guarded(name, "coverage levels", lambda: coverage_dataset(ctx, name, reads, ee, ref, thin), lines, verdicts)
        if res is None:
            continue
        lines += res["lines"] + [""]
        data.append(res["data"])
        verdicts.append(verdict(f"coverage levels: {name}", None, cost_benefit(res)))
    return {"title": "What leaving out 1%, 0.5%, 0.2%, 0.1% and 0% of the longest amplicons costs and buys", "lines": lines, "verdicts": verdicts, "files": [], "data": data}


def landscape_dataset(ctx: Ctx, name: str, reads: Path, ee: str, ref, span: int, step: int, thin) -> dict:
    from amplitrim.config import parse_ee
    res = recommended(ctx, reads, f"{name}_land", max_ee=parse_ee(ee), **({"max_pairs": 20000} if ctx.quick else {}))
    rec = res.recommendation
    if rec is None:
        raise RuntimeError(f"{name}: AmpliTrim made no recommendation")
    f0, r0, floor = rec.trim_f, rec.trim_r, res.floor
    dreads = _thin(ctx, reads, name, thin)
    offs = list(range(-span, span + 1, step))
    rows = []
    for df in offs:
        for dr in offs:
            tf, tr = f0 + df, r0 + dr
            if tf < 100 or tr < 100 or tf > res.sub.read_len_f or tr > res.sub.read_len_r:
                continue
            r = dada_on(ctx, dreads, tf, tr, ee, ref, f"{name}_{tf}_{tr}")
            r.update({"df": df, "dr": dr, "meets": tf + tr >= floor})
            rows.append(r)
    mine = next((r for r in rows if r["df"] == 0 and r["dr"] == 0), None)
    feasible = [r for r in rows if r["meets"]]
    if mine is None or not feasible:
        raise RuntimeError(f"{name}: no grid point to compare")
    best_all, best_ok = max(rows, key=lambda r: r["merged"]), max(feasible, key=lambda r: r["merged"])
    tfs, trs = sorted({r["tf"] for r in rows}), sorted({r["tr"] for r in rows})
    by = {(r["tf"], r["tr"]): r for r in rows}
    lines = [f"--- {name} (maxEE {ee}); AmpliTrim recommends {f0}/{r0} (coverage floor {floor}); share of pairs merged; * = AmpliTrim, . = total below the floor",
             "  trF\\trR " + " ".join(f"{t:>7d}" for t in trs)]
    for a in tfs:
        cells = []
        for b in trs:
            r = by.get((a, b))
            cells.append("      -" if r is None else f"{100 * r['merged']:6.1f}{'*' if (a, b) == (f0, r0) else ('.' if not r['meets'] else ' ')}")
        lines.append(f"  {a:>8d} " + " ".join(cells))
    regret = 100 * (best_ok["merged"] - mine["merged"])
    lines.append(f"  best of the grid: {best_all['tf']}/{best_all['tr']} at {best_all['merged']:.1%}; best among settings that meet the floor: {best_ok['tf']}/{best_ok['tr']} at {best_ok['merged']:.1%}; "
                 f"AmpliTrim {mine['merged']:.1%}, so {regret:.1f} points below the best that covers the amplicons")
    if ref is not None:
        lines.append(f"  recovered: AmpliTrim {mine['variants']}/{mine['present']} variants, {mine['strains_found']}/{mine['strains_present']} strains; "
                     f"the best feasible setting {best_ok['variants']}/{best_ok['present']}, {best_ok['strains_found']}/{best_ok['strains_present']}")
    return {"name": name, "lines": lines, "regret": regret, "n": len(rows)}


def step_landscape(ctx: Ctx) -> dict:
    jobs = []
    if ctx.p("mockA_high") and ctx.p("reference") and ctx.p("asv_table"):
        jobs.append(("mockA_high", ctx.p("mockA_high"), "2,6", high_input_reference(ctx), 10 if ctx.quick else 20, 1500 if ctx.quick else 10000))
    if ctx.p("v1v2"):
        jobs.append(("v1v2", ctx.p("v1v2"), "2", None, 10, 1500 if ctx.quick else 3000))
    for name, (ready, ref) in parse_extra(ctx.args.extra).items():
        jobs.append((name, ready, "2", ref, 10, 1500 if ctx.quick else 3000))
    lines, verdicts, data = [], [], []
    for name, reads, ee, ref, span, thin in jobs:
        res = guarded(name, "landscape", lambda: landscape_dataset(ctx, name, reads, ee, ref, span, 10, thin), lines, verdicts)
        if res is None:
            continue
        data.append({"name": name, "regret": res["regret"], "n": res["n"]})
        lines += res["lines"] + [""]
        verdicts.append(verdict(f"landscape ({name}): AmpliTrim's recommendation is near the best setting that covers the amplicons", res["regret"] <= MARGINS["landscape_regret_points"],
                                f"{res['regret']:.1f} points below the best of {res['n']} DADA2 runs (margin {MARGINS['landscape_regret_points']})"))
    return {"title": "Is the recommendation near the best of a DADA2 grid around it?", "lines": lines, "verdicts": verdicts, "files": [], "data": data}


def step_asv_classify(ctx: Ctx) -> dict:
    import classify_asvs as ca
    grid = ctx.out / "mock_grid"
    fastas = sorted(grid.glob("asvs_*.fasta")) if grid.exists() else []
    if not fastas:
        raise RuntimeError("no ASV files from mock_grid: run that step first")
    present = high_input_reference(ctx)
    rows_out, tsvs = [], []
    for fa in fastas:
        tag = fa.stem.replace("asvs_", "")
        classes = ca.classify(fa, ctx.p("reference"), present, ctx.p("asv_table"))
        s = ca.summarise(classes)
        with open(ctx.cur / f"classes_{tag}.tsv", "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["asv", "reads", "length", "class", "detail"])
            for r in classes:
                w.writerow([r["id"], r["size"], r["length"], r["class"], r["detail"]])
        m = re.match(r"(\d+)_(\d+)_ee(.+)", tag)
        label = f"{m[1]}/{m[2]} maxEE {m[3].replace('_', ',')}" if m else tag
        rows_out.append({"label": label, "s": s})
        tsvs.append(f"classes_{tag}.tsv")
    show = [c for c in ca.CLASSES if c not in ("expected",)]
    cols = [("setting", "label", 22), ("outside", lambda r: sum(r["s"]["classes"][c][0] for c in ca.CLASSES if c not in ("expected", "listed elsewhere")), 7)]
    cols += [(c, (lambda c: (lambda r: f"{r['s']['classes'][c][0]} ({r['s']['classes'][c][1] / r['s']['reads']:.2%})"))(c), 17) for c in show]
    lines = ["ASVs per class (share of all ASV reads in brackets); outside = every class except expected and listed elsewhere"] + table_lines(rows_out, cols)
    verdicts = []
    mine = next((r for r in rows_out if r["label"] == "257/260 maxEE 2,6"), None)
    theirs = next((r for r in rows_out if r["label"] == "279/263 maxEE 2,6"), None)
    if mine and theirs:
        diff = {c: mine["s"]["classes"][c][0] - theirs["s"]["classes"][c][0] for c in ca.CLASSES}
        total = sum(v for c, v in diff.items() if c not in ("expected", "listed elsewhere"))
        verdicts.append(verdict("ASVs outside the reference: where AmpliTrim's 257/260 has more than the authors' 279/263", None,
                                f"{total:+d} in all; " + ", ".join(f"{c} {v:+d}" for c, v in diff.items() if v and c not in ("expected",))))
    return {"title": "What the ASVs outside the expected sequences are", "lines": lines, "verdicts": verdicts, "files": tsvs}


REGIONS = {"mockA_high": "V1-V3", "v1v2": "V1-V2", "its": "ITS1", "v4": "V4", "mockB": "V1-V3", "v4v5": "V4-V5", "v3v4": "V3-V4",
           "mcm_v46": "V4-V6", "mcm_v12": "V1-V2", "sputum_v1v2": "V1-V2", "sputum_v3v4": "V3-V4", "sputum_v5v7": "V5-V7", "sputum_v7v9": "V7-V9"}


def _num(row: dict, key: str):
    try:
        return float(row[key])
    except (KeyError, ValueError, TypeError):
        return None


def battery_row(name: str, marker: str, region: str, samples: int, rows: list, text: str, numbers: dict = None) -> dict:
    """The numbers of one dataset from its compare_methods run (AmpliTrim against no truncation, and against FIGARO when it ran), in points."""
    mine = [r for r in rows if r.get("name", "").startswith("AmpliTrim")]
    if not mine:
        raise RuntimeError(f"{name}: no AmpliTrim row in the comparison")
    head = (next((r for r in mine if str(r.get("maxEE")) in ("2", "2.0")), None) or mine[0]) if marker == "its" else mine[0]

    def points(r: dict) -> dict:
        ep, dp, ey, mg = (_num(r, k) for k in ("est_pass", "pass", "est_yield", "merged"))
        return {"filter_diff": None if ep is None or dp is None else 100 * (ep - dp), "yield_offset": None if ey is None or mg is None else 100 * (ey - mg),
                "ratio": None if not ey or mg is None else mg / ey, "merged": mg}
    other = "trunc200" if marker == "its" else "full"                               # the other method of the paired comparison, as it is named in the output
    pair = next((q for q in parse_pairwise(text) if q["name"] == other), None)
    vs_other = None if pair is None else {"label": other, "mean": -pair["mean"], "lo": -pair["hi"], "hi": -pair["lo"]}      # printed as other minus AmpliTrim: AmpliTrim minus other is the opposite
    row = {"name": name, "marker": marker, "region": region, "samples": samples, "setting": "no truncation" if marker == "its" else f"{head['truncF']}/{head['truncR']}",
           "maxEE": head["maxEE"], **points(head), "vs_other": vs_other, "vs_figaro": (numbers or {}).get("like_for_like"), "vs_figaro_end": (numbers or {}).get("end_to_end")}
    truth = []
    for r in rows:                                                                      # present when the comparison was scored against expected sequences
        if str(r.get("found", "")).strip() not in ("", "None"):
            nt = str(r.get("truncF", "")) in ("0", "")
            truth.append({"method": r.get("name"), "setting": "no truncation" if nt else f"{r.get('truncF')}/{r.get('truncR')}", "maxEE": r.get("maxEE"), "merged": _num(r, "merged"),
                          "found": _num(r, "found"), "strains": r.get("strains") or None, "identical": _num(r, "identical"), "extra": _num(r, "extra")})
    if truth:
        row["truth"] = truth
    cal = [r for r in mine if None not in (_num(r, "est_pass"), _num(r, "pass"), _num(r, "est_yield"), _num(r, "merged"))]      # the settings with every number: the calibration across settings
    if len(cal) >= 2:
        offs = [100 * (_num(r, "est_yield") - _num(r, "merged")) for r in cal]
        ratios = [_num(r, "merged") / _num(r, "est_yield") for r in cal if _num(r, "est_yield")]
        rank = None
        if len(cal) >= 3:
            try:
                from per_sample import spearman
                rank = spearman([_num(r, "est_yield") for r in cal], [_num(r, "merged") for r in cal])
            except (ValueError, ZeroDivisionError):
                rank = None
        row["calibration"] = {"n": len(cal), "filter_max": max(abs(100 * (_num(r, "est_pass") - _num(r, "pass"))) for r in cal), "offset_mean": sum(offs) / len(offs), "offset_lo": min(offs), "offset_hi": max(offs),
                              "ratio_mean": sum(ratios) / len(ratios), "ratio_lo": min(ratios), "ratio_hi": max(ratios), "rank": rank}
    if marker == "its":
        row["its_rows"] = [{"maxEE": r["maxEE"], **{k: v for k, v in points(r).items() if k in ("filter_diff", "yield_offset", "merged", "ratio")}} for r in mine]
    return row


def battery_dataset(ctx: Ctx, name: str, marker: str, region: str, reads: Path, ref, thin, ee: str) -> tuple:
    """compare_methods on one dataset (AmpliTrim, no truncation or a 200/200 truncation for ITS, FIGARO when it can run) -> (row, lines, verdicts)."""
    from amplitrim.subsample import find_pairs
    base = ctx.cur
    ctx.cur = base / name
    ctx.cur.mkdir(parents=True, exist_ok=True)
    try:
        samples = len(find_pairs(reads))
        out = ctx.cur / "compare"
        if marker == "its":
            extra = ["--marker", "its", "--amplitrim-max-ee", "1", "2", "5", "--setting", "trunc200=200,200,2", "--pairs", str(thin)]
        else:
            extra = ["--amplitrim-max-ee", ee, "1", "5", "--setting", f"full=0,0,{ee}", "--pairs", str(thin)] + (["--reference", ref] if ref else [])      # three maxEE settings, as for ITS
        verdicts, numbers = [], None
        if ctx.p("figaro") and marker != "its":                                          # FIGARO needs reads of one length and refuses ITS
            lines, verdicts, numbers = compare_with_figaro(ctx, reads, out, extra, "battery_" + name, with_numbers=True)
        else:
            ctx.sh([sys.executable, HERE / "compare_methods.py", reads, "--out", out] + extra, "compare_" + name)
            lines = (out / "compare_methods.txt").read_text(encoding="utf-8").splitlines()
        row = battery_row(name, marker, region, samples, read_tsv(out / "compare_methods.tsv"), "\n".join(lines), numbers)
        return row, [f"--- {name} ({region}, {marker.upper()}, {samples} samples)"] + lines, verdicts
    finally:
        ctx.cur = base


def step_battery(ctx: Ctx) -> dict:
    """The standard battery on the datasets built in (mock A high input, V1-V2, ITS): the test-stage table of the manuscript."""
    thin = 1500 if ctx.quick else 5000
    jobs = []
    if ctx.p("mockA_high") and ctx.p("reference") and ctx.p("asv_table"):
        jobs.append(("mockA_high", "16s", ctx.p("mockA_high"), high_input_reference(ctx), "2,6"))
    if ctx.p("v1v2"):
        jobs.append(("v1v2", "16s", ctx.p("v1v2"), None, "2"))
    if ctx.p("its"):
        jobs.append(("its", "its", ctx.p("its"), None, "2"))
    data, verdicts, lines = [], [], []
    for name, marker, reads, ref, ee in jobs:
        out = guarded(name, "battery", lambda: battery_dataset(ctx, name, marker, REGIONS.get(name, "?"), reads, ref, thin, ee), lines, verdicts)
        if out is None:
            continue
        row, ls, vs = out
        data.append(row)
        lines += ls + [""]
        verdicts += vs
    return {"title": "The standard battery on each dataset (the test-stage table)", "lines": lines, "verdicts": verdicts, "files": [], "data": data}


def step_extra_datasets(ctx: Ctx) -> dict:
    """The standard battery on every dataset given with --extra (an independent mock, V4, V3-V4, V4-V5); their coverage levels, landscape and per-sample rank are in the steps of those names."""
    data, verdicts, lines = [], [], []
    for name, (ready, ref) in parse_extra(ctx.args.extra).items():
        out = guarded(name, "battery", lambda: battery_dataset(ctx, name, "16s", REGIONS.get(name, "?"), ready, ref, 1500 if ctx.quick else 5000, "2"), lines, verdicts)
        if out is None:
            continue
        row, ls, vs = out
        data.append(row)
        lines += ls + [""]
        verdicts += vs
    return {"title": "New datasets: the standard battery (their coverage levels are in coverage_levels)", "lines": lines, "verdicts": verdicts, "files": [], "data": data}


# ---------------------------------------------------------------------------------------------------------------- registry
def _r4():
    """The steps of round 4 live in round4_steps.py (which does not import this file): they are given the helpers they need."""
    import types
    import round4_steps
    return round4_steps, types.SimpleNamespace(verdict=verdict, parse_extra=parse_extra, table_lines=table_lines, recommended=recommended,
                                               high_input_reference=high_input_reference, guarded=guarded)


def step_split_tiebreak(ctx: Ctx) -> dict:
    mod, H = _r4()
    return mod.step_split_tiebreak(ctx, H)


def step_missing_variants(ctx: Ctx) -> dict:
    mod, H = _r4()
    return mod.step_missing_variants(ctx, H)


def step_figaro_usable_pairs(ctx: Ctx) -> dict:
    mod, H = _r4()
    return mod.step_figaro_usable_pairs(ctx, H)


def step_asv_length_check(ctx: Ctx) -> dict:
    mod, H = _r4()
    return mod.step_asv_length_check(ctx, H)


def step_figaro_rule_check(ctx: Ctx) -> dict:
    mod, H = _r4()
    return mod.step_figaro_rule_check(ctx, H)


def step_richness_sim(ctx: Ctx) -> dict:
    mod, H = _r4()
    return mod.step_richness_sim(ctx, H)


def step_overmerge_sim(ctx: Ctx) -> dict:
    mod, H = _r4()
    return mod.step_overmerge_sim(ctx, H)


STEPS = [
    # id, tier, estimated minutes (full, quick), needs (paths), tools, function
    ("mock_grid", 1, 40, 8, ["mockA_high", "asv_table", "reference"], ["R", "vsearch", "pandas"], step_mock_grid),
    ("asv_classify", 1, 5, 2, ["mockA_high", "asv_table", "reference"], ["vsearch", "pandas", "rapidfuzz"], step_asv_classify),
    ("mock_clip250", 1, 30, 6, ["mockA_high", "asv_table", "reference"], ["R", "vsearch", "pandas"], step_mock_clip250),
    ("figaro_mockA", 1, 25, 6, ["mockA_high", "asv_table", "reference", "figaro"], ["R", "vsearch", "pandas"], step_figaro_mockA),
    ("mock_dilutions", 1, 90, 6, ["mockA_all", "runinfo", "asv_table", "reference"], ["R", "vsearch", "pandas"], step_mock_dilutions),
    ("coverage_levels", 1, 90, 12, ["anydata"], ["R", "vsearch"], step_coverage_levels),
    ("landscape", 1, 150, 12, ["anydata"], ["R", "vsearch"], step_landscape),
    ("mock_miss_trace", 3, 20, 5, ["mockA_all", "runinfo", "asv_table", "reference"], ["R", "vsearch", "pandas", "rapidfuzz"], step_mock_miss_trace),
    ("mock_pooled", 2, 35, 5, ["mockA_high", "asv_table", "reference"], ["R", "vsearch", "pandas"], step_mock_pooled),
    ("stability", 2, 60, 6, ["mockA_high"], ["vsearch"], step_stability),
    ("per_sample_rank", 2, 45, 8, ["mockA_all"], ["R", "vsearch"], step_per_sample_rank),
    ("its_chimeras", 3, 60, 8, ["its"], ["R", "vsearch"], step_its_chimeras),
    ("figaro_v1v2", 1, 60, 8, ["v1v2", "figaro"], ["R", "vsearch"], step_figaro_v1v2),
    ("v1v2_fulldepth", 2, 120, 12, ["v1v2"], ["R", "vsearch"], step_v1v2_fulldepth),
    ("battery", 1, 60, 8, ["anydata"], ["R", "vsearch"], step_battery),
    ("extra_datasets", 1, 60, 10, ["extra"], ["R", "vsearch"], step_extra_datasets),
    # round 4 (TEST_REGISTER.md, section 9): run by `check_build.py --round4`, in phases
    ("split_tiebreak", 4, 1, 1, ["anydata"], [], step_split_tiebreak),
    ("missing_variants", 4, 5, 2, ["anydata"], ["vsearch"], step_missing_variants),
    ("figaro_usable_pairs", 4, 3, 1, ["anydata"], [], step_figaro_usable_pairs),
    ("asv_length_check", 4, 2, 1, ["anydata"], [], step_asv_length_check),
    ("figaro_rule_check", 4, 5, 2, ["anydata"], [], step_figaro_rule_check),
    ("richness_sim", 4, 240, 10, [], ["R", "vsearch"], step_richness_sim),
    ("overmerge_sim", 4, 60, 8, [], ["R", "vsearch"], step_overmerge_sim),
]
NOT_TESTABLE = [
    "An independent mock community (Zymo D6321, PRJEB67827, with Zymo's own reference sequences): no data on hand (python validation\\check_build.py --download mockB).",
    "V3-V4 2x250 reads (tight overlap on real data) and V4-V5 (PRJNA971391): no data on hand (check_build.py --download v3v4 / v4v5); mock_clip250 is an emulation, not a replacement.",
    "Other platforms (binned quality scores): no data on hand beyond the exactness checks already made.",
    "Housekeeping: run the CI workflow, push to GitHub, choose a licence, archive for a DOI, test on Windows and Linux, regenerate Figures 3 and 5 and Table 6 from the saved tables, verify the references.",
]


def tool_missing(tool: str) -> bool:
    if tool == "R":
        import rtools
        return rtools.find_rscript() is None
    if tool == "vsearch":
        return shutil.which("vsearch") is None
    import importlib.util
    return importlib.util.find_spec(tool) is None


def readiness(ctx: Ctx, needs: list, tools: list) -> list:
    problems = []
    for n in needs:
        if n == "extra":
            if not ctx.args.extra:
                problems.append("no --extra dataset given")
            else:
                for name, (ready, ref) in parse_extra(ctx.args.extra).items():
                    problems += [f"extra {name}: {x} does not exist" for x in (ready, ref) if x is not None and not x.exists()]
            continue
        if n == "anydata":
            if not (ctx.p("v1v2") or ctx.args.extra or (ctx.p("mockA_high") and ctx.p("reference") and ctx.p("asv_table"))):
                problems.append("needs mock A (--mockA-high, --asv-table, --reference), --v1v2 or --extra")
            continue
        p = ctx.p(n)
        if p is None:
            problems.append(f"no --{n.replace('_', '-')} given" if n != "figaro" else "no --figaro-script given")
        elif not p.exists():
            problems.append(f"{n}: {p} does not exist")
    for t in tools:
        if tool_missing(t):
            problems.append(f"{t} is not available")
    return problems


def ref_tag(path) -> str:
    """A short fingerprint of an expected-sequence file, so that a step is repeated when the reference it scored against has changed."""
    try:
        import hashlib
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:8]
    except (OSError, TypeError):
        return ""


def step_datasets(ctx: Ctx, sid: str):
    """The datasets a step would cover with the paths given now (None for a step that does not depend on the list of datasets). A saved result that covered other datasets
    is repeated, so that a dataset added later is included and a finished dataset-independent step is never repeated for nothing."""
    scored = sid in ("battery", "extra_datasets", "coverage_levels", "landscape")                  # the steps that score against expected sequences name the reference they used
    extras = sorted(f"{n}#{ref_tag(ref)}" if (scored and ref and ref_tag(ref)) else n for n, (ready, ref) in parse_extra(ctx.args.extra).items())
    mock = ([f"mockA_high#{ref_tag(ctx.p('reference'))}" if scored and ref_tag(ctx.p("reference")) else "mockA_high"]
            if ctx.p("mockA_high") and ctx.p("reference") and ctx.p("asv_table") else [])
    v1v2 = ["v1v2"] if ctx.p("v1v2") else []
    if sid == "battery":
        return sorted(mock + v1v2 + (["its"] if ctx.p("its") else []))
    if sid == "extra_datasets":
        return extras
    if sid in ("coverage_levels", "landscape"):
        return sorted(mock + v1v2 + extras)
    if sid == "per_sample_rank":
        return sorted((["mockA_all"] if ctx.p("mockA_all") else []) + v1v2 + extras)
    if sid == "stability":
        return sorted([k for k in ("mockA_high", "mockA_all") if ctx.p(k)] + v1v2 + extras)
    return None


def step_state(ctx: Ctx, sid: str):
    f = ctx.out / sid / "result.json"
    if f.exists():
        try:
            return json.loads(f.read_text(encoding="utf-8"))
        except ValueError:
            return None
    return None


def run_steps(ctx: Ctx, selected: list) -> list:
    results = []
    for sid, tier, minutes, quick_minutes, needs, tools, fn in STEPS:
        if sid not in selected:
            continue
        prev = step_state(ctx, sid)
        covers = step_datasets(ctx, sid)
        if prev and prev.get("status") == "ok" and not ctx.args.fresh:
            saved = prev.get("datasets")
            if saved is None and prev.get("data"):                       # a result saved before the datasets were recorded: the names in its data
                saved = sorted({d.get("name") for d in prev["data"] if isinstance(d, dict) and d.get("name")})
            if covers is None or (saved is not None and sorted(saved) == covers):
                ctx.log(f"{sid}: already done (kept; use --fresh to repeat)")
                results.append(prev)
                continue
            ctx.log(f"{sid}: repeated because it covered {saved if saved else 'an unknown set of datasets'} and now covers {covers}")
        problems = readiness(ctx, needs, tools)
        stepdir = ctx.out / sid
        stepdir.mkdir(parents=True, exist_ok=True)
        ctx.cur = stepdir
        res = {"id": sid, "tier": tier, "status": "skipped", "seconds": 0, "title": sid, "lines": [], "verdicts": [], "error": "; ".join(problems)}
        if problems:
            ctx.log(f"{sid}: SKIPPED ({res['error']})")
        else:
            ctx.log(f"{sid}: starting (tier {tier}, about {quick_minutes if ctx.quick else minutes} min)")
            t0 = time.time()
            try:
                out = fn(ctx)
                res.update(out)
                res.update({"status": "ok", "error": ""})
                if covers is not None:
                    res["datasets"] = covers
            except KeyboardInterrupt:
                raise
            except BaseException as err:                           # SystemExit from the helper scripts too
                res.update({"status": "failed", "error": f"{type(err).__name__}: {err}", "trace": traceback.format_exc()[-1500:]})
                ctx.log(f"{sid}: FAILED: {res['error']}")
            res["seconds"] = round(time.time() - t0)
            ctx.log(f"{sid}: {res['status']} in {res['seconds'] / 60:.1f} min")
        (stepdir / "result.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
        results.append(res)
        write_summary(ctx, results)
    return results


def write_summary(ctx: Ctx, results: list) -> Path:
    """The summary of every step that has a saved result, in the order of the registry; the results of this run take precedence. (A run of a few steps with --only
    used to overwrite the summary with those steps alone.)"""
    by_id = {r["id"]: r for r in results}
    results = [r for r in ((by_id.get(sid) or step_state(ctx, sid)) for sid, *_rest in STEPS) if r]
    L = [f"# AmpliTrim overnight validation ({'QUICK slices' if ctx.quick else 'full'}), written {time.strftime('%Y-%m-%d %H:%M')}", "",
         f"Margins fixed before the run: {json.dumps(MARGINS)}", "", "## Steps", "", "| step | tier | status | minutes | note |", "| --- | --- | --- | --- | --- |"]
    for r in results:
        L.append(f"| {r['id']} | {r['tier']} | {r['status']} | {r['seconds'] / 60:.1f} | {(r.get('error') or '')[:110]} |")
    L += ["", "## Verdicts against the pre-set criteria", "", "| step | criterion | outcome | detail |", "| --- | --- | --- | --- |"]
    for r in results:
        for v in r.get("verdicts", []):
            L.append(f"| {r['id']} | {v['criterion']} | {v['outcome']} | {v['detail']} |")
    for r in results:
        if r["status"] == "ok":
            L += ["", f"## {r['id']}: {r['title']}", "", "```"] + r["lines"] + ["```"]
        elif r["status"] == "failed":
            L += ["", f"## {r['id']}: FAILED", "", "```", r.get("trace", r.get("error", "")), "```"]
    L += ["", "## Not testable with the data on hand", ""] + [f"- {x}" for x in NOT_TESTABLE]
    path = ctx.out / "OVERNIGHT_SUMMARY.md"
    path.write_text("\n".join(L) + "\n", encoding="utf-8")
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--mockA-high", type=Path, help="primer-free reads of the six 10^8 and 10^9 cell mock A samples (a `ready` folder)")
    ap.add_argument("--mockA-all", type=Path, help="primer-free reads of all 23 mock A samples")
    ap.add_argument("--runinfo", type=Path, help="ENA table of the mock A project (runs_info.tsv)")
    ap.add_argument("--asv-table", type=Path, help="the authors' Mock_staggered-A_ASVs.csv")
    ap.add_argument("--reference", type=Path, help="mockA_reference.fasta from build_mockA_reference.py (all 53 variants)")
    ap.add_argument("--v1v2", type=Path, help="primer-free reads of the V1-V2 study (PRJEB42528)")
    ap.add_argument("--its", type=Path, help="primer-free reads of the ITS study (PRJEB58189)")
    ap.add_argument("--figaro-script", type=Path, help="path to figaro.py (optional; the FIGARO steps are skipped without it)")
    ap.add_argument("--extra", action="append", metavar="NAME=READY|REFERENCE", help="a new dataset (primer-free reads, and optionally the expected sequences after |); repeat for several")
    ap.add_argument("--its-deep-samples", default=ITS_DEEP_SAMPLES, help="the ITS samples used at full depth (default: the eight used so far)")
    ap.add_argument("--near-id", type=float, default=0.99)
    ap.add_argument("--quick", action="store_true", help="small slices of every step, to check the wiring (about 30 minutes)")
    ap.add_argument("--only", help="comma-separated step names to run (the steps of tier 4, round 4, run only when named here: check_build.py --round4 does that)")
    ap.add_argument("--skip", help="comma-separated step names to leave out")
    ap.add_argument("--fresh", action="store_true", help="repeat steps that already finished")
    ap.add_argument("--list", action="store_true", help="show the steps, what each needs and whether it is ready, then stop")
    ap.add_argument("--summary", action="store_true", help="rewrite OVERNIGHT_SUMMARY.md from the saved results of every step, without running anything, then stop")
    args = ap.parse_args(argv)
    ctx = Ctx(args)
    if args.summary:
        print(write_summary(ctx, []))
        return 0
    ids = [s[0] for s in STEPS]
    standard = [s[0] for s in STEPS if s[1] <= 3]                    # tier 4 (round 4) runs only when named with --only, which is what `check_build.py --round4` does
    selected = [x for x in (args.only.split(",") if args.only else standard) if x]
    unknown = [x for x in selected + (args.skip.split(",") if args.skip else []) if x and x not in ids]
    if unknown:
        raise SystemExit(f"unknown step(s): {', '.join(unknown)}. Steps: {', '.join(ids)}")
    selected = [x for x in selected if not (args.skip and x in args.skip.split(","))]
    if args.list:
        total = 0
        for sid, tier, minutes, quick_minutes, needs, tools, _fn in STEPS:
            problems = readiness(ctx, needs, tools)
            prev = step_state(ctx, sid)
            state = "done" if prev and prev.get("status") == "ok" else ("READY" if not problems else "skipped: " + "; ".join(problems))
            if sid in selected and not problems and state != "done":
                total += quick_minutes if args.quick else minutes
            print(f"  {sid:16s} tier {tier}  ~{quick_minutes if args.quick else minutes:>3d} min  {state}")
        print(f"\nAbout {total / 60:.1f} hours for the steps that are ready (rough: depends on this machine).")
        return 0
    ctx.log(f"overnight run: {'quick' if args.quick else 'full'}; steps: {', '.join(selected)}; output {ctx.out}")
    results = run_steps(ctx, selected)
    path = write_summary(ctx, results)
    bad = [r for r in results if r["status"] == "failed"]
    ctx.log(f"finished: {sum(r['status'] == 'ok' for r in results)} ok, {len(bad)} failed, {sum(r['status'] == 'skipped' for r in results)} skipped. Summary: {path}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
