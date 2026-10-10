#!/usr/bin/env python
"""Compare the settings that different methods recommend, by running DADA2 on the same reads with each of them.

One folder of primer-free paired reads (the `ready` folder of real_data.py) goes in; for every method the setting a user following it
would use is run through the same DADA2 steps, and one table comes out: what each method predicted, what AmpliTrim's engine says about
that setting, and what DADA2 actually did (pairs passing, merged, after chimera removal, ASVs, the abundance spectrum and, with a
reference, how many expected sequences were found and how many extra ASVs appeared).

Methods:
  AmpliTrim   one row per --amplitrim-max-ee value (default 2): the truncation AmpliTrim recommends for that maxEE. When FIGARO is
              included, one more row uses FIGARO's own maxEE, so that only the truncation differs (the like-for-like comparison).
  FIGARO      with --figaro-script: FIGARO's own top entry (trim positions AND the maxEE it pairs with them), run on the pairs whose
              reads have the most common length (FIGARO needs reads of one length); one row per --figaro-amplicon (default: the
              longest amplicon AmpliTrim measured).
  ITS         with --marker its, AmpliTrim's rows are DADA2's ITS workflow (no truncation, minLen from --min-len) at each --amplitrim-max-ee value,
              with AmpliTrim's own estimate for each; add --setting rows such as trunc200=200,200,2 to show what a truncation copied from 16S costs.
  fixed       --setting NAME=F,R,EE for any other setting: a paper's, a tutorial's, ...  EE is one value or two. A truncLen of 0 means
              no truncation, as in DADA2 (full=0,0,2); AmpliTrim then has no estimate to show, since it models cut reads.

It also reports, across the settings, how far AmpliTrim's estimate is from DADA2 (the offset and the rank agreement) and, with six or more
samples, the mean per-sample difference in the share of pairs DADA2 merged against the baseline (the first row, or --baseline), with a bootstrap
95% interval over samples. --pairwise FOLDER redoes that last comparison against another baseline from the saved tables, without running DADA2.

  python validation/compare_methods.py READS --out DIR [--amplitrim-max-ee 2 2,6 3,5] [--figaro-script PATH]
         [--setting authors=279,263,2,6 --setting full=280,283,2] [--reference FASTA] [--pairs N]
"""

import argparse
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def parse_setting(text: str):
    """NAME=F,R,EE or NAME=F,R,EEF,EER -> (name, truncLenF, truncLenR, maxEE text)."""
    name, _, rest = text.partition("=")
    nums = [x.strip() for x in rest.split(",") if x.strip()]
    try:
        ok = bool(name) and len(nums) in (3, 4) and int(nums[0]) >= 0 and int(nums[1]) >= 0 and all(float(x) > 0 for x in nums[2:])
    except ValueError:
        ok = False
    if not ok:
        raise SystemExit(f"--setting {text!r}: write it as NAME=truncLenF,truncLenR,maxEE or NAME=truncLenF,truncLenR,maxEEF,maxEER "
                         "(for example authors=279,263,2,6). A truncLen of 0 means no truncation, as in DADA2: full=0,0,2.")
    return name, int(nums[0]), int(nums[1]), ",".join(nums[2:])


def paired_difference(a: list, b: list, n_boot: int = 10000, seed: int = 1):
    """Mean of a - b over samples, with a bootstrap 95% interval over samples, and how many samples a is above / below b."""
    import numpy as np
    d = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    boots = np.random.default_rng(seed).choice(d, size=(n_boot, len(d)), replace=True).mean(axis=1)
    return float(d.mean()), float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)), int((d > 0).sum()), int((d < 0).sum())


def sample_merged_share(track: Path) -> dict:
    """{sample: pairs merged / pairs given to DADA2} from the per-sample table dada2_real.R writes."""
    out = {}
    with open(track, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if float(r["input"]) > 0:
                out[r["sample"]] = float(r["merged"]) / float(r["input"])
    return out


def _read_fastq_gz(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        while True:
            rec = [fh.readline() for _ in range(4)]
            if not rec[0]:
                return
            yield rec


def subset_reads(reads: Path, out_dir: Path, spec: str) -> Path:
    """Copy only the named samples (comma-separated, or a text file with one name per line) into out_dir, so that DADA2 can be given all
    the reads of a few samples instead of a thin slice of many. A name that is not among the samples stops the run."""
    from amplitrim.subsample import find_pairs
    import shutil
    p = Path(spec)
    names = [ln.strip() for ln in p.read_text(encoding="utf-8-sig").splitlines()] if p.is_file() else spec.split(",")
    names = [n.strip() for n in names if n.strip() and not n.strip().startswith("#")]
    have = {s: (r1, r2) for s, r1, r2 in find_pairs(reads)}
    missing = [n for n in names if n not in have]
    if missing:
        raise SystemExit(f"--only names samples that are not in {reads}: {', '.join(missing)}. Samples found: {', '.join(sorted(have)[:6])}{' ...' if len(have) > 6 else ''}")
    out_dir.mkdir(parents=True, exist_ok=True)
    for n in names:
        for src in have[n]:
            shutil.copy2(src, out_dir / src.name)
    return out_dir


def figaro_input(reads: Path, out_dir: Path, max_pairs: int):
    """FIGARO needs reads of one length, so write an even subsample of the pairs whose two reads have the most common lengths.
    Returns (folder, share of the subsample kept, forward length, reverse length)."""
    from amplitrim.subsample import build_subsample
    sub = build_subsample(reads, out_dir / "all", max_pairs=max_pairs)
    pairs = list(zip(_read_fastq_gz(Path(sub.r1)), _read_fastq_gz(Path(sub.r2))))
    lf = Counter(len(a[1].strip()) for a, _ in pairs).most_common(1)[0][0]
    lr = Counter(len(b[1].strip()) for _, b in pairs).most_common(1)[0][0]
    keep = [(a, b) for a, b in pairs if len(a[1].strip()) == lf and len(b[1].strip()) == lr]
    folder = out_dir / "uniform"
    folder.mkdir(parents=True, exist_ok=True)
    with gzip.open(folder / "subsample_R1.fastq.gz", "wt", newline="\n") as f1, gzip.open(folder / "subsample_R2.fastq.gz", "wt", newline="\n") as f2:
        for a, b in keep:
            f1.write("".join(a))
            f2.write("".join(b))
    return folder, len(keep) / max(len(pairs), 1), lf, lr


def min_len_for(marker: str, trunc_f: int, trunc_r: int, min_len: int) -> int:
    """DADA2's minLen for one arm. The ITS workflow needs one, and so does any arm without truncation (a 0 length): there the filter cuts reads only at a low quality, and a
    read cut to almost nothing stops DADA2 ("sequences must be longer than the kmer-size"). An arm with truncation lengths discards the shorter reads itself."""
    return min_len if (marker == "its" or trunc_f == 0 or trunc_r == 0) else 0


def method_key(name: str) -> str:
    """The part of a method's name used in the file names of its ASVs and per-sample table."""
    return "".join(c if c.isalnum() else "_" for c in name)


def pick_baseline(requested, names: list) -> str:
    """The method the others are compared with: the one asked for (its full name, or text that matches exactly one), else the first."""
    if not requested:
        return names[0]
    if requested in names:
        return requested
    low = requested.lower()
    hits = [n for n in names if n.lower().startswith(low)] or [n for n in names if low in n.lower()]       # a start of the name wins over a piece of it
    if len(hits) != 1:
        raise SystemExit(f"--baseline {requested!r} matches {'more than one' if hits else 'no'} method. The methods are: " + "; ".join(names))
    return hits[0]


def paired_lines(shares: dict, names: list, base: str) -> list:
    """Per-sample difference in the merged share of every method against the baseline, with a bootstrap interval over samples."""
    lines = ["", f"Per-sample difference in the share of pairs DADA2 merged, against {base} (the first method unless --baseline says otherwise); "
                 f"mean, bootstrap 95% interval over samples:"]
    for name in names:
        if name == base:
            continue
        common = sorted(set(shares[name]) & set(shares[base]))
        if len(common) < 6:
            lines.append(f"  {name}: only {len(common)} samples, too few for an interval")
            continue
        mean, lo, hi, up, down = paired_difference([shares[name][s] for s in common], [shares[base][s] for s in common])
        lines.append(f"  {name}: {mean * 100:+.1f} points (95% interval {lo * 100:+.1f} to {hi * 100:+.1f}); better in {up} of {len(common)} samples, worse in {down}")
    return lines


def pairwise_report(folder: Path, baseline=None) -> list:
    """The per-sample comparison again, from the tables an earlier run saved, against any method: no DADA2 is run."""
    table = folder / "compare_methods.tsv"
    if not table.is_file():
        raise SystemExit(f"{table} was not found: give the output folder of an earlier compare_methods.py run.")
    with open(table, encoding="utf-8") as fh:
        names = [r["name"] for r in csv.DictReader(fh, delimiter="\t") if not r.get("error")]
    shares = {n: sample_merged_share(folder / f"track_{method_key(n)}.tsv") for n in names if (folder / f"track_{method_key(n)}.tsv").is_file()}
    if len(shares) < 2:
        raise SystemExit(f"{folder}: fewer than two methods have a saved per-sample table.")
    return paired_lines(shares, list(shares), pick_baseline(baseline, list(shares)))


def calibration_lines(rows: list) -> list:
    """How well AmpliTrim's estimates match what DADA2 did across the settings in the table, and FIGARO's predicted retention against
    DADA2's pass rate. The difference is reported as well as the rank agreement, because a constant offset does not change which setting wins."""
    from per_sample import spearman
    ok = [r for r in rows if not r.get("error") and r.get("est_yield") is not None and r.get("merged") is not None]
    out = []
    if len(ok) >= 3:
        diff = [100 * (r["est_yield"] - r["merged"]) for r in ok]
        dpass = [100 * (r["est_pass"] - r["pass"]) for r in ok]
        line = (f"Across {len(ok)} settings, AmpliTrim's estimate minus DADA2: yield against merged share {sum(diff) / len(diff):+.1f} points "
                f"(range {min(diff):+.1f} to {max(diff):+.1f}); filter pass {sum(dpass) / len(dpass):+.1f} points (range {min(dpass):+.1f} to {max(dpass):+.1f}).")
        if len(ok) >= 4:
            line += f" Rank correlation of the estimated yield with DADA2's merged share: {spearman([r['est_yield'] for r in ok], [r['merged'] for r in ok]):+.2f}."
        ratios = [r["merged"] / r["est_yield"] for r in ok if r["est_yield"] > 0]
        if ratios:                                   # a shortfall that scales with the yield looks like a constant ratio, not a constant offset
            line += f" As a fraction of the estimated yield, DADA2's merged share is {sum(ratios) / len(ratios):.2f} (range {min(ratios):.2f} to {max(ratios):.2f})."
        out.append(line)
        if len(ok) >= 4 and all(r.get("est_strict") is not None for r in ok):
            sdiff = [100 * (r["est_strict"] - r["merged"]) for r in ok]
            out.append(f"The strict estimate (no mismatch allowed in the overlap) minus DADA2's merged share: {sum(sdiff) / len(sdiff):+.1f} points "
                       f"(range {min(sdiff):+.1f} to {max(sdiff):+.1f}); rank correlation {spearman([r['est_strict'] for r in ok], [r['merged'] for r in ok]):+.2f}.")
    fig = [r for r in rows if r["name"].startswith("FIGARO") and not r.get("error") and r.get("own") is not None]
    if fig:
        out.append("FIGARO's predicted retention minus DADA2's pass rate: " + "; ".join(f"{100 * (r['own'] - r['pass']):+.1f} points ({r['name']})" for r in fig) +
                   ". It predicts the filter only, so it is not expected to match the merged share.")
    return out


def pct(x, nd=1):
    return "" if x is None else f"{x:.{nd}%}"


def format_table(rows: list, with_reference: bool) -> list:
    head = (f"{'method':<30s} {'trunc':>8s} {'maxEE':>6s} | {'own prediction':>15s} {'AmpliTrim est.':>14s} | {'DADA2 pass':>10s} {'merged':>7s} "
            f"{'no chim':>8s} {'ASVs':>5s} | {'top 15':>7s} {'>=0.1%':>6s} {'tail':>6s}")
    with_strains = with_reference and any(r.get("strains") for r in rows)
    if with_reference:
        head += f" | {'found':>5s} {'identical':>9s} {'extra':>5s}" + (f" {'strains':>7s}" if with_strains else "")
    out = [head]
    for r in rows:
        if r.get("error"):
            out.append(f"{r['name']:<30s} {r['truncF']}/{r['truncR']:<4} {r['maxEE']:>6s} | {r['error']}")
            continue
        est = "" if r.get("est_yield") is None else (f"{r['est_pass'] * 100:.1f}/{r['est_yield'] * 100:.1f}%" if r.get("est_pass") is not None else pct(r["est_yield"]))
        trunc = "none" if r["truncF"] == 0 and r["truncR"] == 0 else f"{r['truncF']}/{r['truncR']}"
        line = (f"{r['name']:<30s} {trunc:>8s} {r['maxEE']:>6s} | "
                f"{(pct(r['own']) + ' ' + r['own_what']).rjust(15)} {est:>14s} | {pct(r['pass']):>10s} {pct(r['merged']):>7s} "
                f"{pct(r['nonchim']):>8s} {r['asvs']:>5d} | {pct(r['top15'], 2):>7s} {r['n01']:>6d} {pct(r['tail'], 2):>6s}")
        if with_reference:
            line += f" | {r['found']:>5d} {pct(r['identical'], 2):>9s} {r['extra']:>5d}" + (f" {r['strains'] or '':>7s}" if with_strains else "")
        out.append(line)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reads", type=Path, nargs="?", help="folder of primer-free paired FASTQs (the `ready` folder of real_data.py)")
    ap.add_argument("--pairwise", type=Path, metavar="FOLDER", help="no new run: redo the per-sample comparison from the tables saved in FOLDER, against --baseline")
    ap.add_argument("--out", type=Path, default=Path("compare_methods"))
    ap.add_argument("--amplitrim-max-ee", nargs="+", default=["2"], metavar="EE", help="one AmpliTrim row per value, for example 2 2,6 3,5 (default 2)")
    ap.add_argument("--figaro-script", type=Path, help="path to figaro.py: add FIGARO's own recommendation")
    ap.add_argument("--figaro-python", default=sys.executable, help="Python that runs FIGARO (default: this one)")
    ap.add_argument("--figaro-amplicon", type=int, nargs="+", metavar="NT", help="amplicon length(s) to give FIGARO (default: the longest AmpliTrim measured)")
    ap.add_argument("--figaro-pairs", type=int, default=20000, metavar="N", help="pairs given to FIGARO (default 20000)")
    ap.add_argument("--setting", action="append", default=[], metavar="NAME=F,R,EE", help="a fixed setting to include; repeat for several")
    ap.add_argument("--baseline", help="the method the others are compared with: its name, or text that matches only one (default: the first row)")
    ap.add_argument("--only", metavar="NAMES", help="use only these samples (comma-separated, or a file with one name per line); with --pairs 0 DADA2 then gets all their reads")
    ap.add_argument("--pairs", type=int, default=0, metavar="N", help="give DADA2 at most N evenly spaced pairs per sample (default all)")
    ap.add_argument("--allow-unmergeable", type=float, nargs="+", metavar="FRAC", help="also one AmpliTrim row per value (at the first --amplitrim-max-ee) that may leave this share of the longest amplicons unmergeable; the default, 0.002, is the first row")
    ap.add_argument("--min-overlap", type=int, default=12)
    ap.add_argument("--marker", choices=["16s", "its"], default="16s", help="its: amplicons vary in length, so AmpliTrim's rows are DADA2's ITS workflow (no truncation) at each maxEE")
    ap.add_argument("--min-len", type=int, default=50, help="minLen given to DADA2 with --marker its (default 50)")
    ap.add_argument("--reference", type=Path, help="FASTA of the expected sequences (a mock community); name variants of one strain >ID|STRAIN|label")
    ap.add_argument("--near-id", type=float, default=0.99, help="identity at or above which a non-identical ASV counts as near a reference (default 0.99, about 4-5 bases)")
    args = ap.parse_args()
    if args.pairwise:
        print("\n".join(pairwise_report(args.pairwise, args.baseline)))
        return 0
    if args.reads is None:
        ap.error("give the folder of reads, or --pairwise FOLDER")
    # check everything that can be checked before any work starts, so that a typo costs a second and not the first run
    fixed = [parse_setting(t) for t in args.setting]
    if args.marker == "its" and args.figaro_script:
        raise SystemExit("--figaro-script cannot be used with --marker its: FIGARO needs reads of one length and searches truncation lengths, and ITS reads are not truncated.")
    if not args.reads.is_dir():
        raise SystemExit(f"{args.reads} is not a folder of reads.")
    if args.figaro_script and not args.figaro_script.is_file():
        raise SystemExit(f"--figaro-script {args.figaro_script}: no such file. Give the full path of figaro.py (find it with: dir /s /b figaro.py).")
    if args.reference and not args.reference.is_file():
        raise SystemExit(f"--reference {args.reference}: no such file.")

    if args.only:
        from amplitrim.subsample import find_pairs
        have = {s for s, _r1, _r2 in find_pairs(args.reads)}
        wanted = [n.strip() for n in (args.only.split(",") if not Path(args.only).is_file() else Path(args.only).read_text(encoding="utf-8-sig").splitlines()) if n.strip() and not n.strip().startswith("#")]
        if [n for n in wanted if n not in have]:
            raise SystemExit(f"--only names samples that are not in {args.reads}: {', '.join(n for n in wanted if n not in have)}")

    from amplitrim.config import parse_ee
    from amplitrim.pipeline import PipelineError, Settings, run
    from max_ee_sweep import abundance_summary
    from per_sample import per_sample
    from real_data import Dada2Error, dada2_run, ee_text, thin_reads
    if args.reference:
        import shutil
        import mock_check
        if not shutil.which("vsearch"):
            raise SystemExit("vsearch was not found on PATH (needed to score the ASVs against the reference).")
    args.out.mkdir(parents=True, exist_ok=True)
    reads = subset_reads(args.reads, args.out / "subset", args.only) if args.only else args.reads
    reads = thin_reads(reads, args.out / "reads", args.pairs) if args.pairs else reads

    methods, notes, state = [], [], {"longest": None}

    def add_amplitrim(ee: str, name: str, max_unmergeable: float = None):
        extra = {} if max_unmergeable is None else {"max_unmergeable": max_unmergeable}
        try:
            res = run(Settings(input_dir=reads, work_dir=args.out / ("amplitrim_" + ee.replace(",", "_") + ("" if max_unmergeable is None else f"_leave{max_unmergeable:g}")), max_ee=parse_ee(ee), min_overlap=args.min_overlap,
                               marker=args.marker, min_len=args.min_len, **extra))
        except PipelineError as err:
            notes.append(f"{name}: {err}")
            return
        state["longest"] = state["longest"] or res.longest_insert
        if res.recommendation is None:
            notes.append(f"{name}: no setting could merge the longest amplicons")
            return
        rec = res.recommendation
        methods.append({"name": name, "tf": rec.trim_f, "tr": rec.trim_r, "ee": ee_text(ee), "own": rec.yield_est, "own_what": "yield",
                        # ITS is not truncated, so there is no setting to re-estimate later: AmpliTrim's own estimate is the one to show
                        "est": (rec.yield_est, rec.ee_pass_rate, None) if args.marker == "its" else None})     # no strict estimate: long overlaps make it meaningless

    for ee in args.amplitrim_max_ee:
        add_amplitrim(ee, f"AmpliTrim (maxEE {ee})")
    for v in (args.allow_unmergeable or []):                       # the coverage setting: how much of the longest tail may stay unmergeable
        add_amplitrim(args.amplitrim_max_ee[0], f"AmpliTrim (maxEE {args.amplitrim_max_ee[0]}, leave out {v:.1%})", max_unmergeable=v)
    if args.figaro_script:
        from amplitrim.figaro_wrapper import parse_figaro_recommendation, run_figaro
        for a in (args.figaro_amplicon or ([state["longest"]] if state["longest"] else [])):
            try:
                folder, kept, lf, lr = figaro_input(reads, args.out / f"figaro_{a}", args.figaro_pairs)
                rec = parse_figaro_recommendation(run_figaro(args.figaro_script, args.figaro_python, folder, args.out / f"figaro_{a}" / "out", a, args.min_overlap))
            except (RuntimeError, ValueError, OSError) as err:
                notes.append(f"FIGARO (amplicon {a}): {err}")
                continue
            notes.append(f"FIGARO (amplicon {a}) was run on the {kept:.0%} of the sampled pairs whose reads are {lf} and {lr} nt long (it needs reads of one length).")
            methods.append({"name": f"FIGARO (amplicon {a})", "tf": rec.trim_f, "tr": rec.trim_r, "ee": ee_text(",".join(str(x) for x in rec.max_ee)),
                            "own": rec.retention, "own_what": "retained"})
    # like for like: AmpliTrim at FIGARO's own maxEE, so that only the choice of truncation differs between the two
    have = {m["ee"] for m in methods if m["name"].startswith("AmpliTrim")}
    for m in [m for m in methods if m["name"].startswith("FIGARO")]:
        if m["ee"] not in have:
            add_amplitrim(m["ee"], f"AmpliTrim (FIGARO's maxEE {m['ee']})")
            have.add(m["ee"])
    for name, tf, tr, ee in fixed:
        methods.append({"name": name, "tf": tf, "tr": tr, "ee": ee_text(ee), "own": None, "own_what": ""})
    if not methods:
        raise SystemExit("No method produced a setting to compare:\n  " + "\n  ".join(notes))

    rows, shares = [], {}
    for m in methods:
        key = method_key(m["name"])
        row = {"name": m["name"], "truncF": m["tf"], "truncR": m["tr"], "maxEE": m["ee"], "own": m["own"], "own_what": m["own_what"]}
        if m.get("est"):                                  # ITS: the estimate AmpliTrim made for this very setting
            row["est_yield"], row["est_pass"], row["est_strict"] = m["est"]
        elif m["tf"] == 0 or m["tr"] == 0:                # no truncation: AmpliTrim models cut reads, so it has no estimate here
            row["est_yield"] = row["est_pass"] = row["est_strict"] = None
        else:
            est = per_sample(reads, m["tf"], m["tr"], 4000, args.min_overlap, parse_ee(m["ee"]))
            n = sum(x["pairs_used"] for x in est) or 1
            row["est_yield"] = sum(x["estimated_yield"] * x["pairs_used"] for x in est) / n
            row["est_pass"] = sum(x["pass_filter"] * x["pairs_used"] for x in est) / n
            row["est_strict"] = sum(x["strict_yield"] * x["pairs_used"] for x in est) / n
        asv, track = args.out / f"asvs_{key}.fasta", args.out / f"track_{key}.tsv"
        try:
            # no truncation (a 0 length) needs a minimum length, or the shortest reads stop DADA2 ("sequences must be longer than the kmer-size")
            d = dada2_run(reads, m["tf"], m["tr"], m["ee"], args.min_overlap, min_len_for(args.marker, m["tf"], m["tr"], args.min_len), track, asv)
        except Dada2Error as err:
            row["error"] = f"DADA2 could not run: {err}"
            rows.append(row)
            print(format_table([row], False)[1])
            continue
        spectrum = abundance_summary(asv)
        row.update({"pass": d["filtered"], "merged": d["merged"], "nonchim": d["nonchim"], "asvs": d["asvs"],
                    "top15": spectrum["top15"], "n01": spectrum["n_01"], "tail": spectrum["tail"]})
        if args.reference:
            s = mock_check.summarise(asv, args.reference, args.near_id, vsearch="vsearch")
            row.update({"found": len(s["recovered"]), "identical": s["categories"]["exact"][1] / max(s["reads"], 1),
                        "extra": s["categories"]["near"][0] + s["categories"]["none"][0], "strains": f"{s['strains_found']}/{s['strains_total']}" if s.get("strains_total") else None})
        shares[m["name"]] = sample_merged_share(track)
        rows.append(row)
        print(format_table([row], bool(args.reference))[1])

    table = format_table(rows, bool(args.reference))
    report = ["own prediction: AmpliTrim's estimated yield (pairs passing the filter and merging), or FIGARO's predicted retention (the filter only).",
              "AmpliTrim est.: AmpliTrim's engine applied to that method's setting, as filter pass / yield, so every row is judged by one estimator.", *notes]
    report += calibration_lines(rows)
    good = [r["name"] for r in rows if not r.get("error") and r["name"] in shares]
    if len(good) >= 2:
        report += paired_lines(shares, good, pick_baseline(args.baseline, good))
    print("\n" + "\n".join(table) + "\n\n" + "\n".join(report))
    (args.out / "compare_methods.txt").write_text("\n".join(table + [""] + report) + "\n", encoding="utf-8")
    with open(args.out / "compare_methods.tsv", "w", newline="", encoding="utf-8") as fh:
        keys = ["name", "truncF", "truncR", "maxEE", "own", "own_what", "est_yield", "est_pass", "est_strict", "strains", "pass", "merged", "nonchim", "asvs", "top15", "n01", "tail", "found", "identical", "extra", "error"]
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved: {args.out / 'compare_methods.txt'} and {args.out / 'compare_methods.tsv'}; the ASVs and per-sample DADA2 tables of each method are in the same folder.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
