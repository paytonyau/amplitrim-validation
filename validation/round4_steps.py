"""The steps of round 4: the open questions of the validation, each with the criterion fixed in TEST_REGISTER.md (section 9) before it was run.

They are registered in overnight.py (tier 4) and run, in phases with breaks, by `python validation\\check_build.py --round4`.

  Phase A, from results already on disk (no DADA2; minutes)
    split_tiebreak       Q28  does a forward-heavy split beat the balanced one that AmpliTrim prefers among near-ties? (re-reads the landscape grids)
    missing_variants     Q29  which expected variant does one tool recover and the other not, and which read covers the bases that tell it from its sibling?
    figaro_usable_pairs  Q30  how many pairs of each dataset a tool that needs reads of one length can use
    asv_length_check     Q23  do the ASVs of the mocks have the lengths of the expected amplicons? (false joins make wrong-length ASVs)
    figaro_rule_check    Q16  does 'mean expected errors at the cut + 1' reproduce the maxEE that FIGARO reports?
  Phase B, simulations with known truth (DADA2; hours)
    richness_sim         Q0-14  does the excess of the yield estimate over DADA2's merged share grow with the number of taxa?
    overmerge_sim        Q21-27 how many merges are at the wrong length (counted over pairs without an indel), and does the 12 nt minimum overlap hold, with and without repeats, chimeras and indels?

The functions that take (ctx, H) are the steps; H carries the helpers of overnight.py (verdict, parse_extra, table_lines, recommended, high_input_reference), so that this module
does not import overnight.py (which imports it).
"""
import difflib
import re
import shutil
import statistics
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
for p in (HERE, HERE / "simulated", HERE.parent):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

# ----------------------------------------------------------------------------------------------------------------------------------------------- the criteria (fixed in advance)
SPLIT_MIN_DATASETS = 4          # of the datasets with a same-total alternative: the best same-total cell has the longer forward read
SPLIT_MIN_MEAN_GAIN = 1.0       # points of merged pairs, mean over those datasets
RULE_TOLERANCE = 1.0            # expected errors: the predicted maxEE is within this of FIGARO's on both reads
RULE_MIN_DATASETS = 5           # of the six datasets with a FIGARO row
RICHNESS_SUPPORTED = 10.0       # points: the excess at ~500 taxa minus the excess at 25 taxa
RICHNESS_NOT_SUPPORTED = 3.0
LENGTH_TOLERANCE = 5            # nt: an ASV of a mock counts as the length of an expected amplicon within this
LENGTH_MAX_OFF_SHARE = 0.01     # mock A at the default setting: at most this share of the ASV reads off every expected length
MISMERGE_MAX_SHARE = 0.01       # realistic arm, 12 nt: at most this share of the merged reads (of pairs without an indel) at the wrong length
MISMERGE_MAX_RECOVERY_LOSS = 1  # DADA2 at 12 nt may recover at most this many taxa fewer than at 20 nt (recovering more is not a failure)


def mean(values):
    values = [v for v in values if v is not None]
    return statistics.fmean(values) if values else float("nan")


def method_key(name: str) -> str:
    """The part of a method's name used in the file names of its ASVs (the same rule as compare_methods.py)."""
    return "".join(c if c.isalnum() else "_" for c in name)


# ----------------------------------------------------------------------------------------------------------------------------------------------- Q28: the split
HEAD = re.compile(r"^--- (?P<name>\S+) \(maxEE (?P<ee>[^)]*)\); AmpliTrim recommends (?P<f0>\d+)/(?P<r0>\d+) \(coverage floor (?P<floor>\d+)\)")
ROW = re.compile(r"^  (?P<a>[ \d]{8}) ")


def parse_landscape_lines(lines: list) -> list:
    """The grids of the landscape step, read back from the text it saved: [{name, ee, f0, r0, floor, cells: [{tf, tr, merged, meets, mine}]}]."""
    grids, cur = [], None
    for ln in lines:
        m = HEAD.match(ln)
        if m:
            cur = {"name": m["name"], "ee": m["ee"], "f0": int(m["f0"]), "r0": int(m["r0"]), "floor": int(m["floor"]), "trs": None, "cells": []}
            grids.append(cur)
            continue
        if cur is None:
            continue
        if ln.startswith("  trF"):
            cur["trs"] = [int(x) for x in ln.split()[1:]]
            continue
        row = ROW.match(ln)
        if cur["trs"] and row and row["a"].strip().isdigit():
            tf, rest = int(row["a"]), ln[11:]
            for i, tr in enumerate(cur["trs"]):
                cell = rest[i * 8:i * 8 + 7]
                if len(cell) < 7 or cell.strip() == "-":
                    continue
                cur["cells"].append({"tf": tf, "tr": tr, "merged": float(cell[:6]) / 100, "mine": cell[6] == "*", "meets": cell[6] != "."})
    return [g for g in grids if any(c["mine"] for c in g["cells"])]


def same_total_analysis(grid: dict):
    """Among the cells with AmpliTrim's total length, how much better is the best, and does it have the longer forward read? None when AmpliTrim's cell has no same-total neighbour."""
    mine = next((c for c in grid["cells"] if c["mine"]), None)
    if mine is None:
        return None
    same = [c for c in grid["cells"] if c["tf"] + c["tr"] == mine["tf"] + mine["tr"]]
    if len(same) < 2:
        return None
    best = max(same, key=lambda c: c["merged"])
    df = best["tf"] - mine["tf"]
    return {"name": grid["name"], "mine": f"{mine['tf']}/{mine['tr']}", "mine_merged": mine["merged"], "best": f"{best['tf']}/{best['tr']}", "best_merged": best["merged"],
            "gain": 100 * (best["merged"] - mine["merged"]), "df": df, "direction": "forward-heavier" if df > 0 else ("AmpliTrim's own" if df == 0 else "reverse-heavier"), "n_same": len(same)}


def split_decision(results: list) -> dict:
    forward = [r for r in results if r["df"] > 0 and r["gain"] > 0]
    gain = mean([r["gain"] for r in forward])
    return {"n": len(results), "forward": len(forward), "mean_gain": gain, "supported": len(forward) >= SPLIT_MIN_DATASETS and gain >= SPLIT_MIN_MEAN_GAIN}


def step_split_tiebreak(ctx, H) -> dict:
    import json
    f = results_folder(ctx) / "landscape" / "result.json"
    if not f.is_file():
        raise RuntimeError("the landscape step has not run: it holds the grids that this step re-reads (run it first)")
    saved = json.loads(f.read_text(encoding="utf-8"))
    grids = parse_landscape_lines(saved.get("lines", []))
    results = [r for r in (same_total_analysis(g) for g in grids) if r]
    if not results:
        raise RuntimeError("no landscape grid with a same-total neighbour was found in the saved result")
    rows = [{"name": r["name"], "mine": r["mine"], "m": f"{r['mine_merged']:.1%}", "best": r["best"], "b": f"{r['best_merged']:.1%}", "gain": f"{r['gain']:+.1f}", "dir": r["direction"]} for r in results]
    lines = ["Among the cells with AmpliTrim's own total length (the landscape grids, DADA2 merged share): is the best one forward-heavier than AmpliTrim's balanced pick?"]
    lines += H.table_lines(rows, [("dataset", "name", 12), ("AmpliTrim", "mine", 9), ("merged", "m", 7), ("best same total", "best", 15), ("merged", "b", 7), ("gain", "gain", 6), ("the best is", "dir", 16)])
    d = split_decision(results)
    lines += ["", f"the best same-total cell has the longer forward read in {d['forward']} of {d['n']} datasets; mean gain over those {d['mean_gain']:.1f} points"]
    ok = d["supported"] if d["n"] >= 6 else None
    verdicts = [H.verdict(f"a forward-heavy split beats the balanced tie-break (best same-total cell forward-heavier in at least {SPLIT_MIN_DATASETS} datasets, mean gain at least {SPLIT_MIN_MEAN_GAIN} point)", ok,
                          f"{d['forward']} of {d['n']} datasets, mean gain {d['mean_gain']:.1f} points" + ("" if d["n"] >= 6 else " (fewer than six datasets: no verdict)"))]
    return {"title": "Q28: is a forward-heavy split better than the balanced one?", "lines": lines, "verdicts": verdicts, "files": [], "data": results}


# ----------------------------------------------------------------------------------------------------------------------------------------------- Q29: the missing variant
def edit_positions(a: str, b: str) -> list:
    """1-based positions of `a` that differ from `b` (a base of `b` that `a` lacks is marked at the base of `a` before it)."""
    pos = []
    for tag, i1, i2, _j1, _j2 in difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes():
        if tag in ("replace", "delete"):
            pos += list(range(i1 + 1, i2 + 1))
        elif tag == "insert":
            pos.append(max(1, i1))
    return sorted(set(pos))


def nearest_sibling(name: str, seqs: dict):
    best, best_ratio = None, -1.0
    for other, s in seqs.items():
        if other == name:
            continue
        r = difflib.SequenceMatcher(None, seqs[name], s, autojunk=False).ratio()
        if r > best_ratio:
            best, best_ratio = other, r
    return best


def coverage_class(pos: int, length: int, f: int, r: int) -> str:
    """Which read(s) of a pair cover a base: the forward read covers 1..f of the amplicon, the reverse read covers length-r+1..length."""
    by_f, by_r = pos <= f, pos >= length - r + 1
    return "both" if (by_f and by_r) else ("forward only" if by_f else ("reverse only" if by_r else "neither"))


def coverage_counts(positions: list, length: int, f: int, r: int) -> dict:
    out = {"both": 0, "forward only": 0, "reverse only": 0, "neither": 0}
    for p in positions:
        out[coverage_class(p, length, f, r)] += 1
    return out


def _asv_file(folder: Path, name: str) -> Path:
    return folder / f"asvs_{method_key(name)}.fasta"


def results_folder(ctx) -> Path:
    """Where the earlier steps saved their results. A --quick run has a folder of its own (<out>_quick), which holds no landscape grids or comparisons, so it reads the full run's folder when there is one."""
    out = Path(ctx.out)
    if out.name.endswith("_quick") and out.with_name(out.name[:-len("_quick")]).is_dir():
        return out.with_name(out.name[:-len("_quick")])
    return out


def _have(path) -> bool:
    return path is not None and Path(path).exists()


def datasets_with_comparison(ctx, H) -> list:
    """(name, reads folder, reference or None, comparison folder) for every dataset on hand whose compare_methods.py output exists."""
    items = []
    if ctx.p("mockA_high"):
        items.append(("mockA_high", ctx.p("mockA_high"), H.high_input_reference(ctx) if (_have(ctx.p("reference")) and _have(ctx.p("asv_table"))) else None))
    if ctx.p("v1v2"):
        items.append(("v1v2", ctx.p("v1v2"), None))
    for name, (ready, ref) in H.parse_extra(ctx.args.extra).items():
        items.append((name, ready, ref))
    out = []
    for name, ready, ref in items:
        folder = results_folder(ctx) / ("battery" if name in ("mockA_high", "v1v2", "its") else "extra_datasets") / name / "compare"
        if (folder / "compare_methods.tsv").is_file():
            out.append((name, ready, ref, folder))
    return out


def step_missing_variants(ctx, H) -> dict:
    import maxee_frontier
    import mock_check
    lines, data, supports, total = [], [], 0, 0
    for name, _ready, ref, folder in datasets_with_comparison(ctx, H):
        if ref is None:
            continue
        methods = maxee_frontier.read_methods(folder)
        fig = next((m for m in methods if m["name"].upper().startswith("FIGARO")), None)
        mine = next((m for m in methods if m["name"].startswith("AmpliTrim (FIGARO")), None)
        if fig is None or mine is None or not _asv_file(folder, fig["name"]).is_file() or not _asv_file(folder, mine["name"]).is_file():
            lines.append(f"--- {name}: the FIGARO row, the AmpliTrim row at FIGARO's maxEE or their ASV files are missing from {folder}")
            continue
        refs = {n.split()[0]: s for n, s in mock_check.read_fasta(ref)}
        got = {}
        for label, m in (("FIGARO", fig), ("AmpliTrim", mine)):
            s = mock_check.summarise(_asv_file(folder, m["name"]), ref, ctx.args.near_id, "vsearch")
            got[label] = {"set": set(s["recovered"]), "n": len(s["recovered"]), "split": (m["tf"], m["tr"]), "ee": m["ee"]}
        lines.append(f"--- {name} (maxEE {fig['ee'][0]:g},{fig['ee'][1]:g}): FIGARO {got['FIGARO']['split'][0]}/{got['FIGARO']['split'][1]} recovered {got['FIGARO']['n']} of {len(refs)} expected variants; "
                     f"AmpliTrim {got['AmpliTrim']['split'][0]}/{got['AmpliTrim']['split'][1]} recovered {got['AmpliTrim']['n']}")
        for who, other in (("FIGARO", "AmpliTrim"), ("AmpliTrim", "FIGARO")):
            for v in sorted(got[who]["set"] - got[other]["set"]):
                sib = nearest_sibling(v, refs)
                pos = edit_positions(refs[v], refs[sib]) if sib else []
                L = len(refs[v])
                rec, oth = got[who]["split"], got[other]["split"]
                c_rec, c_oth = coverage_counts(pos, L, *rec), coverage_counts(pos, L, *oth)
                total += 1
                supports += int(c_rec["both"] > c_oth["both"])
                shown = ", ".join(map(str, pos[:8])) + (" ..." if len(pos) > 8 else "")
                lines += [f"  only {who} recovered {v} ({L} nt); nearest expected variant {sib}, {len(pos)} differing base(s) at {shown}",
                          f"      under {who}'s split ({rec[0]}/{rec[1]}): " + ", ".join(f"{k} {v_}" for k, v_ in c_rec.items()),
                          f"      under the other split ({oth[0]}/{oth[1]}): " + ", ".join(f"{k} {v_}" for k, v_ in c_oth.items())]
                data.append({"dataset": name, "variant": v, "recovered_by": who, "sibling": sib, "positions": pos, "covered_recovering": c_rec, "covered_other": c_oth})
        if not (got["FIGARO"]["set"] ^ got["AmpliTrim"]["set"]):
            lines.append("  the two recovered the same variants")
        lines.append("")
    if not data and not lines:
        raise RuntimeError("no dataset with a reference and a FIGARO comparison was found")
    verdicts = [H.verdict("a variant that only one tool recovers has more of its distinguishing bases read by both reads under the recovering tool's split (a hypothesis, no pass/fail)", None,
                          f"{supports} of {total} variants" if total else "no variant is recovered by only one tool")]
    return {"title": "Q29: which variant does only one tool recover, and which read covers the bases that tell it apart?", "lines": lines, "verdicts": verdicts, "files": [], "data": data}


# ----------------------------------------------------------------------------------------------------------------------------------------------- Q30: reads of one length
def pair_length_shares(reads_dir, max_samples: int = 12, per_sample: int = 2000) -> dict:
    """How many pairs a tool that needs reads of one length could use: the share with exactly the most common lengths, and the share at least that long."""
    from collections import Counter
    from amplitrim.subsample import find_pairs, read_fastq
    pairs = find_pairs(Path(reads_dir))
    chosen = pairs[::max(1, len(pairs) // max_samples)][:max_samples]
    c, n = Counter(), 0
    for _name, r1, r2 in chosen:
        for k, (a, b) in enumerate(zip(read_fastq(r1), read_fastq(r2))):
            if k >= per_sample:
                break
            c[(len(a[1]), len(b[1]))] += 1
            n += 1
    if not n:
        raise ValueError(f"no reads in {reads_dir}")
    (mf, mr), mc = c.most_common(1)[0]
    return {"samples": len(chosen), "pairs": n, "modal_f": mf, "modal_r": mr, "exact": mc / n, "at_least": sum(v for (lf, lr), v in c.items() if lf >= mf and lr >= mr) / n}


def step_figaro_usable_pairs(ctx, H) -> dict:
    sets = [(k, ctx.p(k)) for k in ("mockA_high", "mockA_all", "v1v2", "its") if ctx.p(k)] + [(n, ready) for n, (ready, _r) in H.parse_extra(ctx.args.extra).items()]
    sets = [(n, r) for n, r in sets if _have(r)]
    if not sets:
        raise RuntimeError("no dataset with a folder of reads on hand")
    rows, data, lines, verdicts = [], [], [], []
    for name, reads in sets:
        r = H.guarded(name, "figaro_usable_pairs", lambda: pair_length_shares(reads, per_sample=300 if ctx.quick else 2000), lines, verdicts)
        if r is None:
            continue
        data.append({"name": name, **r})
        rows.append({"name": name, "s": r["samples"], "n": f"{r['pairs']:,}", "len": f"{r['modal_f']}/{r['modal_r']}", "ex": f"{r['exact']:.1%}", "al": f"{r['at_least']:.1%}"})
    if not data:
        raise RuntimeError("no folder of reads could be read: " + "; ".join(lines))
    lines += ["Pairs a tool needing reads of one length could use (the first reads of up to 12 samples):"]
    lines += H.table_lines(rows, [("dataset", "name", 12), ("samples", "s", 7), ("pairs", "n", 8), ("most common F/R", "len", 15), ("exactly that length", "ex", 19), ("at least that long", "al", 18)])
    v1 = next((d for d in data if d["name"] == "v1v2"), None)
    verdicts += [H.verdict("V1-V2: the share of pairs with exactly the most common lengths is within 5 points of the 47% that FIGARO could use (a check of the earlier figure)",
                          None if v1 is None else abs(100 * v1["exact"] - 47) <= 5, "V1-V2 not on hand" if v1 is None else f"{100 * v1['exact']:.1f}% against 47%")]
    return {"title": "Q30: how many pairs can a tool that needs reads of one length use?", "lines": lines, "verdicts": verdicts, "files": [], "data": data}


# ----------------------------------------------------------------------------------------------------------------------------------------------- Q23: ASV lengths
def asv_length_shares(asv_fasta, reference_fasta, tol: int = LENGTH_TOLERANCE) -> dict:
    import mock_check
    refs = mock_check.read_fasta(reference_fasta)
    lengths = sorted({len(s) for _n, s in refs})
    total = within = 0
    off = []
    for name, seq in mock_check.read_asvs(asv_fasta):
        size = mock_check.size_of(name)
        total += size
        if min(abs(len(seq) - L) for L in lengths) <= tol:
            within += size
        else:
            off.append((size, len(seq), name.split(";")[0]))
    return {"reads": total, "within": within / total if total else 0.0, "off_share": 1 - within / total if total else 0.0, "off": sorted(off, reverse=True)[:3]}


def step_asv_length_check(ctx, H) -> dict:
    lines, verdicts, data = ["Share of the ASV reads whose length is within 5 nt of an expected amplicon (a merge at the wrong offset gives an ASV of another length):"], [], []
    seen = False
    for name, _ready, ref, folder in datasets_with_comparison(ctx, H):
        if ref is None:
            continue
        rows = []
        for asv in sorted(folder.glob("asvs_*.fasta")):
            r = asv_length_shares(asv, ref)
            label = asv.stem[len("asvs_"):]
            data.append({"dataset": name, "method": label, **{k: v for k, v in r.items() if k != "off"}})
            rows.append({"m": label, "w": f"{r['within']:.2%}", "o": f"{r['off_share']:.2%}", "big": ", ".join(f"{s}x{L}nt" for s, L, _n in r["off"])})
            if name == "mockA_high" and label.startswith("AmpliTrim__maxEE_2_6"):
                seen = True
                verdicts.append(H.verdict(f"mock A, AmpliTrim's default setting: at most {LENGTH_MAX_OFF_SHARE:.0%} of the ASV reads off every expected length", r["off_share"] <= LENGTH_MAX_OFF_SHARE, f"{r['off_share']:.2%} off"))
        lines += ["", f"--- {name}"] + H.table_lines(rows, [("method", "m", 44), ("within 5 nt", "w", 11), ("off", "o", 8), ("largest off-length ASVs (reads x length)", "big", 40)])
    if not data:
        raise RuntimeError("no mock community with a reference and a comparison was found")
    if not seen:
        verdicts.append(H.verdict(f"mock A, AmpliTrim's default setting: at most {LENGTH_MAX_OFF_SHARE:.0%} of the ASV reads off every expected length", None, "no ASV file for AmpliTrim's default setting on mock A was found"))
    lines += ["", "On mock B about a third of the reads are contaminants whose lengths need not match an expected amplicon, so only mock A carries a criterion."]
    return {"title": "Q23: do the ASVs of the mocks have the lengths of the expected amplicons?", "lines": lines, "verdicts": verdicts, "files": [], "data": data}


# ----------------------------------------------------------------------------------------------------------------------------------------------- Q16: FIGARO's maxEE rule
def mean_ee_at_cut(pairs, trim_f: int, trim_r: int) -> tuple:
    """The mean expected errors of the reads that reach the cut, counted over the first trim_f (forward) and trim_r (reverse) bases."""

    def one(phred, lengths, cut):
        sel = lengths >= cut
        if not sel.any():
            return float("nan")
        return float((10.0 ** (-phred[sel][:, :cut].astype(float) / 10.0)).sum(axis=1).mean())
    return one(pairs._phred1, pairs.len1, trim_f), one(pairs._phred2, pairs.len2, trim_r)


def rule_check_row(pairs, tf: int, tr: int, fig_ee) -> dict:
    ef, er = mean_ee_at_cut(pairs, tf, tr)
    pf, pr = ef + 1.0, er + 1.0
    return {"mean_f": ef, "mean_r": er, "pred_f": pf, "pred_r": pr, "fig_f": fig_ee[0], "fig_r": fig_ee[1], "within": abs(pf - fig_ee[0]) <= RULE_TOLERANCE and abs(pr - fig_ee[1]) <= RULE_TOLERANCE}


def step_figaro_rule_check(ctx, H) -> dict:
    import maxee_frontier
    from amplitrim.subsample import build_subsample
    from amplitrim.trimmer import ReadPairs
    rows, data, skipped = [], [], []
    for name, ready, _ref, folder in datasets_with_comparison(ctx, H):
        fig = next((m for m in maxee_frontier.read_methods(folder) if m["name"].upper().startswith("FIGARO")), None)
        if fig is None:
            continue
        if not _have(ready):
            skipped.append(f"{name}: the reads ({ready}) are not on hand")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            info = build_subsample(Path(ready), Path(tmp), max_pairs=8000 if ctx.quick else 30000)
            r = rule_check_row(ReadPairs.load(info.r1, info.r2), fig["tf"], fig["tr"], fig["ee"])
        data.append({"name": name, "cut": f"{fig['tf']}/{fig['tr']}", **r})
        rows.append({"name": name, "cut": f"{fig['tf']}/{fig['tr']}", "mean": f"{r['mean_f']:.2f} / {r['mean_r']:.2f}", "pred": f"{r['pred_f']:.1f} / {r['pred_r']:.1f}", "fig": f"{r['fig_f']:g} / {r['fig_r']:g}", "ok": "yes" if r["within"] else "no"})
    if not rows:
        raise RuntimeError("no dataset with a FIGARO comparison and its reads was found" + (": " + "; ".join(skipped) if skipped else ""))
    n_ok = sum(1 for d in data if d["within"])
    lines = ["FIGARO's maxEE at its own cut, against 'the mean expected errors of the reads at that cut + 1' (README, Choosing maxEE):"]
    lines += H.table_lines(rows, [("dataset", "name", 12), ("FIGARO's cut F/R", "cut", 16), ("mean EE at the cut", "mean", 18), ("predicted maxEE", "pred", 15), ("FIGARO's maxEE", "fig", 14), (f"within {RULE_TOLERANCE:g}", "ok", 8)])
    lines += [f"left out: {x}" for x in skipped]
    ok = (n_ok >= RULE_MIN_DATASETS) if len(data) >= 6 else None
    verdicts = [H.verdict(f"the rule 'mean expected errors at the cut + 1' reproduces FIGARO's maxEE (within {RULE_TOLERANCE:g} on both reads in at least {RULE_MIN_DATASETS} of 6 datasets)", ok,
                          f"{n_ok} of {len(data)} datasets" + ("" if len(data) >= 6 else " (fewer than six datasets: no verdict)"))]
    return {"title": "Q16: can FIGARO's maxEE rule be reproduced from the reads?", "lines": lines, "verdicts": verdicts, "files": [], "data": data}


# ----------------------------------------------------------------------------------------------------------------------------------------------- the simulations
def spearman(x: list, y: list) -> float:
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for rank, i in enumerate(order):
            r[i] = float(rank)
        return r
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.fmean(rx), statistics.fmean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def score_asvs(asv_fasta, truth_fasta, info: dict, near_id: float = 0.99) -> dict:
    """The ASVs of one DADA2 run against the simulated truth: taxa recovered exactly (by abundance class), the share of ASV reads in exact ASVs, and the number of ASVs."""
    import community_sim
    import mock_check
    s = mock_check.summarise(asv_fasta, truth_fasta, near_id, "vsearch")
    cls = {n: community_sim.abundance_class(info["counts"].get(n.split("|")[0], 0)) for n in s["recovered"] + s["missing"] + s["only_near"]}
    by = {}
    for c in ("abundant", "medium", "rare"):
        present = [n for n, k in cls.items() if k == c]
        by[c] = (sum(1 for n in present if n in set(s["recovered"])), len(present))
    total = max(s["reads"], 1)
    return {"recovered": len(s["recovered"]), "classes": by, "identical": s["categories"]["exact"][1] / total, "off": s["categories"]["none"][1] / total, "asvs": s["asvs"]}


def richness_decision(rows: list, ee: str = "2") -> dict:
    """The excess (estimate minus merged, in points) at the lowest and at the highest richness, for one maxEE."""
    by = {}
    for r in rows:
        if r["ee"] == ee and r.get("excess") is not None:
            by.setdefault(r["n_taxa"], []).append(r["excess"])
    if len(by) < 2:
        return {"low": None, "high": None, "delta": None, "outcome": None}
    levels = sorted(by)
    low, high = levels[0], (500 if 500 in by else levels[-1])
    d = mean(by[high]) - mean(by[low])
    return {"low": low, "high": high, "delta": d, "outcome": "supported" if d >= RICHNESS_SUPPORTED else ("not supported" if d < RICHNESS_NOT_SUPPORTED else "partly supported")}


def step_richness_sim(ctx, H) -> dict:
    import csv
    import community_sim
    import real_data
    from amplitrim.config import parse_ee
    levels = [25, 100] if ctx.quick else [25, 100, 250, 500, 1000]
    seeds = 1 if ctx.quick else 5
    n_samples, pairs = (3, 1500) if ctx.quick else (10, 5000)
    ees = ["1", "2", "5"]
    qual_from = ctx.p("mockA_high")
    rows, lines, verdicts = [], [], []
    for n_taxa in levels:
        for seed in range(1, seeds + 1):
            tag = f"r{n_taxa}s{seed}"

            def one():
                shutil.rmtree(ctx.cur / tag, ignore_errors=True)                      # the step's own scratch folder from an earlier or interrupted run (the generator itself never overwrites)
                info = community_sim.write_community(ctx.cur / tag / "reads", n_taxa, seed, n_samples, pairs, qual_from=qual_from)
                out = []
                for ee in ees:
                    res = H.recommended(ctx, ctx.cur / tag / "reads", f"{tag}e{ee}", max_ee=parse_ee(ee))
                    rec = res.recommendation
                    if rec is None:
                        raise RuntimeError(f"{tag} maxEE {ee}: AmpliTrim made no recommendation")
                    asv, track = ctx.cur / f"asvs_{tag}e{ee}.fasta", ctx.cur / f"track_{tag}e{ee}.tsv"
                    d = real_data.dada2_run(ctx.cur / tag / "reads", rec.trim_f, rec.trim_r, real_data.ee_text(ee), 12, 0, track, asv)
                    sc = score_asvs(asv, info["truth_fasta"], info, ctx.args.near_id)
                    out.append({"n_taxa": n_taxa, "seed": seed, "ee": ee, "tf": rec.trim_f, "tr": rec.trim_r, "est": rec.yield_est, "pass": d["filtered"], "merged": d["merged"],
                                "excess": 100 * (rec.yield_est - d["merged"]), **sc})
                    ctx.log(f"  {tag} maxEE {ee}: {rec.trim_f}/{rec.trim_r}, estimate {rec.yield_est:.1%}, merged {d['merged']:.1%}, recovered {sc['recovered']}")
                return out
            got = H.guarded(tag, "richness_sim", one, lines, verdicts)
            rows += got or []
    if not rows:
        raise RuntimeError("no community was simulated")
    with open(ctx.cur / "richness_sim.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["n_taxa", "seed", "ee", "trunc_f", "trunc_r", "estimate", "pass", "merged", "excess_points", "taxa_recovered", "identical", "asvs"])
        for r in rows:
            w.writerow([r["n_taxa"], r["seed"], r["ee"], r["tf"], r["tr"], f"{r['est']:.4f}", f"{r['pass']:.4f}", f"{r['merged']:.4f}", f"{r['excess']:.2f}", r["recovered"], f"{r['identical']:.4f}", r["asvs"]])
    table = []
    for n_taxa in sorted({r["n_taxa"] for r in rows}):
        for ee in ees:
            sub = [r for r in rows if r["n_taxa"] == n_taxa and r["ee"] == ee]
            if sub:
                table.append({"n": n_taxa, "ee": ee, "k": len(sub), "est": f"{mean([r['est'] for r in sub]):.1%}", "mer": f"{mean([r['merged'] for r in sub]):.1%}", "exc": f"{mean([r['excess'] for r in sub]):+.1f}",
                              "ab": f"{sum(r['classes']['abundant'][0] for r in sub)}/{sum(r['classes']['abundant'][1] for r in sub)}", "me": f"{sum(r['classes']['medium'][0] for r in sub)}/{sum(r['classes']['medium'][1] for r in sub)}",
                              "ra": f"{sum(r['classes']['rare'][0] for r in sub)}/{sum(r['classes']['rare'][1] for r in sub)}", "id": f"{mean([r['identical'] for r in sub]):.2%}"})
    lines += ["Simulated communities (related taxa, log-normal abundances, real quality strings); AmpliTrim's recommendation at each maxEE, then DADA2 at it; truth known.", ""]
    lines += H.table_lines(table, [("taxa", "n", 6), ("maxEE", "ee", 6), ("runs", "k", 5), ("estimate", "est", 9), ("merged", "mer", 8), ("excess", "exc", 7), ("abundant", "ab", 9), ("medium", "me", 9), ("rare", "ra", 9), ("reads in exact ASVs", "id", 19)])
    ranks = []
    for n_taxa in sorted({r["n_taxa"] for r in rows}):
        for seed in sorted({r["seed"] for r in rows if r["n_taxa"] == n_taxa}):
            sub = sorted([r for r in rows if r["n_taxa"] == n_taxa and r["seed"] == seed], key=lambda r: ees.index(r["ee"]))
            if len(sub) == len(ees):
                ranks.append(spearman([r["est"] for r in sub], [r["merged"] for r in sub]))
    d = richness_decision(rows)
    lines += ["", f"rank correlation of the estimate with DADA2's merged share over the three maxEE settings: mean {mean(ranks):.2f} over {len(ranks)} communities"]
    if d["delta"] is not None:
        lines.append(f"excess at {d['high']} taxa minus excess at {d['low']} taxa (maxEE 2): {d['delta']:+.1f} points")
    ok = None if (ctx.quick or d["outcome"] is None) else d["outcome"] != "not supported"
    verdicts.append(H.verdict(f"the excess of the yield estimate grows with richness (supported: at least {RICHNESS_SUPPORTED:g} points more at 500 than at 25 taxa; not supported: under {RICHNESS_NOT_SUPPORTED:g}; between: partly)", ok,
                              (d["outcome"] or "too few levels") + (f" ({d['delta']:+.1f} points)" if d["delta"] is not None else "") + (" - quick run, no verdict" if ctx.quick else "")))
    verdicts.append(H.verdict("the estimate ranks the three maxEE settings as DADA2 does in every simulated community", None if ctx.quick else (bool(ranks) and min(ranks) >= 0.99), f"lowest rank correlation {min(ranks):.2f}" if ranks else "no community scored"))
    return {"title": "Q0-14: does the yield estimate's excess grow with the richness of the community?", "lines": lines, "verdicts": verdicts, "files": ["richness_sim.tsv"], "data": rows}


ARMS = [("clean", {"repeat_len": 0, "chimera_frac": 0.0, "indel_rate": 0.0}),
        ("realistic", {"repeat_len": 10, "chimera_frac": 0.01, "indel_rate": 1e-5}),
        ("stress", {"repeat_len": 32, "chimera_frac": 0.05, "indel_rate": 1e-4})]
OVERLAPS = (12, 16, 20)


def merged_records(path) -> list:
    """(header, sequence) of a merged-reads file, FASTA or FASTQ."""
    path = Path(path)
    if path.suffix == ".fasta":
        recs, name, seq = [], None, []
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
        if name is not None:
            recs.append((name, "".join(seq)))
        return recs
    from amplitrim.subsample import read_fastq
    return [(h[1:].split()[0], s) for h, s, _p, _q in read_fastq(path)]


def wrong_length_share(records: list, lengths: dict) -> dict:
    """Of the merged reads of normal and repeat taxa, the share whose length is not the true template length (chimeras have no single true length and are left out).
    A pair that carries a simulated indel changes length by one base whatever the merge does, so the false-join signal is the share among pairs WITHOUT an indel ('wrong_clean',
    and by kind); the share among pairs with an indel is reported apart."""
    clean, indel = [0, 0], [0, 0]
    by = {"normal": [0, 0], "repeat": [0, 0]}
    for header, seq in records:
        parts = header.split("|")
        if len(parts) < 3 or parts[2] not in by or parts[1] not in lengths:
            continue
        bad = len(seq) != lengths[parts[1]]
        has_indel = len(parts) > 3 and parts[3] == "indel"
        bucket = indel if has_indel else clean
        bucket[0] += bad
        bucket[1] += 1
        if not has_indel:
            by[parts[2]][0] += bad
            by[parts[2]][1] += 1
    n = clean[1] + indel[1]
    share = lambda b: b[0] / b[1] if b[1] else 0.0
    return {"n": n, "wrong": (clean[0] + indel[0]) / n if n else 0.0, "wrong_clean": share(clean), "wrong_indel": share(indel), "n_indel": indel[1], "normal": share(by["normal"]), "repeat": share(by["repeat"])}


def filter_reads(reads_dir, out_dir, trim_f: int, trim_r: int) -> tuple:
    """Truncate the pairs of the first sample to trim_f/trim_r, filter at maxEE 2 and write them: (forward file, reverse file)."""
    from amplitrim.subsample import find_pairs
    from amplitrim.trimmer import ReadPairs
    _name, r1, r2 = find_pairs(Path(reads_dir))[0]
    out1, out2 = Path(out_dir) / "R1.fastq", Path(out_dir) / "R2.fastq"
    ReadPairs.load(r1, r2).filter(trim_f, trim_r, out1, out2, max_ee=2.0)
    return out1, out2


def run_merge(r1, r2, out_dir, min_overlap: int):
    """vsearch's merge at the given minimum overlap, with its default tolerance for mismatches (the setting AmpliTrim estimates the yield with)."""
    from amplitrim.merge_wrapper import run_vsearch
    return run_vsearch(Path(r1), Path(r2), Path(out_dir), "m", min_overlap, None, "tolerant", keep_merged=True)


def overmerge_decision(rows: list) -> dict:
    """The rule for keeping 12 nt: in the realistic arm the wrong-length share (pairs without an indel) at 12 nt is at most 1%, and DADA2 recovers at 12 nt no more than one taxon fewer than at 20 nt."""
    def at(arm, ov, key):
        return mean([r[key] for r in rows if r["arm"] == arm and r["ov"] == ov])
    w12, r12, r20 = at("realistic", 12, "wrong_clean"), at("realistic", 12, "recovered"), at("realistic", 20, "recovered")
    if any(v != v for v in (w12, r12, r20)):
        return {"ok": None, "wrong12": w12, "diff": None}
    return {"ok": w12 <= MISMERGE_MAX_SHARE and (r12 - r20) >= -MISMERGE_MAX_RECOVERY_LOSS, "wrong12": w12, "diff": r12 - r20}


def step_overmerge_sim(ctx, H) -> dict:
    import csv
    import community_sim
    import real_data
    seeds = 1 if ctx.quick else 3
    n_pairs = 2500 if ctx.quick else 12000
    rows, lines, verdicts = [], [], []
    for arm, kw in ARMS:
        for seed in range(1, seeds + 1):
            tag = f"{arm}{seed}"

            def one():
                shutil.rmtree(ctx.cur / tag, ignore_errors=True)                      # the step's own scratch folder from an earlier or interrupted run (the generator itself never overwrites)
                info = community_sim.write_community(ctx.cur / tag / "reads", 25, seed, 1, n_pairs, qual_from=ctx.p("mockA_high"), **kw)
                out = []
                for ov in OVERLAPS:
                    total = info["max_len"] + ov
                    tf = min(info["read_len"], (total + 1) // 2)
                    tr = min(info["read_len"], total - tf)
                    f1, f2 = filter_reads(ctx.cur / tag / "reads", ctx.cur / tag / f"f{ov}", tf, tr)
                    m = run_merge(f1, f2, ctx.cur / tag / f"m{ov}", ov)
                    sc_m = wrong_length_share(merged_records(m.merged_path), info["lengths"]) if m.merged_path else {"n": 0, "wrong": 0.0, "normal": 0.0, "repeat": 0.0}
                    asv, track = ctx.cur / f"asvs_{tag}_{ov}.fasta", ctx.cur / f"track_{tag}_{ov}.tsv"
                    d = real_data.dada2_run(ctx.cur / tag / "reads", tf, tr, "2", ov, 0, track, asv)
                    sc = score_asvs(asv, info["truth_fasta"], info, ctx.args.near_id)
                    out.append({"arm": arm, "seed": seed, "ov": ov, "tf": tf, "tr": tr, "vsearch_merged": m.merge_rate, "n_merged": sc_m["n"], "wrong": sc_m["wrong"], "wrong_clean": sc_m["wrong_clean"],
                                "wrong_indel": sc_m["wrong_indel"], "wrong_repeat": sc_m["repeat"], "wrong_normal": sc_m["normal"], "dada2_merged": d["merged"], "recovered": sc["recovered"],
                                "identical": sc["identical"], "off": sc["off"]})
                    ctx.log(f"  {tag} overlap {ov}: {tf}/{tr}, vsearch merged {m.merge_rate:.1%}, wrong length {sc_m['wrong_clean']:.2%} (pairs with an indel: {sc_m['wrong_indel']:.1%}), DADA2 merged {d['merged']:.1%}, recovered {sc['recovered']}")
                return out
            rows += H.guarded(tag, "overmerge_sim", one, lines, verdicts) or []
    if not rows:
        raise RuntimeError("no arm was simulated")
    with open(ctx.cur / "overmerge_sim.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["arm", "seed", "overlap", "trunc_f", "trunc_r", "vsearch_merged", "wrong_length_all", "wrong_length_no_indel", "wrong_length_with_indel", "wrong_repeat_no_indel", "dada2_merged", "taxa_recovered", "identical", "off"])
        for r in rows:
            w.writerow([r["arm"], r["seed"], r["ov"], r["tf"], r["tr"], f"{r['vsearch_merged']:.4f}", f"{r['wrong']:.4f}", f"{r['wrong_clean']:.4f}", f"{r['wrong_indel']:.4f}", f"{r['wrong_repeat']:.4f}",
                        f"{r['dada2_merged']:.4f}", r["recovered"], f"{r['identical']:.4f}", f"{r['off']:.4f}"])
    table = []
    for arm, _kw in ARMS:
        for ov in OVERLAPS:
            sub = [r for r in rows if r["arm"] == arm and r["ov"] == ov]
            if sub:
                table.append({"arm": arm, "ov": ov, "k": len(sub), "vm": f"{mean([r['vsearch_merged'] for r in sub]):.1%}", "w": f"{mean([r['wrong_clean'] for r in sub]):.2%}", "wi": f"{mean([r['wrong_indel'] for r in sub]):.1%}", "wr": f"{mean([r['wrong_repeat'] for r in sub]):.2%}",
                              "dm": f"{mean([r['dada2_merged'] for r in sub]):.1%}", "rec": f"{mean([r['recovered'] for r in sub]):.1f}", "id": f"{mean([r['identical'] for r in sub]):.2%}", "off": f"{mean([r['off'] for r in sub]):.2%}"})
    lines += ["25 taxa, reads of 250 nt, the truncation set so that the longest amplicon has the overlap shown; maxEE 2; the arms are described in TEST_REGISTER.md (section 9).", ""]
    lines += H.table_lines(table, [("arm", "arm", 10), ("overlap", "ov", 8), ("runs", "k", 5), ("vsearch merged", "vm", 14), ("wrong length, no indel", "w", 22), ("pairs with an indel", "wi", 19), ("of repeat taxa", "wr", 14), ("DADA2 merged", "dm", 12), ("taxa recovered", "rec", 14),
                                   ("reads in exact ASVs", "id", 19), ("reads in no truth ASV", "off", 21)])
    d = overmerge_decision(rows)
    verdicts.append(H.verdict(f"the 12 nt minimum overlap stands (realistic arm: at most {MISMERGE_MAX_SHARE:.0%} of the merged reads without an indel at the wrong length at 12 nt, and at 12 nt DADA2 recovers no more than {MISMERGE_MAX_RECOVERY_LOSS} taxon fewer than at 20 nt)",
                              None if ctx.quick else d["ok"], "no realistic arm scored" if d["diff"] is None else f"wrong length {100 * d['wrong12']:.2f}% at 12 nt (pairs without an indel); taxa recovered 12 minus 20 nt: {d['diff']:+.1f}" + (" - quick run, no verdict" if ctx.quick else "")))
    return {"title": "Q21-27: how many merges are at the wrong length, and does the 12 nt minimum overlap hold?", "lines": lines, "verdicts": verdicts, "files": ["overmerge_sim.tsv"], "data": rows}
