#!/usr/bin/env python
"""Does AmpliTrim 0.3 predict the share of read pairs that pass DADA2's filter after the trim? (TEST_REGISTER.md, the section of the trim check;
criteria P115 and P116.)

  python validation\\trim_check.py --night C:\\data\\runs\\29_bakker_rerun --results results\\29_bakker_rerun --amplitrim C:\\Users\\Admin\\OneDrive\\Documents\\GitHub\\amplitrim
      --out C:\\data\\runs\\trim_check

For each of the eleven ITS datasets of run 29, AmpliTrim 0.3 (imported from the folder given with --amplitrim, which must hold version 0.3.x) analyses the
thinned reads that DADA2 was given in run 29 (<night>\\<dataset>\\reads) exactly as `amplitrim --its` does: it predicts the share of pairs that pass DADA2's
filter at maxEE 2,2 as they are, and, with the trim that `amplitrim trim` applies, the share that pass after it (pipeline._its_after_trim). DADA2's shares are
those of run 29 (<results>\\<dataset>\\dada2\\std and sw15\\samples.tsv: filtered / input, summed over the samples). Run 29's trim (validation/its_night.py)
and AmpliTrim's (amplitrim/trim.py) are the same code.

P115: on every dataset, the predicted pass rate after the trim is within 1 point of DADA2's (`sw15`).
P116: on every dataset, the predicted pass rate at maxEE 2,2 is within 1 point of DADA2's (`std`), as P92 was for 0.2.0.
Described: the predicted yield after the trim against DADA2's merged share (`sw15`).
Writes TRIM_CHECK_SUMMARY.md and trim_check.tsv to --out.
"""
import argparse
import csv
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

DATASETS = ["atcc_its1", "atcc_its2", "arfken_its1", "arfken_its2", "arfken_mix_its1", "arfken_mix_its2", "its2", "its",
            "bakker_pcr_its1", "bakker_pol_its1", "bakker_pol_its2"]
WITHIN = 1.0


def dada2_share(path: Path) -> dict:
    """Summed over the samples: input, filtered, merged; and the shares filtered / input and merged / input."""
    tot = {"input": 0.0, "filtered": 0.0, "merged": 0.0}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in tot:
                tot[k] += float(r[k])
    n = tot["input"] or float("nan")
    return {**tot, "pass": tot["filtered"] / n, "merged_share": tot["merged"] / n}


def load_amplitrim(folder: Path, want: str = "0.3."):
    """AmpliTrim imported from the repository folder; refuses any version not starting with `want` (0.3. for the trim check, 0.4. for the v0.4 night),
    or a copy installed elsewhere."""
    sys.path.insert(0, str(Path(folder).resolve()))
    for name in [m for m in sys.modules if m == "amplitrim" or m.startswith("amplitrim.")]:
        del sys.modules[name]                                                     # an AmpliTrim already imported (the archive holds an early copy) is set aside
    import amplitrim
    from amplitrim import pipeline
    here = Path(amplitrim.__file__).resolve()
    if not str(here).startswith(str(Path(folder).resolve())) or not amplitrim.__version__.startswith(want):
        raise SystemExit(f"AmpliTrim {want.rstrip('.')} is not the one imported: version {amplitrim.__version__} from {here}")
    return amplitrim.__version__, pipeline


def predict(pipeline, reads: Path, work: Path) -> dict:
    """As `amplitrim --its` does: the analysis at maxEE 2,2, then the prediction of the trim (whatever the share that passes)."""
    s = pipeline.Settings(input_dir=reads, work_dir=work, marker="its")
    first = pipeline._run_once(replace(s, max_ee=2.0))
    if not first.candidates:
        return {"pairs": first.sub.n_pairs, "pass_default": float("nan"), "pass_trim": float("nan"), "yield_trim": float("nan")}
    first.pass_rate_default = first.candidates[0].ee_pass_rate                        # the filter's share: it needs no merging
    t = pipeline._its_after_trim(first, s, work)
    return {"pairs": first.sub.n_pairs, "pass_default": first.pass_rate_default, "pass_trim": t["pass_rate"], "yield_trim": t["yield_est"]}


def verdicts(rows: dict) -> list:
    """rows: {dataset: {"pass_trim", "pass_default", "dada2_sw15", "dada2_std", ...}} -> [(id, ok, text)]."""
    out = []
    for cid, pred, obs, what in (("P115", "pass_trim", "dada2_sw15", "after the trim"), ("P116", "pass_default", "dada2_std", "at maxEE 2,2")):
        diffs = [(lab, (r[pred] - r[obs]) * 100) for lab, r in rows.items()]
        missing = [lab for lab, d in diffs if d != d]
        if missing or len(rows) < len(DATASETS):
            out.append((cid, None, f"no prediction or no DADA2 share for {', '.join(missing) or 'some datasets'}"))
            continue
        worst = max(diffs, key=lambda x: abs(x[1]))
        out.append((cid, all(abs(d) <= WITHIN for _, d in diffs), f"the predicted pass rate {what} is within {WITHIN:g} point of DADA2's on every dataset: "
                    f"{len(diffs)} datasets, the largest {worst[1]:+.2f} points ({worst[0]})"))
    return out


def main(argv=None, pipeline=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--night", type=Path, required=True, help="run 29's folder (the thinned reads)")
    ap.add_argument("--results", type=Path, required=True, help="run 29's results in the archive (DADA2's counts)")
    ap.add_argument("--amplitrim", type=Path, help="the AmpliTrim repository folder holding version 0.3")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    if pipeline is None:
        if a.amplitrim is None:
            ap.error("--amplitrim is required")
        version, pipeline = load_amplitrim(a.amplitrim)
    else:
        version = "stand-in"
    rows = {}
    for lab in DATASETS:
        reads = a.night / lab / "reads"
        std, sw = a.results / lab / "dada2" / "std" / "samples.tsv", a.results / lab / "dada2" / "sw15" / "samples.tsv"
        if not reads.is_dir() or not std.is_file() or not sw.is_file():
            print(f"{lab}: not found ({reads} or the DADA2 counts)")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            p = predict(pipeline, reads, Path(tmp))
        d_std, d_sw = dada2_share(std), dada2_share(sw)
        rows[lab] = {**p, "dada2_std": d_std["pass"], "dada2_sw15": d_sw["pass"], "dada2_sw15_merged": d_sw["merged_share"]}
        print(f"{lab:16s} predicted {p['pass_default']:.2%} -> {p['pass_trim']:.2%} after the trim; DADA2 {d_std['pass']:.2%} -> {d_sw['pass']:.2%}")
    V = verdicts(rows)
    a.out.mkdir(parents=True, exist_ok=True)
    pct = lambda x: "-" if x != x else f"{x:.2%}"
    L = ["# AmpliTrim 0.3: the trim check (summary)", "", f"AmpliTrim {version}; the criteria (P115, P116) are in TEST_REGISTER.md and were fixed before this ran.", "",
         "| Dataset | Pairs analysed | Pass at 2,2: predicted / DADA2 | Pass after the trim: predicted / DADA2 | Yield after the trim: predicted / DADA2 merged |",
         "|---|---|---|---|---|"]
    L += [f"| {lab} | {r['pairs']:,} | {pct(r['pass_default'])} / {pct(r['dada2_std'])} | {pct(r['pass_trim'])} / {pct(r['dada2_sw15'])} | "
          f"{pct(r['yield_trim'])} / {pct(r['dada2_sw15_merged'])} |" for lab, r in rows.items()]
    L += [""] + [f"- **{c}** {'MET' if ok else 'NOT MET' if ok is False else 'not decided'}: {t}" for c, ok, t in V] + ["", "Next: send this file."]
    (a.out / "TRIM_CHECK_SUMMARY.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    cols = ["pairs", "pass_default", "dada2_std", "pass_trim", "dada2_sw15", "yield_trim", "dada2_sw15_merged"]
    (a.out / "trim_check.tsv").write_text("dataset\t" + "\t".join(cols) + "\n" + "".join(lab + "\t" + "\t".join(f"{r[c]:.6g}" for c in cols) + "\n" for lab, r in rows.items()),
                                          encoding="utf-8")
    print("\n".join(L))
    return 0 if all(ok is not None for _, ok, _ in V) else 2


if __name__ == "__main__":
    sys.exit(main())
