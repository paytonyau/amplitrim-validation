#!/usr/bin/env python
"""Night 3: arm G (AmpliTrim at maxEE 2,5) and arm F2 (a careful manual cut), unattended, on every dataset registered in amplitrim_paths.json.

  python validation\\night3.py --dry-run         (prints the plan and the commands; runs nothing)
  python validation\\night3.py                   (the run: about 2 hours; safe to leave; a dataset that fails is noted and the run goes on)
  python validation\\night3.py --fresh           (ignore the results of an earlier start of this same run)
  python validation\\night3.py --repo C:\\amplitrim_v2\\amplitrim

The question (TEST_REGISTER.md, section 14): where AmpliTrim's default maxEE of 2 is too strict, does DADA2's own suggestion for the reverse read (maxEE 2,5) recover most of
what FIGARO's looser pick gains, without losing accuracy on the mock communities? And what would a careful manual cut (F2) give?

For each dataset: AmpliTrim once with --explore-max-ee (for the amplicon lengths and its default cut); the F2 cut from the quality table; then compare_methods.py runs
A (maxEE 2), G (maxEE 2,5), F2 and FIGARO's own pick through DADA2 on the same pairs; the noise of each arm is measured from the qualities; each arm is classified
against FIGARO's pick (or against A where FIGARO gave no result), and G and F2 are also given as differences from A. The coverage rows, FIGARO at the typical length and
the stability of the picks were done in night 2 and are not repeated. NIGHT3_PROTOCOL.md is written into the output folder before any analysis, with the predictions.
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent)]
from night2 import DEFAULT_PAIRS, PAIRS, Night2, fraction, pct          # noqa: E402

ARMS3 = [("A", "AmpliTrim, maxEE 2 (the default; the control)"), ("G", "AmpliTrim's best truncation at maxEE (2,5), DADA2's own suggestion for the reverse read"),
         ("R", "a rule, no extra run: G where A passes under 70% of the pairs, otherwise A"),
         ("F2", "a careful manual cut (post hoc): the median-Q30 cut, never longer than the reads, never too short for the amplicon; maxEE 2, relaxed to (2,5) under 60% pass"),
         ("C", "FIGARO's own pick at the longest amplicon length")]
PREDICTIONS3 = [
    "P6. Where A passes under 70% of the pairs (mock A, mock B, inter-laboratory V4-V6, wastewater, D6300 V4-V5), G merges more pairs than A by more than the margin, with noisier reads.",
    "P7. G is within 10 points of FIGARO's pick on mock A, wastewater and inter-laboratory V4-V6, and 15 or more points behind it on mock B (the cap is 5; FIGARO used 6,9).",
    "P8. Where A passes at least 80% of the pairs (V1-V2, soil V4, both V1-V2/V3-V4 mocks, spike-in, swabs), G merges at most 3 points more pairs than A, with a noisier reverse read.",
    "P9. G is as accurate as A, or more accurate, on every mock with expected sequences; on mock B its share of reads in exact ASVs is at least 3 points above A's.",
    "P10. F2 relaxes to (2,5) on mock B and wastewater.",
    "P11. F2 never merges more pairs than the better of A and G by more than the margin on any 16S dataset.",
    "P12. F2 runs on every 16S dataset and merges at least 20 points more pairs than F did in night 2 wherever F merged under 25% or gave no result (a failed F counts as 0%).",
]


class Night3(Night2):
    TAG = "night3"
    HEAD = "night3.py: arms G and F2"

    # ------------------------------------------------------------------ one dataset
    def compare_command3(self, ds, out, f2, amps):
        cmd = [self.py, self.val / "compare_methods.py", ds["reads"], "--out", out, "--amplitrim-max-ee", "2", self.ma.ALT_EE, "--pairs", str(PAIRS.get(ds["name"], DEFAULT_PAIRS))]
        if ds["marker"] == "its":
            cmd += ["--marker", "its"]
        else:
            if f2:
                cmd += ["--setting", f"f2={f2['cut'][0]},{f2['cut'][1]},{f2['ee']}"]
            if self.cfg.get("figaro") and amps:
                cmd += ["--figaro-script", self.cfg["figaro"], "--figaro-amplicon", *[str(a) for a in amps]]
        if ds["reference"]:
            cmd += ["--reference", ds["reference"]]
        return cmd

    def process(self, ds):
        name, ma = ds["name"], self.ma
        out = self.out / name
        done = out / f"{self.TAG}_result.json"
        if done.is_file() and not self.fresh and not self.dry:
            self.results[name] = json.loads(done.read_text(encoding="utf-8"))
            self.log(f"{name}: finished in an earlier start of this run (kept; --fresh repeats it)")
            return
        if not ds["reads"].is_dir():
            return self.problem(name, f"the reads folder {ds['reads']} does not exist")
        try:
            if not self.dry:
                out.mkdir(parents=True, exist_ok=True)
            info = self._explore(ds["reads"], out / "explore", ds["marker"]) if not self.dry else {"picks": {}, "table": [], "trim": None, "longest": None, "typical": None}
            pairs = None if self.dry else self.build_pairs(ds["reads"])
            sixteen = ds["marker"] == "16s"
            f2 = ma.f2_plan(pairs, info.get("longest")) if (pairs is not None and sixteen and info.get("trim")) else None
            amps = [info["longest"]] if (sixteen and info.get("longest")) else []
            cdir = out / "compare"
            rc, _text = self.run_fn(self.compare_command3(ds, cdir, f2, amps))
            if self.dry:
                return
            if rc != 0 or not (cdir / "compare_methods.tsv").is_file():
                return self.problem(name, f"compare_methods.py ended with code {rc} and left no table")
            res = self.analyse3(ds, info, f2, amps, cdir, pairs)
            done.write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
            self.results[name] = res
        except SystemExit as e:
            self.problem(name, f"stopped: {e}")
        except Exception as e:                                                         # nothing here may stop the other datasets
            self.problem(name, f"failed: {type(e).__name__}: {e}")

    def analyse3(self, ds, info, f2, amps, cdir, pairs):
        ma, cm = self.ma, self.cm
        from amplitrim.config import ee_pair, parse_ee
        rows = list(csv.DictReader(open(cdir / "compare_methods.tsv", newline="", encoding="utf-8"), delimiter="\t"))
        ok = [r for r in rows if not r.get("error") and fraction(r.get("merged")) is not None]          # an arm DADA2 could not run has no merged share
        by_name = {r["name"]: r for r in ok}

        def pair_of(r):
            try:
                return ee_pair(parse_ee(str(r["maxEE"])))
            except Exception:
                return None
        longest = info.get("longest")
        arm = {"A": by_name.get("AmpliTrim (maxEE 2)"), "G": by_name.get(f"AmpliTrim (maxEE {ma.ALT_EE})"), "F2": by_name.get("f2"),
               "C": next((r for r in ok if longest and r["name"] == f"FIGARO (amplicon {longest})"), None)}
        shares = {k: self._shares(cdir / f"track_{cm.method_key(r['name'])}.tsv") for k, r in arm.items() if r}
        noise = {}
        for k, r in arm.items():
            tf, tr, ee = (int(float(r.get("truncF") or 0)), int(float(r.get("truncR") or 0)), pair_of(r)) if r else (0, 0, None)
            noise[k] = pairs.error_load(tf, tr, ee) if (r and tf > 0 and tr > 0 and ee and pairs is not None) else None
        against = "C" if arm["C"] else "A"
        verdicts, vs_A = {}, {}
        for k in ("A", "G", "F2"):
            if not arm.get(k) or k == against or k not in shares or against not in shares:
                continue
            d = ma.paired_diff(shares[k], shares[against])
            ret, nz = ma.retention_class(d["mean"]), ma.noise_class(noise.get(k), noise.get(against))
            verdicts[k] = {"diff": d, "retention": ret, "noise": nz, "verdict": ma.verdict(ret, nz), "accuracy": ma.accuracy_verdict(arm[k], arm[against]) if ds["reference"] else "no truth"}
        for k in ("G", "F2"):
            if arm.get(k) and arm.get("A") and k in shares and "A" in shares:
                d = ma.paired_diff(shares[k], shares["A"])
                vs_A[k] = {"diff": d, "retention": ma.retention_class(d["mean"]), "noise": ma.noise_class(noise.get(k), noise.get("A")),
                           "accuracy": ma.accuracy_verdict(arm[k], arm["A"]) if ds["reference"] else "no truth"}
        a_pass = fraction(arm["A"].get("pass")) if arm.get("A") else None
        uses = "G" if (a_pass is not None and a_pass < ma.POOR_BELOW and arm.get("G")) else "A"
        rule = {"uses": uses, "a_pass": a_pass, "against": against, "verdict": (verdicts.get(uses) or {}).get("verdict"), "diff": (verdicts.get(uses) or {}).get("diff")}
        samples = max((len(v) for v in shares.values()), default=0)
        keys = ("name", "truncF", "truncR", "maxEE")
        return {"dataset": ds["name"], "kind": ds["kind"], "marker": ds["marker"], "samples": samples, "against": against, "amplicons": {"longest": longest, "typical": info.get("typical")},
                "f2": f2, "figaro_failed": [f"FIGARO (amplicon {a})" for a in amps if not arm["C"]],
                "arms": {k: {**{x: r.get(x) for x in keys}, "merged": fraction(r.get("merged")), "pass": fraction(r.get("pass")), "found": r.get("found"), "strains": r.get("strains"),
                             "identical": fraction(r.get("identical")), "noise": noise.get(k)} for k, r in arm.items() if r},
                "verdicts": verdicts, "vs_A": vs_A, "rule": rule, "errors": [r["error"] for r in rows if r.get("error")]}

    # ------------------------------------------------------------------ the protocol and the report
    def protocol_text(self, ds_list):
        ma = self.ma
        L = [f"# Night 3 protocol (written {time.strftime('%Y-%m-%d %H:%M')}, before any analysis)", "", "Fixed in TEST_REGISTER.md, section 14, before the run.", "",
             "Question: where AmpliTrim's default maxEE of 2 is too strict, does maxEE (2,5), DADA2's own suggestion for the reverse read, recover most of what FIGARO's looser pick gains, "
             "without losing accuracy on the mocks; and what does a careful manual cut give?", "", "| arm | what |", "| --- | --- |"]
        L += [f"| {k} | {v} |" for k, v in ARMS3] + [""]
        L += ["F2, exactly: (1) each read cut where its median quality first falls below 30; (2) capped at the length that 90% of the reads reach; (3) if forward + reverse is under the longest "
              "amplicon + 12 nt, the missing nucleotides are added half to each read (the forward takes the odd one), up to each cap, and what one read cannot take goes to the other; "
              f"(4) maxEE (2,2), relaxed to ({ma.ALT_EE}) when fewer than {ma.F2_RELAX_BELOW:.0%} of the pairs pass at (2,2). Rule R uses ({ma.ALT_EE}) where A passes under {ma.POOR_BELOW:.0%}.", "",
              f"Margins as in night 2: merged share {ma.MARGIN_MERGED:g} points; expected errors {ma.MARGIN_NOISE:g} per read; accuracy {ma.MARGIN_ACCURACY}. FIGARO is run with -m 12 at the longest amplicon length; "
              "a FIGARO failure is 'no result'. Fewer than 10 samples: point estimates only.", "", "Predictions, written before the run:"]
        L += [f"- {p}" for p in PREDICTIONS3] + ["", "What follows: if P7 and P9 hold, the planned warning offers (2,5) as its alternative; if either fails, it shows the table without naming one. "
                                                  "F2 is a post hoc baseline and changes nothing.", "", "Datasets:", ""] + [f"- {d['name']} ({d['kind']}, {d['marker'].upper()})" for d in ds_list] + [""]
        return "\n".join(L)

    def table(self, r):
        a, f2 = r["arms"], r.get("f2")
        L = [f"### {r['dataset']}  ({r['kind']}, {r['marker'].upper()}, {r['samples']} samples; verdicts against {'FIGARO' if r['against'] == 'C' else 'A (FIGARO gave no result)'})", ""]
        if r["samples"] and r["samples"] < 10:
            L += ["Fewer than 10 samples: point estimates only, no interval.", ""]
        L += ["| arm | cut F/R | maxEE | merged | passing the filter | errors per kept read F / R | vs " + r["against"] + ": merged difference | noise | verdict | accuracy | vs A: merged difference |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for k in ("A", "G", "F2", "C"):
            if k not in a:
                continue
            x, v, w = a[k], r["verdicts"].get(k), r["vs_A"].get(k)
            nz = "-" if not x["noise"] else f"{x['noise'][0]:.2f} / {x['noise'][1]:.2f}"

            def dtxt(d):
                if not d or d["mean"] is None:
                    return "-"
                return f"{d['mean']:+.1f}" + ("" if d["lo"] is None else f" ({d['lo']:+.1f} to {d['hi']:+.1f})")
            label = "C FIGARO (longest amplicon)" if k == "C" else k
            L.append(f"| {label} | {x['truncF']}/{x['truncR']} | {x['maxEE']} | {pct(x['merged'])} | {pct(x['pass'])} | {nz} | {dtxt(v['diff']) if v else '-'} | {v['noise'] if v else '-'} | {v['verdict'] if v else '-'} | {v['accuracy'] if v else '-'} | {dtxt(w['diff']) if w else '-'} |")
        if f2:
            L += ["", f"F2: median-Q30 cut {f2['q30'][0]}/{f2['q30'][1]}, caps {f2['caps'][0]}/{f2['caps'][1]}, amplicon needs {f2['need']}; final cut {f2['cut'][0]}/{f2['cut'][1]}"
                      + ("" if f2["reaches"] else " (DOES NOT REACH THE AMPLICON)") + f"; {pct(f2['pass22'])} pass at (2,2), so maxEE {f2['ee']}."]
        if r["figaro_failed"]:
            L += ["", "FIGARO gave no result at: " + ", ".join(r["figaro_failed"]) + " (not a difference in outcome)."]
        if r["kind"] == "mock":
            L += ["", "Accuracy (variants / strains / reads in exact ASVs): " + "; ".join(f"{k}: {x.get('found')} / {x.get('strains')} / {pct(x.get('identical'))}" for k, x in a.items())]
        return L + [""]

    def write_report(self, ds_list):
        L = [f"# Night 3 report ({time.strftime('%Y-%m-%d %H:%M')})", "", f"Problems: {len(self.problems)}" + ("" if not self.problems else "  (listed at the end)"), ""]
        rs = [self.results[d["name"]] for d in ds_list if d["name"] in self.results]
        L += ["## Scorecard: each arm against FIGARO's pick (datasets where FIGARO gave a result; real and mock apart)", ""]
        for kind in ("mock", "real"):
            sub = [r for r in rs if r["kind"] == kind and r["against"] == "C"]
            if not sub:
                continue
            L += [f"**{kind} communities ({len(sub)} datasets)**", "", "| arm | " + " | ".join(r["dataset"] for r in sub) + " |", "| --- |" + " --- |" * len(sub)]
            for k in ("A", "G", "F2"):
                L.append(f"| {k} | " + " | ".join((r["verdicts"].get(k) or {}).get("verdict", "-") for r in sub) + " |")
            L.append("| R (uses) | " + " | ".join(f"{r['rule']['verdict'] or '-'} ({r['rule']['uses']})" for r in sub) + " |")
            L.append("")
        L += ["## Gain over the default: G and F2 minus A (merged pairs, points)", "", "| dataset | A passes | G minus A | G noise | G accuracy | F2 minus A | F2 maxEE |", "| --- | --- | --- | --- | --- | --- | --- |"]
        for r in rs:
            g, f = r["vs_A"].get("G"), r["vs_A"].get("F2")
            gd = "-" if not g or g["diff"]["mean"] is None else f"{g['diff']['mean']:+.1f}"
            fd = "-" if not f or f["diff"]["mean"] is None else f"{f['diff']['mean']:+.1f}"
            L.append(f"| {r['dataset']} | {pct(r['rule']['a_pass'])} | {gd} | {g['noise'] if g else '-'} | {g['accuracy'] if g else '-'} | {fd} | {(r['f2'] or {}).get('ee', '-')} |")
        L += ["", "## Datasets", ""]
        for r in rs:
            L += self.table(r)
        if self.problems:
            L += ["## Problems", ""] + [f"- {n}: {m}" for n, m in self.problems] + [""]
        L += ["Next: send this file, NIGHT3_PROTOCOL.md and night3.log."]
        text = "\n".join(L) + "\n"
        if self.dry:
            print("\n(dry run: nothing written)")
        else:
            (self.out / "NIGHT3_REPORT.md").write_text(text, encoding="utf-8")
            self.log(f"report: {self.out / 'NIGHT3_REPORT.md'}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", type=Path, default=HERE.parent)
    ap.add_argument("--config", type=Path, help="another amplitrim_paths.json (default: the one next to the package)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--fresh", action="store_true")
    a = ap.parse_args(argv)
    return Night3(a.repo, a.dry_run, False, a.fresh, cfg_path=a.config).go()


if __name__ == "__main__":
    sys.exit(main())
