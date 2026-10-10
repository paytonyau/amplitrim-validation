#!/usr/bin/env python
"""Night 4b: the percentile rule through DADA2 (arm H), against the default (A), rule R's relaxed setting (G) and FIGARO (C), unattended.

  python validation\\night4.py --dry-run                     (prints the plan and the commands; runs nothing)
  python validation\\night4.py                               (the run: about 1 to 2 hours; safe to leave; a dataset that fails is noted and the run goes on)
  python validation\\night4.py --only mockB wastewater_V3-V4  (some datasets)
  python validation\\night4.py --fresh                       (ignore the results of an earlier start of this same run)

The question (TEST_REGISTER.md, section 18): the percentile rule (per read, int(q-th percentile of the expected errors at AmpliTrim's cut) + 1, never below 2) left the
default alone on every clean dataset and gave FIGARO's maxEE within 1 at q = 0.80 (section 17). Chosen for the release, before this run: q = 0.75 and never above 5
(DADA2's tutorial example is c(2, 5)). Does that rule (arm H) merge at least as many pairs as rule R's relaxed setting (2,5), without losing accuracy on the mocks?

For each 16S dataset: AmpliTrim once at maxEE 2 (its cut and the amplicon lengths); the rule's maxEE from the reads at that cut: arm H (q = 0.75, cap 5) and, as a
comparison that shows what the cap costs, arm H80 (q = 0.80, no cap: the rule of section 17). Where both give (2,2), the result is arm A by construction and nothing is
run. Otherwise compare_methods.py runs A, G (2,5), H, H80 and FIGARO's own pick through DADA2 on the same pairs, and each arm is classified against FIGARO (or against A
where FIGARO gave no result) and as a difference from G and from A. ITS is not part of this test: nothing is truncated, so there is no cut to compute the rule at.
NIGHT4_PROTOCOL.md is written into the output folder before any analysis, with the predictions.
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent)]
from night2 import DEFAULT_PAIRS, PAIRS, fraction, pct          # noqa: E402
from night3 import Night3                                       # noqa: E402

Q_RULE, CAP = 0.75, 5                                            # the rule chosen for the release before this run: q = 0.75, never above 5
Q_ALT = 0.80                                                     # the rule of section 17 (q = 0.80, no cap), kept as arm H80
PASSING_ACCURACY = ("as accurate", "more accurate")             # 'mixed' and 'less accurate' do not pass P20
ARMS4 = [("A", "AmpliTrim, maxEE 2 (the default; the control)"),
         ("G", "AmpliTrim's best truncation at maxEE (2,5): rule R's relaxed setting, the one to beat"),
         ("H", f"AmpliTrim's best truncation at the percentile rule's maxEE, q = {Q_RULE:g}, never below 2 and never above {CAP} (the rule chosen for the release)"),
         ("H80", f"the same at q = {Q_ALT:g} with no cap (the rule of section 17: it shows what the cap costs)"),
         ("C", "FIGARO's own pick at the longest amplicon length")]
PREDICTIONS4 = [
    "P19. On every dataset that was run, H merges at least as many pairs as G within the margin (H minus G at least -2 points). Offline, H is expected to FAIL this on mock A "
    "((2,3) passes 71.5% of the pairs against 79.0% for (2,5)) and probably on inter-laboratory V4-V6 ((3,3): 59.7% against 62.5%): q = 0.75 relaxes less than (2,5) there. "
    "It is registered because adoption was fixed before the run.",
    "P20. On every mock with expected sequences that was run (mock A, mock B, D6300 V4-V5), H is as accurate as G or more accurate (variants and strains within 1, reads in identical ASVs within 0.1 points; 'mixed' fails).",
    "P21. H merges more pairs than G by more than the margin (2 points) on wastewater (offline: (4,3) passes 75.8% against 67.8%) and on mock B (offline for the uncapped (4,8): 69.5% against 50.9%; "
    "the capped (4,5) is not measured offline).",
    f"P22 (exploratory). Where the uncapped rule gives a maxEE above {CAP} on a read (expected: mock B), H80 merges more pairs than H by more than the margin: that is what the cap costs.",
]
FOLLOW4 = ("What follows: if P19 and P20 hold, the percentile rule (q = 0.75, cap 5) replaces rule R in the release, with its post hoc status stated. P19 is expected to fail on mock A "
           "and probably on V4-V6: then the rule is NOT shown to be at least as good as R, and it is not adopted on the strength of arm H80 or of any other q chosen after this run "
           "(that would be a new post hoc rule). It may still be chosen for the release as a stated trade, conformity with DADA2's guidance and cleaner reads against fewer pairs, "
           "but not as a validated improvement over R. If P20 fails, rule R stays. Rule R is the rule of the release until this has been read. P21 and P22 are descriptions.")


def margin():
    import maxee_arms as ma
    return ma.MARGIN_MERGED


class Night4(Night3):
    TAG = "night4"
    HEAD = "night4.py: the percentile rule (arm H)"
    only = None

    def datasets(self):
        return [d for d in super().datasets() if not self.only or d["name"] in self.only]

    # ------------------------------------------------------------------ one dataset
    def compare_command4(self, ds, out, ees, amps):
        cmd = [self.py, self.val / "compare_methods.py", ds["reads"], "--out", out, "--amplitrim-max-ee", *ees, "--pairs", str(PAIRS.get(ds["name"], DEFAULT_PAIRS))]
        if self.cfg.get("figaro") and amps:
            cmd += ["--figaro-script", self.cfg["figaro"], "--figaro-amplicon", *[str(a) for a in amps]]
        if ds["reference"]:
            cmd += ["--reference", ds["reference"]]
        return cmd

    def skip(self, ds, reason, plan=None):
        res = {"dataset": ds["name"], "kind": ds["kind"], "marker": ds["marker"], "skipped": reason, "plan": plan}
        self.log(f"{ds['name']}: {reason}")
        if not self.dry:
            (self.out / ds["name"]).mkdir(parents=True, exist_ok=True)
            (self.out / ds["name"] / f"{self.TAG}_result.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
        self.results[ds["name"]] = res

    def process(self, ds):
        import maxee_rules as mr
        name, ma = ds["name"], self.ma
        out = self.out / name
        done = out / f"{self.TAG}_result.json"
        if done.is_file() and not self.fresh and not self.dry:
            self.results[name] = json.loads(done.read_text(encoding="utf-8"))
            self.log(f"{name}: finished in an earlier start of this run (kept; --fresh repeats it)")
            return
        if ds["marker"] == "its":
            return self.skip(ds, "ITS: nothing is truncated, so there is no cut to compute the rule at; not part of this test")
        if not ds["reads"].is_dir():
            return self.problem(name, f"the reads folder {ds['reads']} does not exist")
        try:
            if self.dry:
                self.run_fn(self.compare_command4(ds, out / "compare", ["2", ma.ALT_EE, "<the rule's maxEE at q = 0.75, cap 5>", "<the rule's maxEE at q = 0.80, no cap>"], ["<longest amplicon>"]))
                return
            out.mkdir(parents=True, exist_ok=True)
            info = self._explore(ds["reads"], out / "explore", ds["marker"])
            if not info.get("trim"):
                return self.problem(name, "AmpliTrim gave no recommendation at maxEE 2, so there is no cut to compute the rule at")
            pairs = self.build_pairs(ds["reads"])
            tf, tr = info["trim"]
            ef, er = mr.read_errors(pairs, tf, tr)
            if not len(ef) or not len(er):
                return self.problem(name, f"no reads are long enough for AmpliTrim's cut {tf}/{tr}")
            h = {"H": mr.q_rule(ef, er, Q_RULE, cap=CAP)[1], "H80": mr.q_rule(ef, er, Q_ALT)[1]}
            plan = {"cut": [tf, tr], "H": ma.ee_str(h["H"]), "H80": ma.ee_str(h["H80"]), "mean_ee": [float(ef.mean()), float(er.mean())],
                    "p75": [float(np.percentile(ef, 75)), float(np.percentile(er, 75))], "p80": [float(np.percentile(ef, 80)), float(np.percentile(er, 80))]}
            if h["H"] == (2, 2) and h["H80"] == (2, 2):
                return self.skip(ds, f"the rule gives (2,2) at AmpliTrim's cut {tf}/{tr} for both rules: identical to A by construction, nothing to run", plan)
            ees = ["2", ma.ALT_EE]
            for ee in (plan["H"], plan["H80"]):
                if ee not in ees:
                    ees.append(ee)
            amps = [info["longest"]] if info.get("longest") else []
            cdir = out / "compare"
            rc, _text = self.run_fn(self.compare_command4(ds, cdir, ees, amps))
            if rc != 0 or not (cdir / "compare_methods.tsv").is_file():
                return self.problem(name, f"compare_methods.py ended with code {rc} and left no table")
            res = self.analyse4(ds, info, plan, amps, cdir, pairs)
            done.write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
            self.results[name] = res
        except SystemExit as e:
            self.problem(name, f"stopped: {e}")
        except Exception as e:                                                         # nothing here may stop the other datasets
            self.problem(name, f"failed: {type(e).__name__}: {e}")

    def analyse4(self, ds, info, plan, amps, cdir, pairs):
        ma, cm = self.ma, self.cm
        from amplitrim.config import ee_pair, parse_ee
        rows = list(csv.DictReader(open(cdir / "compare_methods.tsv", newline="", encoding="utf-8"), delimiter="\t"))
        ok = [r for r in rows if not r.get("error") and fraction(r.get("merged")) is not None]          # an arm DADA2 could not run has no merged share
        by_name = {r["name"]: r for r in ok}
        nm = lambda ee: f"AmpliTrim (maxEE {ee})"                                                        # noqa: E731
        longest = info.get("longest")
        arm = {"A": by_name.get(nm("2")), "G": by_name.get(nm(ma.ALT_EE)), "H": by_name.get(nm(plan["H"])), "H80": by_name.get(nm(plan["H80"])),
               "C": next((r for r in ok if longest and r["name"] == f"FIGARO (amplicon {longest})"), None)}
        shares = {k: self._shares(cdir / f"track_{cm.method_key(r['name'])}.tsv") for k, r in arm.items() if r}

        def pair_of(r):
            try:
                return ee_pair(parse_ee(str(r["maxEE"])))
            except Exception:
                return None
        noise = {}
        for k, r in arm.items():
            tf, tr, ee = (int(float(r.get("truncF") or 0)), int(float(r.get("truncR") or 0)), pair_of(r)) if r else (0, 0, None)
            noise[k] = pairs.error_load(tf, tr, ee) if (r and tf > 0 and tr > 0 and ee and pairs is not None) else None
        against = "C" if arm["C"] else "A"
        mock = bool(ds["reference"])
        verdicts, vs_G, vs_A = {}, {}, {}
        for k in ("A", "G", "H", "H80"):
            if arm.get(k) and k != against and k in shares and against in shares:
                d = ma.paired_diff(shares[k], shares[against])
                ret, nz = ma.retention_class(d["mean"]), ma.noise_class(noise.get(k), noise.get(against))
                verdicts[k] = {"diff": d, "retention": ret, "noise": nz, "verdict": ma.verdict(ret, nz), "accuracy": ma.accuracy_verdict(arm[k], arm[against]) if mock else "no truth"}
        for base, store in (("G", vs_G), ("A", vs_A)):
            for k in ("H", "H80") if base == "G" else ("G", "H", "H80"):
                if arm.get(k) and arm.get(base) and k in shares and base in shares:
                    d = ma.paired_diff(shares[k], shares[base])
                    store[k] = {"diff": d, "retention": ma.retention_class(d["mean"]), "noise": ma.noise_class(noise.get(k), noise.get(base)),
                                "accuracy": ma.accuracy_verdict(arm[k], arm[base]) if mock else "no truth"}
        h = vs_G.get("H")
        adopt = {"merged": None if not h or h["diff"]["mean"] is None else h["diff"]["mean"] >= -margin(),
                 "accuracy": None if not (mock and h) else h["accuracy"] in PASSING_ACCURACY}
        cap_cost = ma.paired_diff(shares["H80"], shares["H"]) if ("H" in shares and "H80" in shares) else None            # what the cap costs: H80 minus H
        samples = max((len(v) for v in shares.values()), default=0)
        keys = ("name", "truncF", "truncR", "maxEE")
        return {"dataset": ds["name"], "kind": ds["kind"], "marker": ds["marker"], "samples": samples, "against": against, "plan": plan, "amplicons": {"longest": longest},
                "figaro_failed": [f"FIGARO (amplicon {a})" for a in amps if not arm["C"]],
                "arms": {k: {**{x: r.get(x) for x in keys}, "merged": fraction(r.get("merged")), "pass": fraction(r.get("pass")), "found": r.get("found"), "strains": r.get("strains"),
                             "identical": fraction(r.get("identical")), "noise": noise.get(k)} for k, r in arm.items() if r},
                "verdicts": verdicts, "vs_G": vs_G, "vs_A": vs_A, "adopt": adopt, "cap_cost": cap_cost, "errors": [r["error"] for r in rows if r.get("error")]}

    # ------------------------------------------------------------------ the predictions, read from the results
    def outcomes(self, rs):
        run = [r for r in rs if "arms" in r]
        m = margin()
        mean = lambda r, k, base: ((r.get(base) or {}).get(k) or {}).get("diff", {}).get("mean")                 # noqa: E731
        p19 = [(r["dataset"], mean(r, "H", "vs_G")) for r in run if mean(r, "H", "vs_G") is not None]
        o19 = None if not p19 else all(v >= -m for _, v in p19)
        p20 = [(r["dataset"], r["vs_G"]["H"]["accuracy"]) for r in run if r["kind"] == "mock" and (r["vs_G"].get("H") or {}).get("accuracy") not in (None, "no truth")]
        o20 = None if not p20 else all(a in PASSING_ACCURACY for _, a in p20)
        p21 = [(d, mean(r, "H", "vs_G")) for r in run for d in (r["dataset"],) if d in ("mockB", "wastewater_V3-V4") and mean(r, "H", "vs_G") is not None]
        o21 = None if len(p21) < 2 else all(v > m for _, v in p21)
        above = lambda r: max(float(x) for x in r["plan"]["H80"].split(",")) > CAP                                  # noqa: E731
        p22 = [(r["dataset"], r["cap_cost"]["mean"]) for r in run if r.get("cap_cost") and r["cap_cost"]["mean"] is not None and above(r)]
        o22 = None if not p22 else all(v > m for _, v in p22)
        fmt = lambda xs: ", ".join(f"{d} {v:+.1f}" if isinstance(v, float) else f"{d} {v}" for d, v in xs) or "no dataset"          # noqa: E731
        return [("P19", o19, f"H minus G, merged points: {fmt(p19)}"), ("P20", o20, f"H against G, accuracy: {fmt(p20)}"),
                ("P21", o21, f"H minus G, merged points: {fmt(p21)}"), ("P22", o22, f"H80 minus H (the cost of the cap), merged points: {fmt(p22)}")]

    # ------------------------------------------------------------------ the protocol and the report
    def protocol_text(self, ds_list):
        ma = self.ma
        L = [f"# Night 4b protocol (written {time.strftime('%Y-%m-%d %H:%M')}, before any analysis)", "", "Fixed in TEST_REGISTER.md, section 18, before the run.", "",
             "Question: does the percentile rule chosen for the release (q = 0.75, never above 5) merge at least as many pairs as rule R's relaxed setting (2,5), without losing accuracy on the mocks?", "", "| arm | what |", "| --- | --- |"]
        L += [f"| {k} | {v} |" for k, v in ARMS4] + [""]
        L += [f"The rule, exactly: per read, int(q-th percentile of the expected errors of the reads at AmpliTrim's cut at maxEE 2) + 1, never below 2; H: q = {Q_RULE:g} and never above {CAP}; "
              f"H80: q = {Q_ALT:g}, no cap; the search is then repeated at that maxEE. Where both are (2,2), the arm is A and nothing is run. ITS is not tested.", "",
              f"Margins as in night 2: merged share {ma.MARGIN_MERGED:g} points; expected errors {ma.MARGIN_NOISE:g} per read; accuracy {ma.MARGIN_ACCURACY}. FIGARO is run with -m 12 at the longest "
              "amplicon length; a FIGARO failure is 'no result'. Fewer than 10 samples: point estimates only. For accuracy, only 'as accurate' and 'more accurate' pass.", "", "Predictions, written before the run:"]
        L += [f"- {p}" for p in PREDICTIONS4] + ["", FOLLOW4, "", "Datasets:", ""] + [f"- {d['name']} ({d['kind']}, {d['marker'].upper()})" for d in ds_list] + [""]
        return "\n".join(L)

    def table(self, r):
        a = r["arms"]
        L = [f"### {r['dataset']}  ({r['kind']}, {r['marker'].upper()}, {r['samples']} samples; verdicts against {'FIGARO' if r['against'] == 'C' else 'A (FIGARO gave no result)'})", ""]
        p = r["plan"]
        L += [f"The rule at AmpliTrim's cut {p['cut'][0]}/{p['cut'][1]}: mean expected errors {p['mean_ee'][0]:.2f} / {p['mean_ee'][1]:.2f}; 75th percentile {p['p75'][0]:.2f} / {p['p75'][1]:.2f}; "
              f"80th percentile {p['p80'][0]:.2f} / {p['p80'][1]:.2f}; maxEE {p['H']} for H (q = {Q_RULE:g}, cap {CAP}) and {p['H80']} for H80 (q = {Q_ALT:g}, no cap).", ""]
        if r["samples"] and r["samples"] < 10:
            L += ["Fewer than 10 samples: point estimates only, no interval.", ""]

        def dtxt(d):
            if not d or d["mean"] is None:
                return "-"
            return f"{d['mean']:+.1f}" + ("" if d["lo"] is None else f" ({d['lo']:+.1f} to {d['hi']:+.1f})")
        L += ["| arm | cut F/R | maxEE | merged | passing the filter | errors per kept read F / R | vs " + r["against"] + ": merged difference | verdict | vs G: merged difference | accuracy vs G |",
              "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for k in ("A", "G", "H", "H80", "C"):
            if k not in a:
                continue
            x, v, g = a[k], r["verdicts"].get(k), r["vs_G"].get(k)
            nz = "-" if not x["noise"] else f"{x['noise'][0]:.2f} / {x['noise'][1]:.2f}"
            label = "C FIGARO (longest amplicon)" if k == "C" else k
            L.append(f"| {label} | {x['truncF']}/{x['truncR']} | {x['maxEE']} | {pct(x['merged'])} | {pct(x['pass'])} | {nz} | {dtxt(v['diff']) if v else '-'} | {v['verdict'] if v else '-'} | {dtxt(g['diff']) if g else '-'} | {g['accuracy'] if g else '-'} |")
        ad = r["adopt"]
        L += ["", "Adoption check for H: merged pairs at least G's within the margin: " + ("-" if ad["merged"] is None else ("yes" if ad["merged"] else "**NO**")) + "; accuracy not below G's: "
              + ("no truth" if ad["accuracy"] is None else ("yes" if ad["accuracy"] else "**NO**")) + "."]
        cc = r.get("cap_cost")
        if cc and cc["mean"] is not None:
            L += ["", f"What the cap costs (H80 minus H, merged pairs): {cc['mean']:+.1f} points" + ("" if cc["lo"] is None else f" ({cc['lo']:+.1f} to {cc['hi']:+.1f})") + "."]
        if r["figaro_failed"]:
            L += ["", "FIGARO gave no result at: " + ", ".join(r["figaro_failed"]) + " (not a difference in outcome)."]
        if r["kind"] == "mock":
            L += ["", "Accuracy (variants / strains / reads in exact ASVs): " + "; ".join(f"{k}: {x.get('found')} / {x.get('strains')} / {pct(x.get('identical'))}" for k, x in a.items())]
        return L + [""]

    def write_report(self, ds_list):
        L = [f"# Night 4b report ({time.strftime('%Y-%m-%d %H:%M')})", "", f"Problems: {len(self.problems)}" + ("" if not self.problems else "  (listed at the end)"), ""]
        rs = [self.results[d["name"]] for d in ds_list if d["name"] in self.results]
        L += ["## Predictions", ""]
        for p, ok, txt in self.outcomes(rs):
            L.append(f"- **{p}** {'MET' if ok else 'NOT MET' if ok is not None else 'not decided (no dataset to decide it)'}: {txt}")
        L += ["", "## The rule's maxEE and the gain over rule R's setting (merged pairs, points)", "", "| dataset | rule's maxEE (H / H80) | H minus G | H minus A | H accuracy vs G | merged ok | accuracy ok |", "| --- | --- | --- | --- | --- | --- | --- |"]
        for r in rs:
            if "arms" not in r:
                continue
            g, a, ad = r["vs_G"].get("H"), r["vs_A"].get("H"), r["adopt"]
            f = lambda x: "-" if not x or x["diff"]["mean"] is None else f"{x['diff']['mean']:+.1f}"                 # noqa: E731
            yn = lambda v: "-" if v is None else ("yes" if v else "NO")                                              # noqa: E731
            L.append(f"| {r['dataset']} | {r['plan']['H']} / {r['plan']['H80']} | {f(g)} | {f(a)} | {g['accuracy'] if g else '-'} | {yn(ad['merged'])} | {yn(ad['accuracy'])} |")
        sk = [r for r in rs if "skipped" in r]
        if sk:
            L += ["", "## Not run", ""] + [f"- {r['dataset']}: {r['skipped']}" for r in sk]
        L += ["", "## Datasets", ""]
        for r in rs:
            if "arms" in r:
                L += self.table(r)
        if self.problems:
            L += ["## Problems", ""] + [f"- {n}: {m}" for n, m in self.problems] + [""]
        L += ["Next: send this file, NIGHT4_PROTOCOL.md and night4.log."]
        text = "\n".join(L) + "\n"
        if self.dry:
            print("\n(dry run: nothing written)")
        else:
            (self.out / "NIGHT4_REPORT.md").write_text(text, encoding="utf-8")
            self.log(f"report: {self.out / 'NIGHT4_REPORT.md'}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", type=Path, default=HERE.parent)
    ap.add_argument("--config", type=Path, help="another amplitrim_paths.json (default: the one next to the package)")
    ap.add_argument("--only", nargs="+", help="run only these datasets")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--fresh", action="store_true")
    a = ap.parse_args(argv)
    n = Night4(a.repo, a.dry_run, False, a.fresh, cfg_path=a.config)
    n.only = a.only
    return n.go()


if __name__ == "__main__":
    sys.exit(main())
