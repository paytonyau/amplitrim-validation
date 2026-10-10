#!/usr/bin/env python
"""Night 2: the maxEE head-to-head, unattended, on every dataset registered in amplitrim_paths.json.

  python validation\\night2.py --dry-run         (prints the plan and the commands; runs nothing)
  python validation\\night2.py                   (the run: about 5 to 7 hours; safe to leave; a dataset that fails is noted and the run goes on)
  python validation\\night2.py --no-stability    (skip the stability of the picks, which saves about an hour)
  python validation\\night2.py --fresh           (ignore the results of an earlier start of this same run)
  python validation\\night2.py --repo C:\\amplitrim_v2\\amplitrim

The question: when AmpliTrim chooses the maxEE as well as the truncation, does it end up where FIGARO ends up, with as many merged pairs and no noisier reads?
It needs the build that contains validation\\maxee_arms.py (run python validation\\check_build.py first). The arms, the margins and the verdict words are in that
file's docstring and in TEST_REGISTER.md, section 12; NIGHT2_PROTOCOL.md is written into the output folder before any analysis, with the predictions.

For each dataset: AmpliTrim once with --explore-max-ee (its picks at the floors 50% and 80%, the amplicon lengths); the rule 'mean expected errors at the cut + 1';
the median-Q30 cut; then compare_methods.py runs every arm through DADA2 on the same pairs (FIGARO at the longest amplicon and, when a tail of longer amplicons makes them
differ, at the typical one); the noise of each arm is measured from the qualities; each arm is classified against FIGARO's pick; with a tail, AmpliTrim also runs leaving
out 1% and 5% of the longest amplicons; and, with 9 or more samples, each tool picks again on three disjoint sets of samples to see whether its pick holds.
"""
import argparse
import copy
import csv
import json
import shutil
import sys
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace

PAIRS = {"mockA_high": 0}                     # DADA2 gets all the reads of mock A (six samples); every other dataset 5,000 pairs per sample, as in night 1
DEFAULT_PAIRS = 5000
MOCKS = {"mockA_high", "mockB", "mcm_v12", "mcm_v46", "d6300_V3-V4", "d6300_V4-V5"}
MIN_FREE = 6e9
ARMS = [("A", "AmpliTrim, maxEE 2 (the default)"), ("B", "AmpliTrim at its own quality-first maxEE, floor 50%"), ("B2", "the same, floor 80%"),
        ("B4", "AmpliTrim at 'mean expected errors at the cut + 1' (the rule FIGARO is documented to use)"), ("D", "AmpliTrim's truncation at FIGARO's maxEE"),
        ("E", "no truncation"), ("F", "median-Q30 cut")]
PREDICTIONS = [
    "P1. Where AmpliTrim's default keeps 80% or more of the pairs (V1-V2, soil V4, both inter-laboratory V1-V2 and the NovaSeq V3-V4 sets), the floors are not binding, so B and B2 equal A, "
    "and FIGARO's stricter maxEE gives 'trade-off: more pairs, noisier' for A against FIGARO.",
    "P2. On poor-quality reads (mock B, wastewater, inter-laboratory V4-V6) B merges fewer pairs than FIGARO by more than the margin, with cleaner reads (about 18 points on mock B); "
    "B2 narrows the gap and B4 comes nearer FIGARO than B does. AmpliTrim's grid stops at (5,5), so B and B2 cannot reach FIGARO's (6,9) on mock B.",
    "P3. On mock B the accuracy of B is lower than FIGARO's (fewer variants or fewer reads in exact ASVs).",
    "P4. On D6300 V4-V5 FIGARO runs at the typical amplicon length, chooses much shorter reads than AmpliTrim's default, and leaving out 1% of the longest amplicons raises AmpliTrim's merged share.",
    "P5. The picks of both tools hold across three sets of samples, except where read quality differs strongly between samples (wastewater).",
]


def fraction(x):
    """A value of compare_methods.tsv as a fraction (the table may hold percents)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v / 100.0 if v > 1.5 else v


def pct(x, nd=1):
    return "-" if x is None else f"{100 * x:.{nd}f}%"


def gb(n):
    return f"{n / 1e9:.1f} GB"


class Night2:
    TAG = "night2"                                    # names the output folder, the log, the protocol and the result files (night3.py reuses this class)
    HEAD = "night2.py: the maxEE head-to-head"

    def __init__(self, repo, dry=False, stability=True, fresh=False, run=None, explore=None, pairs=None, shares=None, stable=None, disk_free=None, now=time.strftime, cfg_path=None):
        self.repo, self.val, self.dry, self.do_stability, self.fresh = Path(repo), Path(repo) / "validation", dry, stability, fresh
        sys.path[:0] = [str(self.val), str(self.repo)]
        if not (self.val / "maxee_arms.py").is_file():
            raise SystemExit(f"{self.val / 'maxee_arms.py'} is missing: unzip night2_build.zip over {self.repo} first, then run python validation\\check_build.py")
        import check_build
        import compare_methods
        import maxee_arms
        self.cb, self.cm, self.ma = check_build, compare_methods, maxee_arms
        self.cfg_path = Path(cfg_path) if cfg_path else self.repo / "amplitrim_paths.json"
        self.cfg = check_build.load_config(self.cfg_path)
        self.data = Path(self.cfg["data"])
        self.py = sys.executable
        self._exec = run or self._popen
        self._explore, self._pairs, self._shares, self._stable = explore or maxee_arms.explore, pairs, shares or compare_methods.sample_merged_share, stable
        self.disk_free = disk_free or (lambda p: shutil.disk_usage(p).free)
        self.out = self.data / "runs" / f"{self.TAG}_{now('%Y-%m-%d')}"
        self.report, self.problems, self.results = [], [], {}
        self.logf = None

    # ------------------------------------------------------------------ plumbing
    def log(self, msg=""):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line, flush=True)
        if self.logf:
            self.logf.write(line + "\n")
            self.logf.flush()

    def problem(self, ds, msg):
        self.problems.append((ds, msg))
        self.log(f"PROBLEM ({ds}): {msg}")

    def run_fn(self, cmd):
        self.log("$ " + " ".join(str(c) for c in cmd))
        return (0, "") if self.dry else self._exec(cmd)

    def _popen(self, cmd):
        import subprocess
        p = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", bufsize=1)
        lines = []
        for line in p.stdout:
            lines.append(line.rstrip("\n"))
            print("   " + lines[-1], flush=True)
            if self.logf:
                self.logf.write("   " + lines[-1] + "\n")
        p.wait()
        return p.returncode, "\n".join(lines)

    # ------------------------------------------------------------------ the datasets
    def datasets(self):
        cfg, out = self.cfg, []
        for key, marker in (("mockA_high", "16s"), ("v1v2", "16s"), ("its", "its")):
            if cfg.get(key):
                ref = cfg.get("reference") if key == "mockA_high" else None
                out.append({"name": key, "reads": Path(cfg[key]), "reference": Path(ref) if ref else None, "marker": marker})
        for name, e in (cfg.get("extras") or {}).items():
            if e.get("ready"):
                out.append({"name": name, "reads": Path(e["ready"]), "reference": Path(e["reference"]) if e.get("reference") else None, "marker": "16s"})
        for d in out:
            d["kind"] = "mock" if d["name"] in MOCKS or d["reference"] else "real"
        return out

    def build_pairs(self, reads):
        if self._pairs:
            return self._pairs(reads)
        tmp = tempfile.mkdtemp(prefix="night2_pairs_")
        return self.ma.load_pairs(reads, Path(tmp))

    # ------------------------------------------------------------------ one dataset
    def compare_command(self, ds, out, ees, tf_tr_q30, amps, tail):
        thin = PAIRS.get(ds["name"], DEFAULT_PAIRS)
        cmd = [self.py, self.val / "compare_methods.py", ds["reads"], "--out", out, "--amplitrim-max-ee", *ees, "--pairs", str(thin)]
        if ds["marker"] == "its":
            cmd += ["--marker", "its", "--setting", "trunc200=200,200,2"]
        else:
            cmd += ["--setting", "full=0,0,2", "--setting", f"q30={tf_tr_q30[0]},{tf_tr_q30[1]},2"]
            if self.cfg.get("figaro"):
                cmd += ["--figaro-script", self.cfg["figaro"], "--figaro-amplicon", *[str(a) for a in amps]]
            if tail:
                cmd += ["--allow-unmergeable", "0.01", "0.05"]
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
            rule = ma.rule_ee(pairs, *info["trim"]) if (pairs is not None and ds["marker"] == "16s" and info.get("trim") and all(info["trim"])) else None
            plan = ma.arm_plan(info, rule)
            ees = ma.unique_ees(plan) or ["2"]
            q30 = ma.q30_cut(pairs) if pairs is not None and ds["marker"] == "16s" else (0, 0)
            tail = ds["marker"] == "16s" and ma.has_tail(info.get("longest"), info.get("typical"))
            amps = [a for a in dict.fromkeys([info.get("longest"), info.get("typical") if tail else None]) if a]
            cdir = out / "compare"
            rc, text = self.run_fn(self.compare_command(ds, cdir, ees, q30, amps, tail))
            if self.dry:
                return
            if rc != 0 or not (cdir / "compare_methods.tsv").is_file():
                return self.problem(name, f"compare_methods.py ended with code {rc} and left no table")
            res = self.analyse(ds, info, plan, q30, tail, amps, cdir, pairs)
            res["stability"] = self.stability(ds, info, out) if (self.do_stability and ds["marker"] == "16s") else None
            done.write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
            self.results[name] = res
        except SystemExit as e:
            self.problem(name, f"stopped: {e}")
        except Exception as e:                                                         # nothing here may stop the other datasets
            self.problem(name, f"failed: {type(e).__name__}: {e}")

    def analyse(self, ds, info, plan, q30, tail, amps, cdir, pairs):
        ma, cm = self.ma, self.cm
        from amplitrim.config import ee_pair, parse_ee
        rows = list(csv.DictReader(open(cdir / "compare_methods.tsv", newline="", encoding="utf-8"), delimiter="\t"))
        ok = [r for r in rows if not r.get("error") and fraction(r.get("merged")) is not None]          # an arm DADA2 could not run has no merged share: it is an error, not an arm
        by_name = {r["name"]: r for r in ok}

        def pair_of(r):
            try:
                return ee_pair(parse_ee(str(r["maxEE"])))
            except Exception:
                return None
        amp_rows = [r for r in ok if r["name"].startswith("AmpliTrim (") and "leave out" not in r["name"]]
        fig = {a: next((r for r in ok if r["name"] == f"FIGARO (amplicon {a})"), None) for a in amps}
        c_row = fig.get(info.get("longest")) if info.get("longest") else None
        arm = {"A": by_name.get("AmpliTrim (maxEE 2)")}
        for label in ("B", "B2", "B4"):
            arm[label] = by_name.get(f"AmpliTrim (maxEE {plan[label]})") if plan.get(label) else None
        arm["C"] = c_row
        arm["C_typical"] = fig.get(info.get("typical")) if tail and info.get("typical") != info.get("longest") else None
        arm["D"] = next((r for r in amp_rows if c_row and pair_of(r) == pair_of(c_row)), None)
        arm["E"] = by_name.get("trunc200" if ds["marker"] == "its" else "full")
        arm["F"] = by_name.get("q30")
        shares = {k: self._shares(cdir / f"track_{cm.method_key(r['name'])}.tsv") for k, r in arm.items() if r}
        noise = {}
        for k, r in arm.items():
            tf, tr, ee = int(float(r.get("truncF") or 0)) if r else 0, int(float(r.get("truncR") or 0)) if r else 0, pair_of(r) if r else None
            noise[k] = pairs.error_load(tf, tr, ee) if (r and tf > 0 and tr > 0 and ee and pairs is not None) else None
        against = "C" if c_row else "A"
        verdicts = {}
        for k in ("A", "B", "B2", "B4", "D", "E", "F", "C_typical"):
            if not arm.get(k) or k == against or k not in shares or against not in shares:
                continue
            d = ma.paired_diff(shares[k], shares[against])
            ret, nz = ma.retention_class(d["mean"]), ma.noise_class(noise.get(k), noise.get(against))
            verdicts[k] = {"diff": d, "retention": ret, "noise": nz, "verdict": ma.verdict(ret, nz), "accuracy": ma.accuracy_verdict(arm[k], arm[against]) if ds["reference"] else "no truth"}
        coverage = [{"name": r["name"], "truncF": r.get("truncF"), "truncR": r.get("truncR"), "merged": fraction(r.get("merged")), "pass": fraction(r.get("pass"))} for r in ok if "leave out" in r["name"]]
        if coverage and arm["A"]:
            coverage.insert(0, {"name": arm["A"]["name"] + " (default 0.2%)", "truncF": arm["A"].get("truncF"), "truncR": arm["A"].get("truncR"), "merged": fraction(arm["A"].get("merged")), "pass": fraction(arm["A"].get("pass"))})
        failed = [n for n in (f"FIGARO (amplicon {a})" for a in amps) if n not in by_name]
        samples = max((len(v) for v in shares.values()), default=0)
        return {"dataset": ds["name"], "kind": ds["kind"], "marker": ds["marker"], "samples": samples, "against": against, "tail": tail, "amplicons": {"longest": info.get("longest"), "typical": info.get("typical")},
                "plan": plan, "picks": {str(k): v for k, v in info.get("picks", {}).items()}, "q30": list(q30),
                "arms": {k: {"name": r["name"], "truncF": r.get("truncF"), "truncR": r.get("truncR"), "maxEE": r.get("maxEE"), "merged": fraction(r.get("merged")), "pass": fraction(r.get("pass")),
                             "nonchim": fraction(r.get("nonchim")), "found": r.get("found"), "strains": r.get("strains"), "identical": fraction(r.get("identical")), "noise": noise.get(k)} for k, r in arm.items() if r},
                "verdicts": verdicts, "coverage": coverage, "figaro_failed": failed, "errors": [r["error"] for r in rows if r.get("error")]}

    def stability(self, ds, info, out):
        """Each tool picks again on three disjoint sets of samples. Injectable for tests; the real version needs vsearch and FIGARO."""
        if self._stable:
            return self._stable(ds, info, out)
        from amplitrim.subsample import find_pairs
        names = sorted(s for s, _a, _b in find_pairs(ds["reads"]))
        if len(names) < 9:
            return {"skipped": f"{len(names)} samples (9 are needed)"}
        ma, cm = self.ma, self.cm
        picks = {"AmpliTrim": [], "AmpliTrim B": [], "FIGARO": []}
        for g in range(3):
            sub = cm.subset_reads(ds["reads"], out / "stability" / f"g{g}", ",".join(names[g::3]))
            inf = self._explore(sub, out / "stability" / f"g{g}_work", ds["marker"])
            picks["AmpliTrim"].append({"trim": list(inf["trim"] or ()), "ee": "2"})
            pk = inf["picks"].get(0.5)
            picks["AmpliTrim B"].append({"trim": list((pk and inf["trim"]) or ()), "ee": ma.ee_str(pk["ee"]) if pk else None})
            if self.cfg.get("figaro") and info.get("longest"):
                try:
                    from amplitrim.figaro_wrapper import parse_figaro_recommendation, run_figaro
                    folder, _kept, _lf, _lr = cm.figaro_input(sub, out / "stability" / f"g{g}_figaro", 20000)
                    rec = parse_figaro_recommendation(run_figaro(Path(self.cfg["figaro"]), sys.executable, folder, out / "stability" / f"g{g}_figaro" / "out", info["longest"], 12))
                    picks["FIGARO"].append({"trim": [rec.trim_f, rec.trim_r], "ee": ma.ee_str(rec.max_ee)})
                except Exception as e:
                    picks["FIGARO"].append({"error": str(e)[:80]})
        return {"picks": picks, "stable": {k: self.is_stable(v) for k, v in picks.items() if v}}

    def is_stable(self, picks):
        good = [p for p in picks if p.get("trim")]
        if len(good) < len(picks) or not good:
            return False
        spread = [max(p["trim"][i] for p in good) - min(p["trim"][i] for p in good) for i in (0, 1)]
        return max(spread) <= self.ma.STABILITY_NT and len({p["ee"] for p in good}) == 1

    # ------------------------------------------------------------------ the run and the report
    def protocol_text(self, ds_list):
        ma = self.ma
        L = [f"# Night 2 protocol (written {time.strftime('%Y-%m-%d %H:%M')}, before any analysis)", "",
             "Question: when AmpliTrim chooses the maxEE as well as the truncation, does it end up where FIGARO ends up, with as many merged pairs and no noisier reads?", "",
             "Every arm is a truncation and a maxEE, run through DADA2 on the same pairs with the same 12 nt overlap for both tools. FIGARO is run with -m 12 and given the longest amplicon AmpliTrim measured "
             "(and, when a tail of longer amplicons makes them differ, also the typical one); a FIGARO failure is reported as 'no result'.", "", "| arm | what |", "| --- | --- |"]
        L += [f"| {k} | {v} |" for k, v in ARMS] + ["| C | FIGARO's own pick |", ""]
        L += ["Measures: the paired per-sample difference in the share of pairs DADA2 merges (an interval from 10 samples); the mean expected errors per kept read, forward and reverse, measured from the qualities "
              "at the arm's truncation and maxEE; on mocks the variants and strains recovered and the share of reads in exact ASVs.", "",
              f"Margins, fixed here: merged share {ma.MARGIN_MERGED:g} points; expected errors {ma.MARGIN_NOISE:g} per read; accuracy {ma.MARGIN_ACCURACY} (the battery's margins); stability {ma.STABILITY_NT} nt.",
              "Verdicts against FIGARO's pick: matches / better / worse / trade-off (more pairs but noisier, or fewer pairs but cleaner) / mixed noise / no result. 'Better' needs a gain on one axis and no loss on the other. "
              "On a mock, accuracy is reported beside the verdict and is what decides a disagreement. Where FIGARO has no result the arms are classified against A. No winner is declared on real data unless one pick is better on both axes.", "",
              "AmpliTrim's maxEE grid stops at (5,5), so its picks cannot reach FIGARO's (6,9) on mock B; its quality-first pick is never stricter than 2.", "", "Predictions, written before the run:"]
        L += [f"- {p}" for p in PREDICTIONS] + ["", "Datasets:", ""] + [f"- {d['name']} ({d['kind']}, {d['marker'].upper()})" for d in ds_list] + [""]
        return "\n".join(L)

    def go(self):
        ds_list = self.datasets()
        free = self.disk_free(self.data if self.data.exists() else Path("."))
        self.log(f"{self.HEAD} on {len(ds_list)} datasets; free space {gb(free)}")
        if free < MIN_FREE:
            self.log(f"STOP: less than {gb(MIN_FREE)} free.")
            return 2
        if not self.dry:
            self.out.mkdir(parents=True, exist_ok=True)
            self.logf = open(self.out / f"{self.TAG}.log", "a", encoding="utf-8")
            proto = self.out / f"{self.TAG.upper()}_PROTOCOL.md"
            if not proto.exists():                                                      # a restart keeps the protocol that was written before the first analysis
                proto.write_text(self.protocol_text(ds_list), encoding="utf-8")
        try:
            with self.cb.keep_awake():
                for ds in ds_list:
                    self.log(f"=== {ds['name']}")
                    self.process(ds)
            self.write_report(ds_list)
        finally:                                                                          # the log is always closed (Windows cannot delete or move a file that is still open)
            if self.logf:
                self.logf.close()
                self.logf = None
        return 0

    def table(self, r):
        a = r["arms"]
        L = [f"### {r['dataset']}  ({r['kind']}, {r['marker'].upper()}, {r['samples']} samples; verdicts against {'FIGARO' if r['against'] == 'C' else 'A (FIGARO gave no result)'})", ""]
        if r["samples"] and r["samples"] < 10:
            L += ["Fewer than 10 samples: point estimates only, no interval.", ""]
        L += ["| arm | cut F/R | maxEE | merged | passing the filter | errors per kept read F / R | vs " + r["against"] + ": merged difference | noise | verdict | accuracy |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        for k in ("A", "B", "B2", "B4", "C", "C_typical", "D", "E", "F"):
            if k not in a:
                continue
            x, v = a[k], r["verdicts"].get(k)
            nz = "-" if not x["noise"] else f"{x['noise'][0]:.2f} / {x['noise'][1]:.2f}"
            diff = "-"
            if v and v["diff"]["mean"] is not None:
                ci = "" if v["diff"]["lo"] is None else f" ({v['diff']['lo']:+.1f} to {v['diff']['hi']:+.1f})"
                diff = f"{v['diff']['mean']:+.1f}{ci}"
            label = {"C": "C FIGARO (longest amplicon)", "C_typical": "C FIGARO (typical amplicon)"}.get(k, k)
            L.append(f"| {label} | {x['truncF']}/{x['truncR']} | {x['maxEE']} | {pct(x['merged'])} | {pct(x['pass'])} | {nz} | {diff} | {v['noise'] if v else '-'} | {v['verdict'] if v else '-'} | {v['accuracy'] if v else '-'} |")
        if r["figaro_failed"]:
            L += ["", "FIGARO gave no result at: " + ", ".join(r["figaro_failed"]) + " (not a difference in outcome)."]
        if r["kind"] == "mock":
            L += ["", "Accuracy (variants / strains / reads in exact ASVs): " + "; ".join(f"{k}: {x.get('found')} / {x.get('strains')} / {pct(x.get('identical'))}" for k, x in a.items())]
        return L + [""]

    def write_report(self, ds_list):
        L = [f"# Night 2 report ({time.strftime('%Y-%m-%d %H:%M')})", "", f"Problems: {len(self.problems)}" + ("" if not self.problems else "  (listed at the end)"), ""]
        rs = [self.results[d["name"]] for d in ds_list if d["name"] in self.results]
        L += ["## Scorecard: each AmpliTrim arm against FIGARO's pick (datasets where FIGARO gave a result; real and mock apart)", ""]
        for kind in ("mock", "real"):
            sub = [r for r in rs if r["kind"] == kind and r["against"] == "C"]
            if not sub:
                continue
            L += [f"**{kind} communities ({len(sub)} datasets)**", "", "| arm | " + " | ".join(r["dataset"] for r in sub) + " |", "| --- |" + " --- |" * len(sub)]
            for k in ("A", "B", "B2", "B4", "D", "E", "F"):
                L.append(f"| {k} | " + " | ".join((r["verdicts"].get(k) or {}).get("verdict", "-") for r in sub) + " |")
            L.append("")
        L += ["## Datasets", ""]
        for r in rs:
            L += self.table(r)
        cov = [r for r in rs if r.get("coverage")]
        if cov:
            L += ["## The coverage setting (share of the longest amplicons AmpliTrim may leave unmergeable)", ""]
            for r in cov:
                L += [f"**{r['dataset']}**: longest amplicon {r['amplicons']['longest']}, typical {r['amplicons']['typical']}", "", "| setting | cut F/R | passing the filter | merged |", "| --- | --- | --- | --- |"]
                L += [f"| {c['name']} | {c['truncF']}/{c['truncR']} | {pct(c['pass'])} | {pct(c['merged'])} |" for c in r["coverage"]] + [""]
        stab = [r for r in rs if r.get("stability")]
        if stab:
            L += ["## Stability of the picks across three disjoint sets of samples", "", "| dataset | tool | picks (cut F/R; maxEE) | holds |", "| --- | --- | --- | --- |"]
            for r in stab:
                s = r["stability"]
                if "skipped" in s:
                    L.append(f"| {r['dataset']} | - | {s['skipped']} | - |")
                    continue
                for tool, plist in s["picks"].items():
                    txt = " ; ".join(("error" if p.get("error") else f"{'/'.join(str(i) for i in p['trim'])}; {p['ee']}") for p in plist)
                    L.append(f"| {r['dataset']} | {tool} | {txt} | {'yes' if s['stable'].get(tool) else 'no'} |")
            L.append("")
        if self.problems:
            L += ["## Problems", ""] + [f"- {n}: {m}" for n, m in self.problems] + [""]
        L += ["Next: send this file, NIGHT2_PROTOCOL.md and night2.log."]
        text = "\n".join(L) + "\n"
        if self.dry:
            print("\n(dry run: nothing written)")
        else:
            (self.out / "NIGHT2_REPORT.md").write_text(text, encoding="utf-8")
            self.log(f"report: {self.out / 'NIGHT2_REPORT.md'}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--config", type=Path, help="another amplitrim_paths.json (default: the one next to the package)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-stability", action="store_true")
    ap.add_argument("--fresh", action="store_true")
    a = ap.parse_args(argv)
    return Night2(a.repo, a.dry_run, not a.no_stability, a.fresh, cfg_path=a.config).go()


if __name__ == "__main__":
    sys.exit(main())
