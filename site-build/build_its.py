"""Build the ITS pages of the validation archive's website, and the pages that serve both markers: docs/its/index.html, docs/its/datasets/<dataset>/index.html,
docs/its/trims/, docs/its/settings/, docs/recovery/, docs/criteria/, docs/definitions/, docs/runs/ and docs/datasets/index.html (the catalogue of all datasets).
Every number is read from the files the runs committed (results/..., docs/DATASET_NAMES.tsv, docs/DEFINITIONS.md, TEST_REGISTER.md), at build time; nothing is typed in.
Datasets and settings have readable names (names.py); the archive's own codes appear only inside <code> file paths.
Usage: python build_its.py [path to the validation archive]   (default: the folder above this one)"""
import csv
import json
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import names as N                                                      # noqa: E402
from shell import ARCH, BLOB, ENA, SOFT, esc, page_shell               # noqa: E402

NOTE = "The numbers are read from the files the runs committed to the validation archive (the paths are given on each page); the register fixed every criterion before its run."
ITS_RUNS = [(26, "26_its_night", "its_night_metrics.tsv"), (27, "27_its_rules", "its_rules_metrics.tsv"),
            (29, "29_bakker_rerun", "its_bakker_metrics_min3.tsv"), (31, "31_v04_night", "v04_metrics.tsv")]
CRITERIA_FILES = [(25, "25_its_followup_summary/ITS_FOLLOWUP_SUMMARY.md"), (26, "26_its_night/ITS_NIGHT_SUMMARY.md"), (27, "27_its_rules/ITS_RULES_SUMMARY.md"), (28, "28_its_bakker/ITS_BAKKER_SUMMARY.md"),
                  (29, "29_bakker_rerun/ITS_BAKKER_SUMMARY_min3.md"), (30, "30_trim_check/TRIM_CHECK_SUMMARY.md"), (31, "31_v04_night/V04_SUMMARY.md")]
OLD_PAGE = {"mockA_high": "mock-a", "mockB": "mock-b", "mcm_v12": "interlab-v1v2", "mcm_v46": "interlab-v4v6", "d6300_V3-V4": "d6300-v3v4", "d6300_V4-V5": "d6300-v4v5",
            "v1v2": "v1v2", "v4": "soil-v4", "its": "its", "wastewater_V3-V4": "wastewater", "spikein_V3-V4": "spike-in", "swabs_V3-V4": "swabs"}
TRUTH = ["atcc_its1", "atcc_its2", "arfken_mix_its1", "arfken_mix_its2", "bakker_pcr_its1", "bakker_pol_its1", "bakker_pol_its2", "pauvert_its1"]


# ------------------------------------------------------------------ small helpers
def f(x):
    return None if x in (None, "") else float(x)


def pc(x, d=1):
    v = f(x)
    return "\u2013" if v is None else f"{v * 100:.{d}f}%"


def nn(x, d=1):
    v = f(x)
    return "\u2013" if v is None else f"{v:.{d}f}"


def sg(v, d=1):
    s = f"{v:+.{d}f}"
    return s.replace("+0.0", "0.0").replace("-0.0", "0.0") if d == 1 else s.replace("+0.00", "0.00").replace("-0.00", "0.00")


def points(v, d=1):
    return "\u2013" if v is None else sg(v * 100, d)


def mono(s):
    return f"<code>{esc(s)}</code>"


def cell(text, cls=""):
    return (text, cls)


def table(headers, rows, cls="", num_cols=()):
    h = "".join(f'<th scope="col"{" class=\"n\"" if i in num_cols else ""}>{esc(t)}</th>' for i, t in enumerate(headers))
    body = []
    for r in rows:
        tds = []
        for i, c in enumerate(r):
            text, k = c if isinstance(c, tuple) else (c, "")
            k = (k + (" n" if i in num_cols else "")).strip()
            tds.append(f'<td{f" class=\"{k}\"" if k else ""}>{text}</td>')
        body.append("<tr>" + "".join(tds) + "</tr>")
    return f'<div class="tablewrap"><table{f" class=\"{cls}\"" if cls else ""}><thead><tr>{h}</tr></thead><tbody>' + "".join(body) + "</tbody></table></div>"


def tiles(items):
    return '<ul class="glance">' + "".join(f'<li><span class="lab">{esc(a)}</span><span class="big">{esc(b)}</span><span class="lab">{esc(c)}</span></li>' for a, b, c in items) + "</ul>"


def chip(text, cls):
    return f'<span class="chip {cls}">{esc(text)}</span>'


def figure(svg, caption):
    return f'<figure class="chart">{svg}<figcaption>{esc(caption)}</figcaption></figure>'


def section(title, body, sid):
    return f'<section class="wrap" id="{sid}"><h2>{esc(title)}</h2>\n{body}\n</section>'


def hero(crumb, h1, intro="", extra=""):
    return f'<section class="hero"><div class="wrap">\n<p class="crumb">{crumb}</p>\n<h1>{esc(h1)}</h1>\n{extra}{f"<p>{intro}</p>" if intro else ""}\n</div></section>'


def page(title, desc, body, depth, here):
    return page_shell(title, desc, body, "../" * depth, NOTE, here)


# ------------------------------------------------------------------ the data
class Data:
    def __init__(self, root):
        self.root = pathlib.Path(root)
        self.names = N.load(self.root)
        self.by_code = N.by_code(self.names)
        self.its = [r for r in self.names if N.is_its(r)]
        self.short = {r["archive_code"]: r["short_name"] for r in self.names}
        self.full = {r["archive_code"]: f"{r['id']} {r['full_name']}" for r in self.names}
        self._arm_re = re.compile(r"\b(" + "|".join(sorted(N.SETTINGS, key=len, reverse=True)) + r")\b")
        self._ds_re = re.compile(r"\b(" + "|".join(sorted((c for c in self.short if c not in ("its", "its2", "v1v2", "v4")), key=len, reverse=True)) + r")\b")
        self.metrics = self._metrics()
        self.runs = self._runs()
        self.criteria = self._criteria()
        self.trims = self._trims()
        self.predictions = self._predictions()
        self.trim_check = self._read_tsv("results/30_trim_check/trim_check.tsv")
        self.recovery = self._recovery()
        self.definitions = (self.root / "docs" / "DEFINITIONS.md").read_text(encoding="utf-8")

    def _read_tsv(self, rel):
        with open(self.root / rel, encoding="utf-8", newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def _metrics(self):
        out = []
        for run, folder, fname in ITS_RUNS:
            for r in self._read_tsv(f"results/{folder}/{fname}"):
                out.append({"run": run, "file": f"results/{folder}/{fname}", **r})
        return out

    def _runs(self):
        sys.path.insert(0, str(self.root / "validation"))
        import organise_data as od
        reg = (self.root / "TEST_REGISTER.md").read_text(encoding="utf-8")
        titles = {int(m.group(1)): m.group(2) for m in re.finditer(r"^## (\d+)\. (.+)$", reg, re.M)}
        out = []
        for run, name, sections, scripts, _ in od.RUNS:
            folders = sorted(p.name for p in (self.root / "results").iterdir() if p.is_dir() and re.match(rf"0*{run}_", p.name))
            secs = []
            for part in sections.split(","):
                part = part.strip()
                if re.fullmatch(r"\d+-\d+", part):
                    a, b = map(int, part.split("-"))
                    secs += list(range(a, b + 1))
                elif part.isdigit():
                    secs.append(int(part))
            out.append({"run": run, "name": name, "sections": sections, "scripts": scripts, "folders": folders, "titles": [(s, titles.get(s, "")) for s in secs]})
        return out

    def readable(self, text):
        """The archive's codes in a line of a summary, in words."""
        text = text.replace("**", "").replace("`", "")
        text = self._ds_re.sub(lambda m: self.short[m.group(1)], text)
        text = re.sub(r"\bits2\b", self.short["its2"], text)
        text = re.sub(r"\b(?:[Tt]he )?second ITS stud(?:y|ies)\b", lambda m: ("The " if m.group(0)[0] == "T" else "") + self.short["its2"] if "stud" in m.group(0) else m.group(0), text)
        text = re.sub(r"\b([Tt]he )?first ITS study\b", lambda m: (m.group(1) or "") + self.short["its"], text)
        text = re.sub(r"\bits(?= [+\-\u2212]\d)", self.short["its"], text)
        return self._arm_re.sub(lambda m: N.setting(m.group(1)), text)

    def _criteria(self):
        pat = re.compile(r"^- \*\*(P\d+)(?: ([^*]+?))?\*\* (MET|NOT MET|not decided)(?: \(([^)]*)\))?:? ?(.*)$")
        tab = re.compile(r"^\| (P\d+) \| (MET|NOT MET|not decided):? (.*) \|$")             # run 25's summary lists its criteria in a table
        out = []
        for run, rel in CRITERIA_FILES:
            for line in (self.root / "results" / rel).read_text(encoding="utf-8").splitlines():
                m = pat.match(line)
                t = tab.match(line)
                if t and not m:
                    out.append({"run": run, "id": t.group(1), "arm": "", "outcome": t.group(2), "note": "", "text": self.readable(t.group(3)), "file": f"results/{rel}"})
                    continue
                if m:
                    arm = m.group(2)
                    out.append({"run": run, "id": m.group(1), "arm": N.setting(arm) if arm in N.SETTINGS else (arm or ""), "outcome": m.group(3), "note": m.group(4) or "",
                                "text": self.readable(m.group(5)), "file": f"results/{rel}"})
        return out

    def _trims(self):
        text = (self.root / "results" / "31_v04_night" / "V04_SUMMARY.md").read_text(encoding="utf-8")
        out = []
        for line in text.splitlines():
            m = re.match(r"^\| (\w+)( \(best\))? \| (True|False) \| ([\d.]+)% \| ([\d.]+)% \| ([\d.]+) \| ([\d.]+)% \|$", line)
            if m:
                out.append({"code": m.group(1), "best": bool(m.group(2)), "safe": m.group(3) == "True", "yield": float(m.group(4)), "atcc2": float(m.group(5)),
                            "cut": float(m.group(6)), "merge": float(m.group(7))})
        return out

    def _predictions(self):
        flow = {(r["dataset"], r["arm"]): float(r["pass"]) for r in self.metrics if r["run"] == 31}
        arm_of = lambda kind, v: {"truncQ": "std" if v == "2" else f"tq{v}", "window": "sw15" if v == "4" else f"w{v}", "quality": f"q{v}", "minLen": f"ml{v}"}[kind]
        out = []
        for r in self.its:
            p = self.root / "results" / "31_v04_night" / r["archive_code"] / "predict.json"
            if not p.is_file():
                continue
            d = json.loads(p.read_text(encoding="utf-8"))
            for kind in ("truncQ", "window", "quality", "minLen"):
                for v, pred in d.get(kind, {}).items():
                    a = arm_of(kind, v)
                    if (r["archive_code"], a) in flow and pred is not None:
                        out.append({"dataset": r["archive_code"], "setting": a, "predicted": float(pred), "dada2": flow[(r["archive_code"], a)]})
        return out

    def _recovery(self):
        text = (self.root / "results" / "32_recovery" / "RECOVERY.md").read_text(encoding="utf-8")
        rows = []
        for line in text.splitlines():
            if re.match(r"^\| (ITS|16S) \|", line):
                c = [x.strip() for x in line.strip("|").split("|")]
                rows.append({"marker": c[0], "code": c[1], "arm": c[2], "samples": c[3], "expected": c[4], "per_sample": c[5], "any_read": c[6], "any_sample": c[7], "all_samples": c[8]})
        pub = re.search(r"^Published, for comparison: (.+)$", text, re.M)
        var = (re.search(r"^16S exact variants recovered \(pooled\): (.+)$", text, re.M) or [None, ""])[1]
        return {"rows": rows, "published": self.readable(pub.group(1)) if pub else "", "variants": self.readable(var)}

    # --- lookups
    def rows(self, code, run=None, arm=None):
        return [r for r in self.metrics if r["dataset"] == code and (run is None or r["run"] == run) and (arm is None or r["arm"] == arm)]

    def one(self, code, arm, run=31):
        r = self.rows(code, run, arm)
        return r[0] if r else None


def d_link(D, code, depth_to_root):
    """A link to the page of an ITS dataset."""
    r = D.by_code[code]
    return f'<a href="{"../" * depth_to_root}its/datasets/{N.slug(r)}/">{esc(r["short_name"])}</a>'


# ------------------------------------------------------------------ charts
def dumbbell(D, codes):
    """Reads kept with the DADA2 default (grey) and after the trim (blue), per dataset."""
    rows = [(c, D.one(c, "std"), D.one(c, "sw15")) for c in codes]
    h = 70 + 26 * len(rows)
    x = lambda v: 250 + v * 250
    b = ['<line class="ink" x1="250" y1="46" x2="500" y2="46"/>']
    for t in (0, 0.25, 0.5, 0.75, 1):
        b.append(f'<line class="ink" x1="{x(t):.0f}" y1="42" x2="{x(t):.0f}" y2="50"/><text class="tm" x="{x(t):.0f}" y="38" text-anchor="middle">{int(t * 100)}%</text>')
    for i, (c, s, t) in enumerate(rows):
        y = 70 + 26 * i
        a, z = float(s["yield"]), float(t["yield"])
        b.append(f'<text class="tm" x="244" y="{y + 4}" text-anchor="end">{esc(D.short[c])}</text><line class="ink" x1="{x(a):.1f}" y1="{y}" x2="{x(z):.1f}" y2="{y}"/>'
                 f'<circle class="dflt" data-v="{a:.4f}" cx="{x(a):.1f}" cy="{y}" r="5"/><circle class="acc" data-v="{z:.4f}" cx="{x(z):.1f}" cy="{y}" r="6"/>')
    label = "Reads kept after chimera removal: DADA2 default and after the trim, per ITS dataset"
    return f'<svg viewBox="0 0 520 {h}" role="img" aria-label="{esc(label)}"><title>{esc(label)}</title>{"".join(b)}</svg>'


TRIM_LABEL = {"sw15": ("Window 4, Quality 15 (best)", -12, 4, "end"), "w3": ("Window 3", 10, -2, "start"), "q18": ("Quality 18", 9, 14, "start"), "w5": ("Window 5", 9, 14, "start"),
              "w2": ("Window 2", 0, 22, "middle"), "w6": ("Window 6", -9, 4, "end"), "w8": ("Window 8", -9, -4, "end"), "q12": ("Quality 12", 9, 10, "start"),
              "w10": ("Window 10", 9, 4, "start"), "w12": ("Window 12", 9, 4, "start"), "w16": ("Window 16", 9, 4, "start"), "w20": ("Window 20", 9, 4, "start")}


def trims_chart(D):
    xs = lambda v: 70 + v / 45 * 410
    ys = lambda v: 250 - (v - 70) / 20 * 200
    b = []
    for t in (70, 75, 80, 85, 90):
        b.append(f'<line class="ink" x1="70" y1="{ys(t):.0f}" x2="480" y2="{ys(t):.0f}"/><text class="tm" x="64" y="{ys(t) + 4:.0f}" text-anchor="end">{t}%</text>')
    for t in (0, 10, 20, 30, 40):
        b.append(f'<text class="tm" x="{xs(t):.0f}" y="272" text-anchor="middle">{t}</text>')
    for t in D.trims:
        cx, cy = xs(t["cut"]), ys(t["yield"])
        lab, dx, dy, an = TRIM_LABEL[t["code"]]
        cls = "acc" if t["safe"] else "open"
        b.append(f'<circle class="{cls}" data-v="{t["yield"]}" cx="{cx:.1f}" cy="{cy:.1f}" r="{7 if t["best"] else 5}"/>'
                 f'<text class="tm" x="{cx + dx:.1f}" y="{cy + dy:.1f}" text-anchor="{an}">{esc(lab)}</text>')
    b.append('<text class="tm" x="275" y="292" text-anchor="middle">Mean bases cut per read</text>')
    label = "Twelve trims: mean reads kept against mean bases cut per read; filled points are safe on every dataset with a known answer"
    return f'<svg viewBox="0 0 520 300" role="img" aria-label="{esc(label)}"><title>{esc(label)}</title>{"".join(b)}</svg>'


# ------------------------------------------------------------------ pages
def arms_table(D, code):
    rows = []
    for r in sorted(D.rows(code), key=lambda r: (r["run"], N.setting_rank(r["arm"]))):
        truth = f(r.get("recall")) is not None
        rows.append([f'<a href="../../../runs/#run-{r["run"]}">{r["run"]}</a>', cell(esc(N.setting(r["arm"])), "name"), f"{int(float(r['input'])):,}", pc(r["pass"]), pc(r["merged"]), pc(r["yield"]),
                     nn(r["asv_median"], 0), pc(r.get("spurious"), 2) if truth else "\u2013", pc(r.get("recall")) if truth else "\u2013", nn(r.get("found")) if truth else "\u2013"])
    return table(["Run", "Setting", "Input pairs", "Pass the filter", "Merged", "Reads kept", "ASVs (median)", "False reads", "Recall", "Members found"], rows, num_cols=range(2, 10))


def its_dataset_page(D, row):
    code, truth = row["archive_code"], row["archive_code"] in TRUTH
    s, t = D.one(code, "std"), D.one(code, "sw15")
    gain = (float(t["yield"]) - float(s["yield"])) * 100
    tl = [("Pass the filter at maxEE 2,2", pc(s["pass"]), "of the read pairs pass DADA2's default filter"),
          ("Reads kept", f"{pc(s['yield'])} \u2192 {pc(t['yield'])}", f"DADA2 default \u2192 after the trim ({sg(gain)} points)")]
    if truth:
        tl.append(("Members found", f"{pc(s['recall'])} \u2192 {pc(t['recall'])}", f"DADA2 default \u2192 after the trim ({points(float(t['recall']) - float(s['recall']))} points)"))
        tl.append(("False reads", f"{pc(s['spurious'], 2)} \u2192 {pc(t['spurious'], 2)}", "DADA2 default \u2192 after the trim"))
    else:
        if row["kind"].startswith("Mock"):
            tl.append(("Known answer", "Isolates", "the isolates' own sequences are the reference for the mixes (DS18, DS19); members found are counted on the mixes"))
        else:
            tl.append(("Known answer", "None", "a real study: the members found cannot be counted"))
        tl.append(("ASVs per sample (median)", f"{nn(s['asv_median'], 0)} \u2192 {nn(t['asv_median'], 0)}", "DADA2 default \u2192 after the trim"))
    facts = (f'<dl class="facts"><dt>Marker, platform</dt><dd>{esc(row["marker_region"])}, {esc(row["platform"])}</dd><dt>Kind</dt><dd>{esc(row["kind"])}</dd>'
             f'<dt>Data</dt><dd>{esc(row["accession"])}</dd><dt>Source</dt><dd>{esc(row["source"])}</dd><dt>In the archive</dt><dd>{mono(row["data_folder"])}</dd></dl>')
    old = OLD_PAGE.get(code)
    more = f'<p>The first validation (the 16S-era tests of the paper) also has a page on this dataset: <a href="../../../datasets/{old}/">{esc(row["short_name"])}</a>.</p>' if old else ""
    rec = [r for r in D.recovery["rows"] if r["code"] == code]
    recov = ""
    if rec:
        recov = section("Recovery of the known members", table(["Setting", "Samples", "Expected per sample", "Found per sample (at least 0.1% of the reads / any read)", "In any sample", "In all samples"],
                        [[cell(esc(N.setting(r["arm"])), "name"), r["samples"], r["expected"], esc(r["per_sample"] + " / " + r["any_read"]), esc(r["any_sample"]), esc(r["all_samples"])] for r in rec],
                        num_cols=(1, 2)) + '<p>What the words mean: <a href="../../../definitions/">definitions</a>. Found is a match of sequences, not a taxonomic assignment.</p>', "recovery")
    pred = [p for p in D.predictions if p["dataset"] == code]
    tc = [r for r in D.trim_check if r["dataset"] == code]
    pred_html = ""
    if pred:
        pred_html = section("AmpliTrim's predicted pass rate against DADA2's", table(["Setting", "AmpliTrim predicted", "DADA2 measured", "Difference (points)"],
                            [[cell(esc(N.setting(p["setting"])), "name"), pc(p["predicted"], 2), pc(p["dada2"], 2), points(p["predicted"] - p["dada2"], 2)] for p in sorted(pred, key=lambda p: N.setting_rank(p["setting"]))],
                            num_cols=(1, 2, 3)) + (f'<p>At maxEE 2,2 without the trim, AmpliTrim predicted {pc(tc[0]["pass_default"], 2)} and DADA2 measured {pc(tc[0]["dada2_std"], 2)}; after the trim, '
                                                  f'{pc(tc[0]["pass_trim"], 2)} against {pc(tc[0]["dada2_sw15"], 2)} (run 30).</p>' if tc else ""), "predictions")
    srcs = sorted({r["file"] for r in D.rows(code)})
    body = (hero('<a href="../../">ITS</a> &rsaquo; Dataset', f"{row['id']} {row['full_name']}", "", f'<p class="chips">{chip(row["kind"], "measured")}{chip(row["marker_region"], "set")}</p>{facts}{tiles(tl)}')
            + section("Every setting tried", f'<p class="legend">One row for each run and setting. False reads and members found need a known answer; a dash means there is none.</p>' + arms_table(D, code), "settings")
            + recov + pred_html
            + section("Where the numbers come from", "<p>" + ", ".join(code_(p) for p in srcs) + f", and <code>results/31_v04_night/{code}/predict.json</code> in the <a href=\"{BLOB}results\">results</a> of the validation archive. {more[3:-4] if more else ''}</p>", "sources"))
    return page(f"{row['full_name']}: ITS results of the AmpliTrim validation", f"Every ITS setting tried on {row['full_name']}.", body, 3, "its/")


def code_(s):
    return mono(s)


def its_overview(D):
    truth = [c for c in TRUTH]
    gains = [(c, (float(D.one(c, "sw15")["yield"]) - float(D.one(c, "std")["yield"])) * 100) for c in truth]
    dr = [(c, (float(D.one(c, "sw15")["recall"]) - float(D.one(c, "std")["recall"])) * 100) for c in truth]
    sp = [(c, (float(D.one(c, "sw15")["spurious"]) - float(D.one(c, "std")["spurious"])) * 100) for c in truth]
    best = next(t for t in D.trims if t["best"])
    tl = [("Reads kept", f"{sg(min(g for _, g in gains))} to {sg(max(g for _, g in gains))} points", f"more read pairs kept after the trim, on all {len(truth)} datasets with a known answer"),
          ("Largest loss of members found", f"{abs(min(d for _, d in dr)):.1f} points", "after the trim, on the dataset where it lost most"),
          ("Largest rise of false reads", f"{max(sp_ for _, sp_ in sp):+.2f} points".replace("+0.00", "0.00"), "after the trim (the limit set beforehand was 1 point)"),
          ("The best of twelve trims", N.setting(best["code"]).replace("Trim: ", ""), f"mean reads kept {best['yield']:.1f}% on the {len(truth)} datasets with a known answer")]
    rows = []
    for r in D.its:
        c = r["archive_code"]
        s, t = D.one(c, "std"), D.one(c, "sw15")
        tr = c in TRUTH
        rows.append([cell(f'<a href="datasets/{N.slug(r)}/">{esc(r["id"] + " " + r["short_name"])}</a>', "name"), esc(r["kind"]), pc(s["pass"]), pc(s["yield"]), pc(t["yield"]),
                     points(float(t["yield"]) - float(s["yield"])), pc(s["recall"]) if tr else "\u2013", pc(t["recall"]) if tr else "\u2013",
                     pc(s["spurious"], 2) if tr else "\u2013", pc(t["spurious"], 2) if tr else "\u2013"])
    body = (hero('<a href="../">Results</a> &rsaquo; ITS', "ITS: cutting poor read ends",
                 "Fungal reads often have poor ends. This page shows what happened on every ITS dataset when the poor read ends were cut before DADA2 (window 4, mean quality 15), against DADA2's default settings. "
                 "Datasets with a known answer are mock communities of fungi.", tiles(tl))
            + section("The default against the trim, dataset by dataset", figure(dumbbell(D, [r["archive_code"] for r in D.its]),
                      "Grey: reads kept with DADA2's default settings. Blue: reads kept after the trim. The line between them is the gain. Reads kept is the share of the input pairs left after DADA2 has filtered, denoised, merged and removed chimeras.")
                      + table(["Dataset", "Kind", "Pass at maxEE 2,2", "Reads kept: default", "Reads kept: trim", "Gain (points)", "Members found: default", "Members found: trim", "False reads: default", "False reads: trim"], rows,
                              num_cols=range(2, 10)) + '<p class="legend">Members found and false reads need a known answer, so real studies and the isolate datasets show a dash. Words: <a href="../definitions/">definitions</a>.</p>', "glance")
            + section("More on the ITS runs", '<ul><li><a href="trims/">The twelve trims compared</a></li><li><a href="settings/">The other settings: minLen, truncQ, forward reads only, and the predictions</a></li>'
                      '<li><a href="../recovery/">Recovery of the known members</a></li><li><a href="../criteria/">Every criterion and its outcome</a></li><li><a href="../runs/">Every run</a></li></ul>', "more"))
    return page("ITS: cutting poor read ends: AmpliTrim validation results", "What cutting the poor ends of ITS reads did on every ITS dataset of the AmpliTrim validation.", body, 1, "its/")


def trims_page(D):
    rows = [[cell(esc(N.setting(t["code"])) + (" (best)" if t["best"] else ""), "name"), ("Yes", "ok") if t["safe"] else ("No", "bad"), f"{t['yield']:.2f}%", f"{t['atcc2']:.2f}%", f"{t['cut']:.1f}", f"{t['merge']:.2f}%",
             sg(t["yield"] - next(x for x in D.trims if x["code"] == "sw15")["yield"], 2)] for t in D.trims]
    crit = [c for c in D.criteria if c["run"] == 31 and c["id"] in ("P119", "P120", "P130")]
    body = (hero('<a href="../../">Results</a> &rsaquo; <a href="../">ITS</a> &rsaquo; Trims', "The twelve trims compared",
                 "A trim cuts each read before the first window of bases whose mean quality is below a limit. Twelve were compared on the eight datasets with a known answer: windows of 2 to 20 bases at mean quality 15, and mean quality 12 and 18 at window 4.")
            + section("What each trim kept and cost", figure(trims_chart(D), "Mean reads kept against mean bases cut per read. Filled points: safe on every dataset with a known answer. Open points: not safe. Larger point: the best by the rule fixed beforehand.")
                      + table(["Trim", "Safe on every dataset", "Mean reads kept", "Reads kept: ATCC Standard ITS2", "Mean bases cut per read", "Mean merge share", "Against window 4, quality 15 (points)"], rows, num_cols=range(2, 7))
                      + '<p class="legend">Safe: false reads rose by at most 1 point and members found fell by at most 2 points, against the default, on every dataset with a known answer. The best: the highest mean reads kept among the safe trims; trims within 0.5 points are tied, and the one that cuts the fewest bases wins.</p>', "trims")
            + section("The criteria", table(["Criterion", "Trim", "Outcome", "What was measured"], [[c["id"], cell(esc(c["arm"]), "name"), (c["outcome"], "ok" if c["outcome"] == "MET" else "bad"), cell(esc(c["text"]), "name")] for c in crit]) +
                      '<p>All criteria: <a href="../../criteria/">Criteria</a>.</p>', "criteria")
            + section("Where the numbers come from", f"<p>{mono('results/31_v04_night/V04_SUMMARY.md')} (register sections 43 and 44).</p>", "sources"))
    return page("The twelve trims compared: AmpliTrim validation results", "Twelve ways of cutting poor read ends in ITS reads, compared on the datasets with a known answer.", body, 2, "its/")


def settings_page(D):
    truth = TRUTH
    groups = [("minLen after the trim", ["ml20", "ml100", "ml150"], "sw15"), ("truncQ", ["tq3", "tq5", "tq7"], "std"), ("Forward reads only", ["fwd", "fwd15"], "std")]
    out = []
    for title, arms, base in groups:
        rows = []
        for a in arms:
            for c in truth:
                x, y = D.one(c, a), D.one(c, base)
                if not x or not y:
                    continue
                rows.append([cell(esc(N.setting(a)), "name"), cell(esc(D.short[c]), "name"), points(float(x["yield"]) - float(y["yield"])), points(float(x["recall"]) - float(y["recall"])),
                             points(float(x["spurious"]) - float(y["spurious"]), 2)])
        crit = [c for c in D.criteria if c["run"] == 31 and any(c["arm"] == N.setting(a) for a in arms)]
        out.append(section(title, f'<p class="legend">Change in points against {esc(N.setting(base))}: reads kept, members found and false reads, on each dataset with a known answer.</p>'
                           + table(["Setting", "Dataset", "Reads kept", "Members found", "False reads"], rows, num_cols=(2, 3, 4))
                           + (table(["Criterion", "Setting", "Outcome"], [[c["id"], cell(esc(c["arm"]), "name"), (c["outcome"], "ok" if c["outcome"] == "MET" else "bad")] for c in crit]) if crit else ""), "g" + str(len(out))))
    rows = [[cell(esc(D.short[p["dataset"]]), "name"), cell(esc(N.setting(p["setting"])), "name"), pc(p["predicted"], 2), pc(p["dada2"], 2), points(p["predicted"] - p["dada2"], 2)]
            for p in sorted(D.predictions, key=lambda p: (D.short[p["dataset"]], N.setting_rank(p["setting"])))]
    big = max(abs(p["predicted"] - p["dada2"]) for p in D.predictions) * 100
    out.append(section("AmpliTrim's predictions at each setting", f"<p>AmpliTrim predicted the share of pairs passing DADA2's filter at every setting a user can choose. The largest difference from DADA2's own pass rate was {big:.2f} points, over {len(rows)} predictions.</p>"
                       + '<details class="more"><summary>All predictions</summary>' + table(["Dataset", "Setting", "AmpliTrim predicted", "DADA2 measured", "Difference (points)"], rows, num_cols=(2, 3, 4)) + "</details>", "predictions"))
    body = hero('<a href="../../">Results</a> &rsaquo; <a href="../">ITS</a> &rsaquo; Settings', "The other ITS settings", "What minLen, truncQ and forward reads only did on the datasets with a known answer, and how well AmpliTrim predicted the filter at each setting.") + "".join(out) \
        + section("Where the numbers come from", f"<p>{mono('results/31_v04_night/v04_metrics.tsv')}, the <code>predict.json</code> of each dataset in {mono('results/31_v04_night/')} and {mono('results/31_v04_night/V04_SUMMARY.md')}.</p>", "sources")
    return page("The other ITS settings: AmpliTrim validation results", "minLen, truncQ, forward reads only, and AmpliTrim's predictions at each setting.", body, 2, "its/")


def recovery_page(D):
    def name(r):
        return D.full[r["code"]] if r["code"] in D.full else r["code"]
    rows = [[cell(esc(r["marker"]), ""), cell(esc(name(r)), "name"), cell(esc(N.setting(r["arm"])), "name"), r["samples"], r["expected"], esc(r["per_sample"]), esc(r["any_read"]), esc(r["any_sample"]), esc(r["all_samples"])]
            for r in D.recovery["rows"]]
    body = (hero('<a href="../">Results</a> &rsaquo; Recovery', "Recovery of the known members",
                 "Of the fungi and bacteria known to be in each mock community, how many did DADA2 find, per sample and in total? ITS: the DADA2 default and the trim. 16S: DADA2 at AmpliTrim's cut. "
                 "A member counts as found when it holds at least 0.1% of a sample's reads; the second count accepts any read.")
            + section("The table", '<p class="legend">Found per sample is the mean share of the expected members found in a sample. In any sample and in all samples give counts of members.</p>'
                      + table(["Marker", "Dataset", "Setting", "Samples", "Expected per sample", "Found per sample (at least 0.1%)", "Found per sample (any read)", "Found in any sample (0.1% / any read)", "Found in all samples (0.1% / any read)"], rows, num_cols=(3, 4)), "table")
            + section("Compared with the authors", f"<p>{esc(D.recovery['published'])}</p><p>Exact 16S variants recovered, pooled over the samples of each community: {esc(D.recovery['variants'])}.</p>", "published")
            + section("What \u2018found\u2019 means", '<p>A match of sequences, not a taxonomic assignment: no reference database is used and no rank is involved. See the <a href="../definitions/">definitions</a>.</p>', "meaning")
            + section("Where the numbers come from", f"<p>{mono('results/32_recovery/RECOVERY.md')} (register section 45), made from runs 31 (ITS) and 26 (16S).</p>", "sources"))
    return page("Recovery of the known members: AmpliTrim validation results", "How many known members each mock community gave back, ITS and 16S.", body, 1, "recovery/")


def criteria_page(D):
    runs = {}
    for c in D.criteria:
        runs.setdefault(c["run"], []).append(c)
    names = {r["run"]: r["name"] for r in D.runs}
    secs = []
    for run, items in sorted(runs.items()):
        note = " This run was set aside: its judgements are not used." if run == 28 else ""
        rows = [[c["id"], cell(esc(c["arm"]), "name"), (c["outcome"], "ok" if c["outcome"] == "MET" else "bad" if c["outcome"] == "NOT MET" else "dim"), cell(esc(c["text"] + (f" ({c['note']})" if c["note"] else "")), "name")] for c in items]
        met = sum(c["outcome"] == "MET" for c in items)
        secs.append(section(f"Run {run}: {names.get(run, '')}", f"<p>{len(items)} criteria: {met} met, {sum(c['outcome'] == 'NOT MET' for c in items)} not met, {sum(c['outcome'] == 'not decided' for c in items)} not decided.{note}</p>"
                            + table(["Criterion", "Setting", "Outcome", "What was measured"], rows), f"run-{run}"))
    body = hero('<a href="../">Results</a> &rsaquo; Criteria', "Every criterion and its outcome",
                f"Each test was written down, with what would count as success, before it ran. These are the {len(D.criteria)} criteria of the ITS runs 26 to 31 and what each found. "
                "Run 25's criteria are in the register (sections 31 and 32).") + "".join(secs) + section("Where the numbers come from", "<p>The summary file of each run in <a href=\"" + BLOB + "results\">results</a>, and the <a href=\"" + BLOB + "TEST_REGISTER.md\">test register</a>.</p>", "sources")
    return page("Every criterion and its outcome: AmpliTrim validation results", "The criteria of the ITS runs and what each found.", body, 1, "criteria/")


def md_inline(s):
    s = esc(s)
    s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"`(.+?)`", r"<code>\1</code>", s)
    return re.sub(r"(?<!\*)\*(?!\s)(.+?)\*(?!\*)", r"<em>\1</em>", s)


def definitions_page(D):
    out, table_rows, paras = [], [], []
    for line in D.definitions.splitlines():
        if line.startswith("| ") and not line.startswith("|---"):
            table_rows.append([c.strip() for c in line.strip("|").split("|")])
        elif line.strip() and not line.startswith("#") and not line.startswith("|---"):
            paras.append(line.strip())
    head, rows = table_rows[0], table_rows[1:]
    body = hero('<a href="../">Results</a> &rsaquo; Definitions', "Definitions", " ".join(md_inline(p) for p in paras[:1])) \
        + section("The words", table(head, [[cell(md_inline(r[0]), "name"), cell(md_inline(r[1]), "name")] for r in rows]), "words") \
        + section("What \u2018found\u2019 is not", "".join(f"<p>{md_inline(p)}</p>" for p in paras[1:]), "not") \
        + section("Settings", table(["Setting", "Name"], [[cell(mono(k), ""), cell(esc(v), "name")] for k, v in N.SETTINGS.items()]) + "<p>The codes are the register's; the pages use the names.</p>", "settings")
    return page("Definitions: AmpliTrim validation results", "What member, true, false and found mean, and that they are sequence matches, not taxonomic ranks.", body, 1, "definitions/")


def runs_page(D):
    rows = []
    for r in D.runs:
        rows.append([f'<span id="run-{r["run"]}">{r["run"]}</span>', cell(esc(r["name"]), "name"), cell("; ".join(f"{s}. {esc(t)}" for s, t in r["titles"]) or "none (descriptive work)", "name"),
                     cell(", ".join(mono(x) for x in r["folders"]) or "\u2013", "name")])
    body = hero('<a href="../">Results</a> &rsaquo; Runs', "Every run of the validation", "Which runs were made, under which register sections, and where the archive keeps their results.") \
        + section("The runs", table(["Run", "What", "Register sections", "Results in the archive"], rows), "runs") \
        + section("Where the numbers come from", f"<p>{mono('validation/organise_data.py')} and the <a href=\"{BLOB}TEST_REGISTER.md\">test register</a>.</p>", "sources")
    return page("Every run of the validation: AmpliTrim validation results", "The runs, their register sections and where their results are kept.", body, 1, "runs/")


def catalog_page(D):
    rows = []
    for r in D.names:
        c = r["archive_code"]
        links = []
        if c in OLD_PAGE:
            links.append(f'<a href="{OLD_PAGE[c]}/">first validation</a>')
        if N.is_its(r):
            links.append(f'<a href="../its/datasets/{N.slug(r)}/">ITS runs</a>')
        if any(x["code"] == c for x in D.recovery["rows"]):
            links.append('<a href="../recovery/">recovery</a>')
        rows.append([r["id"], cell(esc(r["full_name"]), "name"), esc(r["kind"]), esc(r["marker_region"]), esc(r["platform"]), esc(r["accession"]), esc(r["source"]), cell(", ".join(links) or "\u2013", "name")])
    body = hero('<a href="../">Results</a> &rsaquo; Datasets', "All datasets", f"The {len(D.names)} public datasets of the validation, each with one ID and one name, the same in the paper, the supplementary workbook and these pages. There is no DS13.") \
        + section("The datasets", table(["ID", "Name", "Kind", "Marker", "Platform", "Accession", "Source", "Results"], rows, cls="catalog"), "datasets") \
        + section("Where the names come from", f"<p><a href=\"{BLOB}docs/DATASET_NAMES.tsv\">docs/DATASET_NAMES.tsv</a> links every ID and name to the folder on the data drive and to the code the archive's scripts use.</p>", "sources")
    return page("All datasets: AmpliTrim validation results", "The datasets of the AmpliTrim validation, with their IDs, names and results.", body, 1, "datasets/")


def build(root):
    D = Data(root)
    out = {"its/index.html": its_overview(D), "its/trims/index.html": trims_page(D), "its/settings/index.html": settings_page(D), "recovery/index.html": recovery_page(D),
           "criteria/index.html": criteria_page(D), "definitions/index.html": definitions_page(D), "runs/index.html": runs_page(D), "datasets/index.html": catalog_page(D)}
    for r in D.its:
        out[f"its/datasets/{N.slug(r)}/index.html"] = its_dataset_page(D, r)
    return out


def main(root=None):
    root = pathlib.Path(root or (sys.argv[1] if len(sys.argv) > 1 else HERE.parent)).resolve()
    pages = build(root)
    for rel, text in pages.items():
        p = root / "docs" / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    print("written:", len(pages), "pages")
    return pages


if __name__ == "__main__":
    main()
