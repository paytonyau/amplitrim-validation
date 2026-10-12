"""Build the results site of the validation archive: docs/index.html, docs/datasets/<id>/index.html and docs/reports/<id>.md.
Written for a reader who has not seen our notes: plain words, one meaning for each colour, and every term explained once on the overview.
Usage: python build_results.py [path to the validation archive]   (default: the folder above this one)"""
import html
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE.parent
sys.path.insert(0, str(HERE))
import results_data as R

ARCH, SOFT = "https://github.com/paytonyau/amplitrim-validation", "https://github.com/paytonyau/amplitrim"
BLOB, ENA = ARCH + "/blob/main/", "https://www.ebi.ac.uk/ena/browser/view/"
esc = lambda s: html.escape(str(s), quote=True)
D, META = R.load(ROOT)
ORDER = R.IDS
CSS = "".join((HERE / n).read_text(encoding="utf-8") for n in ("style.css", "design_css.css", "results_css.css"))
THEMEHEAD, THEME = (HERE / "theme_head.js").read_text(encoding="utf-8"), (HERE / "theme.js").read_text(encoding="utf-8")
FAVICON = "data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 32 32'%3E%3Crect width='32' height='32' fill='%231b201d'/%3E%3Ctext x='16' y='24' font-size='24' text-anchor='middle' fill='%23e8ebe6' font-family='Georgia,serif' font-weight='700'%3EA%3C/text%3E%3C/svg%3E"
SLOT_TABLE = {"t1": "1", "t3": "3", "t4": "4", "t5": "5", "t6": "6", "t7": "7", "t8": "8", "t9": "9", "t11": "11"}


def inline(s):
    return re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", esc(s))

def tick(t): return "0" if t == 0 else f"{int(t):+d}"
def fmt(v, plus=True):
    s = f"{v:+.1f}" if plus else f"{v:.1f}"
    return s.replace("+0.0", "0.0").replace("-0.0", "0.0")
def off_by(vs):
    v = R.num(vs)
    return "exact" if v == 0 else f"{abs(v):.1f} points too {'high' if v > 0 else 'low'}"
def merged_vs(vs):
    v = R.num(vs)
    return "the same" if v == 0 else f"{abs(v):.1f} points {'more' if v > 0 else 'fewer'}"
def role_parts(d):
    cls, label, desc = R.role_plain(d["t1"]["Role"]); extra = re.search(r"\((.*)\)", d["t1"]["Role"])
    return cls, label, desc + (f" ({extra.group(1)})" if extra else "")


# ------------------------------------------------------------------ charts: one meaning per colour. Blue = more than the comparison, orange = fewer; the band is the margin; every mark carries its value in data-v
W = 520
def svg(h, label, body): return f'<svg viewBox="0 0 {W} {h}" role="img" aria-label="{esc(label)}"><title>{esc(label)}</title>{body}</svg>'
def nice(x):
    for n in (1, 2, 3, 4, 5, 10, 20, 30, 40, 60, 100):
        if n >= x: return n
    return 100
def cls_of(v, signed): return "acc" if (not signed or v >= 0) else "warn"


def marker_chart(value_s, lim, band, band_label, title, ends, signed=False):
    v = R.num(value_s); xs = lambda t: 30 + (t + lim) / (2 * lim) * 460
    b = [f'<rect class="band2" x="{xs(-band):.1f}" y="26" width="{xs(band) - xs(-band):.1f}" height="26" rx="4"/>', '<line class="ink" x1="30" y1="39" x2="490" y2="39"/>']
    for t in (-lim, 0, lim): b.append(f'<line class="ink" x1="{xs(t):.1f}" y1="34" x2="{xs(t):.1f}" y2="44"/><text class="tm" x="{xs(t):.1f}" y="62" text-anchor="middle">{tick(t)}</text>')
    b.append(f'<circle class="{cls_of(v, signed)}" data-v="{esc(value_s)}" cx="{xs(v):.1f}" cy="39" r="7"/><text class="{"ta" if cls_of(v, signed) == "acc" else "tw"}" x="{xs(v):.1f}" y="18" text-anchor="middle">{esc(value_s)}</text>')
    b.append(f'<text class="tm" x="30" y="82">\u2190 {esc(ends[0])}</text><text class="tm" x="490" y="82" text-anchor="end">{esc(ends[1])} \u2192</text>')
    b.append(f'<text class="ta" x="{(xs(-band) + xs(band)) / 2:.1f}" y="100" text-anchor="middle">shaded: {esc(band_label)}</text>')
    return svg(108, title, "".join(b))


def strip_chart(items, lim, band, band_label, title, ends):
    x0, x1 = 215, 500; xs = lambda t: x0 + (t + lim) / (2 * lim) * (x1 - x0); rh = 24; top = 20
    h = top + rh * len(items) + 54; bottom = top + rh * len(items)
    b = [f'<rect class="band2" x="{xs(-band):.1f}" y="{top - 4}" width="{xs(band) - xs(-band):.1f}" height="{rh * len(items) + 2}"/>', f'<line class="ink" x1="{xs(0):.1f}" y1="{top - 6}" x2="{xs(0):.1f}" y2="{bottom}"/>']
    for i, (label, vs) in enumerate(items):
        y = top + i * rh + 9; v = R.num(vs)
        b.append(f'<text class="t" x="8" y="{y + 4}" style="font-size:12px">{esc(label)}</text><line class="grid" x1="{x0}" y1="{y}" x2="{x1}" y2="{y}"/><circle class="acc" data-v="{esc(vs)}" cx="{xs(v):.1f}" cy="{y}" r="5"/>')
    for t in (-lim, 0, lim): b.append(f'<text class="tm" x="{xs(t):.1f}" y="{bottom + 14}" text-anchor="middle">{tick(t)}</text>')
    b.append(f'<text class="tm" x="{x0}" y="{bottom + 32}">\u2190 {esc(ends[0])}</text><text class="tm" x="{x1}" y="{bottom + 32}" text-anchor="end">{esc(ends[1])} \u2192</text><text class="ta" x="{(xs(-band) + xs(band)) / 2:.1f}" y="{bottom + 48}" text-anchor="middle">shaded: {esc(band_label)}</text>')
    return svg(h, title, "".join(b))


def diverging(items, lim, band, title, ends):
    """Bars from zero. items: (label, value string) or ('#', group title)."""
    x0, x1 = 215, 500; xs = lambda t: x0 + (t + lim) / (2 * lim) * (x1 - x0); rh = 24; top = 12
    rows = len(items); bottom = top + rh * rows; h = bottom + 50
    b = []
    if band: b.append(f'<rect class="band2" x="{xs(-band):.1f}" y="{top - 4}" width="{xs(band) - xs(-band):.1f}" height="{rh * rows + 2}"/>')
    b.append(f'<line class="ink" x1="{xs(0):.1f}" y1="{top - 6}" x2="{xs(0):.1f}" y2="{bottom}"/>')
    for i, (label, vs) in enumerate(items):
        y = top + i * rh
        if label == "#":
            b.append(f'<text class="ta" x="8" y="{y + 16}">{esc(vs)}</text>'); continue
        v = R.num(vs); w = abs(xs(v) - xs(0)); x = min(xs(v), xs(0))
        b.append(f'<text class="t" x="8" y="{y + 16}" style="font-size:12px">{esc(label)}</text>')
        b.append(f'<rect class="{"bar" if v >= 0 else "barn"}" data-v="{esc(vs)}" x="{x:.1f}" y="{y + 5}" width="{max(w, 1.5):.1f}" height="14" rx="2"/>')
        b.append(f'<text class="tm" x="{xs(v) + (6 if v >= 0 else -6):.1f}" y="{y + 16}" text-anchor="{"start" if v >= 0 else "end"}">{esc(vs)}</text>')
    for t in (-lim, 0, lim): b.append(f'<text class="tm" x="{xs(t):.1f}" y="{bottom + 14}" text-anchor="middle">{tick(t)}</text>')
    b.append(f'<text class="tm" x="{x0}" y="{bottom + 32}">\u2190 {esc(ends[0])}</text><text class="tm" x="{x1}" y="{bottom + 32}" text-anchor="end">{esc(ends[1])} \u2192</text>')
    if band: b.append(f'<text class="ta" x="{(xs(-band) + xs(band)) / 2:.1f}" y="{bottom + 48}" text-anchor="middle">shaded: within {band} points</text>')
    return svg(h, title, "".join(b))


def shares(items, title):
    x0, x1 = 250, 470; rh = 28; top = 8; h = top + rh * len(items) + 24
    b = [f'<line class="ink" x1="{x0}" y1="{top - 2}" x2="{x0}" y2="{top + rh * len(items)}"/>']
    for i, (label, vs) in enumerate(items):
        y = top + i * rh; p = R.pct(vs)
        b.append(f'<text class="t" x="8" y="{y + 18}" style="font-size:12px">{esc(label)}</text><rect class="bar" data-v="{esc(vs)}" x="{x0}" y="{y + 5}" width="{max((x1 - x0) * p / 100, 1.5):.1f}" height="16" rx="2"/><text class="tm" x="{x0 + (x1 - x0) * p / 100 + 6:.1f}" y="{y + 18}">{esc(vs)}</text>')
    for t in (0, 50, 100): b.append(f'<text class="tm" x="{x0 + (x1 - x0) * t / 100:.1f}" y="{top + rh * len(items) + 14}" text-anchor="middle">{t}%</text>')
    return svg(h, title, "".join(b))


# ------------------------------------------------------------------ a section is written once and comes out as a web page and as a plain document
class Sec:
    def __init__(self, sid, title): self.id, self.title, self.h, self.m = sid, title, [], []
    def p(self, text): self.h.append(f"<p>{inline(text)}</p>"); self.m.append(text)
    def h3(self, text): self.h.append(f"<h3>{esc(text)}</h3>"); self.m.append(f"**{text}**")
    def fine(self, text): self.h.append(f'<p class="pending">{inline(text)}</p>'); self.m.append(f"*{text}*")
    def chart(self, svg_, caption): self.h.append(f'<figure class="chart">{svg_}<figcaption>{esc(caption)}</figcaption></figure>'); self.m.append(f"*{caption}*")
    def table(self, header, rows):
        th = "".join(f'<th scope="col">{esc(c)}</th>' for c in header)
        tr = "".join("<tr>" + "".join((f'<th scope="row">{esc(c)}</th>' if i == 0 else f"<td>{esc(c)}</td>") for i, c in enumerate(r)) + "</tr>" for r in rows)
        self.h.append(f'<div class="scroll"><table><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table></div>')
        self.m.append("| " + " | ".join(header) + " |\n|" + "|".join("---" for _ in header) + "|\n" + "\n".join("| " + " | ".join(str(c).replace("|", "/") for c in r) + " |" for r in rows))
    def html(self): return f'<section class="wrap" id="{self.id}"><h2>{esc(self.title)}</h2>{"".join(self.h)}</section>'
    def md(self): return f"## {self.title}\n\n" + "\n\n".join(self.m)


def plain(v): return R.VALUE_PLAIN.get(v, v)


def dataset_sections(d):
    out = []; t3 = d["t3"]; name = d["name"]
    # 1. the estimate
    s = Sec("estimate", "How good was AmpliTrim's estimate?")
    s.p("Before DADA2 runs, AmpliTrim estimates two things for each setting: the share of read pairs that will pass DADA2's quality filter, and the share that will go on to merge. We then ran DADA2 for real and compared.")
    s.chart(marker_chart(t3["Filter error"], 3, 2, "within 2 points, the target we set in advance", f"{name}: how far AmpliTrim's estimate of the pass rate was from the measured pass rate, in points", ("estimate too low", "estimate too high")),
            "How far AmpliTrim's estimate of the share of pairs passing the filter was from what DADA2 measured. Zero is exact; inside the shaded band meets the target.")
    s.table(["Samples", "Pass-rate estimate off by (points)", "Merge estimate off by (points; average and range over the settings tried)", "Same order as DADA2 (1.00 = identical)", "The other tool's pass-rate estimate off by (points)"],
            [[t3["Samples"], t3["Filter error"], t3["Yield excess"], t3["Rank correlation"], t3["the other tool error"]]])
    s.p("The merge estimate is too high because it stops before DADA2 denoises the reads. Use it to compare settings, not as a forecast.")
    out.append(s)
    # 2. the other tool, same maxEE
    s = Sec("versus", "Compared with the other tool, at the same maxEE")
    s.p("The other tool is an existing automatic tool for choosing truncation lengths. Given the same maxEE, we compared how many read pairs DADA2 merges with AmpliTrim's cut and with the other tool's cut. A cut such as 247/199 means the forward reads are cut to 247 nt and the reverse reads to 199 nt.")
    if "t4" in d:
        t4 = d["t4"]; cutk = [k for k in t4 if k.startswith("the other tool")][0]; dv = t4["Difference"].split(" ")[0]
        s.chart(marker_chart(dv, 4, 2, "within 2 points, the margin set in advance: equivalent", f"{name}: pairs merged by DADA2 with AmpliTrim's cut minus with the other tool's, in points", ("fewer pairs merged", "more pairs merged"), signed=True),
                "The difference in the share of read pairs DADA2 merged: AmpliTrim's cut minus the other tool's. Blue means AmpliTrim merged more, orange fewer. Inside the shaded band the two cuts are equivalent.")
        s.table(["Samples", "The other tool's cut and maxEE", "AmpliTrim's cut", "Difference in pairs merged (points)"], [[t4["Samples"], t4[cutk], t4["AmpliTrim: cut"], t4["Difference"]]])
        if "to" in t4["Difference"]: s.fine("The numbers in brackets are a 95% interval from resampling the samples; it is given only when there are 10 or more samples.")
    else:
        s.p(R.NO_COMPARISON[t3["the other tool error"]])
    out.append(s)
    # 3. maxEE
    s = Sec("maxee", "Choosing the maxEE limit")
    s.p("maxEE is the most expected errors DADA2 allows in a read; reads over the limit are dropped. DADA2's default is 2. AmpliTrim's released rule sets the limit for each read direction from your own reads, never below 2, with no upper limit.")
    s.h3("What the released rule does")
    if "rule" in d:
        r = d["rule"]; rk = "the rule's maxEE"
        note = R.CAP_NOTE.get(next(k for k, v in D.items() if v is d))
        s.p(f"On this dataset the rule changes the limit to **maxEE {r[rk]}**. DADA2 merged this share of the read pairs:")
        third = f"Released rule (maxEE {r[rk]})"
        s.chart(shares([("Default maxEE 2", r["maxEE 2"]), ("maxEE 2,5 (looser on the reverse read)", r["2,5"]), (third, r["the rule"])], f"{name}: share of pairs merged by DADA2 under three maxEE settings"),
                "The share of the input read pairs that DADA2 merged with three maxEE settings, from the final build of the software.")
        s.table(["Default maxEE 2", "maxEE 2,5", "Released rule", "The rule's maxEE"], [[r["maxEE 2"], r["2,5"], r["the rule"], r[rk]]])
        if note:
            s.p(f"The rule has no cap. With a cap of 5 it would have set **maxEE {note['capped_max_ee']}** here, and DADA2 merged {note['capped_merged']} of the pairs; without one it sets maxEE {r[rk]} and DADA2 merged {r['the rule']}.")
    else:
        s.p("On this dataset the released rule keeps DADA2's default of 2, as the validation page records for the clean datasets and for ITS.")
    t11 = d["t11"]
    s.h3("What a looser limit would do")
    s.p("We also tried relaxing the limit on the reverse read to 5 (maxEE 2,5) and measured the change in pairs merged, in points of the input pairs.")
    groups = [("#", "Relaxing to maxEE 2,5, compared with the default"), ("Change in pairs merged", t11["G minus A"])]
    other = []
    if "t5" in d and R.num(d["t5"]["A"]) is not None: other.append(("Default maxEE 2", d["t5"]["A"]))
    if R.num(t11["G against the other tool"]) is not None: other.append(("maxEE 2,5", t11["G against the other tool"]))
    if other: groups += [("#", "AmpliTrim against the other tool, which picks its own maxEE")] + other
    lim = nice(max(abs(R.num(v)) for k, v in groups if k != "#"))
    s.chart(diverging(groups, lim, None, f"{name}: effect of a looser maxEE, in points", ("fewer pairs merged", "more pairs merged")),
            "Changes in the share of pairs DADA2 merged, in points. Blue means more pairs merged, orange fewer. The lower group compares AmpliTrim with the other tool using its own choice of cut and maxEE.")
    s.table(["Change in pairs merged (points)", "Quality of the extra reads kept", "Accuracy (needs expected sequences)"], [[t11["G minus A"], plain(t11["G noise"]), plain(t11["G accuracy"])]])
    if other:
        s.p("Against the other tool, which picks its own cut and maxEE, DADA2 merged this many more (+) or fewer (-) pairs, in points:")
        s.table(["AmpliTrim setting", "Difference from the other tool's own choice (points)"], [[{"Default maxEE 2": "Default maxEE 2", "maxEE 2,5": "maxEE 2,5"}[k], v] for k, v in other])
    if "t6" in d:
        hd = ["Setting", "Variants found", "Strains found", "Reads matching an expected sequence exactly"]
        s.h3("Accuracy on this mock community")
        s.p("This dataset has expected sequences, so each setting can be scored.")
        s.table(hd, [[R.ARM_PLAIN[r["Arm"]], r["Variants"], r["Strains"], r["Reads in exact ASVs"]] for r in d["t6"]])
    out.append(s)
    # 4. amplicon length
    if "t7" in d or "t8" in d:
        s = Sec("length", "Amplicon length and the longest amplicons")
        if "t7" in d:
            t7 = d["t7"]; ks = [k for k in t7 if k != "Dataset"]; lens = t7[ks[0]].split(" / ")
            items = [(lab, t7[k].split(", ")[-1]) for lab, k in (("Other tool: longest length", ks[1]), ("Other tool: typical length", ks[2]), ("AmpliTrim, default", ks[3])) if R.pct(t7[k]) is not None]
            s.h3("The amplicon length the other tool needs")
            s.p(f"The other tool needs you to tell it the amplicon length. We tried the longest ({lens[0]} nt) and the typical ({lens[1]} nt); it chose different cuts. AmpliTrim measures the lengths from your reads, so it needs no such number.")
            s.chart(shares(items, f"{name}: share of pairs merged by DADA2"), "The share of the input read pairs DADA2 merged with the cut each tool chose.")
            s.table(["Length given: longest / typical", "Other tool, given the longest", "Other tool, given the typical", "AmpliTrim, default"], [[t7[ks[0]], t7[ks[1]], t7[ks[2]], t7[ks[3]]]])
        if "t8" in d:
            t8 = d["t8"]
            s.h3("Leaving out the longest amplicons")
            s.p("AmpliTrim sets the shortest cut that can merge your amplicons, but it leaves out the very longest few so that a handful of outliers do not force long, noisy reads. By default it leaves out the longest 0.2%. DADA2's merged share, at maxEE 2, when the longest 0.2%, 1% and 5% are left out:")
            s.chart(shares([("Longest 0.2% left out (default)", t8["0.2%"]), ("Longest 1% left out", t8["1%"]), ("Longest 5% left out", t8["5%"])], f"{name}: share of pairs merged by DADA2 at three settings"), "The share of the input read pairs DADA2 merged.")
            s.table(["Longest 0.2% left out (default)", "Longest 1%", "Longest 5%"], [[t8["0.2%"], t8["1%"], t8["5%"]]])
        out.append(s)
    # 5. stability
    if "t9" in d:
        t9 = d["t9"]; s = Sec("stability", "Would the choice change on other samples?")
        s.p("Each tool chose a cut again on three separate sets of samples from this dataset. The choice counts as the same if the cuts differ by at most 10 nt and the maxEE is the same. Cuts are forward/reverse; the maxEE is in brackets.")
        s.table(["AmpliTrim's cuts", "Same each time?", "The other tool's cuts", "Same each time?"], [[t9["AmpliTrim"], t9["AmpliTrim holds"], t9["the other tool"], t9["the other tool holds"]]])
        out.append(s)
    return out


def tests_not_taken(d):
    names = {"t4": "the comparison with the other tool at the same maxEE", "t6": "the accuracy scoring (it needs expected sequences)", "t7": "the test of the amplicon length given to the other tool", "t8": "the test of leaving out the longest amplicons", "t9": "the test of choosing again on other samples"}
    return [v for k, v in names.items() if k not in d]


def sources_line(d):
    tabs = sorted({SLOT_TABLE[k] for k in SLOT_TABLE if k in d}, key=int)
    return f"Tables {', '.join(tabs)} of the findings" + (" and the table of the final build on the validation page" if "rule" in d else "")


def page_shell(title, desc, body, base):
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(title)}</title>
<meta name="description" content="{esc(desc)}">
<meta name="theme-color" content="#1b201d">
<meta property="og:type" content="website">
<meta property="og:title" content="{esc(title)}">
<meta property="og:description" content="{esc(desc)}">
<link rel="icon" href="{FAVICON}">
<script>{THEMEHEAD}</script>
<style>
{CSS}
</style>
</head>
<body>
<a class="skip" href="#main">Skip to content</a>
<header class="site"><div class="wrap bar">
<a class="brand" href="{base}">Ampli<span class="trim">Trim</span> results</a>
<nav aria-label="Main"><ul>
<li><a href="{base}">Overview</a></li>
<li><a href="{base}#datasets">Datasets</a></li>
<li><a href="{base}how-it-works/">How it works</a></li>
<li><a href="{SOFT}">Software</a></li>
<li><a href="{ARCH}">Archive</a></li>
</ul></nav>
</div></header>
<main id="main">
{body}
</main>
<footer><div class="wrap">
<p><a href="{base}">Overview</a> · <a href="{SOFT}">AmpliTrim on GitHub</a> · <a href="{ARCH}">Validation archive</a> · <a href="{BLOB}docs/FINDINGS.md">The findings</a> · <a href="{BLOB}TEST_REGISTER.md">Test register</a></p>
<p>By <a href="https://orcid.org/0000-0002-3283-0370">Payton (Tung On) Yau</a>. The numbers come from the findings of the validation archive. The tool we compare with is not named on these pages; the findings name it.</p>
</div></footer>
<script>
{THEME}
</script>
</body>
</html>
'''


def dataset_page(i):
    d = D[i]; n = ORDER.index(i); prev_i = ORDER[n - 1] if n else None; next_i = ORDER[n + 1] if n < len(ORDER) - 1 else None
    t1, t3 = d["t1"], d["t3"]; cls, label, desc = role_parts(d)
    ref = f"reference [{d['ref']}] in the findings (to be verified)" if d["ref"] else "to add"
    pd = d.get("t5", {}).get("pass_default") or d["t11"]["pass_default"]
    glance = [("Pass-rate estimate", off_by(t3["Filter error"]), "AmpliTrim's estimate of the share of pairs passing DADA2's filter, against what DADA2 measured (target: within 2 points)"),
              ("Pass DADA2's default filter", f"{pd}%", "of the read pairs pass at DADA2's default maxEE of 2")]
    glance.append(("Against the other tool", merged_vs(d["t4"]["Difference"].split(" ")[0]), "pairs merged at the same maxEE") if "t4" in d else ("Against the other tool", "no comparison", esc(R.NO_COMPARISON[t3["the other tool error"]].split(",")[0].rstrip("."))))
    glance.append(("maxEE the released rule sets", d["rule"]["the rule's maxEE"], f"merging {d['rule']['the rule']} of the pairs (default maxEE 2: {d['rule']['maxEE 2']})") if "rule" in d else ("maxEE the released rule sets", "2", "the rule keeps DADA2's default here"))
    tiles = "".join(f'<li><span class="lab">{esc(a)}</span><span class="big">{esc(b)}</span><span class="lab">{esc(c) if "&" not in c else c}</span></li>' for a, b, c in glance)
    secs = dataset_sections(d); missing = tests_not_taken(d)
    miss = f"<p>This dataset was not part of: {esc('; '.join(missing))}.</p>" if missing else "<p>This dataset took part in every test.</p>"
    pager = '<nav class="pager" aria-label="Datasets">' + (f'<a href="../{prev_i}/">&larr; {esc(D[prev_i]["name"])}</a>' if prev_i else "<span></span>") + (f'<a href="../{next_i}/">{esc(D[next_i]["name"])} &rarr;</a>' if next_i else "<span></span>") + "</nav>"
    body = f'''<section class="hero"><div class="wrap">
<p class="crumb"><a href="../../">Results</a> &rsaquo; Dataset {n + 1} of {len(ORDER)}</p>
<h1>{esc(d["name"])}</h1>
<p class="chips"><span class="chip {cls}">{esc(label)}</span><span class="chip measured">{"Mock community" if t3["Kind"] == "mock" else "Real samples"}</span></p>
<dl class="facts">
<dt>Marker, region, platform</dt><dd>{esc(t1["Marker, region, platform"])}</dd>
<dt>Samples</dt><dd>{esc(t1["Samples"])}</dd>
<dt>Known truth</dt><dd>{esc(t1["Truth"])}</dd>
<dt>How we used it</dt><dd>{esc(label)}: {esc(desc)}</dd>
<dt>Data</dt><dd>{", ".join(f'<a href="{ENA}{a}">{a}</a>' for a in d["accessions"])}</dd>
<dt>Source study</dt><dd>{esc(ref)}</dd>
</dl>
<ul class="glance">{tiles}</ul>
<div class="notebox"><p><strong>In short.</strong> {esc(d["notes"])}</p></div>
</div></section>
{"".join(s.html() for s in secs)}
<section class="wrap" id="tests"><h2>Which tests it took part in</h2>
{miss}
</section>
<section class="wrap" id="sources"><h2>Where the numbers come from</h2>
<p>{esc(sources_line(d))}, in the <a href="{BLOB}docs/FINDINGS.md">findings</a> and on the <a href="{BLOB}docs/validation.md">validation page</a>. The criteria were written in the <a href="{BLOB}TEST_REGISTER.md">test register</a> before the runs. The same numbers, as a plain document, are in the <a href="{BLOB}docs/reports/{i}.md">evidence report for this dataset</a>.</p>
<p class="pending">The report that AmpliTrim itself wrote on this dataset is not in the archive yet. It is added when the results of the runs are collected.</p>
{pager}
</section>'''
    return page_shell(f"{d['name']}: AmpliTrim test results", f"How AmpliTrim did on {d['name']}: the estimate, the comparison with the other tool and the choice of maxEE.", body, "../../")


def dataset_report(i):
    d = D[i]; t1 = d["t1"]; cls, label, desc = role_parts(d); pd = d.get("t5", {}).get("pass_default") or d["t11"]["pass_default"]
    lines = [f"# {d['name']}: evidence report", "",
             "*AmpliTrim validation archive. The numbers come from [`docs/FINDINGS.md`](../FINDINGS.md) and [`docs/validation.md`](../validation.md). The tool we compare with is not named here; the findings name it.*", "",
             f"- **Marker, region, platform:** {t1['Marker, region, platform']}", f"- **Samples:** {t1['Samples']}", f"- **Known truth:** {t1['Truth']}", f"- **How we used it:** {label}: {desc}",
             "- **Data:** " + ", ".join(f"[{a}]({ENA}{a})" for a in d["accessions"]), "- **Source study:** " + (f"reference [{d['ref']}] in the findings (to be verified)" if d["ref"] else "to add"),
             f"- **Share of pairs passing DADA2's default filter (maxEE 2):** {pd}%", "", "## In short", "", d["notes"], ""]
    for s in dataset_sections(d): lines += [s.md(), ""]
    miss = tests_not_taken(d); lines += ["## Which tests it took part in", "", ("This dataset was not part of: " + "; ".join(miss) + ".") if miss else "This dataset took part in every test.", ""]
    lines += ["## Where the numbers come from", "", sources_line(d) + ".", "", f"The web page of this dataset is in the results site of the archive: `datasets/{i}/`.", ""]
    return "\n".join(lines)


# ------------------------------------------------------------------ the overview
def overview_page():
    dv = R.derived(D); st = dv["stats"]
    roles = {"dev": 0, "heldout": 0, "case": 0}
    for i in ORDER: roles[role_parts(D[i])[0]] += 1
    sb = [("Held", "met", "AmpliTrim's estimate of the share of pairs that pass DADA2's filter", f"On all {len(ORDER)} datasets the estimate was within {fmt(dv['max_abs_fe'], False)} points of what DADA2 measured. The target we set before testing was 2 points. The average miss was {fmt(dv['mean_abs_fe'], False)} points."),
          ("Partly", "partly", "AmpliTrim's estimate of the share of pairs that merge", f"On {dv['n_rank_one']} of the {dv['n_rank']} datasets where the order could be compared, it put the settings in exactly the order DADA2 did, and it did less well on " + " and ".join(f"{D[i]['name']} ({r:.2f})" for i, r in dv["low_rank"]) + f". The estimate is also off in size, because it stops before DADA2 denoises the reads: from {fmt(abs(dv['ye_min']), False)} points too low to {fmt(dv['ye_max'], False)} points too high, depending on the dataset. Use it to compare settings, not as a forecast."),
          ("Held", "met", "AmpliTrim's cut against the other tool's, at the same maxEE", f"AmpliTrim's cut merged more pairs on {st['all']['ahead']} of {st['all']['n']} datasets, the same on {st['all']['tied']} and fewer on {st['all']['behind']} (average {fmt(st['all']['mean'])} points). Only one dataset differs by more than the 2-point margin we set, so the cuts are equivalent, not better."),
          ("Partly", "partly", "Getting the same cut on other samples", f"Run on three separate sets of samples, AmpliTrim chose the same cut on {dv['hold_a']} of {dv['n_t9']} datasets and the other tool on {dv['hold_o']} of {dv['n_t9']}. AmpliTrim's misses are small moves along a flat stretch where the yield barely changes."),
          ("Sensible default, not independently tested", "partly", "The maxEE rule", f"It follows the advice of DADA2's author: maxEE 2 is a good general-purpose value, and relaxing it somewhat is also OK. The rule keeps 2 unless your reads need more. The numbers are ours: we chose the 75th percentile after seeing our test data, and the author gives no ceiling, so the rule has none; the report notes a limit above 5. So it is a sensible default that has not been tested on independent data, and it did not reach the result we had decided in advance would make us adopt it (that test used an earlier version of the rule, with a cap of 5, which merged 60.7% on average). On the five datasets where it differs from the default, DADA2 merged {META['rule_mean'][1]} of the pairs on average with the default maxEE, {META['rule_mean'][2]} with maxEE 2,5 and {META['rule_mean'][3]} with the released rule, measured afterwards. We did not apply the criteria again to the new numbers, because we had seen them first. We keep the rule and say so openly.")]
    sbh = "".join(f'<tr><th scope="row">{esc(t)}</th><td><span class="chip {c}">{esc(l)}</span><p>{esc(x)}</p></td></tr>' for l, c, t, x in sb)
    fe_items = [(D[i]["name"], D[i]["t3"]["Filter error"]) for i in ORDER]
    g_items = []
    for g, title in (("real", "Real samples"), ("mock", "Mock communities"), ("wastewater", "Wastewater*")):
        g_items.append(("#", title))
        for i in ORDER:
            if "t4" in D[i] and ((i == "wastewater") == (g == "wastewater")) and (g == "wastewater" or D[i]["t3"]["Kind"] == g):
                g_items.append((D[i]["name"], D[i]["t4"]["Difference"].split(" ")[0]))
    cards = []
    for i in ORDER:
        d = D[i]; cls, label, _ = role_parts(d)
        cmp_ = f"same maxEE, against the other tool: {merged_vs(d['t4']['Difference'].split(' ')[0])}" if "t4" in d else "no comparison with the other tool"
        cards.append(f'<li class="dcard"><h3><a href="datasets/{i}/">{esc(d["name"])}</a></h3><p class="chips"><span class="chip {cls}">{esc(label)}</span></p><p>{esc(d["t1"]["Marker, region, platform"].split(";")[0])}</p><p class="nums">Pass-rate estimate: {esc(off_by(d["t3"]["Filter error"]))}. {esc(cmp_[0].upper() + cmp_[1:])}.</p></li>')
    body = f'''<section class="hero" id="top"><div class="wrap">
<h1>How well does AmpliTrim work?</h1>
<p class="tagline">Test results on {len(ORDER)} public datasets, one page each.</p>
<p class="lead">AmpliTrim chooses DADA2 settings for you. To find out how good its choices are, we tested it against real DADA2 runs on {len(ORDER)} public datasets: mock communities, real samples and one ITS case study. We wrote down what we expected before the runs. These pages show what held and what did not.</p>
</div></section>
<section class="wrap" id="how"><h2>How to read these pages</h2>
<dl class="facts">
<dt>Estimate and measurement</dt><dd>AmpliTrim <strong>estimates</strong> how many read pairs will pass DADA2's filter and merge. We then ran DADA2 and <strong>measured</strong> it. The pages compare the two.</dd>
<dt>Points</dt><dd>Differences are in percentage points of the read pairs that went in. A difference of +2 means 2 more pairs in every 100.</dd>
<dt>maxEE</dt><dd>The most expected errors DADA2 allows in a read. Reads over the limit are dropped. DADA2's default is 2.</dd>
<dt>A cut</dt><dd>Where each read is cut: 247/199 means the forward reads at 247 nt and the reverse reads at 199 nt.</dd>
<dt>The other tool</dt><dd>An existing automatic tool for choosing truncation lengths. We compare AmpliTrim with it; the findings name it.</dd>
<dt>The shaded band</dt><dd>In the charts, the margin we set in advance. A result inside the band counts as equivalent.</dd>
<dt>Colours</dt><dd>Blue: more than the comparison. Orange: fewer. Marks drawn in a single colour are plain measurements.</dd>
<dt>Our datasets</dt><dd>{roles['heldout']} datasets were kept aside and not used in the design (independent tests). {roles['dev']} were used while designing AmpliTrim (Mock community A and Soil, V4), so they are shown but are not independent evidence. One is the ITS case study.</dd>
<dt>Rounds of tests</dt><dd>We tested in three rounds. The first round called two cuts equivalent if they merged within 1 point of each other; later rounds used 2 points, which is what these pages show. A dataset's page says when the stricter margin changes the picture.</dd>
</dl>
</section>
<section class="wrap" id="scoreboard"><h2>What held</h2>
<div class="scroll"><table class="decisions"><thead><tr><th scope="col">What we tested</th><th scope="col">What happened</th></tr></thead><tbody>{sbh}</tbody></table></div>
<p>The full list of criteria and outcomes, written before the runs, is in the <a href="{BLOB}TEST_REGISTER.md">test register</a>.</p>
</section>
<section class="wrap" id="charts"><h2>All datasets at a glance</h2>
<h3>How close was the pass-rate estimate?</h3>
<figure class="chart">{strip_chart(fe_items, 3, 2, "within 2 points, the target we set in advance", "How far AmpliTrim's estimate of the pass rate was from the measured pass rate, in points, on each dataset", ("estimate too low", "estimate too high"))}<figcaption>Each mark is one dataset: how far AmpliTrim's estimate of the share of pairs passing DADA2's filter was from what DADA2 measured. Zero is exact. Every mark is inside the band.</figcaption></figure>
<h3>AmpliTrim's cut against the other tool's, at the same maxEE</h3>
<figure class="chart">{diverging(g_items, 3, 2, "Pairs merged by DADA2: AmpliTrim's cut minus the other tool's, in points, by group", ("AmpliTrim merged fewer", "AmpliTrim merged more"))}<figcaption>The difference in the share of read pairs DADA2 merged, AmpliTrim's cut minus the other tool's. Blue: AmpliTrim merged more. Orange: fewer. Inside the shaded band, the cuts are equivalent.</figcaption></figure>
<p>* Wastewater is the only real dataset in these tests with poor read quality: its reads are noisy (an estimated 54% of its pairs pass DADA2's default filter), and the loss is to read quality, not to short reads. Its reads were also trimmed to mixed lengths, and every pair was in reverse orientation. It is the one case where the other tool was ahead. Its page explains why.</p>
<p>By group: real samples ({st['real']['n']} datasets): ahead on {st['real']['ahead']}, tied on {st['real']['tied']} (average {fmt(st['real']['mean'])} points). Mock communities ({st['mock']['n']}): ahead on {st['mock']['ahead']}, tied on {st['mock']['tied']} (average {fmt(st['mock']['mean'])} points). Wastewater ({st['wastewater']['n']}): behind by {fmt(abs(st['wastewater']['mean']), False)} point.</p>
</section>
<section class="wrap" id="datasets"><h2>The {len(ORDER)} datasets</h2>
<ul class="dcards">{"".join(cards)}</ul>
</section>
<section class="wrap" id="more"><h2>Where the full record is</h2>
<ul>
<li><a href="{BLOB}docs/FINDINGS.md">The findings</a>: every table, with the method.</li>
<li><a href="{BLOB}docs/validation.md">The validation page</a>: the evidence behind each recommendation of the released package.</li>
<li><a href="{BLOB}TEST_REGISTER.md">The test register</a>: the criteria and predictions, written before the runs.</li>
<li><a href="{ARCH}">The archive</a>: the scripts that reproduce the runs.</li>
</ul>
</section>'''
    return page_shell("How well does AmpliTrim work? Test results on 12 datasets", "What held and what did not when AmpliTrim was tested against real DADA2 runs on 12 public datasets.", body, "./")


def main():
    docs = ROOT / "docs"; (docs / "datasets").mkdir(parents=True, exist_ok=True); (docs / "reports").mkdir(exist_ok=True)
    (docs / "index.html").write_text(overview_page(), encoding="utf-8"); (docs / ".nojekyll").write_text("", encoding="utf-8")
    for i in ORDER:
        (docs / "datasets" / i).mkdir(exist_ok=True)
        (docs / "datasets" / i / "index.html").write_text(dataset_page(i), encoding="utf-8")
        (docs / "reports" / f"{i}.md").write_text(dataset_report(i), encoding="utf-8")
    print("written: docs/index.html,", len(ORDER), "dataset pages and", len(ORDER), "reports")


if __name__ == "__main__":
    main()
