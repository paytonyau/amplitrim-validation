"""Build the 16S section of the validation archive's website and the pages that sit above both markers: docs/16s/index.html, docs/16s/datasets/<dataset>/index.html,
docs/index.html (the overview of 16S and ITS) and the 'moved' pages at the addresses of the first version (docs/datasets/<old name>/index.html).
The 16S pages have the same layout as the ITS pages (build_its.py): tiles, then a table for each question. Their numbers are those of docs/FINDINGS.md (read by results_data.py);
FIGARO is named; the letters the nights gave to their arms are replaced by words.
Usage: python build_16s.py [path to the validation archive]   (default: the folder above this one)"""
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else HERE.parent
sys.path.insert(0, str(HERE))
import names as N                                                          # noqa: E402
import results_data as R                                                   # noqa: E402
import build_its as I                                                      # noqa: E402
from shell import BLOB, ENA, esc, page_shell                               # noqa: E402

NOTE = "The numbers are read from the findings of the validation archive (docs/FINDINGS.md) and the results folders it lists; the register fixed every criterion before its run."
IDS16 = [i for i in R.IDS if i != "its"]


def fig(s):
    """The findings say 'the other tool'; the pages say FIGARO."""
    return (s.replace("The other tool's", "FIGARO's").replace("the other tool's", "FIGARO's").replace("The other tool", "FIGARO").replace("the other tool", "FIGARO"))


def vs_fig(v):
    """A difference against FIGARO's own cut and maxEE, or the reason there is none."""
    return esc(v) if num(v) is not None else "FIGARO gave no result"


def page(title, desc, body, depth, here, head_extra=""):
    return page_shell(title, desc, body, "./" if depth == 0 else "../" * depth, NOTE, here, head_extra)


class Site:
    def __init__(self, root):
        self.root = pathlib.Path(root)
        self.D, self.META = R.load(self.root)
        self.names = N.load(self.root)
        self.by_code = N.by_code(self.names)
        self.row = {i: self.by_code[R._CODE[i]] for i in R.IDS}
        self.its = I.Data(self.root)

    def slug(self, i):
        return N.slug(self.row[i])


def num(s):
    return R.num(s)


def sign(v, d=1):
    return f"{v:+.{d}f}".replace("+0.0", "0.0").replace("-0.0", "0.0")


# ------------------------------------------------------------------ one 16S dataset
def dataset_page(S, i, k):
    d, row = S.D[i], S.row[i]
    t1, t3 = d["t1"], d["t3"]
    t4, t5, t6, t11, rule = d.get("t4"), d.get("t5"), d.get("t6"), d["t11"], d.get("rule")
    cls, label, desc = R.role_plain(t1["Role"])
    extra = re.search(r"\((.*)\)", t1["Role"])
    role = desc + (f" ({extra.group(1)})" if extra else "")
    chips = I.chip("Mock community" if t3["Kind"] == "mock" else "Real samples", "measured") + I.chip(label, cls) + I.chip(row["marker_region"], "set")
    acc = f'<a href="{ENA}{row["accession"]}">{esc(row["accession"])}</a>' if row["accession"].startswith("PRJ") else esc(row["accession"])         # the master list of names gives the accession
    facts = (f'<dl class="facts"><dt>Marker, region, platform</dt><dd>{esc(fig(t1["Marker, region, platform"]))}</dd><dt>Samples</dt><dd>{esc(fig(t1["Samples"]))}</dd>'
             f'<dt>Known truth</dt><dd>{esc(fig(R.VALUE_PLAIN.get(t1["Truth"], t1["Truth"])))}</dd><dt>How we used it</dt><dd>{esc(fig(role))}</dd><dt>Data</dt><dd>{acc}</dd>'
             f'<dt>Source</dt><dd>{esc(row["source"])}</dd><dt>In the archive</dt><dd>{I.mono(row["data_folder"])}</dd></dl>')
    pass_default = t11["pass_default"]
    rk = "the rule's maxEE"
    tl = [("Pass-rate estimate off by", f"{t3['Filter error']} points", "AmpliTrim's estimate of the share of pairs passing DADA2's filter, minus DADA2's own (target: within 2 points)"),
          ("Pass at maxEE 2", f"{pass_default}%", "of the read pairs pass DADA2's default filter")]
    if t4:
        tl.append(("AmpliTrim against FIGARO", f"{t4['Difference'].split(' ')[0]} points", "pairs merged with AmpliTrim's cut minus FIGARO's cut, at the same maxEE"))
    else:
        tl.append(("AmpliTrim against FIGARO", "no comparison", fig(R.NO_COMPARISON.get(t3["the other tool error"], "FIGARO gave no comparison on this dataset."))))
    tl.append(("maxEE the rule sets", rule[rk] if rule else "2", f"merging {rule['the rule']} of the pairs (default maxEE 2: {rule['maxEE 2']})" if rule else "the rule keeps DADA2's default of 2 here"))
    secs = [I.section("How good were AmpliTrim's estimates?", I.table(["Samples", "Pass-rate estimate off by (points)", "Merge estimate off by (points; mean and range)", "Same order as DADA2 (1.00 = identical)", "FIGARO's pass-rate estimate off by (points)"],
            [[esc(t3["Samples"]), esc(t3["Filter error"]), esc(t3["Yield excess"]), esc(t3["Rank correlation"]), esc(R.VALUE_PLAIN.get(t3["the other tool error"], "not run" if t3["the other tool error"] in R.NO_COMPARISON else t3["the other tool error"]))]], num_cols=range(0, 5))
            + '<p class="legend">The merge estimate stops before DADA2 denoises the reads, so it is too high: use it to compare settings, not as a forecast.</p>', "estimates")]
    if t4:
        cutk = [x for x in t4 if x.startswith("the other tool")][0]
        secs.append(I.section("AmpliTrim against FIGARO at the same maxEE", '<p class="legend">FIGARO is an automatic tool for choosing truncation lengths. A cut such as 247/199 means the forward reads are cut to 247 nt and the reverse reads to 199 nt.</p>'
                              + I.table(["Samples", "FIGARO's cut and maxEE", "AmpliTrim's cut", "Difference in pairs merged (points)"], [[esc(t4["Samples"]), esc(t4[cutk]), esc(t4["AmpliTrim: cut"]), esc(t4["Difference"])]], num_cols=(0, 3)), "figaro"))
    else:
        secs.append(I.section("AmpliTrim against FIGARO at the same maxEE", f"<p>{esc(fig(R.NO_COMPARISON.get(t3['the other tool error'], 'FIGARO gave no comparison on this dataset.')))}</p>", "figaro"))
    rows = []
    if rule:
        rows = [[I.cell("Default maxEE 2", "name"), esc(rule["maxEE 2"]), "\u2013", vs_fig(t5["A"]) if t5 else "\u2013"],
                [I.cell("maxEE 2,5 (looser on the reverse read)", "name"), esc(rule["2,5"]), esc(t11["G minus A"]), vs_fig(t11["G against the other tool"])],
                [I.cell(f"The rule (maxEE {esc(rule[rk])})", "name"), esc(rule["the rule"]), "\u2013", "\u2013"]]
    else:
        rows = [[I.cell("Default maxEE 2", "name"), "\u2013", "\u2013", vs_fig(t5["A"]) if t5 else "\u2013"],
                [I.cell("maxEE 2,5 (looser on the reverse read)", "name"), "\u2013", esc(t11["G minus A"]), vs_fig(t11["G against the other tool"])]]
    maxee = (I.table(["Setting", "Pairs merged by DADA2", "Change against the default (points)", "Against FIGARO's own cut and maxEE (points)"], rows, num_cols=(1, 2, 3))
             + '<p class="legend">maxEE is the most expected errors DADA2 allows in a read; reads over the limit are dropped, and DADA2\'s default is 2. The rule sets the limit for each read direction from your own reads, never below 2.</p>'
             + I.table(["Change in pairs merged by maxEE 2,5 (points)", "Quality of the extra reads kept", "Accuracy (needs expected sequences)"],
                       [[esc(t11["G minus A"]), esc(R.VALUE_PLAIN.get(t11["G noise"], t11["G noise"])), esc(R.VALUE_PLAIN.get(t11["G accuracy"], t11["G accuracy"]))]], num_cols=(0,)))
    secs.append(I.section("The maxEE limit", maxee, "maxee"))
    if t6:
        secs.append(I.section("Accuracy on this mock community", "<p>This dataset has expected sequences, so each setting can be scored.</p>" + I.table(["Setting", "Variants found", "Strains found", "Reads matching an expected sequence exactly"],
                              [[I.cell(esc(fig(R.ARM_PLAIN.get(r["Arm"], r["Arm"]))), "name"), esc(r["Variants"]), esc(r["Strains"]), esc(r["Reads in exact ASVs"])] for r in t6], num_cols=(1, 2, 3)), "accuracy"))
    code = R._CODE[i]
    rec = [r for r in S.its.recovery["rows"] if r["code"] == code]
    if rec:
        secs.append(I.section("Recovery of the known members", I.table(["Setting", "Samples", "Expected per sample", "Found per sample (at least 0.1% of the reads / any read)", "In any sample", "In all samples"],
                              [[I.cell(esc(N.setting(r["arm"])), "name"), r["samples"], r["expected"], esc(r["per_sample"] + " / " + r["any_read"]), esc(r["any_sample"]), esc(r["all_samples"])] for r in rec], num_cols=(1, 2))
                              + '<p>What the words mean: <a href="../../../definitions/">definitions</a>. A strain is found when one of its variants is found exactly.</p>', "recovery"))
    secs.append(I.section("About this dataset", f"<p>{esc(fig(d['notes']))}</p>", "about"))
    secs.append(I.section("Where the numbers come from", f'<p><a href="{BLOB}docs/FINDINGS.md">docs/FINDINGS.md</a> (the tables of sections 4.1 to 4.10) and the results folders it names, in <a href="{BLOB}results">results</a>.</p>', "sources"))
    prev = f'<a href="../{S.slug(IDS16[k - 1])}/">&larr; {esc(S.row[IDS16[k - 1]]["short_name"])}</a>' if k > 0 else ""
    nxt = f'<a href="../{S.slug(IDS16[k + 1])}/">{esc(S.row[IDS16[k + 1]]["short_name"])} &rarr;</a>' if k < len(IDS16) - 1 else ""
    pager = f'<nav class="pager" aria-label="Datasets"><span>{prev}</span><span>{nxt}</span></nav>'
    body = I.hero(f'<a href="../../">16S</a> &rsaquo; Dataset {k + 1} of {len(IDS16)}', f"{row['id']} {row['full_name']}", "", f'<p class="chips">{chips}</p>{facts}{I.tiles(tl)}') + "\n".join(secs) + f'<section class="wrap">{pager}</section>'
    return page(f"{row['full_name']}: 16S results of the AmpliTrim validation", f"How AmpliTrim did on {row['full_name']}.", body, 3, "16s/")


# ------------------------------------------------------------------ the 16S overview
def error_chart(S):
    rows = []
    for i in IDS16:
        t3 = S.D[i]["t3"]
        rows.append((S.row[i]["short_name"], num(t3["Filter error"]), num(t3["the other tool error"])))
    lo, hi, h = -4, 8, 70 + 26 * len(rows)
    x = lambda v: 250 + (v - lo) / (hi - lo) * 250
    b = [f'<rect class="band2" x="{x(-2):.1f}" y="52" width="{x(2) - x(-2):.1f}" height="{h - 70}"/>']
    for t in range(-4, 9, 2):
        b.append(f'<line class="ink" x1="{x(t):.0f}" y1="{h - 22}" x2="{x(t):.0f}" y2="{h - 16}"/><text class="tm" x="{x(t):.0f}" y="{h - 4}" text-anchor="middle">{t:+d}'.replace("+0", "0") + "</text>")
    b.append(f'<line class="ink" x1="{x(0):.0f}" y1="48" x2="{x(0):.0f}" y2="{h - 22}"/>')
    for j, (name, a, f_) in enumerate(rows):
        y = 70 + 26 * j
        b.append(f'<text class="tm" x="244" y="{y + 4}" text-anchor="end">{esc(name)}</text>')
        if f_ is not None:
            b.append(f'<line class="ink" x1="{x(a):.1f}" y1="{y}" x2="{x(f_):.1f}" y2="{y}"/><circle class="warn" data-v="{f_}" cx="{x(f_):.1f}" cy="{y}" r="5"/>')
        b.append(f'<circle class="acc" data-v="{a}" cx="{x(a):.1f}" cy="{y}" r="6"/>')
    label = "How far each tool's estimate of the pass rate was from DADA2's measured pass rate, per 16S dataset: AmpliTrim and FIGARO"
    return f'<svg viewBox="0 0 520 {h}" role="img" aria-label="{esc(label)}"><title>{esc(label)}</title>{"".join(b)}</svg>'


def overview16(S):
    errs = [abs(num(S.D[i]["t3"]["Filter error"])) for i in IDS16]
    diffs = [num(S.D[i]["t4"]["Difference"]) for i in IDS16 if S.D[i].get("t4")]
    rules = [S.D[i]["rule"]["the rule's maxEE"] for i in IDS16 if S.D[i].get("rule")]
    low = [num(S.D[i]["t11"]["G minus A"]) for i in IDS16 if S.D[i]["t11"]["pass_default"] < 70]
    worst_v, worst_i = min((num(S.D[i]["t11"]["G minus A"]), i) for i in IDS16)
    tl = [("Pass-rate estimate", f"within {max(errs):.1f} points", f"of DADA2's measured pass rate on all {len(IDS16)} 16S datasets (target: within 2 points)"),
          ("AmpliTrim against FIGARO", f"{sign(min(diffs))} to {sign(max(diffs))} points", f"pairs merged with AmpliTrim's cut minus FIGARO's, at the same maxEE, on the {len(diffs)} datasets where FIGARO ran"),
          ("The maxEE rule", f"raises it on {len(rules)}", f"of {len(IDS16)} datasets (to {', '.join(sorted(set(rules)))}); it keeps DADA2's default of 2 elsewhere"),
          ("A looser limit (maxEE 2,5)", f"{sign(min(low))} to {sign(max(low))} points", f"more pairs merged on the {len(low)} datasets where the default passed under 70%; on {S.row[worst_i]['short_name']}, where {S.D[worst_i]['t11']['pass_default']}% passed, it merged {abs(worst_v):.1f} points fewer")]
    rows = []
    for k, i in enumerate(IDS16):
        d, row = S.D[i], S.row[i]
        t3, t4, rule = d["t3"], d.get("t4"), d.get("rule")
        rows.append([I.cell(f'<a href="datasets/{S.slug(i)}/">{esc(row["id"] + " " + row["short_name"])}</a>', "name"), esc("Mock" if t3["Kind"] == "mock" else "Real"), esc(t3["Samples"]), f'{d["t11"]["pass_default"]}%',
                     esc(t3["Filter error"]), esc(t3["Yield excess"].split(" ")[0]), esc(t3["Rank correlation"]), esc(t4["Difference"].split(" ")[0]) if t4 else "\u2013", esc(rule["the rule's maxEE"]) if rule else "2"])
    body = (I.hero('<a href="../">Results</a> &rsaquo; 16S', "16S: the cut and the maxEE limit",
                   "For 16S reads AmpliTrim recommends where to cut each read and how many expected errors to allow (maxEE). This page shows how its estimates compared with what DADA2 measured, and how its cut compared with FIGARO's, on every 16S dataset.", I.tiles(tl))
            + I.section("Dataset by dataset", I.figure(error_chart(S), "Blue: how far AmpliTrim's estimate of the pass rate was from DADA2's measured pass rate, in points. Orange: how far FIGARO's was, where it ran. The shaded band is the target set in advance, within 2 points.")
                        + I.table(["Dataset", "Kind", "Samples", "Pass at maxEE 2", "Pass-rate estimate off by (points)", "Merge estimate too high by (points, mean)", "Same order as DADA2", "AmpliTrim minus FIGARO (points)", "maxEE the rule sets"], rows, num_cols=range(2, 9))
                        + '<p class="legend">Merge estimate too high: it stops before DADA2 denoises the reads, so use it to compare settings, not as a forecast. Same order: 1.00 means AmpliTrim ranks the settings as DADA2 does. A dash: FIGARO could not be run, or gave no result. Words: <a href="../definitions/">definitions</a>.</p>', "glance")
            + I.section("More on the 16S runs", '<ul><li><a href="../recovery/">Recovery of the known members</a></li><li><a href="../criteria/">Every criterion and its outcome</a></li><li><a href="../runs/">Every run</a></li>'
                        f'<li><a href="{BLOB}docs/FINDINGS.md">The full findings</a></li></ul>', "more"))
    return page("16S: the cut and the maxEE limit: AmpliTrim validation results", "How AmpliTrim's estimates and cuts compared with DADA2 and FIGARO on every 16S dataset.", body, 1, "16s/")


# ------------------------------------------------------------------ the overview of both markers
def overview(S):
    errs = [abs(num(S.D[i]["t3"]["Filter error"])) for i in IDS16]
    its = S.its
    gains = [(float(its.one(c, "sw15")["yield"]) - float(its.one(c, "std")["yield"])) * 100 for c in I.TRUTH]
    drops = [(float(its.one(c, "sw15")["recall"]) - float(its.one(c, "std")["recall"])) * 100 for c in I.TRUTH]
    best = next(t for t in its.trims if t["best"])
    t16 = [("Pass-rate estimate", f"within {max(errs):.1f} points", f"of DADA2's measured pass rate on all {len(IDS16)} 16S datasets")]
    t_its = [("Reads kept by the trim", f"{sign(min(gains))} to {sign(max(gains))} points", f"on all {len(I.TRUTH)} ITS datasets with a known answer"),
             ("Largest loss of members found", f"{abs(min(drops)):.1f} points", f"the best of twelve trims: {N.setting(best['code']).replace('Trim: ', '').lower()}")]
    body = (I.hero("Results", "Test results of AmpliTrim",
                   "AmpliTrim recommends the settings of DADA2's filter from your own reads. It was tested against real DADA2 runs on public datasets: 16S amplicon data, and ITS data from fungi, including mock communities with a known answer. "
                   "Each marker has its own section with the same layout: an overview, a page for every dataset, and the settings compared.")
            + I.section("16S", f'<p>The cut and the maxEE limit, on {len(IDS16)} datasets.</p>' + I.tiles(t16) + '<p><a class="btn primary" href="16s/">16S results</a></p>', "s16")
            + I.section("ITS", f'<p>Cutting the poor ends of the reads, on {len(its.its)} datasets, {len(I.TRUTH)} of them with a known answer.</p>' + I.tiles(t_its) + '<p><a class="btn primary" href="its/">ITS results</a></p>', "its")
            + I.section("Across both markers", '<ul><li><a href="datasets/">All 24 datasets, with one ID and one name each</a></li><li><a href="recovery/">Recovery of the known members</a></li><li><a href="criteria/">Every criterion and its outcome</a></li>'
                        '<li><a href="definitions/">Definitions</a></li><li><a href="runs/">Every run</a></li><li><a href="how-it-works/">How AmpliTrim works</a></li></ul>', "both"))
    return page("Test results of AmpliTrim", "How AmpliTrim did against real DADA2 runs on 16S and ITS datasets.", body, 0, "")


# ------------------------------------------------------------------ the addresses of the first version
def moved(S, old, target_rel, new_name):
    meta = f'<meta http-equiv="refresh" content="0; url=../../{target_rel}">'
    body = I.hero('<a href="../../">Results</a> &rsaquo; Moved', "This page has moved", f'The results for {esc(new_name)} are now at <a href="../../{target_rel}">{esc(target_rel)}</a>.')
    return page("This page has moved: AmpliTrim validation results", "The page of this dataset has moved.", body, 2, "datasets/", meta)


def build(root):
    S = Site(root)
    out = {"index.html": overview(S), "16s/index.html": overview16(S)}
    for k, i in enumerate(IDS16):
        out[f"16s/datasets/{S.slug(i)}/index.html"] = dataset_page(S, i, k)
    for i in R.IDS:
        target = f"its/datasets/{S.slug(i)}/" if i == "its" else f"16s/datasets/{S.slug(i)}/"
        out[f"datasets/{i}/index.html"] = moved(S, i, target, S.row[i]["full_name"])
    return out


def main(root=None):
    root = pathlib.Path(root or ROOT).resolve()
    pages = build(root)
    for rel, text in pages.items():
        p = root / "docs" / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    (root / "docs" / ".nojekyll").write_text("", encoding="utf-8")
    print("written:", len(pages), "pages")
    return pages


if __name__ == "__main__":
    main()
