"""The results site (docs/index.html, docs/datasets/*/index.html, docs/reports/*.md) is true to the findings: every row it shows is a row of the findings, every number is in the sources, the charts carry the values of their tables, and the pages are accessible and self-contained. No network."""
import re
import shutil
import subprocess
import sys
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "site-build"))
import results_data as R                                                                           # noqa: E402

DOCS = ROOT / "docs"
FINDINGS = (DOCS / "FINDINGS.md").read_text(encoding="utf-8")
VALIDATION = (DOCS / "validation.md").read_text(encoding="utf-8")
D, META = R.load(ROOT)
ORDER = R.IDS
norm = lambda s: " ".join(s.split())
VOID = {"meta", "link", "img", "br", "hr", "input"}
BLOCKS = {"p", "li", "ul", "ol", "h1", "h2", "h3", "h4", "tr", "th", "td", "div", "section", "header", "footer", "nav", "main", "pre", "summary", "details", "table", "svg", "text", "title", "figure", "figcaption", "dd", "dt", "dl"}
PAGES = [("overview", DOCS / "index.html")] + [(i, DOCS / "datasets" / i / "index.html") for i in ORDER]


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.stack, self.errors, self.headings, self.ids, self.links, self.imgs, self.scripts, self.sheets = [], [], [], [], [], [], [], []
        self.text, self.svgs, self.figures, self.datav, self.lang, self.title, self._in_body, self._skip, self._link = [], [], [], [], None, "", False, 0, None
        self._svg = None; self._fig = None
        self.feed(text); self.close()

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag not in VOID: self.stack.append(tag)
        if "id" in a: self.ids.append(a["id"])
        if tag == "html": self.lang = a.get("lang")
        if tag == "body": self._in_body = True
        if tag in BLOCKS and self._in_body: self.text.append(" ")
        if tag in ("script", "style"): self._skip += 1
        if tag == "script": self.scripts.append(a)
        if tag == "link": self.sheets.append(a)
        if tag == "img": self.imgs.append(a)
        if re.fullmatch(r"h[1-6]", tag): self.headings.append(int(tag[1]))
        if tag == "a": self.links.append({"href": a.get("href", ""), "text": ""}); self._link = self.links[-1]
        if tag == "svg": self._svg = {"role": a.get("role"), "label": a.get("aria-label"), "title": ""}; self.svgs.append(self._svg)
        if tag == "figure": self._fig = {"svg": False, "caption": ""}; self.figures.append(self._fig)
        if "data-v" in a: self.datav.append(a["data-v"])

    def handle_endtag(self, tag):
        if tag in VOID: return
        if tag in BLOCKS and self._in_body: self.text.append(" ")
        if not self.stack or self.stack[-1] != tag: self.errors.append(f"</{tag}> closes {self.stack[-1:] or 'nothing'}")
        if self.stack and tag in self.stack:
            while self.stack and self.stack.pop() != tag: pass
        if tag in ("script", "style"): self._skip -= 1
        if tag == "a": self._link = None
        if tag == "svg": self._svg = None
        if tag == "figure": self._fig = None

    def handle_data(self, data):
        if self._skip: return
        if not self._in_body:
            if self.stack and self.stack[-1] == "title": self.title += data
            return
        self.text.append(data)
        if self._link is not None: self._link["text"] += data
        if self._svg is not None and self.stack and self.stack[-1] == "title": self._svg["title"] += data
        if self._fig is not None:
            if self._svg is not None: self._fig["svg"] = True
            if self.stack and self.stack[-1] == "figcaption": self._fig["caption"] += data

    @property
    def body_text(self): return norm("".join(self.text))


PARSED = {name: Page(path.read_text(encoding="utf-8")) for name, path in PAGES}
strip_commas = lambda s: re.sub(r"(?<=\d),(?=\d{3}\b)", "", s)
NUMBER = re.compile(r"\d+(?:\.\d+)*")
clean = lambda s: R.anon(re.sub(r"\\\[(\d+)\\\]", r"[\1]", s).replace("\\_", "_").replace("&#91;", "["))


def test_the_data_covers_the_datasets_and_tables_as_the_findings_do():
    assert ORDER == [d[0] for d in R.DATASETS] and len(ORDER) == 12
    want = {"t1": 12, "t3": 12, "t4": 10, "t5": 10, "t6": 4, "t7": 6, "t8": 6, "t9": 6, "t11": 12, "rule": 5}
    for slot, n in want.items():
        assert sum(1 for i in ORDER if slot in D[i]) == n, (slot, n)
    assert sum(len(D[i]["t6"]) for i in ORDER if "t6" in D[i]) == 12
    st = R.groups(D)                                                                                # the grouping the README shows
    assert (st["real"]["n"], st["real"]["ahead"], st["real"]["tied"], st["real"]["behind"]) == (4, 2, 2, 0) and round(st["real"]["mean"], 1) == 0.6
    assert (st["mock"]["n"], st["mock"]["ahead"], st["mock"]["tied"], st["mock"]["behind"]) == (5, 2, 3, 0) and round(st["mock"]["mean"], 1) == 0.6
    assert (st["wastewater"]["n"], st["wastewater"]["behind"]) == (1, 1) and st["wastewater"]["mean"] == -1.0
    assert (st["all"]["ahead"], st["all"]["tied"], st["all"]["behind"]) == (4, 5, 1) and round(st["all"]["mean"], 1) == 0.4


def test_every_row_of_the_data_is_a_row_of_the_findings_word_for_word():
    text = clean(FINDINGS); vtext = VALIDATION.replace("`", "")
    lines = {norm(l.strip().strip("|")) for l in text.splitlines() if l.startswith("|")}
    vlines = {norm(l.strip().strip("|")) for l in vtext.splitlines() if l.startswith("|")}
    checked = 0
    for i in ORDER:
        for slot in ("t1", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t11"):
            if slot not in D[i]: continue
            for rec in (D[i][slot] if slot == "t6" else [D[i][slot]]):
                row = norm(" | ".join(str(v) for k, v in rec.items() if k != "pass_default"))
                assert row in lines, (i, slot, row)                                                 # the extractor changed nothing
                checked += 1
        if "rule" in D[i]:
            row = norm(" | ".join(D[i]["rule"].values()).replace("2,5 |", "2,5 |"))
            assert any(row.split(" | ")[1:] == l.split(" | ")[1:] for l in vlines), (i, row)
            checked += 1
    assert checked == 12 + 12 + 10 + 10 + 12 + 6 + 6 + 6 + 12 + 5


def test_the_pages_and_reports_all_exist():
    assert (DOCS / ".nojekyll").exists()
    for i in ORDER:
        assert (DOCS / "datasets" / i / "index.html").exists() and (DOCS / "reports" / (i + ".md")).exists(), i


def test_the_pages_are_self_contained_well_formed_and_accessible():
    for name, path in PAGES:
        p = PARSED[name]; raw = path.read_text(encoding="utf-8")
        assert not p.errors and not p.stack, (name, p.errors, p.stack)
        assert p.lang == "en" and p.title.strip() and '<meta name="viewport"' in raw, name
        scripts = re.findall(r"<script>(.*?)</script>", raw, flags=re.S)
        assert not p.imgs and all("src" not in x for x in p.scripts) and len(scripts) == 2 and all("amplitrim-theme" in x for x in scripts), name   # only the two inline theme scripts
        assert all(l.get("rel") == "icon" and l["href"].startswith("data:") for l in p.sheets), name                                              # no external file
        assert p.headings.count(1) == 1 and all(b - a <= 1 for a, b in zip(p.headings, p.headings[1:])), (name, p.headings)
        assert len(p.ids) == len(set(p.ids)) and p.links[0]["href"] == "#main", name
        assert all(norm(l["text"]) for l in p.links) and all(l["href"][1:] in p.ids for l in p.links if l["href"].startswith("#")), name
        assert "http://" not in raw.replace("http://www.w3.org/2000/svg", ""), name
        for s in p.svgs: assert s["role"] == "img" and s["label"] and s["title"].strip(), (name, s)                                          # every chart is labelled
        assert len(p.figures) == len(p.svgs) and all(f["svg"] and len(f["caption"].strip()) > 40 for f in p.figures), name                        # and says how to read it


def test_every_link_goes_somewhere_that_exists():
    for name, path in PAGES:
        for l in PARSED[name].links:
            h = l["href"]
            if h.startswith("#"): continue
            if h.startswith("https://github.com/paytonyau/amplitrim-validation/blob/main/"):
                assert (ROOT / h.split("/blob/main/")[1].split("#")[0]).exists(), (name, h)
            elif not h.startswith("https://"):
                target = (path.parent / h.split("#")[0]).resolve()
                if h.endswith("/") or target.is_dir(): target = target / "index.html"
                assert target.exists() and DOCS.resolve() in target.parents, (name, h)
    ov = PARSED["overview"]
    assert all(any(l["href"] == f"datasets/{i}/" for l in ov.links) for i in ORDER)                      # the overview reaches every dataset
    for n, i in enumerate(ORDER):
        hrefs = [l["href"] for l in PARSED[i].links]
        prev_, next_ = (ORDER[n - 1] if n > 0 else None), (ORDER[n + 1] if n < 11 else None)
        assert [h for h in hrefs if re.fullmatch(r"\.\./[a-z0-9-]+/", h)] == [f"../{x}/" for x in (prev_, next_) if x], i      # previous and next, in order
        assert "../../" in hrefs


def test_no_page_or_report_names_the_other_tool_or_makes_a_claim_we_decided_not_to_make():
    for name, path in PAGES:
        raw = path.read_text(encoding="utf-8").lower().split("<style>")[0] + PARSED[name].body_text.lower()
        for word in ("figaro", "guarantee", "resolution", "fastest", "faster", "state-of-the-art", "outperform"):
            assert word not in raw, (name, word)
    for i in ORDER:
        assert "figaro" not in (DOCS / "reports" / (i + ".md")).read_text(encoding="utf-8").lower(), i


def test_every_number_on_the_pages_is_in_the_findings_or_the_validation_page():
    known = set(NUMBER.findall(strip_commas(FINDINGS + VALIDATION)))
    ticks = {"1", "2", "3", "4", "5", "10", "20", "30", "40", "60", "100"}                              # the labels of chart axes
    accessions = set(NUMBER.findall(" ".join(a for d in R.DATASETS for a in d[2])))                     # run accessions such as SRR16961582
    for name, path in PAGES:
        found = set(NUMBER.findall(strip_commas(PARSED[name].body_text)))
        extra = found - known - ticks - accessions
        assert not extra, (name, extra)
    for i in ORDER:
        for tok in NUMBER.findall(strip_commas(D[i]["notes"])): assert tok in known, (i, tok)           # and every number in the short notes


def test_the_charts_carry_the_values_of_their_tables():
    allowed = {}
    for i in ORDER:
        vals = set()
        for slot in ("t3", "t4", "t5", "t7", "t8", "t11", "rule"):
            if slot in D[i]:
                for v in D[i][slot].values():
                    vals.add(str(v)); vals.add(str(v).split(" ")[0]); vals.update(re.findall(r"[+-]?\d+(?:\.\d+)?%", str(v))); vals.update(re.findall(r"[+-]\d+\.\d", str(v)))
        allowed[i] = vals
    for i in ORDER:
        marks = PARSED[i].datav
        assert marks and all(m in allowed[i] for m in marks), (i, [m for m in marks if m not in allowed[i]])
    every = set().union(*[allowed[i] for i in ORDER])
    ov = PARSED["overview"].datav
    assert len(ov) == 12 + 10 and all(m in every for m in ov)                                           # twelve filter errors and ten differences
    assert Counter(ov[:12]) == Counter(D[i]["t3"]["Filter error"] for i in ORDER)                       # the strip chart shows Table 3, row by row
    assert sorted(ov[12:]) == sorted(D[i]["t4"]["Difference"].split(" ")[0] for i in ORDER if "t4" in D[i])


def test_each_dataset_page_says_what_the_data_say():
    for n, i in enumerate(ORDER):
        d, p = D[i], PARSED[i]; t = p.body_text
        label = R.role_plain(d["t1"]["Role"])[1]
        assert d["name"] in t and f"Dataset {n + 1} of 12" in t and d["t1"]["Samples"] in t and label in t, i
        assert d["notes"] in t and any(a in " ".join(l["href"] for l in p.links) for a in d["accessions"]), i
        assert d["t3"]["Filter error"] in t and d["t3"]["Yield excess"] in t and d["t3"]["Rank correlation"] in t
        if "t4" in d: assert d["t4"]["Difference"] in t and "Compared with the other tool, at the same maxEE" in t and d["t4"]["AmpliTrim: cut"] in t
        else: assert R.NO_COMPARISON[d["t3"]["the other tool error"]] in t                            # the reason there is no comparison, in words
        assert ("What the released rule does" in t) and (("The rule's maxEE" in t) == ("rule" in d)), i
        assert ("Would the choice change on other samples?" in t) == ("t9" in d) and ("Leaving out the longest amplicons" in t) == ("t8" in d) and ("Accuracy on this mock community" in t) == ("t6" in d), i
        assert d["t11"]["G minus A"] in t
        assert ("to add" in t) != bool(d["ref"]) or d["ref"]                                          # a missing reference says "to add"

def test_the_overview_states_the_outcomes_as_the_data_give_them():
    t = PARSED["overview"].body_text; dv = R.derived(D); st = dv["stats"]
    fmt1 = lambda v: f"{v:.1f}"
    assert f"within {fmt1(dv['max_abs_fe'])} points of what DADA2 measured" in t and f"The average miss was {fmt1(dv['mean_abs_fe'])} points" in t
    assert f"On {dv['n_rank_one']} of the {dv['n_rank']} datasets where the order could be compared" in t and "Mock community A (0.90)" in t and "Wastewater, V3-V4 (0.80)" in t
    assert f"merged more pairs on {st['all']['ahead']} of {st['all']['n']} datasets, the same on {st['all']['tied']} and fewer on {st['all']['behind']}" in t
    assert f"chose the same cut on {dv['hold_a']} of {dv['n_t9']} datasets and the other tool on {dv['hold_o']} of {dv['n_t9']}" in t
    assert "Sensible default, not independently tested" in t and "after seeing our test data" in t and "the author gives no ceiling" in t    # what the maxEE rule is, and what it is not
    assert "did not reach the result we had decided in advance would make us adopt it" in t                                         # the facts stay, in plain words
    assert "Wastewater*" in t and "the only real dataset in these tests with poor read quality" in t
    kinds = [R.role_plain(D[i]["t1"]["Role"])[0] for i in ORDER]
    assert f"{kinds.count('heldout')} datasets were kept aside" in t and f"{kinds.count('dev')} were used while designing AmpliTrim" in t
    for term in ("Estimate and measurement", "Points", "maxEE", "A cut", "The other tool", "The shaded band", "Colours", "Our datasets", "Rounds of tests"):
        assert term in t, term                                                                          # every term is explained once, up front

def test_the_reports_repeat_the_tables_of_the_data():
    for i in ORDER:
        rep = (DOCS / "reports" / (i + ".md")).read_text(encoding="utf-8"); d = D[i]
        assert rep.startswith("# " + d["name"]) and d["notes"] in rep and d["t3"]["Filter error"] in rep, i
        for slot in ("t4", "t7", "t8", "t9"):
            if slot in d: assert all(str(v) in rep for k, v in d[slot].items() if k not in ("Dataset", "pass_default") and not k.startswith("Dataset")), (i, slot)
        assert d["t11"]["G minus A"] in rep and ("t5" not in d or d["t5"]["A"] in rep), i
        if R.num(d["t11"]["G against the other tool"]) is not None: assert d["t11"]["G against the other tool"] in rep, i
        if "t6" in d: assert all(r["Reads in exact ASVs"] in rep for r in d["t6"]), i

def test_the_site_is_up_to_date_with_the_findings(tmp_path):
    """Rebuild into a copy: the pages and reports on disk must be exactly what the builder makes from the findings."""
    work = tmp_path / "archive"; (work / "docs").mkdir(parents=True)
    for f in ("FINDINGS.md", "validation.md"): shutil.copy2(DOCS / f, work / "docs" / f)
    shutil.copytree(ROOT / "site-build", work / "site-build")
    done = subprocess.run([sys.executable, str(work / "site-build" / "build_results.py"), str(work)], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-500:]
    for name, path in PAGES:
        assert (work / path.relative_to(ROOT)).read_text(encoding="utf-8") == path.read_text(encoding="utf-8"), name + " is out of date: run site-build/build_results.py"
    for i in ORDER:
        rp = Path("docs") / "reports" / (i + ".md")
        assert (work / rp).read_text(encoding="utf-8") == (ROOT / rp).read_text(encoding="utf-8"), i


BANNED_PHRASES = [r"\bnights?\b", r"\barms?\b", r"\bB2\b", r"\bB4\b", r"\bF2\b", r"\bG minus\b", r"\bheld[- ]out\b", r"\bdevelopment\b", r"pre-?registered", r"adoption criteri", r"rank correlation", r"yield excess", r"filter error", r"Table \d+:", r"role in the tests"]


def test_the_pages_and_reports_use_plain_words_and_not_our_internal_ones():
    texts = {name: PARSED[name].body_text for name, _ in PAGES}
    texts.update({i + ".md": (DOCS / "reports" / (i + ".md")).read_text(encoding="utf-8") for i in ORDER})
    texts.update({"note:" + i: D[i]["notes"] for i in ORDER})
    for name, text in texts.items():
        for pat in BANNED_PHRASES:
            assert not re.search(pat, text, flags=re.I), (name, pat, re.search(pat, text, flags=re.I).group(0))


def _charts(name):
    return re.findall(r"<svg\b[^>]*aria-label=\"([^\"]*)\"[^>]*>(.*?)</svg>", (DOCS / ("index.html" if name == "overview" else f"datasets/{name}/index.html")).read_text(encoding="utf-8"), flags=re.S)


def test_every_colour_means_one_thing_in_every_chart():
    """Orange only ever means fewer (a negative value), blue never does; a chart of errors, which is neither good nor bad, uses one colour."""
    seen = 0
    for name, _ in PAGES:
        for label, body in _charts(name):
            for cls, v in re.findall(r'<rect class="(bar|barn)" data-v="([^"]+)"', body):
                assert (cls == "barn") == (R.num(v) is not None and R.num(v) < 0 and "%" not in v), (name, label, cls, v); seen += 1
            for cls, v in re.findall(r'<circle class="(acc|warn)" data-v="([^"]+)"', body):
                signed = "pairs merged by DADA2 with AmpliTrim" in label
                assert cls == ("warn" if signed and R.num(v) < 0 else "acc"), (name, label, cls, v); seen += 1
    assert seen > 100                                                                                   # the patterns really found the marks


def test_every_chart_says_which_way_is_which_and_what_the_band_is():
    for name, _ in PAGES:
        for label, body in _charts(name):
            texts = " ".join(re.findall(r"<text[^>]*>([^<]*)</text>", body))
            if "share of pairs merged by DADA2" in label: continue                                       # bars from zero: nothing to point at
            assert "\u2190" in texts and "\u2192" in texts, (name, label)                                # an arrow at each end of the axis, with words
            if "effect of a looser maxEE" not in label: assert "shaded:" in texts, (name, label)         # and the shaded band is named


def test_the_themes_are_soft_readable_and_can_be_switched():
    def pairs(block): return dict(re.findall(r"--([a-z-]+):(#[0-9a-fA-F]{3,6})", block))
    def lum(h):
        h = h.lstrip("#"); h = h if len(h) == 6 else "".join(c * 2 for c in h)
        f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
        r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
        return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)
    def contrast(a, b):
        la, lb = sorted((lum(a), lum(b)), reverse=True)
        return (la + 0.05) / (lb + 0.05)
    for name, path in PAGES:
        raw = path.read_text(encoding="utf-8"); css = raw[raw.index("<style>"):raw.index("</style>")]
        light, auto, chosen = {}, {}, {}
        for b in re.findall(r"(?m)^:root\{([^}]*)\}", css): light.update(pairs(b))
        for b in re.findall(r'@media \(prefers-color-scheme:dark\)\{:root:not\(\[data-theme="light"\]\)\{([^}]*)\}\}', css): auto.update(pairs(b))
        for b in re.findall(r':root\[data-theme="dark"\]\{([^}]*)\}', css): chosen.update(pairs(b))
        assert auto and auto == chosen, name
        assert lum(light["bg"]) < 0.97 and lum(chosen["bg"]) > 0.012, name                             # neither pure white nor near-black
        for label, theme in (("light", light), ("dark", chosen)):
            for fg in ("fg", "muted", "accent", "warn", "ok", "rev"):
                for bg in ("bg", "card"): assert contrast(theme[fg], theme[bg]) >= 4.5, (name, label, fg, bg)
        head = re.findall(r"<script>(.*?)</script>", raw[:raw.index("<style>")], flags=re.S)
        assert any("localStorage.getItem('amplitrim-theme')" in x and "setAttribute('data-theme'" in x and x.count("{") == x.count("}") for x in head), name
        assert "theme-toggle" in raw and "prefers-color-scheme" in raw, name


def test_mock_b_gives_the_measured_result_of_the_released_rule_and_what_the_cap_would_have_cost_and_no_page_still_states_a_cap():
    """The released rule has no cap; on mock B, the one dataset where it matters, DADA2 was run at its setting (TEST_REGISTER.md, section 24)."""
    note, rule = R.CAP_NOTE["mock-b"], D["mock-b"]["rule"]
    page, report = PARSED["mock-b"].body_text, (DOCS / "reports" / "mock-b.md").read_text(encoding="utf-8")
    assert rule["the rule's maxEE"] == "4,8" and rule["the rule"] == "70.4%"                                                    # the measured result of the released rule is what the table holds
    for text in (page, report.replace("**", "")):                                                                                 # the report keeps the Markdown bold markers
        assert f"The rule has no cap. With a cap of 5 it would have set maxEE {note['capped_max_ee']} here, and DADA2 merged {note['capped_merged']} of the pairs; without one it sets maxEE 4,8 and DADA2 merged 70.4%." in text
        assert "has not been run" not in text and "not yet run" not in text and "Rule with a cap of 5" not in text
    for i in ORDER:
        if i != "mock-b": assert "The rule has no cap" not in PARSED[i].body_text, i                                              # the cap changed nothing anywhere else
    ov = PARSED["overview"].body_text
    assert "measured afterwards" in ov and "that test used an earlier version of the rule, with a cap of 5" in ov and "We did not apply the criteria again" in ov
    assert "has not been run on yet" not in ov and "has since been removed" not in ov
    for name, _ in PAGES:
        assert "never above 5" not in PARSED[name].body_text, name


def test_every_page_of_the_results_site_links_to_the_how_it_works_page_that_lives_here():
    target = (DOCS / "how-it-works" / "index.html").resolve()
    assert target.exists()
    for name, path in PAGES:
        hrefs = [l["href"] for l in PARSED[name].links if l["href"].endswith("how-it-works/")]
        assert hrefs and all(((path.parent / h) / "index.html").resolve() == target for h in hrefs), name       # the menu of every page points at the page of this site
