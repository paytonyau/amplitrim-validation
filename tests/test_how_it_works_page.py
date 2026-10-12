"""The how-it-works page of this archive's site (docs/how-it-works/index.html, moved here from the software repository): self-contained and accessible, every figure a labelled
schematic, the status of every choice stated (the maxEE rule included), soft light and dark themes with a switch. Its numbers and its module table are made from the code of the
software repository by site-build/build_design.py, and test_site_build.py rebuilds the page and compares it byte for byte, so they cannot drift. No network."""
import re
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HOW_FILE = ROOT / "docs" / "how-it-works" / "index.html"
SOFTWARE = ROOT.parent / "amplitrim"                                                                   # the software repository next to this archive (optional)
HOME, REPO = "https://paytonyau.github.io/amplitrim/", "https://github.com/paytonyau/amplitrim"
BLOB, ARCH_BLOB = REPO + "/blob/main/", "https://github.com/paytonyau/amplitrim-validation/blob/main/"
VOID = {"meta", "link", "img", "br", "hr", "input"}
BLOCKS = {"p", "li", "ul", "ol", "h1", "h2", "h3", "h4", "tr", "th", "td", "div", "section", "header", "footer", "nav", "main", "pre", "summary", "details", "table", "svg", "text", "title", "figure", "figcaption"}      # picture labels are separate pieces of text
norm = lambda s: " ".join(s.split())


class Page(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.stack, self.errors, self.headings, self.ids, self.links, self.imgs = [], [], [], [], [], []
        self.scripts, self.sheets, self.pres, self.tables, self.flow = [], [], [], [], []
        self.text, self.attrs_html, self.html_lang, self.title = [], 0, None, ""
        self._skip, self._cell, self._row, self._table, self._pre, self._flow = 0, None, None, None, None, False
        self._in_body, self._sr, self._link, self._step = False, 0, None, False
        self.feed(text)
        self.close()

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag not in VOID: self.stack.append(tag)
        if "id" in a: self.ids.append(a["id"])
        if tag == "html": self.html_lang = a.get("lang")
        if tag == "body": self._in_body = True
        if tag in BLOCKS and self._in_body: self.text.append(" ")                                          # block boundaries are spaces; inline tags add nothing
        if tag in ("script", "style"): self._skip += 1
        if tag == "script": self.scripts.append(a)
        if tag == "link": self.sheets.append(a)
        if tag == "img": self.imgs.append(a)
        if re.fullmatch(r"h[1-6]", tag): self.headings.append(int(tag[1]))
        if tag == "a":
            self.links.append({"href": a.get("href", ""), "text": ""}); self._link = self.links[-1]
        if tag == "img" and self._link is not None: self._link["text"] += a.get("alt", "")
        if {"sr", "num"} & set(a.get("class", "").split()): self._sr += 1                                 # screen-reader-only text and step badges are not prose
        if tag == "ol" and "flow" in a.get("class", "").split(): self._flow = True
        if tag == "span" and "step" in a.get("class", "").split() and self._flow: self.flow.append(""); self._step = True
        if tag == "table": self._table = []
        if tag == "tr" and self._table is not None: self._row = []
        if tag in ("th", "td") and self._row is not None: self._cell = []
        if tag == "br" and self._cell is not None: self._cell.append("\n")
        if tag == "pre": self._pre = []

    def handle_endtag(self, tag):
        if tag in VOID: return
        if tag in BLOCKS and self._in_body: self.text.append(" ")
        if not self.stack or self.stack[-1] != tag: self.errors.append(f"</{tag}> closes {self.stack[-1:] or 'nothing'}")
        if self.stack and tag in self.stack:
            while self.stack and self.stack.pop() != tag: pass
        if tag in ("script", "style"): self._skip -= 1
        if tag == "a": self._link = None
        if tag == "ol": self._flow = False
        if tag == "span" and self._sr: self._sr -= 1
        if tag == "span": self._step = False
        if tag in ("th", "td") and self._cell is not None:
            lines = tuple(x for x in (norm(p) for p in "".join(self._cell).split("\n")) if x); self._row.append(lines); self._cell = None
        if tag == "tr" and self._row is not None: self._table.append(tuple(self._row)); self._row = None
        if tag == "table": self.tables.extend(self._table); self._table = None
        if tag == "pre": self.pres.append("".join(self._pre)); self._pre = None

    def handle_data(self, data):
        if self._skip: return
        if not self._in_body:
            if self.stack and self.stack[-1] == "title": self.title += data
            return
        if self._sr: return
        self.text.append(data)
        if self._cell is not None: self._cell.append(data)
        if self._pre is not None: self._pre.append(data)
        if self._link is not None: self._link["text"] += data
        if self._step: self.flow[-1] += data

    @property
    def body_text(self): return norm("".join(self.text))



HOW = Page(HOW_FILE.read_text(encoding="utf-8"))


def test_the_page_is_self_contained_and_well_formed():
    raw = HOW_FILE.read_text(encoding="utf-8")
    assert not HOW.errors and not HOW.stack, (HOW.errors, HOW.stack)                                       # every tag opened is closed
    assert HOW.html_lang == "en" and HOW.title.strip() and '<meta name="viewport"' in raw
    assert all("src" not in x for x in HOW.scripts)                                                         # only inline scripts
    assert all(l.get("rel") == "icon" and l["href"].startswith("data:") for l in HOW.sheets)               # no external stylesheet or font
    assert not HOW.imgs                                                                                     # the figures are inline SVG, there is no image file
    assert not re.search(r"<(iframe|form|object|embed)\b", raw) and "http://" not in raw.replace("http://www.w3.org/2000/svg", "")


def test_the_page_is_accessible_in_the_basics():
    p = HOW
    assert p.headings.count(1) == 1                                                                         # one h1
    assert all(b - a <= 1 for a, b in zip(p.headings, p.headings[1:])), p.headings                          # no skipped heading level
    assert all(norm(l["text"]) for l in p.links), [l["href"] for l in p.links if not norm(l["text"])]
    assert len(p.ids) == len(set(p.ids)) and p.links[0]["href"] == "#main"                                  # unique ids; the skip link comes first
    assert all(l["href"][1:] in p.ids for l in p.links if l["href"].startswith("#"))                        # every in-page link has a target


def test_every_link_goes_somewhere_that_exists():
    for l in HOW.links:
        href = l["href"]
        if href.startswith("#"):
            continue
        if not href.startswith("https://"):                                                                 # a link inside this site: the results overview
            target = (HOW_FILE.parent / href.split("#")[0]).resolve()
            if href.endswith("/") or target.is_dir():
                target = target / "index.html"
            assert target.exists() and (ROOT / "docs").resolve() in target.parents, href
        elif href.startswith(ARCH_BLOB):
            assert (ROOT / href[len(ARCH_BLOB):].partition("#")[0]).exists(), href                           # a file of this archive
        elif href.startswith(BLOB):
            if SOFTWARE.exists():                                                                           # a file of the software repository, when it is next to this archive
                assert (SOFTWARE / href[len(BLOB):].partition("#")[0]).exists(), href
        else:
            assert href in (HOME, REPO, REPO + "/issues") or href.startswith("https://orcid.org/"), href    # nothing else leaves the site
    hrefs = [l["href"] for l in HOW.links]
    assert hrefs.count(HOME) == 3 and hrefs.count("../") == 2                                               # back to the software's home page and to the results: the logo, the menu and the footer; the menu and the footer


def test_the_page_makes_no_claim_we_decided_not_to_make():
    raw = HOW_FILE.read_text(encoding="utf-8").lower().split("<style>")[0]
    for word in ("FIGARO", "guarantee", "resolution", "fastest", "faster", "state-of-the-art", "outperform", "best-in-class"):
        assert word.lower() not in HOW.body_text.lower() and word.lower() not in raw, word


def test_every_figure_is_a_labelled_schematic():
    raw = HOW_FILE.read_text(encoding="utf-8")
    svgs = re.findall(r"<svg\b[^>]*>.*?</svg>", raw, flags=re.S)
    assert len(svgs) == 11, len(svgs)                                                                   # three cut cases, seven steps, and the ITS branch
    for svg in svgs:
        assert 'role="img"' in svg and 'aria-label="' in svg and re.search(r"<title>[^<]+</title>", svg)
    figures = re.findall(r"<figure>.*?</figure>", raw, flags=re.S)
    assert len(figures) == len(svgs) and all("<figcaption>Schematic, not your data.</figcaption>" in f for f in figures)   # never mistaken for results


def test_the_page_says_how_well_each_choice_held_up_including_the_maxee_rule():
    raw = HOW_FILE.read_text(encoding="utf-8")
    assert set(re.findall(r'class="chip (\w+)"', raw)) >= {"held", "design", "optimistic", "timing", "unvalidated", "exact", "none", "measured"}
    assert raw.count('<li class="step">') == 8                                                          # seven steps and the ITS branch
    assert "Sensible default, not independently tested" in HOW.body_text and "the author gives no ceiling" in HOW.body_text         # what the rule is, and what it is not
    assert "made after seeing the validation data" in HOW.body_text and "did not reach the result we had decided in advance would make us adopt it" in HOW.body_text      # the facts stay on the page


def theme_vars(raw):
    """The colour variables of the light theme, of the dark theme chosen by the system, and of the dark theme chosen with the switch."""
    css = raw[raw.index("<style>"):raw.index("</style>")]
    pairs = lambda block: dict(re.findall(r"--([a-z-]+):(#[0-9a-fA-F]{3,6})", block))
    light = {}
    for b in re.findall(r"(?m)^:root\{([^}]*)\}", css): light.update(pairs(b))
    auto, chosen = {}, {}
    for b in re.findall(r'@media \(prefers-color-scheme:dark\)\{:root:not\(\[data-theme="light"\]\)\{([^}]*)\}\}', css): auto.update(pairs(b))
    for b in re.findall(r':root\[data-theme="dark"\]\{([^}]*)\}', css): chosen.update(pairs(b))
    return light, auto, chosen

def luminance(h):
    h = h.lstrip("#"); h = h if len(h) == 6 else "".join(c * 2 for c in h)
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)

def contrast(a, b):
    la, lb = sorted((luminance(a), luminance(b)), reverse=True)
    return (la + 0.05) / (lb + 0.05)


def test_the_themes_are_soft_readable_and_can_be_switched():
    raw = HOW_FILE.read_text(encoding="utf-8")
    light, auto, chosen = theme_vars(raw)
    assert auto and auto == chosen                                                                          # the system's dark theme and the chosen one cannot differ
    assert luminance(light["bg"]) < 0.97 and luminance(chosen["bg"]) > 0.012                                # neither pure white nor near-black
    for label, theme in (("light", light), ("dark", chosen)):
        for fg in ("fg", "muted", "accent", "warn", "ok", "rev"):
            for bg in ("bg", "card"):
                assert contrast(theme[fg], theme[bg]) >= 4.5, (label, fg, bg, round(contrast(theme[fg], theme[bg]), 2))      # readable text and marks
        assert contrast(theme["accent-fg"], theme["accent"]) >= 4.5, label                                   # the text on a filled button
    assert "localStorage.getItem('amplitrim-theme')" in raw and "theme-toggle" in raw and "prefers-color-scheme" in raw
    head_scripts = re.findall(r"<script>(.*?)</script>", raw[:raw.index("<style>")], flags=re.S)
    assert any("localStorage.getItem('amplitrim-theme')" in x and "setAttribute('data-theme'" in x and x.count("{") == x.count("}") for x in head_scripts)     # a whole script before the first paint


def test_the_page_says_that_the_maxee_rule_has_no_upper_limit_and_notes_a_limit_above_5():
    found = re.search(r"WARN_ABOVE = (\d+)", (SOFTWARE / "amplitrim" / "pipeline.py").read_text(encoding="utf-8")) if SOFTWARE.exists() else None
    warn = found.group(1) if found else "5"                                                                 # read from the software repository when it is next to this archive
    assert "There is no upper limit" in HOW.body_text and "never above" not in HOW.body_text and f"the report notes it when a limit goes above {warn}" in HOW.body_text
