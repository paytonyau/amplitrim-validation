"""The ITS pages and the pages that serve both markers (docs/its/, docs/recovery/, docs/criteria/, docs/definitions/, docs/runs/, docs/datasets/index.html) are true to the files the
runs committed: they are exactly what site-build/build_its.py makes from those files, they are well formed and self-contained, every link goes somewhere, no archive code stands
in for a name, and the values and counts agree with the sources."""
import re
import sys
from html.parser import HTMLParser
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "site-build"))
import build_its as B          # noqa: E402
import names as N              # noqa: E402

DOCS = ROOT / "docs"
D = B.Data(ROOT)
PAGES = B.build(ROOT)
VOID = {"meta", "link", "img", "br", "hr", "input"}
INTERNAL = re.compile(r"\b(sw10|sw15|sw20|sw25|sweer|fwd15|mfirst|ee23|ee24|ee25|ee33|ee55|tq\d+|ml\d+|mo\d+|atcc_its\d|arfken_\w+|bakker_\w+|pauvert_its1|mockA_high|mockB|mcm_v\d+|d6300_V\d-V\d|miseqsop_mock|spikein_V3-V4|swabs_V3-V4|wastewater_V3-V4)\b")


class P(HTMLParser):
    def __init__(self, text):
        super().__init__()
        self.stack, self.errors, self.headings, self.ids, self.links, self.scripts, self.svgs, self.figs, self.lang = [], [], [], [], [], [], [], [], None
        self.in_code = self.in_style = self.in_svg = 0
        self.text, self._cap, self._in_fig, self.title = [], None, False, ""
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "html":
            self.lang = a.get("lang")
        if "id" in a:
            self.ids.append(a["id"])
        if tag == "a":
            self.links.append(a.get("href", ""))
        if tag in ("h1", "h2", "h3", "h4"):
            self.headings.append(int(tag[1]))
        if tag == "script":
            self.scripts.append(a)
        if tag == "svg":
            self.svgs.append({"role": a.get("role"), "label": a.get("aria-label")})
            self.in_svg += 1
        if tag == "figure":
            self._in_fig = True
            self.figs.append({"svg": False, "caption": ""})
        if tag == "figcaption":
            self._cap = []
        if tag == "code":
            self.in_code += 1
        if tag in ("style", "script"):
            self.in_style += 1
        if self._in_fig and tag == "svg":
            self.figs[-1]["svg"] = True
        if tag not in VOID:
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if not self.stack or self.stack[-1] != tag:
            self.errors.append((tag, list(self.stack[-3:])))
            return
        self.stack.pop()
        if tag == "svg":
            self.in_svg -= 1
        if tag == "code":
            self.in_code -= 1
        if tag in ("style", "script"):
            self.in_style -= 1
        if tag == "figcaption":
            self.figs[-1]["caption"] = "".join(self._cap)
            self._cap = None
        if tag == "figure":
            self._in_fig = False

    def handle_data(self, data):
        if self._cap is not None:
            self._cap.append(data)
        if not (self.in_code or self.in_style or self.in_svg):
            self.text.append(data)


PARSED = {rel: P(text) for rel, text in PAGES.items()}


def body(rel):
    return " ".join(" ".join(PARSED[rel].text).split())


def rows_of(html_text, section_id):
    seg = html_text.split(f'id="{section_id}"')[1].split("</section>")[0]
    return re.findall(r"<tr>(.*?)</tr>", seg.split("<tbody>")[1].split("</tbody>")[0], flags=re.S)


# ---------------------------------------------------------------- the pages are what the builder makes
def test_the_pages_on_disk_are_exactly_what_the_builder_makes_from_the_committed_files():
    assert len(PAGES) == 8 + len(D.its)
    for rel, text in PAGES.items():
        assert (DOCS / rel).read_text(encoding="utf-8") == text, f"docs/{rel} is out of date: run site-build/build_its.py"


def test_every_its_dataset_has_a_page_and_the_catalogue_lists_all_datasets_once():
    assert len(D.its) == 12 and len(D.names) == 24
    ids = [r["id"] for r in D.names]
    assert ids == ["DS%02d" % n if n < 23 else s for n, s in zip(range(1, 25), ids)] or len(set(ids)) == 24
    assert "DS13" not in ids and {"DS23a", "DS23b"} <= set(ids)
    cat = PAGES["datasets/index.html"]
    for r in D.names:
        assert r["full_name"].replace("&", "&amp;").replace("'", "&#x27;") in cat, r["id"]
        if N.is_its(r):
            assert f'href="../its/datasets/{N.slug(r)}/"' in cat


# ---------------------------------------------------------------- structure, links, accessibility
@pytest.mark.parametrize("rel", sorted(PAGES))
def test_each_page_is_well_formed_self_contained_and_accessible(rel):
    p, raw = PARSED[rel], PAGES[rel]
    assert not p.errors and not p.stack, (rel, p.errors, p.stack)
    assert p.lang == "en" and '<meta name="viewport"' in raw and "<title>" in raw
    scripts = re.findall(r"<script>(.*?)</script>", raw, flags=re.S)
    assert len(scripts) == 2 and all("amplitrim-theme" in x for x in scripts) and not any("src" in s for s in p.scripts)       # only the two inline theme scripts
    assert p.headings.count(1) == 1 and all(b - a <= 1 for a, b in zip(p.headings, p.headings[1:])), (rel, p.headings)
    assert len(p.ids) == len(set(p.ids)) and p.links[0] == "#main"
    assert "http://" not in raw.replace("http://www.w3.org/2000/svg", "")
    for s in p.svgs:
        assert s["role"] == "img" and s["label"], (rel, s)
    assert len(p.figs) == len(p.svgs) and all(f["svg"] and len(f["caption"].strip()) > 40 for f in p.figs), rel


def test_every_link_goes_somewhere_that_exists():
    for rel, p in PARSED.items():
        for h in p.links:
            if h.startswith("#"):
                assert h[1:] in p.ids, (rel, h)
            elif h.startswith("https://github.com/paytonyau/amplitrim-validation/blob/main/"):
                assert (ROOT / h.split("/blob/main/")[1].split("#")[0]).exists(), (rel, h)
            elif not h.startswith("https://"):
                target = ((DOCS / rel).parent / h.split("#")[0]).resolve()
                if h.endswith("/") or target.is_dir():
                    target = target / "index.html"
                assert target.exists() and DOCS.resolve() in target.parents, (rel, h)
            if "#run-" in h:
                assert f'id="run-{h.split("#run-")[1]}"' in PAGES["runs/index.html"], (rel, h)


def test_the_old_dataset_pages_and_the_overview_reach_the_new_pages():
    for name in ("index.html", "datasets/mock-a/index.html", "datasets/its/index.html"):
        raw = (DOCS / name).read_text(encoding="utf-8")
        depth = name.count("/")
        base = "./" if depth == 0 else "../" * depth
        for tail in ("its/", "recovery/", "criteria/", "definitions/", "datasets/"):
            assert f'href="{base}{tail}"' in raw, (name, tail)


# ---------------------------------------------------------------- names, not codes
@pytest.mark.parametrize("rel", sorted(PAGES))
def test_no_archive_code_stands_for_a_name(rel):
    found = set(INTERNAL.findall(" ".join(PARSED[rel].text)))
    assert not found, (rel, found)          # codes may appear inside <code> file paths and the definitions' settings table, never as a name in the text


def test_the_pages_use_the_new_names_and_ids():
    t = body("datasets/index.html")
    for word in ("DS01", "H\u00fclp\u00fcsch", "Tedersoo", "Pauvert", "DS23a", "DS23b"):
        assert word in t, word
    assert "mock community A" not in t.lower() and "second ITS" not in t


@pytest.mark.parametrize("rel", sorted(PAGES))
def test_the_old_study_labels_do_not_appear(rel):
    t = body(rel)
    for phrase in ("second ITS", "first ITS", "the ITS dataset", "Mock community A", "Mock community B"):
        assert phrase not in t, (rel, phrase)


# ---------------------------------------------------------------- values and counts agree with the sources
def test_known_values_of_the_its_overview_and_a_dataset_page():
    t = body("its/index.html")
    assert "DS21 ATCC Standard ITS2 | Mock | 26.9% | 25.7% | 83.0% | +57.3 | 66.7% | 94.4% | 0.00% | 0.00%".replace(" | ", " ") in t
    assert "DS24 Pauvert Mock Mock 78.6% 73.1% 83.6% +10.5 57.8% 56.7% 2.36% 2.35%" in t
    assert "+5.1 to +57.3 points" in t and "Window 4, Quality 15" in t
    page = body("its/datasets/ds21-atcc-standard-its2/index.html")
    assert "25.7% \u2192 83.0%" in page and "66.7% \u2192 94.4%" in page


def test_each_dataset_page_lists_every_row_of_its_sources():
    for r in D.its:
        c = r["archive_code"]
        rows = rows_of(PAGES[f"its/datasets/{N.slug(r)}/index.html"], "settings")
        assert len(rows) == len(D.rows(c)), (c, len(rows), len(D.rows(c)))
        for src, html_row in zip(sorted(D.rows(c), key=lambda x: (x["run"], N.setting_rank(x["arm"]))), rows):
            cells = [re.sub(r"<[^>]+>", "", x) for x in re.findall(r"<td[^>]*>(.*?)</td>", html_row, flags=re.S)]
            assert cells[0] == str(src["run"]) and cells[2] == f"{int(float(src['input'])):,}" and cells[3] == B.pc(src["pass"]) and cells[5] == B.pc(src["yield"]), (c, src["arm"])


def test_the_catalogue_counts_the_runs_and_the_criteria_agree_with_the_files():
    n = 0
    for run, rel in B.CRITERIA_FILES:
        n += sum(1 for line in (ROOT / "results" / rel).read_text(encoding="utf-8").splitlines() if re.match(r"^- \*\*P\d+", line) or re.match(r"^\| P\d+ \|", line))
    assert n == len(D.criteria) == 163                                                                                       # the same 163 as the workbook's All criteria sheet
    rows = sum(len(rows_of(PAGES["criteria/index.html"], f"run-{r}")) for r in {c["run"] for c in D.criteria})
    assert rows == n
    assert len(rows_of(PAGES["recovery/index.html"], "table")) == len(D.recovery["rows"]) == 21
    assert len(rows_of(PAGES["runs/index.html"], "runs")) == len(D.runs)
    assert len(D.trims) == 12 and sum(t["best"] for t in D.trims) == 1


def test_the_best_trim_and_the_safe_trims_are_those_of_the_summary():
    best = [t for t in D.trims if t["best"]][0]
    assert best["code"] == "sw15" and {t["code"] for t in D.trims if t["safe"]} >= {"sw15", "w5", "w8", "w10", "w12", "w16", "w20", "q12"}
    assert not {t["code"] for t in D.trims if t["safe"]} & {"w2", "w3", "w6", "q18"}


def test_the_predictions_are_those_of_the_files():
    assert len(D.predictions) > 200
    assert max(abs(p["predicted"] - p["dada2"]) for p in D.predictions) * 100 < 1.0                      # P126: within 1 point at every setting


def test_run_28_is_marked_as_set_aside_and_the_not_decided_criteria_are_kept():
    t = body("criteria/index.html")
    assert "This run was set aside" in t and "not decided" in t
