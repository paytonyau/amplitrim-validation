"""The header is the same on every page of both websites (the software's home page, the how-it-works page and the validation results): the same brand, the same three tabs with the
current one marked, a box for the theme button, and a second row of links that fits on one row. It is written by site-build/header.py only; this test fails if a builder writes its own."""
import os
import re
from pathlib import Path
from urllib.parse import urlparse

import pytest

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"
SOFTWARE_SITE, RESULTS_SITE = "https://paytonyau.github.io/amplitrim/", "https://paytonyau.github.io/amplitrim-validation/"
PAGES = {str(p.relative_to(DOCS).as_posix()): p for p in sorted(DOCS.rglob("index.html"))}
SOFT = Path(os.environ["AMPLITRIM_REPO"]) / "docs" / "index.html" if os.environ.get("AMPLITRIM_REPO") else None
if SOFT is not None and SOFT.is_file():
    PAGES["(software) index.html"] = SOFT


def parse(path):
    raw = path.read_text(encoding="utf-8")
    h = re.search(r'<header class="site">(.*?)</header>', raw, re.S).group(1)
    sites = re.search(r'<nav aria-label="Sites"><ul>(.*?)</ul></nav>', h, re.S).group(1)
    main = re.search(r'<nav aria-label="Main"><ul>(.*?)</ul></nav>', h, re.S).group(1)
    link = lambda block: [(re.sub(r"<[^>]+>", "", m.group(3)), m.group(1), "aria-current" in m.group(2)) for m in re.finditer(r'<a href="([^"]*)"([^>]*)>(.*?)</a>', block)]
    return {"brand": re.sub(r"<[^>]+>", "", re.search(r'<a class="brand"[^>]*>(.*?)</a>', h, re.S).group(1)), "tabs": link(sites), "row2": link(main), "tools": '<div class="tools"></div>' in h,
            "subbar": 'class="wrap subbar"' in h}


def family(name):
    return "software" if name.startswith("(software)") else "how" if name.startswith("how-it-works/") else "results"


@pytest.mark.parametrize("name", sorted(PAGES))
def test_the_first_row_is_the_same_on_every_page(name):
    p = parse(PAGES[name])
    assert p["brand"] == "AmpliTrim" and p["tools"] and p["subbar"], name
    assert [t[0] for t in p["tabs"]] == ["Software", "How it works", "Validation"], name
    assert [t[0] for t in p["tabs"] if t[2]] == [{"software": "Software", "how": "How it works", "results": "Validation"}[family(name)]], name       # exactly one, the right one


@pytest.mark.parametrize("name", sorted(PAGES))
def test_the_tabs_lead_to_the_same_three_places(name):
    p, tabs = parse(PAGES[name]), dict((t[0], t[1]) for t in parse(PAGES[name])["tabs"])
    def where(href):
        if href.startswith("https://") or href.startswith("#"):
            return RESULTS_SITE + "how-it-works/" if False else href
        return (PAGES[name].parent / href.split("#")[0]).resolve()
    if family(name) == "software":
        assert tabs["Software"] in ("#top", SOFTWARE_SITE) and tabs["How it works"] == RESULTS_SITE + "how-it-works/" and tabs["Validation"] == RESULTS_SITE
        return
    assert tabs["Software"] == SOFTWARE_SITE, name
    how = where(tabs["How it works"])
    if tabs["How it works"] == "#top":
        how = PAGES[name].resolve()
    else:
        how = (how / "index.html") if how.is_dir() else how
    assert how == (DOCS / "how-it-works" / "index.html").resolve(), name
    res = (PAGES[name].parent / tabs["Validation"].split("#")[0]).resolve()
    res = res / "index.html" if res.is_dir() else res
    assert res == (DOCS / "index.html").resolve(), name


@pytest.mark.parametrize("name", sorted(PAGES))
def test_the_second_row_has_a_few_short_links_and_marks_the_section(name):
    p = parse(PAGES[name])
    assert 3 <= len(p["row2"]) <= 9 and all(len(t[0]) <= 26 for t in p["row2"]), name           # a row that cannot overflow
    if family(name) == "results" and not name.startswith("how-it-works"):
        assert sum(t[2] for t in p["row2"]) == 1, name                                           # the section of this page is marked
        want = "Datasets" if name.startswith("datasets/") else "Overview" if name == "index.html" else {"16s": "16S", "its": "ITS", "recovery": "Recovery", "criteria": "Criteria", "definitions": "Definitions", "runs": "Runs"}[name.split("/")[0]]
        assert [t[0] for t in p["row2"] if t[2]] == [want], name


def test_the_theme_button_goes_in_the_box_of_row_one_and_the_script_is_in_every_page():
    js = (ROOT / "site-build" / "theme.js").read_text(encoding="utf-8")
    assert "header.site .tools" in js and "amplitrim-theme" in js
    for name, path in PAGES.items():
        assert js.strip() in path.read_text(encoding="utf-8"), name


def test_only_header_py_writes_a_header():
    for f in (ROOT / "site-build").iterdir():
        if f.suffix in (".py", ".html") and f.name != "header.py":
            assert '<header class="site">' not in f.read_text(encoding="utf-8"), f.name
