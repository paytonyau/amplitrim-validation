"""The 16S section of the website (docs/16s/), the overview of both markers (docs/index.html) and the 'moved' pages at the first version's addresses are true to the findings: they are
exactly what site-build/build_16s.py makes from docs/FINDINGS.md, every number on them is a number of the findings (or of the names and recovery tables), FIGARO is named, none of the
letters the nights gave their arms shows, and every link goes somewhere."""
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "site-build"))
import build_16s as B          # noqa: E402
import results_data as R       # noqa: E402
import test_its_site as T      # noqa: E402  (the page parser and the ITS pages)

DOCS = ROOT / "docs"
PAGES = B.build(ROOT)
PARSED = {rel: T.P(text) for rel, text in PAGES.items()}
S = B.Site(ROOT)
NUMBER = re.compile(r"[+\-\u2212]?\d+(?:\.\d+)?")
PAGES16 = [r for r in PAGES if r.startswith("16s/datasets/")]


def text(rel):
    return " ".join(" ".join(PARSED[rel].text).split())


def numbers(s):
    return {n.lstrip("+-\u2212") for n in NUMBER.findall(s.replace(",", ""))}


def test_the_pages_on_disk_are_exactly_what_the_builder_makes_from_the_findings():
    assert len(PAGES) == 2 + len(B.IDS16) + len(R.IDS) and len(PAGES16) == 11
    for rel, html in PAGES.items():
        assert (DOCS / rel).read_text(encoding="utf-8") == html, f"docs/{rel} is out of date: run site-build/build_16s.py"


@pytest.mark.parametrize("rel", sorted(PAGES))
def test_each_page_is_well_formed_self_contained_and_accessible(rel):
    p, raw = PARSED[rel], PAGES[rel]
    assert not p.errors and not p.stack, (rel, p.errors, p.stack)
    assert p.lang == "en" and '<meta name="viewport"' in raw and "<title>" in raw
    assert len(re.findall(r"<script>(.*?)</script>", raw, flags=re.S)) == 2 and not any("src" in s for s in p.scripts)
    assert p.headings.count(1) == 1 and all(b - a <= 1 for a, b in zip(p.headings, p.headings[1:])), (rel, p.headings)
    assert len(p.ids) == len(set(p.ids)) and p.links[0] == "#main"
    assert all(s["role"] == "img" and s["label"] for s in p.svgs) and len(p.figs) == len(p.svgs) and all(len(f["caption"].strip()) > 40 for f in p.figs), rel


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


def test_every_number_on_the_pages_is_a_number_of_the_findings_the_names_or_the_recovery_table():
    known = numbers((DOCS / "FINDINGS.md").read_text(encoding="utf-8") + (DOCS / "validation.md").read_text(encoding="utf-8") + (DOCS / "DATASET_NAMES.tsv").read_text(encoding="utf-8")
                    + (ROOT / "results" / "32_recovery" / "RECOVERY.md").read_text(encoding="utf-8"))
    for rel in PAGES:
        if rel.startswith("datasets/") or rel == "index.html":                                    # the root overview also quotes the ITS runs (tested with the ITS pages)
            continue
        extra = numbers(text(rel)) - known
        extra = {n for n in extra if not (n.isdigit() and int(n) <= 31)}                       # counts of datasets and rows
        assert not extra, (rel, sorted(extra))


@pytest.mark.parametrize("rel", sorted(PAGES))
def test_figaro_is_named_and_the_nights_letters_and_placeholders_are_gone(rel):
    t = text(rel)
    assert "the other tool" not in t.lower() and "to be verified" not in t and "in the findings" not in t.replace("docs/FINDINGS.md", ""), rel
    assert not re.search(r"\b(B2|B4|F2)\b|\bArm [A-G]\b", t), rel
    if rel.startswith("16s/datasets/"):
        assert "FIGARO" in t, rel


def test_the_overview_of_a_dataset_agrees_with_the_findings():
    t = text("16s/datasets/ds01-hulpusch-mock-a/index.html")
    assert "-0.8 points" in t and "65%" in t and "0.0 points" in t and "2,3" in t and "12 of 14" in t
    assert "DS01 Mock A (H\u00fclp\u00fcsch 2023), 16S V1-V3" in t and "Hülpüsch et al., 2023" in t
    o = text("16s/index.html")
    assert "within 0.8 points" in o and all(S.row[i]["short_name"] in o for i in B.IDS16)
    assert len(re.findall(r"<tr>", PAGES["16s/index.html"])) == 1 + len(B.IDS16)


def test_every_old_address_has_a_moved_page_that_leads_to_a_page_that_exists():
    for i in R.IDS:
        raw = PAGES[f"datasets/{i}/index.html"]
        target = re.search(r'content="0; url=([^"]+)"', raw).group(1)
        assert (DOCS / "datasets" / i / target).resolve().joinpath("index.html").exists() or (DOCS / "datasets" / i / target / "index.html").exists(), (i, target)


def test_the_root_overview_leads_to_both_markers_and_states_the_headline_numbers():
    t = PAGES["index.html"]
    assert 'href="16s/"' in t and 'href="its/"' in t and "AmpliTrim" in t
    o = text("index.html")
    assert "within 0.8 points" in o and "+5.1 to +57.3 points" in o and "1.2 points" in o
