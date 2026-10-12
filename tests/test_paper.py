"""The JOSS paper in paper/: its form, that every citation resolves, and that its numbers are the ones of docs/validation.md. No network, no vsearch.
Two checks (the citation file, the command-line options) need the software repository and are skipped here; they run when the paper is moved there for submission."""
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
PAPER = (ROOT / "paper" / "paper.md").read_text(encoding="utf-8")
BIB = (ROOT / "paper" / "paper.bib").read_text(encoding="utf-8")
VALIDATION = (ROOT / "docs" / "validation.md").read_text(encoding="utf-8")
FRONT, BODY = PAPER.split("---", 2)[1], PAPER.split("---", 2)[2]
TEXT = BODY.split("# References")[0]


def test_the_front_matter_has_what_joss_asks_for():
    for key in ("title:", "tags:", "authors:", "affiliations:", "date:", "bibliography: paper.bib"):
        assert key in FRONT, key
    assert re.search(r"orcid:\s*\d{4}-\d{4}-\d{4}-\d{3}[\dX]", FRONT) and "affiliation: 1" in FRONT and "index: 1" in FRONT       # a valid ORCID, and the affiliation is declared
    title = re.search(r"title:\s*'(.+)'", FRONT).group(1)
    assert title.split(":")[0] == "AmpliTrim" and "maxEE" in title


def test_the_front_matter_agrees_with_the_citation_file_of_the_software_repository():
    """Runs when the paper sits next to the software's CITATION.cff, which is the case at submission; here, in the archive, there is none."""
    cff_path = ROOT / "CITATION.cff"
    if not cff_path.exists():
        pytest.skip("no CITATION.cff next to the paper: this check runs when the paper is in the software repository")
    cff = cff_path.read_text(encoding="utf-8")
    assert re.search(r"orcid:\s*(\S+)", FRONT).group(1) in cff                                                                    # the same person
    assert "maxEE" in re.search(r'title:\s*"(.+)"', cff).group(1)                                                                  # title and citation file name the same things


def test_the_paper_has_the_sections_joss_expects_in_order_and_is_not_too_long():
    heads = re.findall(r"^# (.+)$", BODY, flags=re.M)
    assert heads == ["Summary", "Statement of need", "State of the field", "Software design", "Research impact statement", "AI usage disclosure", "Acknowledgements", "References"]
    words = len(re.findall(r"\b[\w'`-]+\b", TEXT))
    assert 750 <= words <= 1750, words                                                              # JOSS: 750 to 1,750 words


def test_every_citation_is_in_the_bibliography_and_every_entry_is_cited():
    cited = set(re.findall(r"@([A-Za-z0-9_]+)", " ".join(re.findall(r"\[(@[^\]]+)\]", TEXT))))
    cited |= set(re.findall(r"(?<![\w.])@([A-Za-z0-9_]+)(?=[\s;,.\]])", " ".join(re.findall(r"\[[^\]]*@[^\]]*\]", TEXT))))
    entries = set(re.findall(r"^@\w+\{([^,]+),", BIB, flags=re.M))
    assert cited and cited <= entries, cited - entries
    assert entries <= cited, entries - cited                                                         # nothing in the bibliography is left uncited
    for key in entries:                                                                              # each entry has what makes it findable
        block = BIB[BIB.index("{" + key + ","):].split("\n}\n")[0]
        assert "title" in block and ("doi" in block or "url" in block), key


def test_every_figure_of_the_paper_is_on_the_validation_page():
    """Each figure with a unit (%, points, nt) that the paper states must be printed on the validation page; a number changed in one place is caught."""
    flat_paper, flat_page = re.sub(r"\s+", " ", TEXT), re.sub(r"\s+", " ", VALIDATION)
    allowed = {"0.2%": "a design constant: the longest 0.2% of amplicons are set aside", "10 points": "5.1 + 6.4 = 11.5 against 1.9: about 10",
               "7 points": "up to +7.0", "32 points": "+32.1 rounded"}
    tokens = set(re.findall(r"\d[\d,]*(?:\.\d+)?(?:%| points| nt)", flat_paper))
    assert len(tokens) >= 15, tokens                                                                 # the pattern finds the figures at all
    missing = sorted(t for t in tokens if t not in allowed and t not in flat_page)
    assert not missing, f"in the paper but not on the validation page: {missing}"
    for rounded in ("+7.0", "+21.1", "+32.1", "11.5 points", "1.9%"):                               # what the allowed roundings stand for
        assert rounded in flat_page, rounded
    for phrase in ("12 public datasets", "9 of the 11", "74 to 76%", "6.1 to 50.4", "8% and 98%", "49,720", "0.0 to 0.3 points"):
        assert phrase in flat_page, phrase
    assert "9 of 11" in flat_paper and "6.1 to 50.4" in flat_paper and "12 public datasets" in flat_paper


def test_what_the_author_has_to_confirm_is_marked_and_nothing_else_is():
    marks = re.findall(r"\[AUTHOR[^\]]*\]", TEXT)
    assert len(marks) == 3, marks                                                                    # the impact sentence, the AI disclosure, the acknowledgements
    assert not re.search(r"\b(TODO|FIXME|XXX|TBD)\b", PAPER)


def test_the_paper_states_the_design_it_describes():
    assert "two options, `--its` and `--verbose`" in TEXT
    for claim in ("0.2%", "12 nt", "75th percentile", "never below 2. There is no upper limit"):
        assert claim in TEXT, claim


def test_the_options_in_the_paper_are_those_of_the_released_command_line():
    """Runs against the released package. The development build in this archive has more options, so the check is skipped here."""
    from amplitrim.cli import build_parser
    options = {a.option_strings[0] for a in build_parser()._actions if a.option_strings and a.option_strings[0] not in ("-h", "--help", "--version")}
    if "--marker" in options:
        pytest.skip("the amplitrim package found here is the development build: this check runs in the software repository")
    assert options == {"--input-dir", "--work-dir", "--out-dir", "--its", "--verbose"}                  # three required and the two named in the paper
    pyproject = ROOT / "pyproject.toml"
    if pyproject.exists():
        assert "3.9" in pyproject.read_text(encoding="utf-8")
