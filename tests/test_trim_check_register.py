"""The trim check in the archive: the criteria, the threshold and the datasets of the code are the ones written in its register section."""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))
import trim_check as T  # noqa: E402
import organise_data as od  # noqa: E402

TITLE = "The trim check: does AmpliTrim 0.3 predict the share of pairs that pass DADA2's filter after the trim?"


def section() -> str:
    text = (ROOT / "TEST_REGISTER.md").read_text(encoding="utf-8")
    m = re.search(rf"^## \d+\. {re.escape(TITLE)}.*?(?=^## \d+\. |\Z)", text, re.M | re.S)
    assert m, f"the register has no section '{TITLE}'"
    return m.group(0)


def test_the_criteria_threshold_and_datasets_are_the_ones_of_the_register():
    s = section()
    crit = {m.group(1): m.group(2) for m in re.finditer(r"^- \*\*(P\d+)\.\*\* (.*)$", s, re.M)}
    assert sorted(crit) == ["P115", "P116"] and all("within 1 point" in t and "eleven datasets" in t for t in crit.values())
    assert T.WITHIN == 1.0 and len(T.DATASETS) == 11 and "`42d795e`" in s


def test_the_section_is_listed_for_the_run_and_its_folder_is_recognised():
    n = od.run_number("trim_check")
    row = [r for r in od.RUNS if r[0] == n][0]
    assert re.match(r"## (\d+)\.", section()).group(1) in re.findall(r"\d+", row[2]) and row[1] == "Trim check"
