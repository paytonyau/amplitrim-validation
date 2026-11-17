"""The v0.4 night in the archive: the arms, grids and thresholds of the code are the ones written in its register section."""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))
import its_night as N  # noqa: E402
import v04_report as V  # noqa: E402
import pauvert_prep as P  # noqa: E402
import organise_data as od  # noqa: E402

TITLE = "The v0.4 night: a fourth mock community, the trim's window, minLen, truncQ 3 to 7, forward reads, a better ITS yield estimate, the user's settings, and 16S unchanged"


def section() -> str:
    text = (ROOT / "TEST_REGISTER.md").read_text(encoding="utf-8")
    m = re.search(rf"^## \d+\. {re.escape(TITLE)}.*?(?=^## \d+\. |\Z)", text, re.M | re.S)
    assert m, "the register has no section of the v0.4 night"
    return m.group(0)


def test_the_criteria_and_the_grids_are_the_registers():
    s = section()
    assert sorted(set(re.findall(r"\*\*(P1[12]\d|P130)", s))) == [f"P{n}" for n in range(117, 131)]
    assert N.V04_WINDOWS == (2, 3, 5, 6, 8, 10, 12, 16, 20) and N.V04_QUALITIES == (12, 18) and N.V04_MINLENS == (20, 100, 150) and N.V04_TRUNCQ == (3, 5, 7)
    assert "windows of 2 to 20 bases" in s.replace("Windows of 2 to 20 bases", "windows of 2 to 20 bases") and "mean quality 12 or 18" in s
    assert (V.SP_UP, V.REC_DOWN, V.TIE, V.FLAG) == (0.01, 0.02, 0.005, 0.25) and "within 0.5 points" in s and "at most 25%" in s
    assert "Twenty-two arms" in s and len(N.v04_rows.__doc__) and V.PAUVERT == P.LABEL and P.MIN_PARTIAL == 150 and "shorter than 150 bases" in s


def test_the_run_is_listed():
    n = od.run_number("v04_night")
    row = [r for r in od.RUNS if r[0] == n][0]
    assert re.match(r"## (\d+)\.", section()).group(1) in re.findall(r"\d+", row[2]) and row[1] == "v0.4 night"
