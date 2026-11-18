"""The recovery report in the archive: its section and its threshold."""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))
import recovery_report as R  # noqa: E402
import organise_data as od  # noqa: E402

TITLE = "The recovery report: how many known members each mock community gave back"


def section() -> str:
    text = (ROOT / "TEST_REGISTER.md").read_text(encoding="utf-8")
    m = re.search(rf"^## \d+\. {re.escape(TITLE)}.*?(?=^## \d+\. |\Z)", text, re.M | re.S)
    assert m, "the register has no section of the recovery report"
    return m.group(0)


def test_the_method_is_the_registers():
    s = section()
    assert R.SEEN == 0.001 and "at least 0.1% of the sample's reads" in s and "at least one read" in s
    assert set(R.SIXTEEN) == {"mockA_high", "mockB", "d6300_V3-V4", "d6300_V4-V5", "miseqsop_mock"} and "75.8%" in s


def test_the_run_is_listed():
    n = od.run_number("recovery")
    row = [r for r in od.RUNS if r[0] == n][0]
    assert re.match(r"## (\d+)\.", section()).group(1) in re.findall(r"\d+", row[2]) and row[1] == "Recovery report"
