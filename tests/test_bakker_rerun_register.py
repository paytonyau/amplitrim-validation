"""The Bakker re-run in the archive: the far-primer cut and the overhang gate of the code are the ones written in its register section; section 37's criteria are kept."""
import inspect
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))
import bakker_prep as BP  # noqa: E402
import its_bakker_report as B  # noqa: E402
import organise_data as od  # noqa: E402

TITLE = "The Bakker re-run: Bakker's reads with the far primer cut"


def section() -> str:
    text = (ROOT / "TEST_REGISTER.md").read_text(encoding="utf-8")
    m = re.search(rf"^## \d+\. {re.escape(TITLE)}.*?(?=^## \d+\. |\Z)", text, re.M | re.S)
    assert m, f"the register has no section '{TITLE}'"
    return m.group(0)


def test_the_far_primer_cut_of_the_code_is_the_one_of_the_register():
    s = section()
    sig = inspect.signature(BP.trim3).parameters
    assert (sig["rate"].default, sig["min_overlap"].default) == (0.1, 3)
    assert "up to 10% mismatches" in s and "at least 3 bases" in s and "`bakker_prep.trim3`" in s


def test_the_overhang_gate_and_section_37s_criteria():
    s = section()
    assert "at most 1% of the reads of `std`" in s and "--overhang-max 0.01" in s and "--overhang-max" in inspect.getsource(B)
    assert "P105 to P114 of section 37" in s and (B.SAFE_SPURIOUS, B.SAFE_RECALL, B.USEFUL) == (1.0, 2.0, 5.0)


def test_the_section_is_listed_for_the_run_and_its_folder_is_recognised():
    n = od.run_number("bakker_rerun")
    row = [r for r in od.RUNS if r[0] == n][0]
    assert re.match(r"## (\d+)\.", section()).group(1) in re.findall(r"\d+", row[2]) and row[1] == "Bakker re-run"
    assert od.run_number("its_bakker") != n and (ROOT / "validation" / "references_bakker" / "PREP_REPORT_rerun.md").is_file()
