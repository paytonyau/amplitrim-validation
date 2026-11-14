"""The ITS rules night in the archive: the thresholds and settings in the code are the ones written in its register section, and its folders are recognised."""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))
import its_night as N  # noqa: E402
import its_rules_report as P  # noqa: E402
import organise_data as od  # noqa: E402

TITLE = "ITS rules night: can AmpliTrim compute truncQ, the trim strength and maxEE for ITS?"


def section() -> str:
    text = (ROOT / "TEST_REGISTER.md").read_text(encoding="utf-8")
    m = re.search(rf"^## \d+\. {re.escape(TITLE)}.*?(?=^## \d+\. |\Z)", text, re.M | re.S)
    assert m, f"the register has no section '{TITLE}'"
    return m.group(0)


def criteria() -> dict:
    return {m.group(1): m.group(2) for m in re.finditer(r"^- \*\*(P\d+)[^*]*\*\* (.*)$", section(), re.M)}


def test_every_criterion_is_written_in_the_register():
    assert sorted(criteria(), key=lambda c: int(c[1:])) == [f"P{n}" for n in range(94, 105)]


def test_the_thresholds_of_the_code_are_the_ones_of_the_register():
    c = criteria()
    assert "rises by at most 1 point" in c["P94"] and "recall falls by at most 2 points" in c["P94"] and "both truths" in c["P94"]
    assert (P.SAFE_SPURIOUS, P.SAFE_RECALL, P.USEFUL, P.MECHANISM, P.FLAG) == (1.0, 2.0, 5.0, 2.0, 0.25)
    assert "at least 5 points" in c["P95"] and "rises by at least 2 points" in c["P96"] and "at most 25%" in c["P97"] and "at least 5 points" in c["P103"]


def test_the_arms_and_the_rule_of_the_code_are_the_ones_of_the_register():
    s = section()
    assert P.TQ_ARMS == ["tq3", "tq4", "tq5", "tq6", "tq8"] and "`tq3`, `tq4`, `tq5`, `tq6`, `tq8`" in s
    assert N.SW_STRENGTHS == (10, 15, 20, 25) and "`sw10`, `sw15`, `sw20`, `sw25`" in s and P.SW_ARMS == ["sw10", "sw20", "sw25"]
    arms = {a[0]: a for a in N.RULES_ARMS}
    assert [arms[f"tq{q}"][7] for q in (3, 4, 5, 6, 8, 10)] == [3, 4, 5, 6, 8, 10] and arms["std"][5:8] == (2, 2, 2)
    assert (N.RULE_Q, N.TRUNC_Q, N.MIN_LEN) == (0.75, 2, 50)
    for phrase in ("75th percentile", "at least 50 bases", "max(2, int(P75) + 1)", "`truncQ = 2`", "window 4"):
        assert phrase in s, phrase


def test_the_section_is_listed_for_the_run_and_its_folders_are_recognised():
    n = od.run_number("its_rules")
    row = [r for r in od.RUNS if r[0] == n][0]
    sec = re.match(r"## (\d+)\.", section()).group(1)
    assert sec in re.findall(r"\d+", row[2]) and row[1].lower() == "its rules night"
    assert od.run_number("its_night") != n and od.run_number("its_rules_smoke") == n            # the smoke folder must be moved away before numbering
