"""The ITS overnight run in the archive: the thresholds and settings in the code are the ones written in its register section, the run's folders are recognised,
and the reference sequences are there. The register is written before the run; these tests fail if a number is changed in the code or in the register only."""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))
import its_night as N  # noqa: E402
import its_night_report as P  # noqa: E402
import organise_data as od  # noqa: E402

TITLE = "ITS overnight run: what AmpliTrim can recommend for ITS"


def section() -> str:
    text = (ROOT / "TEST_REGISTER.md").read_text(encoding="utf-8")
    m = re.search(rf"^## \d+\. {re.escape(TITLE)}.*?(?=^## \d+\. |\Z)", text, re.M | re.S)
    assert m, f"the register has no section '{TITLE}'"
    return m.group(0)


def criteria() -> dict:
    return {m.group(1): m.group(2) for m in re.finditer(r"^- \*\*(P\d+)[^*]*\*\* (.*)$", section(), re.M)}


def test_every_criterion_is_written_in_the_register():
    assert sorted(criteria()) == [f"P{n}" for n in range(83, 94)]


def test_the_thresholds_of_the_code_are_the_ones_of_the_register():
    c = criteria()
    checks = [("P83", "within 10 points", P.P83_POINTS, 10.0), ("P84", "by at most 2 points", P.P84_POINTS, 2.0), ("P85", "rises by at most 1 point", P.P85_SPURIOUS, 1.0),
              ("P85", "recall falls by at most 2 points", P.P85_RECALL, 2.0), ("P86", "at least 5 points", P.P86_GAIN, 5.0), ("P87", "at most 25%", P.P87_RISE, 0.25),
              ("P90", "at least 10 points", P.P90_GAIN, 10.0), ("P91", "at least 0.8", P.P91_RHO, 0.8), ("P92", "within 1 point", P.P92_POINTS, 1.0),
              ("P93", "at most 0.5 point", P.P93_POINTS, 0.5)]
    for cid, phrase, value, expected in checks:
        assert phrase in c[cid], (cid, phrase)
        assert value == expected, (cid, value)


def test_the_definitions_and_the_arms_of_the_code_are_the_ones_of_the_register():
    s = section()
    assert (P.EDITS, P.UNJUDGED_SHARE, P.REF_SHARE, P.SEEN_SHARE) == (2, 0.03, 0.01, 0.001)
    for phrase in ("within 2 edits", "within 3% of the length", "at least 1% of an isolate's reads", "0.1% or more of the sample's reads"):
        assert phrase in s, phrase
    assert N.DIFFS == P.DIFFS == (10, 5, 2, 0) and "allowing 10 (today's tolerant default), 5, 2 and 0 differences" in s
    assert (N.WINDOW, N.WINDOW_Q) == (4, 15) and "window 4, mean quality 15" in s
    arms = {a[0]: a for a in N.ITS_ARMS}
    for name in ("std", "ee23", "ee24", "ee25", "ee33", "ee55", "tq10", "tq20", "sw15", "fwd", "mfirst", "mo20", "mo30"):
        assert f"`{name}`" in s, name
    assert [arms[a][5:7] for a in ("ee23", "ee24", "ee25", "ee33", "ee55")] == [(2, 3), (2, 4), (2, 5), (3, 3), (5, 5)] and "2,3 / 2,4 / 2,5 / 3,3 / 5,5" in s
    assert (arms["tq10"][7], arms["tq20"][7]) == (10, 20) and "`truncQ = 10` / `20`" in s
    assert (arms["mo20"][8], arms["mo30"][8]) == (20, 30) and "`minOverlap = 20` / `30`" in s
    assert P.MAXEE_ARMS == {a: arms[a][5:7] for a in P.MAXEE_ARMS} and P.Q2_ARMS == ["ee23", "ee24", "ee25", "ee33", "ee55", "tq10", "tq20", "sw15"]
    assert "5,000 read pairs per sample" in s and "chimera removal" in s


def test_the_section_is_listed_for_the_run_and_its_folders_are_recognised():
    n = od.run_number("its_night")
    row = [r for r in od.RUNS if r[0] == n][0]
    sec = re.match(r"## (\d+)\.", section()).group(1)
    assert sec in re.findall(r"\d+", row[2]) and row[1].lower() == "its overnight run"
    assert od.run_number("its_followup") != n


def test_the_reference_sequences_of_the_section_are_in_the_archive():
    folder = ROOT / "validation" / "references_its"
    for m in ("ITS1", "ITS2"):
        refs = P.read_references(folder / f"references_{m}.fasta")
        exact = {member for member, kind, _ in refs if kind.startswith("strain")}
        assert exact == {"Af", "Cn", "Ca", "Cg", "Mg", "Sc"}, (m, exact)              # the six strain-exact ATCC members of the section's review
    assert (folder / "REFERENCES_REPORT.md").is_file()
