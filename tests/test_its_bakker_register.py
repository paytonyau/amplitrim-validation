"""The Bakker night in the archive: the thresholds, arms and truths of the code are the ones written in its register section; the references are committed."""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))
import its_bakker_report as B  # noqa: E402
import its_night as N  # noqa: E402
import bakker_prep as BP  # noqa: E402
import organise_data as od  # noqa: E402

TITLE = "The Bakker night: do the ITS rules hold on a third mock community?"


def section() -> str:
    text = (ROOT / "TEST_REGISTER.md").read_text(encoding="utf-8")
    m = re.search(rf"^## \d+\. {re.escape(TITLE)}.*?(?=^## \d+\. |\Z)", text, re.M | re.S)
    assert m, f"the register has no section '{TITLE}'"
    return m.group(0)


def criteria() -> dict:
    return {m.group(1): m.group(2) for m in re.finditer(r"^- \*\*(P\d+)[^*]*\*\* (.*)$", section(), re.M)}


def test_every_criterion_is_written_in_the_register():
    assert sorted(criteria(), key=lambda c: int(c[1:])) == [f"P{n}" for n in range(105, 115)]


def test_the_thresholds_of_the_code_are_the_ones_of_the_register():
    s, c = section(), criteria()
    assert "rises by at most 1 point and recall falls by at most 2 points" in s
    assert (B.SAFE_SPURIOUS, B.SAFE_RECALL, B.USEFUL, B.LOW_PASS, B.TRIGGER_LOSS, B.EER_LOSS, B.COMBO_GAIN, B.COMBO_LOSS, B.FLAG) == (1.0, 2.0, 5.0, 0.5, 2.0, 1.0, 2.0, 1.0, 0.25)
    assert "fewer than 50%" in c["P106"] and "at least 5 points" in c["P106"] and "at most 2 points of yield" in c["P107"]
    assert "at most 1 point of yield" in c["P109"] and "by at least 2 points over `sw15`" in c["P111"] and "at most 1 point against `sw15`" in c["P111"] and "at most 25%" in c["P112"]


def test_the_arms_datasets_and_primers_of_the_code_are_the_ones_of_the_register():
    s = section()
    arms = {a[0]: a for a in N.BAKKER_ARMS}
    assert sorted(arms) == ["std", "sw10", "sw15", "tq10"] and arms["tq10"][7] == 10 and N.BAKKER_STRENGTHS == (10, 15)
    for arm in ("`std`", "`sw15`", "`sw10`", "`eer`", "`sweer`", "`tq10`"):
        assert f"| {arm} |" in s
    assert sorted(N.BAKKER_FOLDERS) == sorted(B.BAKKER) == ["bakker_pcr_its1", "bakker_pol_its1", "bakker_pol_its2"] and all(f"`{k}`" in s for k in B.BAKKER)
    assert BP.AMPLICONS == {"bakker_pcr_its1": ("ITS1", "BITS", "B58S3"), "bakker_pol_its1": ("ITS1", "ITS1F", "ITS2"), "bakker_pol_its2": ("ITS2", "ITS3_KYO2", "ITS4_KYO3")}
    for phrase in ("BITS and B58S3", "ITS1F and ITS2", "ITS3-KYO2 and ITS4-KYO3"):
        assert phrase in s, phrase


def test_the_committed_references_hold_all_19_isolates_on_both_loci():
    for m, n in (("bakker_pcr_its1", 47), ("bakker_pol_its1", 47), ("bakker_pol_its2", 23)):    # ITS1: 22 amplicons, three with ambiguity codes written as their variants
        recs = [l for l in (ROOT / "validation" / "references_bakker" / f"references_{m}.fasta").read_text().splitlines() if l.startswith(">")]
        members = 18 if m in BP.ONE_MEMBER else 19                                # on ITS1 F. oxysporum and F. verticillioides count as one
        assert len(recs) == n and len({r[1:].split("|")[0] for r in recs}) == members and all("|strain_record|" in r for r in recs)


def test_the_section_is_listed_for_the_run_and_its_folders_are_recognised():
    n = od.run_number("its_bakker")
    row = [r for r in od.RUNS if r[0] == n][0]
    assert re.match(r"## (\d+)\.", section()).group(1) in re.findall(r"\d+", row[2]) and row[1] == "Bakker night"
    assert od.run_number("its_bakker_smoke") == n and od.run_number("its_rules") != n
