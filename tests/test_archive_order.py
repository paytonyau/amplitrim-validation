"""The order of the archive: the index lists every file of validation/, the run numbers agree between the code, the validation README and the register, and every path the
root README names exists. These fail when a script is added without being listed, or when a number is changed in one place only."""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))
import organise_data as od  # noqa: E402


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


def sections_of(text: str) -> set:
    lo, _, hi = text.strip().partition("-")
    return set(range(int(lo), int(hi or lo) + 1))


def test_the_index_lists_every_file_of_validation_once_in_a_table():
    readme = read("validation/README.md")
    index = readme[readme.index("## Which script is which"):readme.index("## Notes on the tools")]
    files = sorted(p.name for p in (ROOT / "validation").iterdir() if p.is_file() and p.name != "README.md")
    missing = [f for f in files if f"`{f}`" not in index]
    assert not missing, f"not in the index of validation/README.md: {missing}"
    assert (ROOT / "validation" / "simulated").is_dir() and "`simulated/`" in index


def test_the_run_numbers_agree_between_the_code_the_validation_readme_and_the_register():
    rows = re.findall(r"^\| (\d\d) \| (.+?) \| (.+?) \| ([\d-]+) \|$", read("validation/README.md"), re.M)
    assert [int(n) for n, *_ in rows] == [r[0] for r in od.RUNS] == list(range(1, 15))
    for (n, what, scripts, sections), run in zip(rows, od.RUNS):
        assert sections == run[2] and what.lower() == run[1].lower(), (n, what)
        assert scripts.replace("`", "") == run[3], (n, scripts)                      # the scripts column says what the code says
    register = read("TEST_REGISTER.md")
    contents = re.findall(r"^\| (\d+) \| (\d\d) \| (.+) \|$", register, re.M)
    headings = re.findall(r"^## (\d+)\. ", register, re.M)
    assert [c[0] for c in contents] == headings and len(headings) == 26                    # the table lists every section once, in order
    for section, run, _ in contents:
        owner = [r[0] for r in od.RUNS if int(section) in sections_of(r[2])]
        assert owner == [int(run)], (section, run, owner)


def test_every_path_of_the_root_readme_table_exists():
    readme = read("README.md")
    table = readme[readme.index("**1. The answer**"):readme.index("## Which software was tested")]
    paths = [p for cell in re.findall(r"^\| (.+?) \| ", table, re.M) for p in re.findall(r"`([^`\s]+)`", cell)]
    assert len(paths) >= 12
    missing = [p for p in paths if not (ROOT / p.rstrip("/")).exists() and "/" in p + "/"]
    assert not missing, f"named in README.md but not there: {missing}"


def test_the_scripts_of_the_run_table_exist_except_the_one_flagged_as_not_in_this_archive():
    for number, name, sections, scripts, _ in od.RUNS:
        for item in re.split(r",\s*(?![^()]*\))", scripts):
            m = re.match(r"([\w]+\.(?:py|R))", item.strip())
            assert m, (number, item)
            exists = (ROOT / "validation" / m.group(1)).is_file()
            flagged = "not in this archive" in item
            assert exists != flagged, f"run {number:02d}: {m.group(1)} " + ("is flagged as missing but exists" if exists else "is missing and not flagged")
    assert "not in this archive" in read("validation/README.md")                                # the gap is stated where the run table is
