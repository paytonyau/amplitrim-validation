"""validation/release_check.py: the released package's recommendation against the settings that nights 3 and 4b ran through DADA2. The package is a stand-in."""
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "validation"))


def load():
    spec = importlib.util.spec_from_file_location("release_check", HERE.parent / "validation" / "release_check.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# what night 3 (arm A) and night 4b (arm H, where the rule relaxed maxEE) ran through DADA2; None = the dataset has no truncation (ITS)
A = {"mockA_high": (277, 243, "2"), "v1v2": (195, 168, "2"), "its": (0, 0, "2"), "v4": (124, 115, "2"), "mcm_v12": (160, 165, "2"), "mockB": (279, 248, "2"), "mcm_v46": (273, 265, "2"),
     "d6300_V3-V4": (224, 218, "2"), "d6300_V4-V5": (282, 280, "2"), "wastewater_V3-V4": (254, 186, "2"), "spikein_V3-V4": (229, 210, "2"), "swabs_V3-V4": (224, 218, "2")}
H = {"mockA_high": (273, 247, "2,3"), "mockB": (277, 250, "4,5"), "mcm_v46": (271, 267, "3"), "d6300_V4-V5": (282, 280, "2,3"), "wastewater_V3-V4": (261, 183, "4,3")}


def write_nights(tmp_path):
    n3, n4 = tmp_path / "n3", tmp_path / "n4"
    for name, (f, r, ee) in A.items():
        (n3 / name).mkdir(parents=True)
        (n3 / name / "night3_result.json").write_text(json.dumps({"arms": {"A": {"truncF": str(f), "truncR": str(r), "maxEE": ee}, "G": {"truncF": "1", "truncR": "1", "maxEE": "2,5"}}}))
    for name in A:
        (n4 / name).mkdir(parents=True)
        body = {"arms": {"H": {"truncF": str(H[name][0]), "truncR": str(H[name][1]), "maxEE": H[name][2]}}} if name in H else {"skipped": "the rule gives (2,2)"}
        (n4 / name / "night4_result.json").write_text(json.dumps(body))
    return n3, n4


def write_config(tmp_path):
    cfg = {"extras": {}}
    for name in A:
        d = tmp_path / "data" / name
        d.mkdir(parents=True, exist_ok=True)
        if name in ("mockA_high", "v1v2", "its"):
            cfg[name] = str(d).replace("\\", "/")
        else:
            cfg["extras"][name] = {"ready": str(d).replace("\\", "/")}
    p = tmp_path / "paths.json"
    p.write_text(json.dumps(cfg))
    return p


def package(released, seen=None):
    """A stand-in for the package: what it recommends on a folder, by the folder's name; `seen` records the marker each dataset was run with."""
    def run(settings):
        name = Path(settings.input_dir).name
        if seen is not None:
            seen[name] = settings.marker
        got = released.get(name)
        return SimpleNamespace(recommendation=None if got is None else SimpleNamespace(trim_f=got[0], trim_r=got[1]), cfg=SimpleNamespace(max_ee=got[2] if got else 2.0))
    return run


def what_was_tested():
    return {k: (H[k] if k in H else A[k]) for k in A}


def as_released(settings):
    return {k: (v[0], v[1], tuple(float(x) for x in v[2].split(",")) if "," in v[2] else float(v[2])) for k, v in settings.items()}


def test_maxee_is_read_as_a_number_a_pair_or_text():
    rc = load()
    assert rc.parse_ee("2") == (2.0, 2.0) and rc.parse_ee("2,3") == (2.0, 3.0) and rc.parse_ee(3.0) == (3.0, 3.0) and rc.parse_ee((2.0, 5.0)) == (2.0, 5.0) and rc.parse_ee("4, 3") == (4.0, 3.0)
    assert rc.ee_text((3.0, 3.0)) == "3" and rc.ee_text((2.0, 4.0)) == "2,4"


def test_the_tested_settings_are_arm_h_where_the_rule_relaxed_and_arm_a_otherwise(tmp_path):
    rc = load()
    n3, n4 = write_nights(tmp_path)
    assert rc.tested_settings(n3, n4, "mockB") == ("night 4b, arm H", 277, 250, (4.0, 5.0))
    assert rc.tested_settings(n3, n4, "v4") == ("night 3, arm A", 124, 115, (2.0, 2.0))                      # night 4b skipped it: the rule gave (2,2)
    assert rc.tested_settings(n3, n4, "its") == ("night 3, arm A", 0, 0, (2.0, 2.0))
    assert rc.tested_settings(n3, None, "mockB") == ("night 3, arm A", 279, 248, (2.0, 2.0))                  # without night 4b the arm A is what there is
    assert rc.tested_settings(tmp_path / "nothing", n4, "v4") is None and rc.tested_settings(None, None, "v4") is None


def test_the_settings_are_the_same_only_if_the_cut_and_the_maxEE_are_equal():
    rc = load()
    tested = ("night 4b, arm H", 277, 250, (4.0, 5.0))
    assert rc.compare("mockB", (277, 250, (4.0, 5.0)), tested)["same"] is True
    assert rc.compare("mockB", (277, 250, (4.0, 4.0)), tested)["same"] is False                              # a different maxEE
    assert rc.compare("mockB", (276, 250, (4.0, 5.0)), tested)["same"] is False                              # one nucleotide off
    assert rc.compare("mockB", (277, 251, (4.0, 5.0)), tested)["same"] is False
    assert rc.compare("mockB", None, tested)["same"] is False and rc.compare("mockB", (277, 250, (4.0, 5.0)), None)["same"] is None


def test_the_verdict_needs_every_dataset_and_every_one_the_same(tmp_path):
    rc = load()
    tested = lambda k: ("t", *( (A | H)[k][:2]), rc.parse_ee((A | H)[k][2]))                               # noqa: E731
    rows = [rc.compare(k, (tested(k)[1], tested(k)[2], tested(k)[3]), tested(k)) for k in rc.DATASETS]
    assert rc.verdict(rows) is True
    rows[3] = rc.compare("v4", (125, 115, (2.0, 2.0)), tested("v4"))
    assert rc.verdict(rows) is False
    assert rc.verdict(rows[:11]) is None and rc.verdict([dict(r, same=None) if i == 0 else r for i, r in enumerate(rows)]) is None


def test_the_check_runs_end_to_end_and_the_report_says_what_follows(tmp_path):
    rc = load()
    n3, n4 = write_nights(tmp_path)
    seen = {}
    rc_rc = rc.main(["--config", str(write_config(tmp_path)), "--night3", str(n3), "--night4", str(n4), "--out", str(tmp_path / "o")], pipeline_run=package(as_released(what_was_tested()), seen))
    assert len(seen) == 12 and seen["its"] == "its" and all(v == "16s" for k, v in seen.items() if k != "its")             # ITS is run as ITS, everything else as 16S
    rep = (tmp_path / "o" / "RELEASE_CHECK_REPORT.md").read_text()
    assert rc_rc == 0 and "- **P30** MET" in rep and "12 of 12 the same" in rep and "are the released package's results" in rep
    assert "| mockB | night 4b, arm H | 277/250, maxEE 4,5 | 277/250, maxEE 4,5 | yes |" in rep and "| its | night 3, arm A | 0/0, maxEE 2 | 0/0, maxEE 2 | yes |" in rep
    rows = (tmp_path / "o" / "release_check.tsv").read_text().splitlines()
    assert len(rows) == 13 and rows[0].startswith("dataset\ttested") and "mockB\tnight 4b, arm H\t277\t250\t4,5\t277\t250\t4,5\tTrue" in rows


def test_a_difference_is_named_and_the_nights_are_not_claimed_for_it(tmp_path):
    rc = load()
    n3, n4 = write_nights(tmp_path)
    released = as_released(what_was_tested())
    released["wastewater_V3-V4"] = (261, 183, (4.0, 4.0))                                                    # the package gives 4 where the night ran 4,3
    rc.main(["--config", str(write_config(tmp_path)), "--night3", str(n3), "--night4", str(n4), "--out", str(tmp_path / "o")], pipeline_run=package(released))
    rep = (tmp_path / "o" / "RELEASE_CHECK_REPORT.md").read_text()
    assert "- **P30** NOT MET" in rep and "11 of 12 the same, 1 different" in rep and "| wastewater_V3-V4 | night 4b, arm H | 261/183, maxEE 4,3 | 261/183, maxEE 4 | **NO** |" in rep
    assert "are not claimed for the released package" in rep and "are the released package's results" not in rep


def test_a_part_run_is_not_decided(tmp_path):
    rc = load()
    n3, n4 = write_nights(tmp_path)
    rc.main(["--config", str(write_config(tmp_path)), "--night3", str(n3), "--night4", str(n4), "--only", "mockB", "v4", "--out", str(tmp_path / "o")], pipeline_run=package(as_released(what_was_tested())))
    rep = (tmp_path / "o" / "RELEASE_CHECK_REPORT.md").read_text()
    assert "- **P30** not decided" in rep and "2 of 2 the same" in rep


def test_missing_things_are_said_plainly(tmp_path, capsys):
    rc = load()
    assert rc.main(["--config", str(write_config(tmp_path)), "--out", str(tmp_path / "o")], pipeline_run=package({})) == 2 and "Give --night3 and --night4" in capsys.readouterr().out
    empty = tmp_path / "wrong"
    empty.mkdir()
    assert rc.main(["--config", str(write_config(tmp_path)), "--night3", str(empty), "--package", str(empty), "--out", str(tmp_path / "o")]) == 2
    assert "has no 'amplitrim' package in it" in capsys.readouterr().out
    # the repository's own package is the development build: it has no percentile rule, and the check says so instead of comparing the wrong thing
    assert rc.main(["--config", str(write_config(tmp_path)), "--night3", str(empty), "--out", str(tmp_path / "o")]) == 2
    out = capsys.readouterr().out
    assert "Using the amplitrim package at" in out and "has no percentile rule" in out
