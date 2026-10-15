"""validation/organise_data.py on synthetic data trees: the numbering of runs and datasets, the rewriting of amplitrim_paths.json, the refusals, the rollback, the undo,
and the finding and setting aside of duplicate files."""
import csv
import json
import os
import re
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "validation"))
import organise_data as od  # noqa: E402

DATASET_DIRS = ["mock_V1-V3_mockA_PRJEB67797", "mock_V1-V3_mockB_PRJEB67827", "mock_V4-V6_interlab_PRJEB34919", "mock_V1-V2_interlab_PRJEB34919",
                "real_V3-V4_wastewater_PRJEB75155", "real_V3-V4_swabs_PRJNA1261676", "mock_V3-V4_D6300_PRJNA1261676", "real_V3-V4_buetas_PRJNA933120"]
RUN_DIRS = ["night3_2026-10-09", "night4_2026-10-10", "night4_offline", "night5_mockB", "nocap_check_mockB", "night6_pilot", "its_check", "its_merge_check", "release_check", "round4"]


@pytest.fixture
def tree(tmp_path):
    base = tmp_path / "data"
    for name in DATASET_DIRS + ["mock_V9_mystery_PRJXX1"]:
        ready = base / name / "prepared" / "all" / "ready"
        ready.mkdir(parents=True)
        (ready / "a_R1.fastq.gz").write_bytes(name.encode())
    (base / "backup").mkdir()
    for name in RUN_DIRS + ["mystery_run", "07_already"]:
        (base / "runs" / name).mkdir(parents=True)
        (base / "runs" / name / "result.tsv").write_text(name)
    return base


def test_every_run_is_numbered_by_its_name_and_the_dated_night_4_folder_is_night_4b():
    expected = {"night3_2026-10-09": 5, "night4_2026-10-10": 7, "night4_offline": 6, "night5_mockB": 13, "nocap_check_mockB": 13, "nocap_check": 13, "night6_pilot": 14,
                "night6_PRJNA933120": 14, "its_check": 8, "its_merge_check": 9, "release_check": 10, "round4": 2, "round4_phase1": 2, "night1": 3, "night2_2026-10-09": 4,
                "soil_v4_overmerge": 11, "estimate_check": 12, "night1_x": 3}
    for name, number in expected.items():
        assert od.run_number(name) == number, name
    assert all(od.run_number(n) is None for n in ("mystery_run", "night10", "round40", "night40", "night50", "night60"))
    assert [r[0] for r in od.RUNS] == list(range(1, 15)) and len({r[1] for r in od.RUNS}) == 14


def test_datasets_are_recognised_by_accession_and_by_region_where_one_accession_holds_two():
    expected = {"mock_V1-V3_mockA_PRJEB67797": 1, "mock_V1-V3_mockB_PRJEB67827": 2, "mock_V1-V2_interlab_PRJEB34919": 3, "mock_V4-V6_interlab_PRJEB34919": 4,
                "mock_V3-V4_D6300_PRJNA1261676": 5, "real_V3-V4_swabs_PRJNA1261676": 12, "real_V1-V2_PRJEB42528": 7, "real_V4_soil_PRJNA574030": 8, "real_ITS_wheat_PRJEB58189": 9,
                "real_V3-V4_wastewater_PRJEB75155": 10, "real_V3-V4_spikein_PRJNA1243306": 11, "real_V3-V4_buetas_PRJNA933120": 13, "mock_V4-V5_D6300_SRR16961581": 6}
    for name, number in expected.items():
        assert od.dataset_number(name) == number, name
    assert od.dataset_number("mock_V9_mystery_PRJXX1") is None and od.dataset_number("mock_PRJEB34919") is None            # an ambiguous accession with no region is not guessed
    assert sorted(d[0] for d in od.DATASETS) == list(range(1, 14))


def test_the_plan_names_what_it_renames_and_leaves_or_lists_the_rest(tree):
    s = od.scan(tree)
    new = {old.name: n.name for _, _, old, n in s["plan"]}
    assert new["night3_2026-10-09"] == "05_night3_2026-10-09" and new["night4_2026-10-10"] == "07_night4_2026-10-10" and new["mock_V1-V3_mockB_PRJEB67827"] == "02_mock_V1-V3_mockB_PRJEB67827"
    assert new["mock_V4-V6_interlab_PRJEB34919"] == "04_mock_V4-V6_interlab_PRJEB34919" and new["mock_V1-V2_interlab_PRJEB34919"] == "03_mock_V1-V2_interlab_PRJEB34919"
    assert len(s["plan"]) == len(RUN_DIRS) + len(DATASET_DIRS)
    why = {p.name: w for p, w in s["skipped"]}
    assert why["07_already"] == "already numbered" and "--number mystery_run=NN" in why["mystery_run"] and "looks like a dataset" in why["mock_V9_mystery_PRJXX1"]
    assert [p.name for p in s["left"]] == ["backup"]
    only = od.scan(tree, "runs")
    assert {k for k, *_ in only["plan"]} == {"run"}
    manual = od.scan(tree, "runs", {"mystery_run": 3})
    assert ("run", 3, tree / "runs" / "mystery_run", tree / "runs" / "03_mystery_run") in manual["plan"]


def config_for(base: Path) -> dict:
    return {"extras": {"mockB": {"ready": str(base / "mock_V1-V3_mockB_PRJEB67827" / "prepared" / "all" / "ready"), "reference": str(base / "mock_V1-V3_mockB_PRJEB67827" / "reference" / "reference.fasta")},
                       "wastewater": {"ready": str(base / "real_V3-V4_wastewater_PRJEB75155" / "prepared" / "all" / "ready")}},
            "runs": [str(base / "runs" / "night3_2026-10-09"), str(base / "runs" / "night4_2026-10-10")], "elsewhere": "/somewhere/else/mock_V1-V3_mockB_PRJEB67827",
            "trap": str(base / "mock_V1-V3_mockB_PRJEB67827_extra" / "x"), "number": 5}


def test_apply_renames_rewrites_the_config_keeps_a_copy_and_the_undo_puts_everything_back(tree, capsys):
    cfg = tree.parent / "amplitrim_paths.json"
    original = json.dumps(config_for(tree), indent=2) + "\n"
    cfg.write_text(original)
    assert od.main(["--base", str(tree), "--config", str(cfg)]) == 0
    assert "Nothing was changed" in capsys.readouterr().out and (tree / "runs" / "night3_2026-10-09").is_dir() and cfg.read_text() == original      # a plan changes nothing
    assert od.main(["--base", str(tree), "--config", str(cfg), "--apply"]) == 0
    out = capsys.readouterr().out
    assert "Done. To put everything back: --undo" in out
    assert (tree / "runs" / "05_night3_2026-10-09" / "result.tsv").read_text() == "night3_2026-10-09" and not (tree / "runs" / "night3_2026-10-09").exists()
    assert (tree / "02_mock_V1-V3_mockB_PRJEB67827" / "prepared" / "all" / "ready" / "a_R1.fastq.gz").is_file() and (tree / "mystery_run" if False else tree / "runs" / "mystery_run").is_dir()
    c = json.loads(cfg.read_text())
    assert c["extras"]["mockB"]["ready"] == str(tree / "02_mock_V1-V3_mockB_PRJEB67827" / "prepared" / "all" / "ready")
    assert c["extras"]["mockB"]["reference"].endswith(str(Path("02_mock_V1-V3_mockB_PRJEB67827") / "reference" / "reference.fasta"))
    assert c["runs"] == [str(tree / "runs" / "05_night3_2026-10-09"), str(tree / "runs" / "07_night4_2026-10-10")]
    assert c["elsewhere"] == "/somewhere/else/mock_V1-V3_mockB_PRJEB67827" and c["trap"] == str(tree / "mock_V1-V3_mockB_PRJEB67827_extra" / "x") and c["number"] == 5   # only the folders that moved
    assert list(tree.parent.glob("amplitrim_paths.json.bak_*")) and (tree.parent / next(p.name for p in tree.parent.glob("amplitrim_paths.json.bak_*"))).read_text() == original
    record = next(tree.glob("organise_data_*.csv"))
    rows = list(csv.DictReader(open(record)))
    assert len(rows) == len(RUN_DIRS) + len(DATASET_DIRS) and {r["kind"] for r in rows} == {"run", "dataset"}
    assert od.main(["--base", str(tree), "--config", str(cfg), "--apply"]) == 0
    assert "Nothing to number" in capsys.readouterr().out                                    # the second time there is nothing left to do
    assert od.main(["--base", str(tree), "--config", str(cfg), "--undo", str(record)]) == 0
    assert (tree / "runs" / "night3_2026-10-09" / "result.tsv").is_file() and (tree / "mock_V1-V3_mockB_PRJEB67827").is_dir() and not list((tree / "runs").glob("0[0-9]_night*"))
    assert json.loads(cfg.read_text()) == json.loads(original)


def test_a_path_is_rewritten_only_when_the_folder_is_its_start_in_any_slash_style_and_case():
    m = [(r"C:\data\mock_V1-V3_mockB_PRJEB67827", r"C:\data\02_mock_V1-V3_mockB_PRJEB67827")]
    assert od.rewrite_value(r"C:\data\mock_V1-V3_mockB_PRJEB67827\prepared\all\ready", m) == (r"C:\data\02_mock_V1-V3_mockB_PRJEB67827\prepared\all\ready", True)
    assert od.rewrite_value("C:/data/mock_V1-V3_mockB_PRJEB67827/prepared", m) == ("C:/data/02_mock_V1-V3_mockB_PRJEB67827/prepared", True)
    assert od.rewrite_value(r"c:\DATA\MOCK_V1-V3_MOCKB_PRJEB67827", m)[1] is True
    assert od.rewrite_value(r"C:\data\mock_V1-V3_mockB_PRJEB67827_extra\x", m) == (r"C:\data\mock_V1-V3_mockB_PRJEB67827_extra\x", False)
    assert od.rewrite_value(r"D:\other\mock_V1-V3_mockB_PRJEB67827", m)[1] is False
    obj, n = od.rewrite_json({"a": [r"C:\data\mock_V1-V3_mockB_PRJEB67827\x", 3, None], "b": {"c": r"C:\data\mock_V1-V3_mockB_PRJEB67827"}}, m)
    assert n == 2 and obj["a"][1:] == [3, None]


def test_it_refuses_while_a_run_seems_to_be_going_unless_forced(tree, capsys):
    log = tree / "runs" / "night6_pilot" / "night6.log"
    log.write_text("running")
    assert od.main(["--base", str(tree), "--apply"]) == 1
    assert "A run seems to be going" in capsys.readouterr().out and (tree / "runs" / "night3_2026-10-09").is_dir()
    old = time.time() - 3 * 3600
    os.utime(log, (old, old))
    log2 = tree / "runs" / "night6_pilot" / "night6_state.json"
    log2.write_text("{}")
    assert od.main(["--base", str(tree), "--apply"]) == 1                                    # the state file is fresh
    assert od.main(["--base", str(tree), "--apply", "--force"]) == 0 and (tree / "runs" / "05_night3_2026-10-09").is_dir()


def test_a_failed_rename_puts_back_what_was_renamed_and_leaves_the_config_alone(tree, monkeypatch):
    cfg = tree.parent / "amplitrim_paths.json"
    cfg.write_text(json.dumps(config_for(tree)))
    before = sorted(p.name for p in (tree / "runs").iterdir()) + sorted(p.name for p in tree.iterdir())
    plan = od.scan(tree)["plan"]
    real, calls = Path.rename, []

    def flaky(self, target):
        calls.append(self.name)
        if len(calls) == 4:
            raise PermissionError("in use by another program")
        return real(self, target)
    monkeypatch.setattr(Path, "rename", flaky)
    with pytest.raises(SystemExit, match="3 folder\\(s\\) already renamed were put back"):
        od.apply_plan(plan, tree, [cfg], "t")
    monkeypatch.setattr(Path, "rename", real)
    assert before == sorted(p.name for p in (tree / "runs").iterdir()) + sorted(p.name for p in tree.iterdir()) and json.loads(cfg.read_text()) == config_for(tree)
    assert not list(tree.glob("organise_data_*.csv"))


def test_an_existing_target_stops_everything_before_any_rename(tree):
    (tree / "runs" / "05_night3_2026-10-09").mkdir()
    with pytest.raises(SystemExit, match="already exists"):
        od.apply_plan(od.scan(tree)["plan"], tree, [], "t")
    assert (tree / "runs" / "night3_2026-10-09").is_dir() and (tree / "runs" / "night4_2026-10-10").is_dir()


def test_the_next_run_number_follows_the_highest_one_known_or_found(tmp_path):
    assert od.next_run_number(tmp_path / "none") == 15
    runs = tmp_path / "runs"
    (runs / "15_new").mkdir(parents=True)
    (runs / "night6_old").mkdir()
    assert od.next_run_number(runs) == 16
    assert od.parse_numbers(["mystery run=7"]) == {"mystery run": 7}
    with pytest.raises(SystemExit):
        od.parse_numbers(["nonumber"])
    with pytest.raises(SystemExit):
        od.parse_numbers(["x=100"])


def write(path: Path, content: bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


def test_duplicates_are_found_by_content_not_by_name_hard_links_and_small_files_are_not_duplicates(tmp_path):
    big = b"A" * 1_500_000
    a = write(tmp_path / "a" / "reads.fastq.gz", big)
    b = write(tmp_path / "b" / "deep" / "copy_of_reads.fq.gz", big)
    write(tmp_path / "c" / "same_size_other_content.bin", b"B" * 1_500_000)
    write(tmp_path / "d" / "small1.bin", b"S" * 1000)
    write(tmp_path / "d" / "small2.bin", b"S" * 1000)
    link_ok = True
    try:
        os.link(a, tmp_path / "a" / "hardlink.fastq.gz")
    except OSError:
        link_ok = False
    other = write(tmp_path / "e" / "only_the_start_equal.bin", b"A" * 1_499_999 + b"Z")
    groups = od.find_duplicates(tmp_path, log=lambda *_: None)
    assert len(groups) == 1
    size, h, paths = groups[0]
    assert size == 1_500_000 and paths[0] in (a, tmp_path / "a" / "hardlink.fastq.gz") and paths[1:] == [b]      # the shortest path is kept; the other name of the same file is not an extra copy
    assert other not in paths and (not link_ok or len(paths) == 2)
    text = "\n".join(od.duplicate_lines(groups))
    assert "1 group(s)" in text and "1.5 MB would be freed" in text and "keep" in text and str(b) in text
    assert od.human(2_500_000_000) == "2.50 GB" and od.human(740_000) == "0.7 MB"
    assert od.find_duplicates(tmp_path, min_size=2_000_000, log=lambda *_: None) == []


def test_extra_copies_are_set_aside_never_deleted_and_can_be_restored_and_protected_copies_stay(tmp_path, capsys):
    base = tmp_path / "data"
    data = b"X" * 1_200_000
    keep = write(base / "run" / "x.bin", data)
    extra = write(base / "other" / "deeper" / "x_copy.bin", data)
    assert od.main(["--base", str(base), "--duplicates"]) == 0
    assert "Nothing was changed" in capsys.readouterr().out and extra.exists()
    assert od.main(["--base", str(base), "--duplicates", "--apply"]) == 0
    out = capsys.readouterr().out
    assert not extra.exists() and keep.exists() and "delete that folder" in out
    aside = list((base / "_duplicates").rglob("x_copy.bin"))
    assert len(aside) == 1 and aside[0].read_bytes() == data and aside[0].parent.name == "deeper"
    record = next(base.glob("duplicates_moved_*.csv"))
    row = next(csv.DictReader(open(record)))
    assert row["original"] == str(extra) and int(row["bytes"]) == 1_200_000 and row["kept_copy"] == str(keep)
    assert od.find_duplicates(base, log=lambda *_: None) == []                                   # the set-aside copy is not found again as a duplicate
    assert od.main(["--base", str(base), "--duplicates", "--restore", str(record)]) == 0
    assert extra.read_bytes() == data and "1 file(s) put back" in capsys.readouterr().out
    shutil_target = base / "run2" / "y.bin"
    write(shutil_target, data)
    assert od.main(["--base", str(base), "--duplicates", "--apply", "--protect", "deeper"]) == 0
    assert extra.exists() and not keep.exists() and not shutil_target.exists()                  # the protected copy is the one kept; the unprotected ones are set aside
    assert len(list((base / "_duplicates").rglob("x.bin"))) == 1 and len(list((base / "_duplicates").rglob("y.bin"))) == 1


def test_the_table_prints_every_run_and_dataset(capsys):
    assert od.main(["--table"]) == 0
    text = capsys.readouterr().out
    assert all(f"{n:02d}" in text for n in range(1, 15)) and "night 6" in text and "buetas-2024" in text and "night1_v2.py (not in this archive)" in text


def test_a_missing_base_folder_is_reported(tmp_path, capsys):
    assert od.main(["--base", str(tmp_path / "nowhere")]) == 2 and "does not exist" in capsys.readouterr().out


def test_when_two_copies_are_protected_neither_is_moved_and_only_the_unprotected_copy_is(tmp_path, capsys):
    base = tmp_path / "data"
    data = b"P" * 1_100_000
    first = write(base / "p1" / "protected_a.bin", data)
    second = write(base / "p2" / "protected_b.bin", data)
    third = write(base / "p3" / "ordinary.bin", data)
    assert od.main(["--base", str(base), "--duplicates", "--apply", "--protect", "protected_"]) == 0
    assert first.exists() and second.exists() and not third.exists()
    assert len(list((base / "_duplicates").rglob("ordinary.bin"))) == 1 and not list((base / "_duplicates").rglob("protected_*"))


def test_a_copy_the_config_names_is_kept_even_when_it_has_the_longer_path_and_the_protect_text_ignores_slashes_and_case(tmp_path, capsys):
    base = tmp_path / "data"
    data = b"C" * 1_200_000
    shorter = write(base / "other" / "x.bin", data)                                          # the shortest path: the default keeper
    named = write(base / "ds" / "prepared" / "all" / "ready" / "x.bin", data)                # but the config points here
    cfg = tmp_path / "amplitrim_paths.json"
    cfg.write_text(json.dumps({"extras": {"ds": {"ready": str(named.parent)}}, "note": "not a path", "gone": str(base / "missing" / "folder")}))
    assert od.config_paths([cfg]) == [str(named.parent)]                                      # only strings that are paths that exist
    assert od.main(["--base", str(base), "--duplicates", "--apply", "--config", str(cfg)]) == 0
    out = capsys.readouterr().out
    assert "Protected, because amplitrim_paths.json names them: 1 path(s)" in out
    assert named.exists() and not shorter.exists() and len(list((base / "_duplicates").rglob("x.bin"))) == 1
    assert od.is_protected(Path("C:/Data/Ds/READY/x.bin"), ["ds\\ready"]) and od.is_protected("a/b/ready_old/x", ["Ready"]) and not od.is_protected("a/b/x", ["ready"])


def write_config(path: Path, data) -> Path:
    path.write_text(json.dumps(data))
    return path


def test_check_paths_lists_every_path_of_the_config_and_exits_with_1_when_one_is_missing(tmp_path, capsys):
    here = tmp_path / "ready"
    here.mkdir()
    cfg = write_config(tmp_path / "amplitrim_paths.json", {"datasets": {"a": {"ready": str(here), "reference": str(tmp_path / "gone" / "reference.fasta")}}, "runs": [str(here)],
                                                          "url": "https://example.org/a/b", "label": "mock B", "windows": "C:\\data\\runs\\night3"})
    assert od.main(["--base", str(tmp_path), "--config", str(cfg), "--check-paths"]) == 1
    out = capsys.readouterr().out
    assert out.count("  ok       ") == 2 and out.count("MISSING  ") == 2                         # the Windows path is recognised as a path (it does not exist here); the URL and the label are not paths
    assert "reference.fasta" in out and "C:\\data\\runs\\night3" in out and "example.org" not in out
    assert "2 of 4 path(s) exist; 2 MISSING" in out


def test_check_paths_exits_with_0_when_everything_exists_and_changes_nothing(tmp_path, capsys):
    here = tmp_path / "ready"
    here.mkdir()
    cfg = write_config(tmp_path / "amplitrim_paths.json", {"ready": str(here)})
    before = sorted(p.name for p in tmp_path.iterdir())
    assert od.main(["--base", str(tmp_path), "--config", str(cfg), "--check-paths"]) == 0
    assert "1 of 1 path(s) exist." in capsys.readouterr().out and sorted(p.name for p in tmp_path.iterdir()) == before


def test_check_paths_after_the_numbering_finds_the_renamed_folders(tree, tmp_path, capsys):
    ready = tree / "mock_V1-V3_mockB_PRJEB67827" / "prepared" / "all" / "ready"
    cfg = write_config(tmp_path / "amplitrim_paths.json", {"ready": str(ready)})
    assert od.main(["--base", str(tree), "--config", str(cfg), "--check-paths"]) == 0
    assert od.main(["--base", str(tree), "--config", str(cfg), "--apply", "--force"]) == 0
    capsys.readouterr()
    assert od.main(["--base", str(tree), "--config", str(cfg), "--check-paths"]) == 0                # the config was rewritten with the folder, so it still exists
    assert "02_mock_V1-V3_mockB_PRJEB67827" in capsys.readouterr().out


def test_check_paths_without_a_config_or_with_an_unreadable_one(tmp_path, capsys, monkeypatch):
    assert od.main(["--base", str(tmp_path), "--config", str(tmp_path / "missing.json"), "--check-paths"]) == 1
    assert "cannot be read" in capsys.readouterr().out
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.setattr(od, "default_configs", lambda base: [])
    assert od.main(["--base", str(empty), "--check-paths"]) == 2                                     # no config anywhere: a message and exit code 2, not a crash
    assert "No amplitrim_paths.json found" in capsys.readouterr().out
