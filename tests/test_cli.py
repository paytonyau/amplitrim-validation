"""The command-line tool: what it prints, what it writes, and how it fails."""
import json
import re

import pandas as pd
import pytest

from amplitrim import cli
from conftest import need_tools, need_vsearch

REQUIRED = ["--input-dir", "in", "--work-dir", "work", "--out-dir", "out"]


@pytest.mark.parametrize("option", ["--region", "--primer-f-len", "--splits", "--a-values"])
def test_removed_options_are_rejected(option):
    with pytest.raises(SystemExit):
        cli.main(REQUIRED + [option, "x"])


def test_exhaustive_gets_a_smaller_default_subsample_unless_max_pairs_is_given():
    parse = lambda *extra: cli.resolve_max_pairs(cli.build_parser().parse_args(REQUIRED + list(extra)))
    assert parse() == 50_000
    assert parse("--exhaustive") == 10_000
    assert parse("--exhaustive", "--max-pairs", "0") == 0 and parse("--max-pairs", "20000") == 20_000


def test_missing_vsearch_is_reported_plainly(monkeypatch, capsys):
    monkeypatch.setattr("amplitrim.pipeline.shutil.which", lambda name: None)
    assert cli.main(REQUIRED) == 3
    assert "vsearch was not found" in capsys.readouterr().err


@need_vsearch
def test_a_missing_figaro_script_is_a_plain_error(capsys):
    assert cli.main(REQUIRED + ["--figaro-script", "/nope/figaro.py"]) == 3
    assert "FIGARO script not found" in capsys.readouterr().err


@need_vsearch
def test_writes_exactly_four_outputs(cli_run):
    assert cli_run.rc == 0
    assert {p.name for p in cli_run.out.iterdir()} == {"REPORT.md", "summary.png", "candidates.csv", "summary.json"}


@need_vsearch
def test_screen_shows_the_settings_in_a_few_plain_lines(cli_run):
    lines = [ln for ln in cli_run.stdout.splitlines() if ln.strip()]
    assert any(ln.startswith("Recommended DADA2 settings") for ln in lines)
    assert any(re.fullmatch(r"\s+truncLenF = \d+\s+truncLenR = \d+\s+maxEE = 2\s+minOverlap = 12", ln) for ln in lines)
    assert any("read pairs analysed should pass filtering and merge" in ln for ln in lines)
    assert any(ln.startswith("Report: ") for ln in lines)
    assert len(lines) <= 9
    assert not re.search(r"\b(INFO|DEBUG)\b|20\d\d-\d\d-\d\d|FIGARO|nt: ", cli_run.stdout)   # no log noise or progress


@need_vsearch
def test_screen_and_summary_json_agree_and_cover_the_longest_amplicon(cli_run):
    rec = json.loads((cli_run.out / "summary.json").read_text())["recommended"]
    assert f"truncLenF = {rec['truncLenF']}" in cli_run.stdout
    assert rec["longest_mergeable_insert"] >= 434 and rec["estimated_yield"] > 0.7   # simulated max insert: 434


@need_vsearch
def test_default_mode_never_calls_figaro(cli_run):
    assert not (cli_run.work / "figaro").exists()
    assert set(pd.read_csv(cli_run.out / "candidates.csv").suggested_by) <= {"quality table", "full reads"}


@need_vsearch
def test_min_overlap_option_changes_the_rule_the_tool_tests_under(reads, tmp_path, capsys):
    args = ["--input-dir", str(reads), "--work-dir", str(tmp_path / "w"), "--out-dir", str(tmp_path / "o"),
            "--min-overlap", "20", "--max-ee", "1"]
    assert cli.main(args) == 0
    assert "maxEE = 1   minOverlap = 20" in capsys.readouterr().out
    rec = json.loads((tmp_path / "o" / "summary.json").read_text())["recommended"]
    assert rec["longest_mergeable_insert"] == rec["truncLenF"] + rec["truncLenR"] - 20


@need_vsearch
def test_its_marker_recommends_no_truncation_and_writes_the_same_four_files(reads, tmp_path, capsys):
    args = ["--input-dir", str(reads), "--work-dir", str(tmp_path / "w"), "--out-dir", str(tmp_path / "o"),
            "--marker", "its", "--min-len", "60"]
    assert cli.main(args) == 0
    out = capsys.readouterr().out
    assert "no truncation (truncLen = 0)   maxEE = 2   minLen = 60   minOverlap = 12" in out
    assert {p.name for p in (tmp_path / "o").iterdir()} == {"REPORT.md", "summary.png", "candidates.csv", "summary.json"}
    rec = json.loads((tmp_path / "o" / "summary.json").read_text())["recommended"]
    assert (rec["truncLenF"], rec["truncLenR"], rec["minLen"]) == (0, 0, 60)


@need_tools
def test_figaro_cross_check_is_opt_in_and_adds_its_own_candidates(cli_run_figaro):
    assert cli_run_figaro.rc == 0
    table = pd.read_csv(cli_run_figaro.out / "candidates.csv")
    assert table.suggested_by.str.contains("FIGARO").any()
    runs = [d for d in (cli_run_figaro.work / "figaro").iterdir() if d.is_dir()]
    assert len(runs) == table.total_length.nunique() <= 4
    assert table.covers_longest_amplicon.all()
