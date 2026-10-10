"""validation/collect_results.py: the text outputs of the runs go to the archive, the reads and the work files do not, and a Windows user folder is redacted."""
import hashlib
import importlib.util
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent


def load():
    spec = importlib.util.spec_from_file_location("collect_results", HERE.parent / "validation" / "collect_results.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make_runs(root: Path):
    files = {"night3_2026-10-09/NIGHT3_REPORT.md": "# report\nC:\\Users\\Admin\\miniconda3\\envs\\amplitrim\\python.exe\n", "night3_2026-10-09/night3.log": "log at C:/Users/jsmith/data\n",
             "night3_2026-10-09/its/night3_result.json": '{"arms": {"A": {"truncF": "0"}}}', "night3_2026-10-09/its/compare/track_a.tsv": "sample\tinput\ns1\t10\n",
             "night3_2026-10-09/its/compare/asvs_a.fasta": ">ASV1\nACGT\n", "night3_2026-10-09/its/compare/reads/s1_R1.fastq.gz": "reads", "night3_2026-10-09/its/work/x.md": "work",
             "night3_2026-10-09/its/compare/big.tsv": "x" * 2_500_000, "night3_2026-10-09/its/compare/picture.png": "png", "night3_2026-10-09/its/compare/loose_R1.fastq.gz": "reads", "night3_2026-10-09/its/compare/expected.fasta": ">x\nACGT\n", "night3_2026-10-09/its/compare/blob": "no suffix",
             "its_check/ITS_CHECK_REPORT.md": "# its\n", "its_check/its_check.tsv": "quantity\tamplitrim\n"}
    for rel, text in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    (root / "stray_file.md").write_text("not in a run folder")                                          # a file directly in --runs is not a run
    return files


def test_text_outputs_are_collected_and_reads_work_files_and_large_files_are_not(tmp_path):
    cr = load()
    make_runs(tmp_path / "runs")
    r = cr.collect(tmp_path / "runs", tmp_path / "out")
    got = {p for p, _, _ in r["copied"]}
    assert got == {"night3_2026-10-09/NIGHT3_REPORT.md", "night3_2026-10-09/night3.log", "night3_2026-10-09/its/night3_result.json", "night3_2026-10-09/its/compare/track_a.tsv",
                   "its_check/ITS_CHECK_REPORT.md", "its_check/its_check.tsv"}
    assert {k: len(v) for k, v in r["left_out"].items()} == {"ASV FASTA": 1, "work folder": 2, "over 2 MB": 1, "not a text output": 4}                       # the png, the loose reads, a FASTA that is not an ASV file, a file without a suffix
    assert (tmp_path / "out" / "night3_2026-10-09" / "its" / "compare" / "track_a.tsv").read_text() == "sample\tinput\ns1\t10\n"          # the structure is kept
    assert not (tmp_path / "out" / "stray_file.md").exists() and not list((tmp_path / "out").rglob("*.fastq*")) and not list((tmp_path / "out").rglob("*.fasta"))


def test_a_windows_user_folder_is_redacted_in_the_copies_and_not_in_the_originals(tmp_path):
    cr = load()
    files = make_runs(tmp_path / "runs")
    cr.collect(tmp_path / "runs", tmp_path / "out")
    rep = (tmp_path / "out" / "night3_2026-10-09" / "NIGHT3_REPORT.md").read_text()
    log = (tmp_path / "out" / "night3_2026-10-09" / "night3.log").read_text()
    assert "C:\\Users\\<user>\\miniconda3\\envs\\amplitrim\\python.exe" in rep and "Admin" not in rep
    assert "C:/Users/<user>/data" in log and "jsmith" not in log                                          # with forward slashes too
    assert (tmp_path / "runs" / "night3_2026-10-09" / "NIGHT3_REPORT.md").read_text() == files["night3_2026-10-09/NIGHT3_REPORT.md"]   # the originals are untouched
    assert cr.redact("no user folder here, C:\\data\\runs stays") == "no user folder here, C:\\data\\runs stays"


def test_the_manifest_lists_every_copied_file_with_its_checksum(tmp_path):
    cr = load()
    make_runs(tmp_path / "runs")
    cr.collect(tmp_path / "runs", tmp_path / "out")
    rows = [ln.split("\t") for ln in (tmp_path / "out" / "MANIFEST.tsv").read_text().splitlines()]
    assert rows[0] == ["path", "bytes", "sha256"] and len(rows) == 7
    for path, size, digest in rows[1:]:
        data = (tmp_path / "out" / path).read_bytes()
        assert int(size) == len(data) and digest == hashlib.sha256(data).hexdigest(), path               # the manifest describes the copies, so the record can be checked


def test_a_dry_run_writes_nothing_and_a_bad_out_folder_is_refused(tmp_path, capsys):
    cr = load()
    make_runs(tmp_path / "runs")
    assert cr.main(["--runs", str(tmp_path / "runs"), "--out", str(tmp_path / "out"), "--dry-run"]) == 0 and not (tmp_path / "out").exists()
    assert "Would copy 6 files" in capsys.readouterr().out
    assert cr.main(["--runs", str(tmp_path / "runs"), "--out", str(tmp_path / "runs")]) == 2                  # into itself
    assert cr.main(["--runs", str(tmp_path / "runs"), "--out", str(tmp_path / "runs" / "results")]) == 2     # inside --runs
    assert cr.main(["--runs", str(tmp_path / "absent"), "--out", str(tmp_path / "out")]) == 2
    (tmp_path / "empty").mkdir()
    assert cr.main(["--runs", str(tmp_path / "empty"), "--out", str(tmp_path / "out2")]) == 2 and "nothing to copy" in capsys.readouterr().out


def test_a_real_run_prints_what_was_copied_and_what_was_left_out(tmp_path, capsys):
    cr = load()
    make_runs(tmp_path / "runs")
    assert cr.main(["--runs", str(tmp_path / "runs"), "--out", str(tmp_path / "out")]) == 0
    out = capsys.readouterr().out
    assert "Copied 6 files" in out and "from 2 run folders" in out and "left out, ASV FASTA: 1 files" in out and "left out, over 2 MB: 1 files" in out
