"""validation/its_check.py: AmpliTrim's ITS numbers against DADA2's filter on the same reads, and DADA2's loss after the filter. R and the pipeline are stand-ins."""
import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "validation"))


def load():
    spec = importlib.util.spec_from_file_location("its_check", HERE.parent / "validation" / "its_check.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def plain_pass(path, min_len=50, max_ee=2.0):
    """An independent, per-read version of DADA2's filter with no truncLen: cut at the first base of Q <= 2, then minLen, no N, expected errors."""
    out = []
    with gzip.open(path, "rt") as fh:
        lines = fh.read().splitlines()
    for i in range(0, len(lines), 4):
        seq, qual = lines[i + 1], [ord(c) - 33 for c in lines[i + 3]]
        cut = next((k for k, q in enumerate(qual) if q <= 2), len(qual))
        out.append(cut >= min_len and "N" not in seq[:cut] and sum(10 ** (-q / 10) for q in qual[:cut]) <= max_ee)
    return np.array(out)


def stand_ins(reads, tmp_path, r_shift=0, pairs_in_shift=0):
    """AmpliTrim's side from the real quality table; DADA2's side from the plain reference; `r_shift` moves DADA2's pair count, `pairs_in_shift` its input count."""
    from amplitrim.subsample import build_subsample
    from amplitrim.trimmer import ReadPairs
    info = build_subsample(reads, tmp_path / "sub", max_pairs=3000)
    pairs = ReadPairs.load(info.r1, info.r2)
    l1, l2, short, has_n, ee_fail = pairs.untruncated(2.0, 50)
    kept = ~(short | has_n | ee_fail)
    e1, _ = pairs._ee_limits(2.0)
    fwd = (np.minimum(pairs.len1, pairs.low1) >= 50) & (pairs.n1 >= np.minimum(pairs.len1, pairs.low1)) & (e1 >= np.minimum(pairs.len1, pairs.low1))
    n = pairs.n
    loss = {"pairs": int(n), "too_short": float(short.sum()) / n, "has_n": float(has_n.sum()) / n, "over_maxEE": float(ee_fail.sum()) / n, "pass": float(kept.mean()),
            "pass_not_merged": 0.1, "pass_and_merged": float(kept.mean()) - 0.1, "overlap_at_full_length": 0.9, "forward_only_pass": float(fwd.mean())}
    f1, f2 = plain_pass(info.r1), plain_pass(info.r2)
    calls = []

    def run(cmd, **kw):
        calls.append(cmd)
        out = f"{n + pairs_in_shift}\t{int((f1 & f2).sum()) + r_shift}\t{n + pairs_in_shift}\t{int(f1.sum())}\n"
        return SimpleNamespace(returncode=0, stdout="Loading dada2\n" + out, stderr="")

    return loss, info, run, calls


def write_night3(root, rows, name="AmpliTrim (maxEE 2)"):
    from compare_methods import method_key
    folder = root / "its" / "compare"
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / f"track_{method_key(name)}.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "input", "filtered", "denoised_f", "denoised_r", "merged", "nonchim"])
        w.writerows(rows)
    return root


def paths_json(tmp_path, reads):
    cfg = tmp_path / "paths.json"
    cfg.write_text(json.dumps({"its": str(reads)}))
    return cfg


def test_the_r_output_is_parsed_and_failures_are_reported():
    ic = load()
    ok = lambda cmd, **kw: SimpleNamespace(returncode=0, stdout="Loading dada2\nwarning\n1000\t900\t1000\t950\n", stderr="")      # noqa: E731
    assert ic.dada2_filter_counts("Rscript", Path("a"), Path("b"), run=ok) == {"pairs_in": 1000, "pairs_out": 900, "forward_in": 1000, "forward_out": 950}
    for bad in (SimpleNamespace(returncode=1, stdout="", stderr="Error: no package called 'dada2'"), SimpleNamespace(returncode=0, stdout="nothing useful\n", stderr="")):
        with pytest.raises(RuntimeError):
            ic.dada2_filter_counts("Rscript", Path("a"), Path("b"), run=lambda cmd, **kw: bad)


def test_the_filters_are_compared_on_the_same_reads_only():
    ic = load()
    loss = {"pairs": 1000, "pass": 0.742, "forward_only_pass": 0.861}
    same = ic.compare_filters(loss, {"pairs_in": 1000, "pairs_out": 740, "forward_in": 1000, "forward_out": 860})
    assert same["same_reads"] and same["pass_gap"] == pytest.approx(0.2) and same["forward_gap"] == pytest.approx(0.1)
    assert not ic.compare_filters(loss, {"pairs_in": 1200, "pairs_out": 890, "forward_in": 1200, "forward_out": 1030})["same_reads"]


def test_dadas_loss_after_the_filter_is_split_into_denoising_and_merging():
    ic = load()
    rows = [{"sample": "a", "input": 1000.0, "filtered": 800.0, "denoised_f": 700.0, "denoised_r": 650.0, "merged": 500.0, "nonchim": 480.0},
            {"sample": "b", "input": 1000.0, "filtered": 700.0, "denoised_f": 600.0, "denoised_r": 650.0, "merged": 400.0, "nonchim": 380.0}]
    x = ic.post_filter_loss(rows)
    assert x["filtered"] == pytest.approx(0.75) and x["merged"] == pytest.approx(0.45) and x["nonchim"] == pytest.approx(0.43) and x["samples"] == 2
    assert x["lost_denoising"] == pytest.approx((1500 - 1250) / 2000) and x["lost_merging"] == pytest.approx((1250 - 900) / 2000)           # min(F, R) per sample: 650 + 600
    assert x["lost_denoising"] + x["lost_merging"] == pytest.approx(x["filtered"] - x["merged"])
    assert ic.post_filter_loss([]) == {}
    heavy = ic.post_filter_loss([{"sample": "a", "input": 1000.0, "filtered": 800.0, "denoised_f": 400.0, "denoised_r": 450.0, "merged": 350.0, "nonchim": 340.0}])
    assert heavy["lost_denoising"] > heavy["lost_merging"]


def test_night3_tracks_are_found_by_the_method_key_and_the_newest_night_is_chosen(tmp_path):
    ic = load()
    root = write_night3(tmp_path / "runs" / "night3_2026-10-09", [("s1", 10, 8, 7, 7, 5, 5)])
    assert ic.find_night3_track(root) is not None and ic.find_night3_track(root, "AmpliTrim (maxEE 2,5)") is None and ic.find_night3_track(tmp_path / "nothing") is None
    (tmp_path / "runs" / "night3_2026-10-11").mkdir()
    (tmp_path / "runs" / "night4_2026-10-12").mkdir()
    assert ic.newest_night3(tmp_path / "runs").name == "night3_2026-10-11" and ic.newest_night3(tmp_path / "missing") is None


def test_the_newest_night_3_is_found_after_the_runs_have_been_numbered(tmp_path):
    ic = load()
    runs = tmp_path / "runs"
    for name in ("05_night3_2026-10-09", "07_night4_2026-10-10", "13_night5_mockB", "notes"):
        (runs / name).mkdir(parents=True)
    assert ic.newest_night3(runs).name == "05_night3_2026-10-09"
    (runs / "05_night3_2026-10-12").mkdir()
    assert ic.newest_night3(runs).name == "05_night3_2026-10-12"                           # the date decides, not the number
    (runs / "night3_2026-10-14").mkdir()
    assert ic.newest_night3(runs).name == "night3_2026-10-14"                              # numbered and un-numbered names compare by their dates
    (runs / "night33_2026-10-20").mkdir()
    (runs / "mynight3_2026-10-20").mkdir()
    assert ic.newest_night3(runs).name == "night3_2026-10-14"                              # a name that only contains night3 is not a night 3 folder
    assert ic.newest_night3(tmp_path / "empty") is None


def test_the_verdicts_follow_the_registered_criteria():
    ic = load()
    base = {"same_reads": True, "pass_amplitrim": 0.742, "pass_dada2": 0.743, "pass_gap": -0.1, "forward_amplitrim": 0.861, "forward_dada2": 0.86, "forward_gap": 0.1}
    v = {p: ok for p, _, ok in ic.verdicts(base, {"lost_denoising": 0.2, "lost_merging": 0.05})}
    assert v == {"P23": True, "P24": True, "P25": True}
    v = {p: ok for p, _, ok in ic.verdicts(dict(base, pass_gap=0.5, forward_gap=-0.21), {"lost_denoising": 0.05, "lost_merging": 0.2})}
    assert v == {"P23": False, "P24": False, "P25": False}                                                  # 0.2 points is the limit, 0.21 is over it
    v = {p: ok for p, _, ok in ic.verdicts(dict(base, same_reads=False), {})}
    assert v == {"P23": None, "P24": None, "P25": None}                                                      # different reads: nothing is decided


def test_the_check_runs_end_to_end_and_writes_its_report(reads, tmp_path):
    ic = load()
    loss, info, run, calls = stand_ins(reads, tmp_path)
    night3 = write_night3(tmp_path / "night3", [("s1", 1000, 800, 500, 520, 400, 390), ("s2", 1000, 700, 300, 310, 280, 270)])
    class Res: its_loss, sub = loss, info
    rc = ic.main(["--config", str(paths_json(tmp_path, reads)), "--night3", str(night3), "--out", str(tmp_path / "o")], run=run, find_rscript=lambda: "Rscript", pipeline_run=lambda s: Res)
    assert rc == 0 and len(calls) == 1 and calls[0][0] == "Rscript" and calls[0][1].endswith("its_filter_check.R") and calls[0][4:] == ["50", "2"]
    rep = (tmp_path / "o" / "ITS_CHECK_REPORT.md").read_text()
    assert "- **P23** MET" in rep and "- **P24** MET" in rep and "- **P25** MET" in rep and "## B. Where DADA2 loses the pairs" in rep and "from the 2 samples" not in rep
    assert "at least 35.0 points to denoising" in rep and "at most 6.0 points to merging" in rep and "loses 82.0 points" not in rep       # (1500 - 800) / 2000 and (800 - 680) / 2000 of the input pairs
    assert "DADA2 loses 41.0 points of the input pairs" in rep                                                                          # filtered 75.0% - merged 34.0%
    rows = list(csv.reader(open(tmp_path / "o" / "its_check.tsv"), delimiter="\t"))
    assert rows[0] == ["quantity", "amplitrim", "dada2", "gap_points"] and abs(float(rows[1][3])) < 1e-6 and abs(float(rows[2][3])) < 1e-6        # the real quality table equals the plain reference


def test_a_disagreement_is_reported_as_not_met_and_not_hidden(reads, tmp_path):
    ic = load()
    loss, info, run, _ = stand_ins(reads, tmp_path, r_shift=-30)                                              # DADA2 passes 30 fewer pairs of 3,000: 1 point
    class Res: its_loss, sub = loss, info
    assert ic.main(["--config", str(paths_json(tmp_path, reads)), "--night3", str(tmp_path / "none"), "--out", str(tmp_path / "o")], run=run, find_rscript=lambda: "Rscript", pipeline_run=lambda s: Res) == 0
    rep = (tmp_path / "o" / "ITS_CHECK_REPORT.md").read_text()
    assert "- **P23** NOT MET" in rep and "+1.00 points" in rep and "- **P24** MET" in rep and "- **P25** not decided" in rep and "was not found" in rep


def test_different_reads_make_the_comparison_invalid(reads, tmp_path):
    ic = load()
    loss, info, run, _ = stand_ins(reads, tmp_path, pairs_in_shift=10)
    class Res: its_loss, sub = loss, info
    ic.main(["--config", str(paths_json(tmp_path, reads)), "--night3", str(tmp_path / "none"), "--out", str(tmp_path / "o")], run=run, find_rscript=lambda: "Rscript", pipeline_run=lambda s: Res)
    rep = (tmp_path / "o" / "ITS_CHECK_REPORT.md").read_text()
    assert "THE READS DIFFER" in rep and "the comparison below is not valid" in rep and "- **P23** not decided" in rep


def test_missing_things_are_said_plainly(reads, tmp_path, capsys):
    ic = load()
    assert ic.main(["--config", str(tmp_path / "absent.json"), "--out", str(tmp_path / "o")]) == 2 and "No reads folder" in capsys.readouterr().out
    assert ic.main(["--config", str(paths_json(tmp_path, reads)), "--out", str(tmp_path / "o")], find_rscript=lambda: None) == 2 and "Rscript was not found" in capsys.readouterr().out
    loss, info, run, _ = stand_ins(reads, tmp_path)
    class Res: its_loss, sub = loss, info
    boom = lambda cmd, **kw: SimpleNamespace(returncode=1, stdout="", stderr="Error in library(dada2): there is no package called 'dada2'")      # noqa: E731
    assert ic.main(["--config", str(paths_json(tmp_path, reads)), "--out", str(tmp_path / "o")], run=boom, find_rscript=lambda: "Rscript", pipeline_run=lambda s: Res) == 2
    assert "DADA2 could not run" in capsys.readouterr().out
    # the package of this repository is the validated build, which has no 'where the pairs go': the check says so instead of crashing
    assert ic.main(["--config", str(paths_json(tmp_path, reads)), "--out", str(tmp_path / "o")], find_rscript=lambda: "Rscript") == 2
    out = capsys.readouterr().out
    assert "Using the amplitrim package at" in out and "has no ITS 'where the pairs go'" in out and "(it was not given)" in out


def test_the_package_folder_is_found_from_what_the_user_gave(tmp_path):
    ic = load()
    repo = tmp_path / "release" / "amplitrim"
    (repo / "amplitrim").mkdir(parents=True)
    (repo / "amplitrim" / "__init__.py").write_text("")
    assert ic.resolve_package(repo) == (repo, None)                                     # the repository: the folder that contains the package
    assert ic.resolve_package(repo / "amplitrim") == (repo, None)                       # the package folder itself
    assert ic.resolve_package(tmp_path / "release") == (repo, None)                     # the folder above the repository (unzipped one level deeper)
    root, tried = ic.resolve_package(tmp_path / "elsewhere")
    assert root is None and any(str(t).endswith("__init__.py") for t in tried)


def test_a_package_folder_without_the_package_is_refused_and_never_replaced_by_another(reads, tmp_path, capsys):
    """The first failure on the author's machine: a wrong --package silently fell back to the repository's own, older package."""
    ic = load()
    empty = tmp_path / "wrong"
    empty.mkdir()
    rc = ic.main(["--config", str(paths_json(tmp_path, reads)), "--package", str(empty), "--out", str(tmp_path / "o")], find_rscript=lambda: "Rscript")
    out = capsys.readouterr().out
    assert rc == 2 and str(empty) in out and "has no 'amplitrim' package in it" in out and "README.md, pyproject.toml and a folder named amplitrim" in out
    assert "Using the amplitrim package at" not in out                                    # it stopped before any package was used


def test_a_given_package_is_the_one_that_is_imported(reads, tmp_path, capsys):
    """A package folder that is given must win over whatever is already imported: here a stand-in package with an ITS result class."""
    ic = load()
    repo = tmp_path / "clean"
    (repo / "amplitrim").mkdir(parents=True)
    (repo / "amplitrim" / "__init__.py").write_text("")
    (repo / "amplitrim" / "pipeline.py").write_text("from dataclasses import dataclass\n@dataclass\nclass Result:\n    its_loss: object = None\nclass Settings:\n    def __init__(self, **kw): self.kw = kw\n"
                                                    "def run(settings): raise SystemExit('stand-in package reached')\n")
    import amplitrim.pipeline as already_imported                                          # the repository's own (older) package is imported, as it is on the author's machine
    old = {m: sys.modules[m] for m in sys.modules if m == "amplitrim" or m.startswith("amplitrim.")}
    assert sys.modules["amplitrim.pipeline"] is already_imported
    try:
        with pytest.raises(SystemExit, match="stand-in package reached"):
            ic.main(["--config", str(paths_json(tmp_path, reads)), "--package", str(repo), "--out", str(tmp_path / "o")], find_rscript=lambda: "Rscript")
        assert f"Using the amplitrim package at {repo.resolve() / 'amplitrim'}" in capsys.readouterr().out
    finally:
        for m in [m for m in sys.modules if m == "amplitrim" or m.startswith("amplitrim.")]:
            del sys.modules[m]
        sys.modules.update(old)
        if str(repo.resolve()) in sys.path:
            sys.path.remove(str(repo.resolve()))

