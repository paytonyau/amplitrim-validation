"""FIGARO wrapper: compatibility shim and candidate generation."""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from amplitrim import figaro_wrapper
from amplitrim.subsample import build_subsample
from conftest import FIGARO, need_tools


@pytest.mark.skipif(sys.version_info < (3, 10), reason="the aliases still exist before Python 3.10")
def test_shim_restores_the_removed_collections_aliases():
    code = "import collections; collections.Iterable; collections.Mapping; print('ok')"
    assert subprocess.run([sys.executable, "-c", code], capture_output=True).returncode != 0
    shimmed = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                             env=figaro_wrapper.figaro_env())
    assert shimmed.returncode == 0 and "ok" in shimmed.stdout


def test_env_puts_the_shim_first_and_keeps_an_existing_pythonpath(monkeypatch):
    monkeypatch.setenv("PYTHONPATH", "/some/where")
    paths = figaro_wrapper.figaro_env(True)["PYTHONPATH"].split(os.pathsep)
    assert paths[0] == str(figaro_wrapper.COMPAT_DIR) and "/some/where" in paths
    assert figaro_wrapper.figaro_env(False)["PYTHONPATH"] == "/some/where"


def test_candidates_for_a_total_use_figaros_best_pairs(tmp_path, monkeypatch):
    seen = {}

    def fake_run(script, py, inp, out, a, m, compat=True):
        seen.update(a=a, m=m)
        out.mkdir(parents=True, exist_ok=True)
        path = out / "trimParameters.json"
        path.write_text(json.dumps([{"trimPosition": [222, 226], "score": 100.0},
                                    {"trimPosition": [223, 225], "score": 100.0},
                                    {"trimPosition": [224, 224], "score": 99.0}]))
        return path

    monkeypatch.setattr(figaro_wrapper, "run_figaro", fake_run)
    cands = figaro_wrapper.candidates_for_total(Path("f.py"), "python", tmp_path, tmp_path, total=448,
                                                min_overlap=12, top_n=2)
    assert (seen["a"], seen["m"]) == (436, 12)                  # a + m is the requested total
    assert [(c.trim_f, c.trim_r) for c in cands] == [(222, 226), (223, 225)]


@need_tools
def test_figaro_output_depends_only_on_the_sum_of_a_and_m(reads, tmp_path):
    sub = build_subsample(reads, tmp_path / "s")
    one = figaro_wrapper.run_figaro(Path(FIGARO), sys.executable, sub.dir, tmp_path / "x1", 426, 18)
    two = figaro_wrapper.run_figaro(Path(FIGARO), sys.executable, sub.dir, tmp_path / "x2", 430, 14)
    assert json.loads(one.read_text()) == json.loads(two.read_text())


def test_missing_scipy_gets_an_install_hint(tmp_path, monkeypatch):
    monkeypatch.setattr(figaro_wrapper.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(
        a, 1, "", "ModuleNotFoundError: No module named 'scipy'\n"))
    with pytest.raises(RuntimeError, match=r"amplitrim\[figaro\]"):
        figaro_wrapper.run_figaro(Path("f.py"), "python", tmp_path, tmp_path / "o", 436, 12)


def test_figaros_own_recommendation_includes_the_maxEE_it_pairs_with_the_trim_positions(tmp_path):
    import json
    from amplitrim.figaro_wrapper import parse_figaro_recommendation
    path = tmp_path / "trimParameters.json"
    path.write_text(json.dumps([{"trimPosition": [255, 217], "maxExpectedError": [2, 3], "readRetentionPercent": 87.04, "score": 82.0},
                                {"trimPosition": [250, 222], "maxExpectedError": [1, 2], "readRetentionPercent": 70.0, "score": 60.0}]))
    rec = parse_figaro_recommendation(path)
    assert (rec.trim_f, rec.trim_r, rec.max_ee, rec.score) == (255, 217, (2, 3), 82.0) and abs(rec.retention - 0.8704) < 1e-12      # the first entry, retention as a share
    path.write_text("[]")
    with pytest.raises(RuntimeError, match="no candidates"):
        parse_figaro_recommendation(path)
