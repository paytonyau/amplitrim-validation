"""The builders in site-build/ still reproduce the home page of the software repository next to this archive and the how-it-works page of this archive, byte for byte. Skipped when the software repository is not there."""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def test_the_builders_reproduce_the_shipped_pages(tmp_path):
    repo = ROOT.parent / "amplitrim"
    if not (repo / "docs" / "index.html").exists():
        pytest.skip("the software repository is not next to this archive")
    work = tmp_path / "amplitrim"
    shutil.copytree(repo, work, ignore=shutil.ignore_patterns(".git", "__pycache__", "*.pyc", "*.egg-info", "build", "dist", ".pytest_cache"))
    out = tmp_path / "docs"
    for script, extra in (("build_site.py", []), ("build_design.py", [str(out)])):
        done = subprocess.run([sys.executable, str(ROOT / "site-build" / script), str(work)] + extra, capture_output=True, text=True)
        assert done.returncode == 0, done.stderr[-500:]
    assert (work / "docs/index.html").read_text(encoding="utf-8") == (repo / "docs/index.html").read_text(encoding="utf-8"), "docs/index.html of the software repository is out of date: run the builders"
    assert (out / "how-it-works/index.html").read_text(encoding="utf-8") == (ROOT / "docs/how-it-works/index.html").read_text(encoding="utf-8"), "docs/how-it-works/index.html of this archive is out of date: run the builders"
    assert not (repo / "docs" / "how-it-works").exists()                                                    # the page lives in this archive only
