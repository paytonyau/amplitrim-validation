"""Shared fixtures. Tests that need vsearch are skipped without it; those that also need FIGARO are
skipped unless AMPLITRIM_FIGARO_SCRIPT points at figaro.py."""
import contextlib
import io
import importlib.util
import os
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from simulate import write_dataset  # noqa: E402

# The tests of this archive test its development build (amplitrim/ next to tests/), not the released package of the software repository. Several scripts replace the
# package for the whole process when they are given --package (night6.load_release, release_check.main, its_check): right for a script run on its own, wrong in
# one test process, where every test after such a call would run against the wrong package. So the development build is imported first and kept, and it is put back
# after every test (test_isolation.py checks this).
sys.path.insert(0, str(ROOT))
import amplitrim  # noqa: E402,F401


def _our_modules():
    return {k: v for k, v in sys.modules.items() if k == "amplitrim" or k.startswith("amplitrim.")}


@pytest.fixture(autouse=True)
def keep_the_development_build():
    path, ours = list(sys.path), _our_modules()
    yield
    sys.path[:] = path
    for name in list(_our_modules()):
        del sys.modules[name]
    sys.modules.update(ours)

FIGARO = os.environ.get("AMPLITRIM_FIGARO_SCRIPT")
if FIGARO and importlib.util.find_spec("scipy") is None:      # FIGARO itself needs scipy: skip its tests rather than fail them
    FIGARO = None
need_vsearch = pytest.mark.skipif(shutil.which("vsearch") is None, reason="vsearch not installed")
need_tools = pytest.mark.skipif(not (FIGARO and shutil.which("vsearch")),
                                reason="needs vsearch and AMPLITRIM_FIGARO_SCRIPT")


@pytest.fixture(scope="session")
def sim(tmp_path_factory):
    """Simulated mock community with ground truth (true insert lengths 423-434)."""
    d = tmp_path_factory.mktemp("sim")
    write_dataset(d, n_samples=2, pairs_per_sample=1500, seed=11)
    return d


@pytest.fixture(scope="session")
def reads(sim):
    """The simulated reads as the tool expects them: adapters and primers already removed."""
    return sim / "primer_trimmed"


def _run_cli(reads, root, extra):
    from amplitrim.cli import main
    inp = root / "input"
    inp.mkdir()
    for i, r1 in enumerate(sorted(reads.glob("*_R1.fastq.gz")), 1):          # Illumina-style names
        base = r1.name.replace("_R1.fastq.gz", "")
        shutil.copy(r1, inp / f"{base}_S{i}_L001_R1_001.fastq.gz")
        shutil.copy(r1.with_name(f"{base}_R2.fastq.gz"), inp / f"{base}_S{i}_L001_R2_001.fastq.gz")
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = main(["--input-dir", str(inp), "--work-dir", str(root / "work"), "--out-dir", str(root / "out")] + extra)
    return SimpleNamespace(rc=rc, out=root / "out", work=root / "work", stdout=buf.getvalue())


@pytest.fixture(scope="session")
def cli_run(reads, tmp_path_factory):
    """One real run of the command-line tool in its default mode (vsearch only)."""
    return _run_cli(reads, tmp_path_factory.mktemp("run"), [])


@pytest.fixture(scope="session")
def cli_run_figaro(reads, tmp_path_factory):
    """One real run with the optional FIGARO cross-check."""
    return _run_cli(reads, tmp_path_factory.mktemp("runf"), ["--figaro-script", FIGARO])


@pytest.fixture(scope="session")
def cli_run_ee25(reads, tmp_path_factory):
    """One real run with different forward and reverse maxEE (--max-ee 2,5)."""
    return _run_cli(reads, tmp_path_factory.mktemp("runee"), ["--max-ee", "2,5"])


@pytest.fixture(scope="session")
def cli_run_explore(reads, tmp_path_factory):
    """One real run that also finds the best truncation for several maxEE values."""
    return _run_cli(reads, tmp_path_factory.mktemp("runsug"), ["--explore-max-ee"])


@pytest.fixture(scope="session")
def cli_run_coverage(reads, tmp_path_factory):
    """One real run that also finds the best setting for several coverage levels."""
    return _run_cli(reads, tmp_path_factory.mktemp("runcov"), ["--explore-coverage"])
