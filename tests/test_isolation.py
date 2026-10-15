"""The development build of this archive stays the package that `import amplitrim` finds, whatever a test does to the import system (tests/conftest.py).

Scripts such as night6.load_release and release_check.main replace the package for the whole process when they are given --package, as they should when run on their own;
in a test process that would make every later test run against the released package of the software repository (whose command line has fewer options), which is how a
whole-suite run once failed on a machine where that repository sits next to this archive.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "validation"))
import night6  # noqa: E402


def fake_release(folder: Path) -> Path:
    pkg = folder / "released" / "amplitrim"
    pkg.mkdir(parents=True)
    (pkg / "__init__.py").write_text('__version__ = "0.1.0"\n')
    (pkg / "pipeline.py").write_text("WARN_ABOVE = 5\n")
    return folder / "released"


def test_a_test_that_loads_the_released_package_replaces_it_for_that_test_only(tmp_path):
    import amplitrim
    assert Path(amplitrim.__file__).resolve().parent == ROOT / "amplitrim"                    # the development build to begin with
    night6.load_release(fake_release(tmp_path))
    pipeline = sys.modules["amplitrim.pipeline"]
    assert Path(pipeline.__file__).resolve().parent == (tmp_path / "released" / "amplitrim").resolve()      # inside the test the released package is the one in use
    assert sys.path[0] == str((tmp_path / "released").resolve())


def test_the_next_test_finds_the_development_build_again():
    import amplitrim
    assert Path(amplitrim.__file__).resolve().parent == ROOT / "amplitrim"
    assert sys.path[0] != "" and not any("released" in p for p in sys.path[:3])
    from amplitrim.cli import build_parser                                                    # a module of the development build, absent from the stand-in release
    assert "--explore-max-ee" in build_parser().format_help()


def test_release_check_without_a_package_sees_the_development_build(tmp_path, capsys):
    sys.path.insert(0, str(ROOT / "validation"))
    import release_check as rc
    assert rc.main(["--night3", str(tmp_path), "--out", str(tmp_path / "o")]) == 2
    assert "has no percentile rule" in capsys.readouterr().out

