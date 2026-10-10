#!/usr/bin/env python3
"""
Check that everything needed to run AmpliTrim, its tests and the real-data comparison is installed, and say
how to fix what is missing.

    python validation/check_environment.py
"""
import importlib
import os
import platform
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rtools import find_rscript  # noqa: E402
import subprocess
import sys

WINDOWS = platform.system() == "Windows"
results = []


def report(name, ok, detail, hint="", required=True):
    results.append((name, ok, required))
    mark = "ok     " if ok else ("MISSING" if required else "absent ")
    print(f"[{mark}] {name:22s} {detail}")
    if not ok and hint:
        print(f"{'':11s}-> {hint}")


def run(cmd):
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return res.returncode, (res.stdout + res.stderr).strip()
    except (OSError, subprocess.TimeoutExpired) as err:
        return 1, str(err)


def main() -> int:
    print(f"{platform.system()} {platform.release()}, Python {platform.python_version()} ({sys.prefix})\n")
    report("Python >= 3.9", sys.version_info >= (3, 9), platform.python_version(), "install Python 3.9 or later")
    for mod in ("amplitrim", "numpy", "pandas", "matplotlib"):
        try:
            m = importlib.import_module(mod)
            report(mod, True, getattr(m, "__version__", ""))
        except ImportError:
            report(mod, False, "not installed", 'run: pip install -e .   (from the AmpliTrim folder)')
    try:
        importlib.import_module("pytest")
        report("pytest", True, importlib.import_module("pytest").__version__)
    except ImportError:
        report("pytest", False, "not installed", "run: pip install pytest", required=False)

    exe = shutil.which("vsearch")
    if exe:
        _, out = run([exe, "--version"])
        line = next((l for l in out.splitlines() if l.startswith("vsearch v")), out.splitlines()[0] if out else exe)
        report("vsearch", True, line[:60])
    else:
        hint = ("run: python validation/install_vsearch_windows.py" if WINDOWS else
                "run: conda install -c bioconda vsearch")
        report("vsearch", False, "not on PATH", hint)

    code, out = run([sys.executable, "-m", "cutadapt", "--version"])
    report("cutadapt", code == 0, out.splitlines()[0] if code == 0 else "not installed",
           "run: pip install cutadapt   (needed to remove primers and adapters)", required=False)

    rscript = find_rscript()
    if rscript:
        code, out = run([rscript, "-e", "cat(as.character(packageVersion('dada2')))"])
        ok = code == 0 and out[:1].isdigit()
        off_path = shutil.which("Rscript") is None
        report("R + DADA2", ok, (f"DADA2 {out}" + (f" (R found at {rscript}, not on PATH; the tools use it anyway)" if off_path else "")) if ok else "R found, DADA2 not installed",
               'in R: install.packages("BiocManager"); BiocManager::install("dada2")', required=False)
    else:
        report("R + DADA2", False, "Rscript not found (not on PATH, not in the usual install folders)",
               "install R (https://cran.r-project.org), then install DADA2 in R",
               required=False)

    figaro = os.environ.get("AMPLITRIM_FIGARO_SCRIPT")
    report("FIGARO (optional)", bool(figaro and os.path.exists(figaro)),
           figaro or "AMPLITRIM_FIGARO_SCRIPT not set", "only for the optional cross-check", required=False)

    need = [ok for _n, ok, req in results if req]
    tests_ok = all(need)
    runner_ok = tests_ok and all(ok for n, ok, _r in results if n in ("cutadapt", "R + DADA2"))
    print(f"\nunit tests and AmpliTrim: {'ready' if tests_ok else 'NOT ready'}")
    print(f"real-data comparison (cutadapt + DADA2): {'ready' if runner_ok else 'not ready'}")
    return 0 if tests_ok else 1


if __name__ == "__main__":
    sys.exit(main())
