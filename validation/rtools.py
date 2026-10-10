"""Find Rscript even when R's folder is not on PATH, which is the usual Windows problem: R installs itself under
C:\\Program Files\\R\\R-x.y.z and is only on PATH if someone added it, and then only for the terminals opened afterwards."""

import os
import re
import shutil
from pathlib import Path
from typing import Optional


def _version(path: Path) -> tuple:
    m = re.search(r"R-(\d+)\.(\d+)\.(\d+)", str(path))
    return tuple(int(x) for x in m.groups()) if m else (0, 0, 0)


def find_rscript() -> Optional[str]:
    """Rscript on PATH, else under R_HOME, else the newest R in the usual Windows install folders. None when there is none."""
    found = shutil.which("Rscript")
    if found:
        return found
    candidates = []
    home = os.environ.get("R_HOME")
    if home:
        candidates += [Path(home) / "bin" / "x64" / "Rscript.exe", Path(home) / "bin" / "Rscript.exe", Path(home) / "bin" / "Rscript"]
    for var, sub in (("ProgramFiles", "R"), ("ProgramW6432", "R"), ("ProgramFiles(x86)", "R"), ("LOCALAPPDATA", "Programs/R")):
        root = os.environ.get(var)
        if root:
            for pattern in ("R-*/bin/x64/Rscript.exe", "R-*/bin/Rscript.exe"):
                candidates += sorted((Path(root) / sub).glob(pattern))
    candidates = [c for c in candidates if c.is_file()]
    return str(max(candidates, key=_version)) if candidates else None


NOT_FOUND = ("Rscript was not found, neither on PATH nor in the usual install folders. Install R (https://cran.r-project.org) and DADA2 "
             "(see the Windows steps in validation/README.md); python validation/check_environment.py shows what is missing.")
