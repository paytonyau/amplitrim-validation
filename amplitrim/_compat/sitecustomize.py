"""
Compatibility shim, loaded ONLY inside the FIGARO subprocess (via PYTHONPATH).

Stock FIGARO (Zymo-Research/figaro, last upstream commit 2024-07) calls
`collections.Iterable`, which was removed in Python 3.10, so it crashes on any
current Python. Re-exposing the removed aliases makes the unmodified FIGARO
script run, without editing or redistributing FIGARO's (GPLv3) source.

A no-op on Python <= 3.9, where the aliases still exist. If the environment
already has its own sitecustomize.py, this shim chains to it so it is not
silently shadowed.
"""
import collections
import collections.abc

for _name in ("Iterable", "Iterator", "Mapping", "MutableMapping", "Sequence",
              "MutableSequence", "Set", "MutableSet", "Callable", "Hashable",
              "Sized", "Container", "Generator"):
    if not hasattr(collections, _name) and hasattr(collections.abc, _name):
        setattr(collections, _name, getattr(collections.abc, _name))

import importlib.util as _ilu
import os as _os
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))
for _p in _sys.path:
    if _os.path.abspath(_p or ".") == _here:
        continue
    _candidate = _os.path.join(_p or ".", "sitecustomize.py")
    if _os.path.isfile(_candidate):
        _spec = _ilu.spec_from_file_location("_chained_sitecustomize", _candidate)
        _mod = _ilu.module_from_spec(_spec)
        _spec.loader.exec_module(_mod)
        break
