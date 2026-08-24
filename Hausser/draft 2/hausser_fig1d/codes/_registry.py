"""Load 00_cancer_types.py by file path (its name starts with a digit, so it
can't be imported with a normal dotted import)."""

from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec

_HERE = Path(__file__).resolve().parent
_spec = spec_from_file_location("cancer_types_registry", _HERE / "00_cancer_types.py")
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)

CANCER_TYPES = _mod.CANCER_TYPES
fig1d_types = _mod.fig1d_types
