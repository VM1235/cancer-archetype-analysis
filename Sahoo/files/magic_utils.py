"""Deprecated copy — use Breast Cancer/codes/magic_utils.py.

That module implements the Sahoo et al. 2024 / Rmagic protocol
(lib-size normalize + sqrt + solver='approximate' on the full gene matrix,
returning only KS genes). Kept here so older ADAPT_*.md notes still resolve.
"""

from pathlib import Path
import sys

_CODES = Path(__file__).resolve().parents[2] / "Breast Cancer" / "codes"
sys.path.insert(0, str(_CODES))
from magic_utils import stream_full_sparse_mtx, run_magic_on_ks_genes  # noqa: E402,F401
