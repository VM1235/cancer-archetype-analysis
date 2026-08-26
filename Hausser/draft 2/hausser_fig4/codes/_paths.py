"""Shared paths for Hausser Fig. 4 reproduction."""

from pathlib import Path

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent


def project_root() -> Path:
    for p in Path(__file__).resolve().parents:
        if (p / "src" / "pca.py").is_file():
            return p
    raise RuntimeError("Cannot find repo root (expected src/pca.py)")


ROOT = project_root()
ORIGINAL = ROOT / (
    "Hausser_Original /Universal cancer tasks, evolutionary tradeoffs, "
    " and the functions of driver mutations"
)
SC_DIR = ORIGINAL / "Karaayvaz2018"
SC_COUNTS = SC_DIR / "GSE118389_counts_rsem.txt"
SC_HQ_HEADERS = SC_DIR / "GSE118389_norm_data_headers.txt"
METABRIC_ORIG = SC_DIR / "brca_metabric"
