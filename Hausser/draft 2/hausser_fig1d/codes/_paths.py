"""Shared paths for this analysis folder."""

from pathlib import Path

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent


def project_root() -> Path:
    for p in Path(__file__).resolve().parents:
        if (p / "src" / "pca.py").is_file():
            return p
    raise RuntimeError("Cannot find repo root (expected src/pca.py)")


ROOT = project_root()
