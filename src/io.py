"""Load expression matrices and published subtype labels."""

from pathlib import Path

import pandas as pd

from src.paths import REFERENCE, SCLC


def _first_existing(*candidates):
    for path in candidates:
        if path is not None and Path(path).is_file():
            return Path(path)
    return Path(candidates[0])


DEFAULT_CELL_LINE_MATRIX = _first_existing(
    SCLC / "data" / "SCLC_combined_Minna_CCLE_batch_corrected_wo_lowgenes.csv",
    REFERENCE
    / "data"
    / "bulk-rna-seq"
    / "SCLC_combined_Minna_CCLE_batch_corrected_wo_lowgenes.csv",
)
DEFAULT_CLUSTER_LABELS = _first_existing(
    SCLC / "data" / "combined_clusters_2020-05-27-MC copy.csv",
    REFERENCE
    / "data"
    / "bulk-rna-seq"
    / "combined_clusters_2020-05-27-MC copy.csv",
)
DEFAULT_COMBINED_MATRIX = _first_existing(
    SCLC / "data" / "CCLE_Minna_Thomas_COMBAT.csv",
    REFERENCE / "data" / "bulk-rna-seq" / "CCLE_Minna_Thomas_COMBAT.csv",
)
DEFAULT_COMBINED_METADATA = _first_existing(
    SCLC / "data" / "Metadata_CCLE_Minna_Thomas_COMBAT.csv",
    REFERENCE
    / "data"
    / "bulk-rna-seq"
    / "parti-input"
    / "Metadata_CCLE_Minna_Thomas_COMBAT.csv",
)


def load_expression_csv(path=None):
    """Load a genes × samples expression table."""
    path = Path(path) if path is not None else DEFAULT_CELL_LINE_MATRIX
    expr = pd.read_csv(path, index_col=0)
    expr.columns = expr.columns.astype(str)
    return expr


def sample_source_from_name(name):
    prefix = str(name).split(".")[0]
    if prefix == "m":
        return "Minna"
    if prefix == "c":
        return "CCLE"
    if prefix == "t":
        return "Tumor"
    return "Unknown"


def load_combined_combat(path=None):
    """Load genes × samples combined cell-line + tumor ComBat matrix."""
    path = Path(path) if path is not None else DEFAULT_COMBINED_MATRIX
    expr = pd.read_csv(path, index_col=0)
    expr.columns = expr.columns.astype(str)
    return expr


def load_combined_metadata(path=None):
    path = Path(path) if path is not None else DEFAULT_COMBINED_METADATA
    meta = pd.read_csv(path)
    meta["Sample"] = meta["Sample"].astype(str)
    return meta.set_index("Sample")


def load_author_subtypes(path=None, column="NEW_10_2020"):
    """Published hierarchical-clustering subtype labels (Figure 1B)."""
    path = Path(path) if path is not None else DEFAULT_CLUSTER_LABELS
    labels = pd.read_csv(path, index_col=0)
    labels.index = labels.index.astype(str)
    return labels[column].astype(str)
