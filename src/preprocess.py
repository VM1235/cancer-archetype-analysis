"""Turn a DepMap all-cancer expression dump into a genes × samples Panel A matrix."""

from __future__ import annotations

import re

import pandas as pd

GENE_WITH_ENTREZ = re.compile(r"^(.*) \((\d+)\)$")


def strip_entrez(name):
    name = str(name)
    m = GENE_WITH_ENTREZ.match(name)
    return m.group(1) if m else name


def pick_breast_models(models, prefer_invasive=True):
    """Return a boolean mask of breast / BRCA models.

    Prefers OncoTree invasive breast carcinoma (BRCA) when that column exists
    and has matches; otherwise falls back to lineage == Breast.
    """
    cols = {c.lower(): c for c in models.columns}

    def col(*names):
        for n in names:
            if n.lower() in cols:
                return cols[n.lower()]
        return None

    code = col("OncotreeCode", "OncotreeSubtypeCode")
    disease = col("OncotreePrimaryDisease", "PrimaryDisease", "OncotreeSubtype")
    lineage = col("OncotreeLineage", "Lineage")

    invasive = pd.Series(False, index=models.index)
    if code is not None:
        invasive |= models[code].astype(str).str.upper().eq("BRCA")
    if disease is not None:
        text = models[disease].astype(str)
        invasive |= text.str.contains("invasive breast", case=False, na=False)
        invasive |= text.str.contains(r"breast invasive", case=False, na=False)

    breast = pd.Series(False, index=models.index)
    if lineage is not None:
        breast |= models[lineage].astype(str).str.contains("breast", case=False, na=False)
    breast |= invasive

    if prefer_invasive and invasive.any():
        return invasive, "invasive / OncoTree BRCA"
    return breast, "Breast lineage"


def expression_to_genes_by_samples(expr):
    """DepMap files are usually samples × genes; Panel A wants genes × samples."""
    expr = expr.copy()
    expr.index = expr.index.astype(str)
    expr.columns = [strip_entrez(c) for c in expr.columns]
    # If columns look like ModelIDs (ACH-…) the file is genes × samples already.
    sample_like = sum(str(c).startswith("ACH-") for c in expr.columns)
    gene_like = sum(str(i).startswith("ACH-") for i in expr.index)
    if gene_like > sample_like:
        expr = expr.T
        expr.index = [strip_entrez(i) for i in expr.index]
        expr.columns = expr.columns.astype(str)
    return expr


def drop_low_genes(expr, max_log_threshold=1.0):
    """Drop genes with log2(TPM+1) < threshold in every sample.

    Groves dropped genes with all values < 2 before their log2; DepMap is already
    log2(TPM+1), so threshold=1 means TPM < 1 in every line.
    """
    keep = expr.max(axis=1) >= max_log_threshold
    return expr.loc[keep]
