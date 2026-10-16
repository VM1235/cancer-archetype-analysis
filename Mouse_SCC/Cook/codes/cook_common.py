"""Shared pieces for the Cook (GSE147405) archetype analysis.

Same algorithms as the mouse pipeline: the PCA, PCHA fit, t-ratio and shuffle
helpers are imported from Mouse_SCC/codes/scc_common.py. Only the paths, the
cell selection and the gene set (human KS genes, no case conversion) differ.
"Tumour" in the mouse code corresponds to cell line here, and "EMT state"
corresponds to the time point after TGF-B1.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

COOK = Path(__file__).resolve().parent.parent            # Mouse_SCC/Cook/
sys.path.insert(0, str(COOK.parent / "codes"))
import scc_common as S                                    # noqa: E402
from scc_common import (fit_pca, project, fit_simplex, eta_squared, barycentric,  # noqa: E402,F401
                        composition, purity_and_nearest, match_vertices, save_json,
                        _shuffle_columns, fit_pcha_best, simplex_volume, ROOT, SEED,
                        N_INIT, N_INIT_NULL, PURITY_CUT)

RAW = COOK / "data" / "raw"
PROC = COOK / "data" / "processed"
EXPR = PROC / "expr_log_filtered.h5ad"
RESULTS = COOK / "results"
KS_FILE = S.KS_FILE
K_VALUES = (3, 4, 5, 6, 7)
LINES = ["A549", "DU145", "MCF7", "OVCA420"]
TIME_ORDER = ["0d", "8h", "1d", "3d", "7d", "8h_rm", "1d_rm", "3d_rm"]
MAX_PER_LINE = 550            # pooled fit: smallest line has 558 forward-arm cells


def load_adata():
    import anndata as ad
    import scipy.sparse as sp
    a = ad.read_h5ad(EXPR)
    a.X = sp.csr_matrix(a.X, dtype=np.float32)
    return a


def ks_genes(var_names):
    s = pd.read_csv(KS_FILE)
    present = set(var_names)
    cat = s["category"].str.lower().str[:3]
    epi = [g for g, c in zip(s["gene_symbol"], cat) if c == "epi" and g in present]
    mes = [g for g, c in zip(s["gene_symbol"], cat) if c == "mes" and g in present]
    return epi, mes


def gene_index(adata, geneset):
    if geneset == "all":
        return np.arange(adata.n_vars)
    if geneset == "ks":
        epi, mes = ks_genes(adata.var_names)
        names = set(epi + mes)
        return np.array([i for i, g in enumerate(adata.var_names) if g in names])
    raise ValueError(geneset)


def fit_cells(obs, lines=None, max_per_line=None, include_rm=False, seed=SEED):
    """Row indices of the fit cells: forward time course only (removal arms are projected, not fitted)."""
    df = pd.DataFrame({"i": np.arange(len(obs)), "l": obs["cell_line"].values,
                       "rm": obs["time"].astype(str).str.endswith("_rm").values})
    if lines:
        df = df[df["l"].isin(lines)]
    if not include_rm:
        df = df[~df["rm"]]
    parts = [g.sample(min(len(g), max_per_line), random_state=seed) if max_per_line else g
             for _, g in df.groupby("l", sort=True)]
    return np.sort(pd.concat(parts)["i"].values)


def run_dir(geneset, tag):
    d = RESULTS / f"{geneset}_{tag}"
    d.mkdir(parents=True, exist_ok=True)
    return d
