"""Shared pieces for the mouse PDAC (GSE173958) archetype pipeline, adapted from Mouse_SCC/codes/scc_common.py.

Batch unit ("sample") = 10x library (M1-PT, M2-Met, ...); "assay" = in_vivo / in_vitro.

Everything that more than one script needs lives here, so each step is
done one way only: loading the matrix, choosing genes, choosing fit cells,
PCA, and turning PC coordinates into archetype weights.

The PCHA fit itself comes from the repo engine (src/archetypes.py), the same
code used for the breast work, so numbers are computed the same way.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

SCC = Path(__file__).resolve().parent.parent          # Mouse_PDAC_GSE173958/
ROOT = SCC.parent.parent                                      # Project_1/
sys.path.insert(0, str(ROOT))
from src.archetypes import fit_pcha_best, simplex_volume, _shuffle_columns  # noqa: E402,F401

PROC = SCC / "data" / "processed"
import os
EXPR = PROC / os.environ.get("EXPR_NAME", "expr_tumour.h5ad")   # expr_tumour_clean.h5ad for the cleaned pooled fit
KS_FILE = ROOT / "Breast Cancer" / "data" / "processed" / "ks_cellline_signature_tan2014.csv"
RESULTS = SCC / "results"

SEED = 0
MAX_PER_TUMOUR = 1500     # sample-balanced fit set (4 samples: 4,241 to 11,908 cells each)
N_PCS_KEEP = 50           # PCs stored per cell; the fit uses only the first k-1
K_VALUES = (3, 4, 5, 6, 7)
N_INIT = 15               # PCHA starts for the real fit (same as the Pal fits)
N_INIT_NULL = 5           # PCHA starts per shuffled dataset
PURITY_CUT = 0.35         # specialist rule from the Pal work

STATE_ORDER = ["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"]
CLUSTERU_ORDER = ["K0", "K1", "K2", "K3", "K4"]


# ------------------------------------------------------------------ data
def load_adata():
    import anndata as ad
    a = ad.read_h5ad(EXPR)
    if "counts" in a.layers:          # raw counts are not needed here; free the memory
        del a.layers["counts"]
    a.X = sp.csr_matrix(a.X, dtype=np.float32)
    return a


def ks_genes_mouse(var_names):
    """Tan 2014 KS genes; GSE114687 is human, so the symbols are used unchanged (CDH1 -> CDH1)."""
    s = pd.read_csv(KS_FILE)
    mouse = s["gene_symbol"]   # human data: symbols used as they are
    present = set(var_names)
    cat = s["category"].str.lower().str[:3]
    epi = [g for g, c in zip(mouse, cat) if c == "epi" and g in present]
    mes = [g for g, c in zip(mouse, cat) if c == "mes" and g in present]
    return epi, mes


def hallmark_emt_genes(var_names, human=True):
    """MSigDB Hallmark_2020 'Epithelial Mesenchymal Transition' genes present in the data
    (mouse symbols by capitalisation, human as they are). Added for the Hallmark-EMT pooled fits."""
    from pathlib import Path as _P
    root = next(p for p in _P(__file__).resolve().parents if (p / "Mouse_SCC" / "hallmark_data").exists())
    for line in open(root / "Mouse_SCC" / "hallmark_data" / "MSigDB_Hallmark_2020.txt"):
        parts = line.rstrip("\n").split("\t")
        if parts[0] == "Epithelial Mesenchymal Transition":
            genes = [g for g in parts[2:] if g]
            break
    if not human:
        genes = [g[0] + g[1:].lower() for g in genes]
    present = set(var_names)
    return [g for g in genes if g in present]


def gene_index(adata, geneset):
    """Column indices of the genes used for PCA."""
    if geneset == "all":
        return np.arange(adata.n_vars)
    if geneset == "ks":
        epi, mes = ks_genes_mouse(adata.var_names)
        names = set(epi + mes)
        return np.array([i for i, g in enumerate(adata.var_names) if g in names])
    if geneset == "hallmark":
        names = set(hallmark_emt_genes(adata.var_names))
        return np.array([i for i, g in enumerate(adata.var_names) if g in names])
    raise ValueError("geneset must be 'all' or 'ks'")


def balanced_fit_cells(obs, max_per_tumour=MAX_PER_TUMOUR, exclude_tumours=(),
                       only_assay=None, seed=SEED):
    """Row indices of the fit cells: at most max_per_tumour random cells per tumour."""
    df = pd.DataFrame({"i": np.arange(len(obs)), "t": obs["sample"].values,
                       "a": obs["assay"].values})
    df = df[~df["t"].isin(list(exclude_tumours))]
    if only_assay is not None:
        df = df[df["a"] == only_assay]
    parts = [g.sample(min(len(g), max_per_tumour), random_state=seed)
             for _, g in df.groupby("t", sort=True, observed=True)]
    return np.sort(pd.concat(parts)["i"].values)


# ------------------------------------------------------------------ PCA
def phase_means(Xfit_dense, phase_fit):
    """Per-phase gene means minus the overall mean (used to remove cell cycle)."""
    overall = Xfit_dense.mean(0)
    return {p: Xfit_dense[phase_fit == p].mean(0) - overall for p in np.unique(phase_fit)}


def fit_pca(Xfit_dense, n_pcs=N_PCS_KEEP, method="gram", n_top=200):
    """Exact PCA on the fit cells. Genes centred, not scaled.

    method="svd":  full SVD (slow for thousands of genes; computes every component).
    method="gram": exact top-n_top components from the eigenvectors of the
                   cells x cells matrix Xc Xc^T. Same PCs as the SVD, much faster
                   because only the leading components are computed. Total
                   variance (for % explained) is the exact sum of squares.
    Sign rule: the largest-|loading| gene of each PC is positive, so reruns
    give the same orientation.
    """
    from scipy.linalg import eigh
    mean = Xfit_dense.mean(0)
    Xc = Xfit_dense - mean
    m = Xc.shape[0]
    if method == "svd" or min(Xc.shape) <= n_top:
        U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        total = (S ** 2).sum()
    else:
        G = Xc @ Xc.T
        w, U = eigh(G, subset_by_index=[m - n_top, m - 1])
        w, U = w[::-1], U[:, ::-1]
        S = np.sqrt(np.clip(w, 0, None))
        Vt = (U.T @ Xc) / S[:, None]
        total = float((Xc ** 2).sum())
        del G
    for j in range(Vt.shape[0]):
        if Vt[j, np.argmax(np.abs(Vt[j]))] < 0:
            Vt[j] *= -1
            U[:, j] *= -1
    evr = S ** 2 / total
    n = min(n_pcs, Vt.shape[0])
    return dict(mean=mean, components=Vt[:n].copy(), evr=evr, scores_fit=U[:, :n] * S[:n])


def project(X_sparse, pca, shift=None):
    """PC scores for any cells: (x - mean - shift) @ components.T, without densifying X."""
    V = pca["components"]
    sc = np.asarray(X_sparse @ V.T) - pca["mean"] @ V.T
    if shift is not None:
        sc -= shift @ V.T
    return sc


# ------------------------------------------------------------------ simplex
def barycentric(points, vertices):
    """Barycentric weights of points (n x d) w.r.t. a simplex with d+1 vertices."""
    a = np.vstack([vertices.T, np.ones((1, vertices.shape[0]))])
    b = np.vstack([np.asarray(points).T, np.ones((1, len(points)))])
    return np.linalg.solve(a, b).T


def composition(w):
    """Clip negative weights to 0 and renormalise (cells just outside the simplex)."""
    w = np.clip(np.asarray(w, float), 0, None)
    tot = w.sum(1, keepdims=True)
    tot[tot == 0] = 1.0
    return w / tot


def purity_and_nearest(w):
    c = composition(w)
    s = np.sort(c, 1)
    return s[:, -1] - s[:, -2], np.argmax(c, 1)


def fit_simplex(scores, k, n_init=N_INIT, seed=SEED):
    """Best-of-n_init PCHA on the first k-1 PCs (repo engine). Returns vertices and t-ratio."""
    from scipy.spatial import ConvexHull
    np.random.seed(seed)                       # py_pcha uses the global RNG
    arcs, _, esv, vol, n_ok = fit_pcha_best(scores, k, n_init=n_init, delta=0.0)
    hull = ConvexHull(np.asarray(scores)[:, : k - 1]).volume
    return dict(k=k, archetypes=np.asarray(arcs), esv=float(esv), volume=float(vol),
                t_ratio=float(vol / hull), n_ok=int(n_ok))


def match_vertices(ref_arcs, new_arcs):
    """Order of new vertices that best matches ref vertices (min total distance)."""
    from scipy.optimize import linear_sum_assignment
    d = np.linalg.norm(ref_arcs[:, None, :] - new_arcs[None, :, :], axis=2)
    r, c = linear_sum_assignment(d)
    return c[np.argsort(r)]


# ------------------------------------------------------------------ io
def run_dir(geneset, tag="main"):
    d = RESULTS / f"{geneset}_{tag}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def save_json(obj, path):
    Path(path).write_text(json.dumps(obj, indent=2, default=lambda o: o.tolist()
                                     if isinstance(o, np.ndarray) else str(o)))


def eta_squared(values, groups):
    """Fraction of the variance of `values` explained by a categorical grouping."""
    v = np.asarray(values, float)
    g = pd.Series(groups).astype(str).values
    total = ((v - v.mean()) ** 2).sum()
    between = sum(((v[g == x].mean() - v.mean()) ** 2) * (g == x).sum() for x in np.unique(g))
    return float(between / total) if total > 0 else np.nan
