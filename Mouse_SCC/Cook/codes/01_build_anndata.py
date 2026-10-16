#!/usr/bin/env python3
"""Step 1: read the four TGF-B1 time-course matrices, filter, normalise, save one h5ad.

Genes: present in all four lines, detected in >=1% of cells overall, mito (MT-)
and ribosomal (RPL/RPS) genes dropped. Normalise to 10,000 UMI per cell, log1p
(same as the mouse build). Usage (from Project_1/):
  .venv/bin/python -u Mouse_SCC/Cook/codes/01_build_anndata.py
"""
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc
from pathlib import Path
RAW = Path("Mouse_SCC/Cook/data/raw"); PROC = Path("Mouse_SCC/Cook/data/processed"); PROC.mkdir(parents=True, exist_ok=True)

parts = []
for L in ["A549", "DU145", "MCF7", "OVCA420"]:
    m = pd.read_csv(RAW / f"GSE147405_{L}_TGFB1_TimeCourse_UMI_matrix.csv.gz", index_col=0)
    meta = pd.read_csv(RAW / f"GSE147405_{L}_TGFB1_TimeCourse_metadata.csv.gz", index_col=0)
    assert list(m.columns) == list(meta.index), L
    print(L, m.shape, flush=True)
    parts.append((m, meta))
genes = sorted(set.intersection(*[set(m.index) for m, _ in parts]))
print("genes shared by all four lines:", len(genes))
X = sp.vstack([sp.csr_matrix(m.loc[genes].values.T.astype(np.float32)) for m, _ in parts]).tocsr()
obs = pd.concat([meta for _, meta in parts])
obs = obs.rename(columns={"CellLine": "cell_line", "Time": "time", "Phase": "phase", "Mix": "mix",
                          "percent.mito": "pct_mito", "Pseudotime": "pseudotime"})
obs = obs[["cell_line", "time", "Sample", "mix", "phase", "nCount_RNA", "nFeature_RNA", "pct_mito",
           "S.Score", "G2M.Score", "pseudotime"]]
obs["removal_arm"] = obs["time"].str.endswith("_rm")
a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
det = np.asarray((a.X > 0).sum(0)).ravel() >= 0.01 * a.n_obs
bad = a.var_names.str.match(r"^MT-") | a.var_names.str.match(r"^(RPL|RPS|MRPL|MRPS)\d")
print("genes >=1% cells:", int(det.sum()), "; mito/ribo removed among them:", int((det & bad).sum()))
a = a[:, det & ~bad].copy()
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
print("final:", a.shape)
print(pd.crosstab(a.obs["cell_line"], a.obs["time"])[["0d","8h","1d","3d","7d","8h_rm","1d_rm","3d_rm"]])
print(pd.crosstab(a.obs["time"], a.obs["mix"]))
a.write_h5ad(PROC / "expr_log_filtered.h5ad")
