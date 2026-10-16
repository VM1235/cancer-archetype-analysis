#!/usr/bin/env python3
"""Step 1b: keep only tumour cells.

The 10x libraries hold everything that was captured: tumour cells plus immune,
endothelial, fibroblast, hepatocyte and red blood cells (M1-CTC is ~all blood).
Archetype analysis of "EMT" only makes sense on tumour cells, so:
  1. HVG(3000) -> scale -> PCA(30) on all QC-passing cells
  2. KMeans (k=40, seed 0) on the PCs
  3. per cluster, mean log-expression of marker panels:
       tumour epithelial  Krt8 Krt18 Krt19 Epcam Cdh1 Krt7
       immune             Ptprc Lyz2 Cd3e Cd79a S100a8 Cd74
       endothelial        Pecam1 Cdh5
       fibroblast         Dcn Lum Pdgfra Col1a1
  4. a cluster is tumour if epithelial >= 0.85 and immune, endothelial,
     fibroblast are all < 0.3
The cluster table is written so the rule can be audited.
Output: data/processed/expr_tumour.h5ad (+ cluster_table.csv), obs gets 'km' and 'is_tumour'.
"""
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
from sklearn.cluster import KMeans

P = Path(__file__).resolve().parent.parent / "data" / "processed"
PANEL = {"epi": "Krt8 Krt18 Krt19 Epcam Cdh1 Krt7", "immune": "Ptprc Lyz2 Cd3e Cd79a S100a8 Cd74",
         "endo": "Pecam1 Cdh5", "fibro": "Dcn Lum Pdgfra Col1a1"}
a = ad.read_h5ad(P / "expr_log_filtered.h5ad")
b = a.copy()
sc.pp.highly_variable_genes(b, n_top_genes=3000)
b = b[:, b.var.highly_variable].copy()
sc.pp.scale(b, max_value=10)
sc.tl.pca(b, n_comps=30, random_state=0)
a.obs["km"] = KMeans(40, n_init=3, random_state=0).fit(b.obsm["X_pca"]).labels_
del b
rows = []
for c in range(40):
    m = a.obs["km"].values == c
    r = {"cluster": c, "n": int(m.sum())}
    for k, gs in PANEL.items():
        r[k] = round(float(a[m][:, [g for g in gs.split() if g in a.var_names]].X.mean()), 3)
    r["top_libraries"] = ", ".join(f"{i}:{100 * v:.0f}%" for i, v in
                                   a.obs["sample"][m].value_counts(normalize=True).head(3).items())
    rows.append(r)
tab = pd.DataFrame(rows)
tab["is_tumour"] = (tab.epi >= 0.85) & (tab.immune < 0.3) & (tab.endo < 0.3) & (tab.fibro < 0.3)
tab.to_csv(P / "cluster_table.csv", index=False)
print(tab.to_string(index=False))
a.obs["is_tumour"] = a.obs["km"].isin(tab.loc[tab.is_tumour, "cluster"]).values
print(pd.crosstab(a.obs["sample"], a.obs["is_tumour"], margins=True))
t = a[a.obs["is_tumour"].values].copy()
t.obs["km"] = t.obs["km"].astype(int)
t.write_h5ad(P / "expr_tumour.h5ad")
