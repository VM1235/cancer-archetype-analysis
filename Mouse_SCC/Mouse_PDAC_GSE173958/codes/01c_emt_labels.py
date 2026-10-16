#!/usr/bin/env python3
"""Step 1c: EMT labels for the tumour cells (GEO has none, so they are derived here).

These are NOT the paper's clusters. They are:
  KS_EMT            mean z of the mouse Tan-2014 KS mesenchymal genes minus the epithelial
                    genes (126 epi + 42 mes genes present), z-scored over all tumour cells.
  cluster_annotation  4 ordered states from KS_EMT quantiles over IN-VIVO tumour cells:
                    Epithelial (lowest 40%), Early-Hybrid (40-75%), Late-Hybrid (75-95%),
                    Mesenchymal (top 5%). In-vitro cells are binned with the same cut-offs.
  clusterU          KMeans (k=5, seed 0) on the first 5 PCs of the KS genes: K0..K4.
The fits never use these labels; they are only read off afterwards, as in the SCC work.
"""
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import scc_common as C

P = C.PROC
a = ad.read_h5ad(P / "expr_tumour.h5ad")
epi, mes = C.ks_genes_mouse(a.var_names)
def zm(gs):
    M = a[:, gs].X.toarray()
    return ((M - M.mean(0)) / np.where(M.std(0) > 0, M.std(0), 1)).mean(1)
a.obs["KS_EMT"] = zm(mes) - zm(epi)
iv = (a.obs["assay"] == "in_vivo").values
cuts = np.quantile(a.obs["KS_EMT"][iv], [0.40, 0.75, 0.95])
a.obs["cluster_annotation"] = pd.Categorical(
    np.array(C.STATE_ORDER)[np.searchsorted(cuts, a.obs["KS_EMT"].values)],
    categories=C.STATE_ORDER, ordered=True)
M = a[:, epi + mes].X.toarray(); M -= M.mean(0)
Z = PCA(5, random_state=0).fit_transform(M)
a.obs["clusterU"] = ["K%d" % i for i in KMeans(5, n_init=5, random_state=0).fit(Z).labels_]
a.write_h5ad(P / "expr_tumour.h5ad")
print("cut-offs", cuts.round(3))
print(pd.crosstab(a.obs["sample"], a.obs["cluster_annotation"], margins=True))
print(pd.crosstab(a.obs["sample"], a.obs["clusterU"], margins=True))
a.obs.rename_axis("barcode").to_csv(P / "cell_metadata.csv")   # read by 02b and 06
