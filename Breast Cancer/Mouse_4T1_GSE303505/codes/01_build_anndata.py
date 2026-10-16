#!/usr/bin/env python3
"""Step 1: read the six 'Total tumour' libraries of GSE303505 (Pattabiraman lab, 4T1 clones),
keep the tumour cells, QC, normalise, save.

Libraries: one tumour per clone, E, EM1, EM2, EM3, M (single-cell-derived 4T1 clones along an
epithelial -> mesenchymal spectrum) and the parental 4T1. The 'Total' samples hold tumour
cells plus immune/stromal cells, so the authors' own `Annotation` column in the GEO metadata
is used to keep only 'Tumor' cells (no re-clustering needed). CD45 libraries are not used.
GEO counts are already filtered by the authors (cells x genes), so only a light QC is added:
  nFeature >= MIN_GENES and mitochondrial % <= MAX_MT.
Genes: detected in >=1% of the kept cells; mito and ribosomal genes dropped.
Expression: CP10k then log1p, as in the SCC/PDAC pipelines.
Batch unit ("sample") = clone/tumour, the analogue of the SCC tumour. There is ONE tumour per
clone, so clone and tumour identity cannot be separated.
"""
import re
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad, scanpy as sc

HERE = Path(__file__).resolve().parent.parent
RAW = HERE / "data" / "raw"
OUT = HERE / "data" / "processed"
MIN_GENES, MAX_MT = 500, 15.0
CLONES = ["E", "EM1", "EM2", "EM3", "M", "4T1"]          # authors' epithelial -> mesenchymal order (4T1 = parental)

S_GENES = "Mcm5 Pcna Tyms Fen1 Mcm2 Mcm4 Rrm1 Ung Gins2 Mcm6 Cdca7 Dtl Prim1 Uhrf1 Mlf1ip Hells Rfc2 Rap2c Rad51 Rrm2 Cdc45 Cdc6 Exo1 Tipin Dscc1 Blm Casp8ap2 Usp1 Clspn Pola1 Chaf1b Brip1 E2f8".split()
G2M_GENES = "Hmgb2 Ddx11 Cdk1 Nusap1 Ube2c Birc5 Tpx2 Top2a Ndc80 Cks2 Nuf2 Cks1b Mki67 Tmpo Cenpk Tacc3 Smc4 Ccnb2 Ckap2l Ckap2 Aurkb Bub1 Kif11 Anp32e Tubb4b Gtse1 Kif20b Hjurp Cdca3 Cdc20 Ttk Cdc25c Kif2c Rangap1 Ncapd2 Dlgap5 Cdca2 Cdca8 Ect2 Kif23 Hmmr Aurka Psrc1 Anln Lbr Ckap5 Cenpe Ctcf Nek2 G2e3 Gas2l3 Cbx5 Cenpa".split()

parts, allann = [], []
for clone in CLONES:
    fm = next(RAW.glob(f"*_{clone}_Total_metadata.tsv.gz"))
    fc = next(RAW.glob(f"*_{clone}_Total_filtered_counts.tsv.gz"))
    meta = pd.read_csv(fm, sep="\t", index_col=0)
    allann.append(meta["Annotation"].value_counts().rename(clone))
    tum = meta.index[meta["Annotation"] == "Tumor"]
    # read counts for tumour columns only
    hdr = pd.read_csv(fc, sep="\t", nrows=0, index_col=0).columns
    cols = [h for h in hdr if h in set(tum)]
    m = pd.read_csv(fc, sep="\t", index_col=0, usecols=[0] + [i + 1 for i, h in enumerate(hdr) if h in set(tum)],
                    dtype={c: np.float32 for c in cols})
    X = sp.csr_matrix(m.values.T.astype(np.float32))
    genes = m.index.astype(str)
    mt = np.asarray([g.startswith("mt-") for g in genes])
    n = np.asarray(X.sum(1)).ravel()
    ng = np.diff(X.indptr)
    pm = 100 * np.asarray(X[:, np.where(mt)[0]].sum(1)).ravel() / np.maximum(n, 1)
    keep = (ng >= MIN_GENES) & (pm <= MAX_MT)
    obs = pd.DataFrame({"sample": clone, "mouse": clone, "assay": "in_vivo", "nCount_RNA": n[keep],
                        "nGene": ng[keep], "pct_counts_mt": pm[keep]},
                       index=[f"{clone}_{c}" for c, k in zip(m.columns, keep) if k])
    parts.append(ad.AnnData(X=X[keep], obs=obs, var=pd.DataFrame(index=genes)))
    print(f"{clone}: {len(meta)} cells in file, {len(tum)} Tumor, {keep.sum()} pass QC", flush=True)

pd.concat(allann, axis=1).fillna(0).astype(int).to_csv(OUT / "annotation_counts_by_clone.csv")
a = ad.concat(parts, join="inner")
del parts
det = np.asarray((a.X > 0).sum(0)).ravel() >= 0.01 * a.n_obs
bad = a.var_names.str.match(r"^(mt-|Mt-)") | a.var_names.str.match(r"^(Rpl|Rps|Mrpl|Mrps)\d")
a = a[:, det & ~np.asarray(bad)].copy()
print("cells x genes:", a.shape)
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
sc.tl.score_genes_cell_cycle(a, s_genes=[g for g in S_GENES if g in a.var_names],
                             g2m_genes=[g for g in G2M_GENES if g in a.var_names])
a.obs["sample"] = pd.Categorical(a.obs["sample"], categories=CLONES, ordered=True)
a.write_h5ad(OUT / "expr_tumour.h5ad")
print(pd.crosstab(a.obs["sample"], a.obs["phase"], margins=True))
