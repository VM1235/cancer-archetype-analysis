#!/usr/bin/env python3
"""Step 1: read the 12 transcriptome libraries of GSE173958, QC, normalise, save.

Libraries: M1-PT/SS/Met/Liver/Lung/CTC, M2-PT/Met/Liver/Lung (two mice, several
harvest sites) and iv-5cell / iv-100cell (cultured cells, in vitro).

QC per cell (the raw 10x matrices still hold many empty droplets/dead cells):
  >= MIN_GENES genes detected, mitochondrial fraction <= MAX_MT, <= MAX_GENES genes
  (the last one is a crude doublet cut).
Genes: detected in >=1% of the kept cells; mito and ribosomal genes dropped.
Expression: CP10k then log1p, as in the SCC pipeline.
Batch unit ("sample") = library, the analogue of the SCC tumour.
"""
import gzip, re, glob
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp
import anndata as ad, scanpy as sc

HERE = Path(__file__).resolve().parent.parent
RAW = HERE / "data" / "raw" / "tx"
OUT = HERE / "data" / "processed"
MIN_GENES, MAX_GENES, MAX_MT = 1000, 7000, 0.15

S_GENES = "Mcm5 Pcna Tyms Fen1 Mcm2 Mcm4 Rrm1 Ung Gins2 Mcm6 Cdca7 Dtl Prim1 Uhrf1 Mlf1ip Hells Rfc2 Rap2c Rad51 Rrm2 Cdc45 Cdc6 Exo1 Tipin Dscc1 Blm Casp8ap2 Usp1 Clspn Pola1 Chaf1b Brip1 E2f8".split()
G2M_GENES = "Hmgb2 Ddx11 Cdk1 Nusap1 Ube2c Birc5 Tpx2 Top2a Ndc80 Cks2 Nuf2 Cks1b Mki67 Tmpo Cenpk Tacc3 Smc4 Ccnb2 Ckap2l Ckap2 Aurkb Bub1 Kif11 Anp32e Tubb4b Gtse1 Kif20b Hjurp Cdca3 Cdc20 Ttk Cdc25c Kif2c Rangap1 Ncapd2 Dlgap5 Cdca2 Cdca8 Ect2 Kif23 Hmmr Aurka Psrc1 Anln Lbr Ckap5 Cenpe Ctcf Nek2 G2e3 Gas2l3 Cbx5 Cenpa".split()

feat = [l.rstrip("\n").split("\t")[1] for l in gzip.open(next(RAW.glob("*M1-PT-features*")), "rt")]
dup = pd.Index(feat).duplicated(keep=False)          # 35 symbols are shared by several Ensembl ids
feat = np.array([f"{g}__{i}" if d else g for i, (g, d) in enumerate(zip(feat, dup))])
mt = np.array([g.startswith("mt-") for g in feat])

parts = []
for f in sorted(RAW.glob("*-matrix.mtx.gz")):
    name = re.search(r"GSM\d+_(.*)-matrix", f.name).group(1)
    bc = [l.strip() for l in gzip.open(str(f).replace("matrix.mtx", "barcodes.tsv"), "rt")]
    m = sio.mmread(f).tocsc()
    n = np.asarray(m.sum(0)).ravel()
    ng = np.diff(m.indptr)
    pm = np.asarray(m[mt].sum(0)).ravel() / np.maximum(n, 1)
    keep = (ng >= MIN_GENES) & (ng <= MAX_GENES) & (pm <= MAX_MT)
    mouse = "invitro" if name.startswith("iv") else name.split("-")[0]
    site = "culture" if name.startswith("iv") else name.split("-", 1)[1]
    obs = pd.DataFrame({"sample": name, "mouse": mouse, "site": site,
                        "assay": "in_vitro" if mouse == "invitro" else "in_vivo",
                        "nCount_RNA": n[keep], "nGene": ng[keep], "pct_counts_mt": 100 * pm[keep]},
                       index=[f"{name}_{b}" for b, k in zip(bc, keep) if k])
    parts.append(ad.AnnData(X=m[:, keep].T.tocsr().astype(np.float32), obs=obs,
                            var=pd.DataFrame(index=feat)))
    print(f"{name}: {m.shape[1]} droplets -> {keep.sum()} cells", flush=True)

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
a.write_h5ad(OUT / "expr_log_filtered.h5ad")
print(pd.crosstab(a.obs["sample"], a.obs["phase"], margins=True))
