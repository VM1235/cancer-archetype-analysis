#!/usr/bin/env python3
"""Cleaned copies of the pooled expression files, so that pooled cleaned fits can be run for Cook, PDAC and GSE114687
(4T1 already has expr_tumour_clean.h5ad from 01d_clean_cells.py). Same rules as each dataset's within-unit cleaning arm:
  Cook     >= 10,000 UMI and mito fraction <= 0.07; drop human technical genes (MALAT1 NEAT1 HBB HBA1 HBA2)
  PDAC     >= 10,000 UMI; drop the mouse technical list (hepatocyte ambient RNA, haemoglobin, Malat1, Neat1, Gm26870, Gm42418, AY036118)
  GSE      >= 10,000 UMI; drop human technical genes
Run from Project_1/:  .venv/bin/python docs/emt_single_cell/best_simplex/make_clean_pooled.py
"""
from pathlib import Path
import anndata as ad, numpy as np
ROOT = Path(__file__).resolve().parents[3]
HUMAN = "MALAT1 NEAT1 HBB HBA1 HBA2".split()
MOUSE = ("Alb Apoa1 Apoa2 Apoh Fga Fgb Fgg Hp Saa1 Saa2 Serpina1a Serpina1b Serpina1c Serpina1d "
         "Serpina1e Hbb-bs Hbb-bt Hba-a1 Hba-a2 Malat1 Neat1 Gm26870 Gm42418 AY036118").split()
JOBS = [("Mouse_SCC/Cook/data/processed", "expr_log_filtered.h5ad", "expr_log_filtered_clean.h5ad", HUMAN, 0.07, "pct_mito"),
        ("Mouse_SCC/Mouse_PDAC_GSE173958/data/processed", "expr_tumour.h5ad", "expr_tumour_clean.h5ad", MOUSE, None, None),
        ("Breast Cancer/GSE114687_pseudospace/data/processed", "expr_tumour.h5ad", "expr_tumour_clean.h5ad", HUMAN, None, None)]
for d, src, dst, tech, mito, mcol in JOBS:
    a = ad.read_h5ad(ROOT / d / src)
    keep = (a.obs["nCount_RNA"].values >= 10000)
    if mito is not None: keep &= (a.obs[mcol].values <= mito)
    genes = ~a.var_names.isin(tech)
    b = a[keep, genes].copy()
    b.write_h5ad(ROOT / d / dst)
    print(d.split("/")[-3], f"cells {a.n_obs} -> {b.n_obs} ({100 * b.n_obs / a.n_obs:.0f}%), genes {a.n_vars} -> {b.n_vars}", flush=True)
