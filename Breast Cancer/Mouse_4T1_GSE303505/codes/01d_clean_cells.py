#!/usr/bin/env python3
"""Step 1d: cleaned copy of the tumour cells, for the cleaned pooled fit.

Same rule as the cleaned within-clone fits (08_within_unit.py --clean --decontam), applied once to all clones:
  - drop cells with < 10,000 UMI (low-depth cells carried the immune/stromal and stress corners)
  - drop cells whose immune/stromal panel score (mean log-expression of Ptprc Lyz2 C1qa C1qb Cd3e Cd79a
    Itgam Adgre1 Col1a1 Col3a1 Dcn Lum Pecam1) is > 0.3
  - drop technical genes (Malat1, Neat1, Gm26870, Gm42418, AY036118, haemoglobin)
Expression values are left as they were (log1p CP10k computed before the cell filter).
Output: data/processed/expr_tumour_clean.h5ad  (use with EXPR_NAME=expr_tumour_clean.h5ad)
"""
from pathlib import Path
import numpy as np, anndata as ad
P = Path(__file__).resolve().parent.parent / "data" / "processed"
PANEL = "Ptprc Lyz2 C1qa C1qb Cd3e Cd79a Itgam Adgre1 Col1a1 Col3a1 Dcn Lum Pecam1".split()
TECH = "Hbb-bs Hbb-bt Hba-a1 Hba-a2 Malat1 Neat1 Gm26870 Gm42418 AY036118".split()
a = ad.read_h5ad(P / "expr_tumour.h5ad")
pi = [a.var_names.get_loc(g) for g in PANEL if g in a.var_names]
score = np.asarray(a.X[:, pi].mean(1)).ravel()
a.obs["panel_score"] = score
keep = (a.obs["nCount_RNA"].values >= 10000) & (score <= 0.3)
print(a.obs.assign(keep=keep).groupby("sample", observed=True)["keep"].agg(["size", "sum"]).rename(columns={"size": "before", "sum": "after"}))
a = a[keep, [g not in TECH for g in a.var_names]].copy()
a.write_h5ad(P / "expr_tumour_clean.h5ad")
print(a.shape)
