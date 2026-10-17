#!/usr/bin/env python3
"""Hallmark check of corner names for the two datasets that never had one: mouse SCC and 4T1.
Same method as Mouse_SCC/hallmark_check and Breast Cancer/GSE114687_pseudospace/codes/09_hallmark_corners.py:
the 10% of cells nearest each corner are scored (mean z of the set's genes, z within the fit's cells) on MSigDB Hallmark
sets; the top set per corner is reported. Mouse symbols = capitalised human symbols.
From Project_1/:  .venv/bin/python docs/cleaning_audit/hallmark_scc_4t1.py   -> docs/cleaning_audit/hallmark_scc_4t1.csv
"""
import os, numpy as np, pandas as pd, anndata as ad
from pathlib import Path
R = Path(__file__).resolve().parents[2]
H = {}
for line in open(R / "Mouse_SCC/hallmark_data/MSigDB_Hallmark_2020.txt"):
    f = line.rstrip("\n").split("\t"); H[f[0].strip()] = [g for g in f[2:] if g]
USE = {"Hypoxia": "Hypoxia", "EMT": "Epithelial Mesenchymal Transition", "G2-M": "G2-M Checkpoint", "E2F": "E2F Targets",
       "TGF-b": "TGF-beta Signaling", "ApicalJunction": "Apical Junction", "KRAS-up": "KRAS Signaling Up"}
SCC = R / "Mouse_SCC"; T4 = R / "Breast Cancer/Mouse_4T1_GSE303505"
# (label, adata, result folder, k)
FITS = [("SCC pooled all genes k3", SCC / "data/processed/expr_log_filtered.h5ad", SCC / "results/all_main", 3),
        ("SCC scRNA-only all genes k3", SCC / "data/processed/expr_log_filtered.h5ad", SCC / "results/all_scrna", 3),
        ("SCC six clean scRNA tumours k3", SCC / "data/processed/expr_log_filtered.h5ad", SCC / "results/all_scrna6", 3),
        ("4T1 pooled light-QC all genes k4", T4 / "data/processed/expr_tumour.h5ad", T4 / "results/all_main", 4),
        ("4T1 pooled cleaned all genes k3", T4 / "data/processed/expr_tumour_clean.h5ad", T4 / "results/all_clean", 3)]
for u in ("EM1", "EM2", "EM3", "M", "4T1"):
    FITS.append((f"4T1 within {u} cleaned k3", T4 / "data/processed/expr_tumour.h5ad", T4 / f"results/within_unit/all_{u}_clean_dc", 3))
cache = {}; rows = []
for label, ap, d, k in FITS:
    if not (d / "pca.npz").exists() or not (d / f"archetypes_k{k}.npy").exists():
        print("skip (missing files):", label); continue
    if ap not in cache: cache[ap] = ad.read_h5ad(ap)
    a = cache[ap]; z = np.load(d / "pca.npz"); arcs = np.load(d / f"archetypes_k{k}.npy"); rws = z["rows"] if "rows" in z.files else None
    if "rows" in z.files: sc = z["scores"]
    else:                                   # pooled fit: scores for every cell, matched to the adata by barcode
        sc = z["scores_all"]; pos = pd.Series(np.arange(a.n_obs), index=a.obs_names)
        rws = pos[[str(b) for b in z["barcodes"]]].values
    pts = sc[:, : k - 1]
    dist = np.linalg.norm(pts[:, None, :] - arcs[None, :, : k - 1], axis=2)
    X = a.X[rws]; phase = a.obs["phase"].astype(str).values[rws] if "phase" in a.obs else None
    S = {}
    for short, full in USE.items():
        idx = [a.var_names.get_loc(g.capitalize()) for g in H[full] if g.capitalize() in a.var_names]
        M = X[:, idx].toarray(); M = (M - M.mean(0)) / np.where(M.std(0) > 0, M.std(0), 1); S[short] = M.mean(1)
    for v in range(k):
        close = dist[:, v] <= np.quantile(dist[:, v], 0.10)
        r = dict(fit=label, corner=f"V{v + 1}", n_cells=int(len(rws)))
        if phase is not None: r["cycling_pct"] = round(100 * float((phase[close] != "G1").mean()))
        for s in USE: r[s] = round(float(S[s][close].mean() - S[s][~close].mean()), 2)
        r["top_set"] = max(USE, key=lambda s_: r[s_]); rows.append(r)
t = pd.DataFrame(rows); t.to_csv(R / "docs/cleaning_audit/hallmark_scc_4t1.csv", index=False)
pd.set_option("display.width", 250); print(t.to_string(index=False))
