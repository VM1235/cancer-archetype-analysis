#!/usr/bin/env python3
"""Step 8b: what are the corners of the within-library simplices?

For each fitted unit (08_within_unit.py) and k: place its cells in the simplex, call
specialists (purity gap >= 0.35), and describe each corner by
  - cells, EMT state mix, mouse/site, cell-cycle phase (% S or G2M), depth (median UMI),
    mito %, KS_EMT score
  - signature = mean log-expression in the 10% of the unit's cells closest to the corner
    minus the mean in the others (all genes)
Then correlate every corner signature with every other (different units too) to see whether
corners recur across libraries, and print each corner's top genes.
Also: how much of PC1-3 of each unit is explained by depth / phase / mito / KS_EMT.

Usage (from Project_1/):
  .venv/bin/python Mouse_SCC/Mouse_PDAC_GSE173958/codes/08b_describe_corners.py
"""
import json
import numpy as np, pandas as pd
import scc_common as C

import sys, glob
CLEAN = "--clean" in sys.argv
SFX = "_clean" if CLEAN else ""
if CLEAN:    # every clean all-genes fit with p < 0.05 in a single library (mice M1/M2 skipped)
    JOBS = []
    for f in sorted(glob.glob(str(C.RESULTS / "within_unit" / "all_*_clean" / "fit_k*.json"))):
        r = json.load(open(f))
        if r["p_value"] < 0.05 and r["unit"] not in ("M1", "M2"):
            JOBS.append(("all", r["k"], r["unit"]))
else:
  JOBS = [("all", 4, u) for u in ["M1-PT", "M1-SS", "M1-Liver", "M1-Met", "M2-PT", "M2-Liver", "M2-Met"]] + \
       [("all", 3, u) for u in ["M1-PT", "M1-SS", "M2-Liver"]] + [("ks", 4, "M1-PT")]
OUT = C.RESULTS / "within_unit"

a = C.load_adata()
obs = a.obs.copy()
X = a.X
genes = a.var_names
rows_out, sigs, names, info = [], [], [], []
for gs, k, unit in JOBS:
    d = OUT / f"{gs}_{unit}{SFX}"
    z = np.load(d / "pca.npz"); arcs = np.load(d / f"archetypes_k{k}.npy")
    rows, sc = z["rows"], z["scores"]
    o = obs.iloc[rows]
    # diagnostics: what do the first PCs of this unit track?
    for j in range(3):
        info.append(dict(unit=unit, geneset=gs, PC=j + 1, pct_var=100 * z["evr"][j],
                         r2_log_depth=np.corrcoef(sc[:, j], np.log10(o["nCount_RNA"]))[0, 1] ** 2,
                         r2_mito=np.corrcoef(sc[:, j], o["pct_counts_mt"])[0, 1] ** 2,
                         r2_KS_EMT=np.corrcoef(sc[:, j], o["KS_EMT"])[0, 1] ** 2,
                         eta_phase=C.eta_squared(sc[:, j], o["phase"]),
                         eta_state=C.eta_squared(sc[:, j], o["cluster_annotation"])))
    pts = sc[:, : k - 1]
    w = C.barycentric(pts, arcs)
    pur, near = C.purity_and_nearest(w)
    dist = np.linalg.norm(pts[:, None, :] - arcs[None, :, : k - 1], axis=2)
    for v in range(k):
        spec = (pur >= C.PURITY_CUT) & (near == v)
        close = dist[:, v] <= np.quantile(dist[:, v], 0.10)
        sig = np.asarray(X[rows[close]].mean(0)).ravel() - np.asarray(X[rows[~close]].mean(0)).ravel()
        ph = o["phase"].astype(str).values
        st = o["cluster_annotation"].astype(str).values
        cl = f"{unit}|{gs}|k{k}|V{v + 1}"
        rows_out.append(dict(corner=cl, n_close10=int(close.sum()), n_specialists=int(spec.sum()),
                             pct_cycling=round(100 * np.isin(ph[close], ["S", "G2M"]).mean(), 1),
                             mean_UMI=int(o["nCount_RNA"].values[close].mean()),
                             mean_mito=round(float(o["pct_counts_mt"].values[close].mean()), 1),
                             KS_EMT=round(float(o["KS_EMT"].values[close].mean()), 2),
                             **{s: round(100 * (st[close] == s).mean(), 0) for s in C.STATE_ORDER},
                             top_genes=", ".join(genes[np.argsort(sig)[::-1][:10]])))
        sigs.append(sig); names.append(cl)
tab = pd.DataFrame(rows_out)
tab.to_csv(OUT / f"corner_description{SFX}.csv", index=False)
R = pd.DataFrame(np.corrcoef(np.vstack(sigs)), index=names, columns=names).round(2)
R.to_csv(OUT / f"corner_signature_correlation{SFX}.csv")
pd.DataFrame(info).round(3).to_csv(OUT / f"unit_pc_drivers{SFX}.csv", index=False)
pd.set_option("display.width", 300); pd.set_option("display.max_colwidth", 90)
print(pd.DataFrame(info).round(2).to_string(index=False))
print(tab.drop(columns=["top_genes"]).to_string(index=False))
print(tab[["corner", "top_genes"]].to_string(index=False))
k4 = [n for n in names if "|all|k4|" in n]
print("\nCorrelation of corner signatures across libraries (all genes, k=4):")
print(R.loc[k4, k4].to_string())
