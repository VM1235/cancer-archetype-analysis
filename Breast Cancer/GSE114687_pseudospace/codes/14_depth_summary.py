#!/usr/bin/env python3
"""Step 14: summarise the depth control for the cleaned k=3 within-sample fits.

For each fit, three versions: `clean` (>= 10,000 UMI, as in the main results), `rd` (the same cells with
log10(UMI) regressed out of every gene) and `dm` (every cell thinned to the same depth, then renormalised).
Per version: shuffle-test p, and the corner with the highest cycling share (% of its 10%-nearest cells in S or G2M)
and its mean UMI. Output: results/within_unit/depth_control_summary.csv and results/depth_for_deck.json.
Usage (from Project_1/): .venv/bin/python "Breast Cancer/GSE114687_pseudospace/codes/14_depth_summary.py"
"""
import json
import numpy as np, pandas as pd, anndata as ad
import scc_common as C

FITS = [("all", "MCF10A_Mock"), ("all", "HuMEC_Mock"), ("all", "MCF10A_TGFB"), ("all", "HuMEC_TGFB"),
        ("ks", "MCF10A_Mock"), ("ks", "HuMEC_Mock"), ("ks", "HuMEC_TGFB"), ("all", "MCF10A")]
VAR = {"clean": ("_clean", "expr_tumour.h5ad"), "rd": ("_clean_rd", "expr_tumour.h5ad"), "dm": ("_dm", "expr_tumour_depthmatched.h5ad")}
obs = {f: ad.read_h5ad(C.PROC / f, backed="r").obs for f in {"expr_tumour.h5ad", "expr_tumour_depthmatched.h5ad"}}
rows = []
for gs, u in FITS:
    r = dict(fit=f"{gs} {u}")
    for v, (sfx, f) in VAR.items():
        d = C.RESULTS / "within_unit" / f"{gs}_{u}{sfx}"
        if not (d / "fit_k3.json").exists():
            r[f"{v}_p"] = np.nan; continue
        z = np.load(d / "pca.npz"); arcs = np.load(d / "archetypes_k3.npy")
        o = obs[f].iloc[z["rows"]]
        pts = z["scores"][:, :2]; dist = np.linalg.norm(pts[:, None, :] - arcs[None, :, :2], axis=2)
        ph = o["phase"].astype(str).values; umi = o["nCount_RNA"].values
        best = (-1, None, None)
        for k in range(3):
            close = dist[:, k] <= np.quantile(dist[:, k], 0.10)
            cyc = 100 * (ph[close] != "G1").mean()
            if cyc > best[0]: best = (cyc, int(umi[close].mean()), k + 1)
        j = json.load(open(d / "fit_k3.json"))
        r.update({f"{v}_n": j["n_cells"], f"{v}_p": j["p_value"], f"{v}_t": round(j["t_ratio"], 3),
                  f"{v}_cyc_corner_pct": round(best[0]), f"{v}_cyc_corner_umi": best[1]})
    rows.append(r)
t = pd.DataFrame(rows)
t.to_csv(C.RESULTS / "within_unit" / "depth_control_summary.csv", index=False)
pd.set_option("display.width", 250); print(t.to_string(index=False))
fmt = lambda r, v: ("n/a" if pd.isna(r.get(f"{v}_p")) else f"{r[f'{v}_p']:.3f}")
out = dict(rows=[[r["fit"], fmt(r, "clean"), fmt(r, "rd"), fmt(r, "dm")] for r in rows])
json.dump(out, open(C.RESULTS / "depth_for_deck.json", "w"), indent=1)
