#!/usr/bin/env python3
"""Step 12 (Cook): do the per-line k=3 corners have distinct gene programmes (the PDAC check)?

For each cleaned per-line fit (all genes, k=3), the 10% of cells nearest each corner are scored for
gene sets (mean z within the line): KS epithelial, KS mesenchymal, proliferation, hypoxia,
TGF-B response (PMEPA1 SERPINE1 INHBA TGFBI FN1 SERPINE2 IGFBP7 TAGLN). Plus cycling %, time mix.
This is run for every line whether or not its simplex is significant, so non-significant corners
can be read as "functional but unsupported" rather than assumed empty.
Output: results/within_line/corner_programmes.csv and results/mail_figures/cook_corner_programmes.png
Usage (from Project_1/):  .venv/bin/python Mouse_SCC/Cook/codes/12_corner_programmes.py
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cook_common as C
a = C.load_adata(); o = a.obs
epi, mes = C.ks_genes(a.var_names)
sets = {"KS epithelial": epi, "KS mesenchymal": mes,
        "proliferation": "HMGB2 BIRC5 CENPA CKS1B TOP2A UBE2C MKI67".split(),
        "hypoxia": "VEGFA BNIP3 NDRG1 ERO1L HILPDA SLC2A1".split(),
        "TGF-β response": "PMEPA1 SERPINE1 INHBA TGFBI FN1 SERPINE2 IGFBP7 TAGLN".split()}
rows, labels = [], []
for L in C.LINES:
    d = C.RESULTS / "within_line" / f"all_{L}_clean"; z = np.load(d / "pca.npz"); arcs = np.load(d / "archetypes_k3.npy")
    r = z["rows"]; P = z["scores"][:, :2]; oo = o.iloc[r]
    pv = json.load(open(d / "fit_k3.json"))["p_value"]
    Z = {}
    for n, gs in sets.items():
        idx = [a.var_names.get_loc(g) for g in gs if g in a.var_names]
        M = a.X[r][:, idx].toarray(); M = (M - M.mean(0)) / np.where(M.std(0) > 0, M.std(0), 1); Z[n] = M.mean(1)
    for v in range(3):
        dd = np.linalg.norm(P - arcs[v, :2], axis=1); c = dd <= np.quantile(dd, .10)
        tm = oo["time"].astype(str).iloc[np.where(c)[0]].value_counts(normalize=True)
        rows.append({"line": L, "corner": f"V{v + 1}", "p_k3": round(pv, 3), **{n: round(float(Z[n][c].mean()), 2) for n in sets},
                     "cycling": round(float(np.isin(oo.phase.astype(str).values[c], ["S", "G2M"]).mean()), 2),
                     "main_time": " ".join(f"{k}:{100 * x:.0f}%" for k, x in tm.head(2).items())})
        labels.append(f"{L} V{v + 1}")
df = pd.DataFrame(rows); df.to_csv(C.RESULTS / "within_line" / "corner_programmes.csv", index=False)
print(df.to_string(index=False))
H = df[list(sets)].values
fig, ax = plt.subplots(figsize=(9.5, 8.2)); im = ax.imshow(H, cmap="RdBu_r", vmin=-1.6, vmax=1.6, aspect="auto")
ax.set_xticks(range(len(sets))); ax.set_xticklabels([n.replace(" ", "\n") for n in sets], fontsize=10); ax.xaxis.tick_top()
ax.set_yticks(range(len(labels))); ax.set_yticklabels([f"{l}  ({df.main_time[i]})" for i, l in enumerate(labels)], fontsize=9)
for i in range(H.shape[0]):
    for j in range(H.shape[1]): ax.text(j, i, f"{H[i, j]:+.1f}", ha="center", va="center", fontsize=9, color="white" if abs(H[i, j]) > 0.9 else "#1B2A3A")
for b in (2.5, 5.5, 8.5): ax.axhline(b, color="k", lw=1.5)
ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
plt.colorbar(im, ax=ax, shrink=.5, pad=.03).set_label("gene-set score (z, within line)")
fig.suptitle("Cook: gene programme at each per-line corner (cleaned, k=3; fits not all significant)", y=0.995, fontsize=12, fontweight="bold")
fig.savefig(C.RESULTS / "mail_figures" / "cook_corner_programmes.png", dpi=170, bbox_inches="tight")
