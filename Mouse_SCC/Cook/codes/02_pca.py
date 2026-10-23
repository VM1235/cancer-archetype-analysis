#!/usr/bin/env python3
"""Step 2 (Cook): PCA on the fit cells, then diagnostics. Mirrors Mouse_SCC/codes/02_pca.py.

Fit cells: forward TGF-B1 time course only (0d, 8h, 1d, 3d, 7d). The three
removal arms (_rm) are projected afterwards as guests.
  --lines all        pooled, at most 550 cells per line (balanced, like 800 per tumour in mouse)
  --lines A549       one line, all its forward-arm cells
Usage (from Project_1/):
  .venv/bin/python -u Mouse_SCC/Cook/codes/02_pca.py --geneset all --lines all --tag pooled
"""
import argparse, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
import cook_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--geneset", choices=["all", "ks", "hallmark"], default="all")
ap.add_argument("--lines", default="all")
ap.add_argument("--tag", required=True)
args = ap.parse_args()
out = C.run_dir(args.geneset, args.tag)
t0 = time.time()
a = C.load_adata()
gi = C.gene_index(a, args.geneset)
genes = a.var_names[gi].tolist()
obs = a.obs.copy()
X = a.X[:, gi].tocsr()
lines = None if args.lines == "all" else args.lines.split(",")
fit_idx = C.fit_cells(obs, lines=lines, max_per_line=C.MAX_PER_LINE if lines is None else None)
print(f"geneset {args.geneset}: {len(genes)} genes; fit cells {len(fit_idx)} of {X.shape[0]}", flush=True)
Xfit = X[fit_idx].toarray().astype(np.float64)
pca = C.fit_pca(Xfit)
scores_all = C.project(X, pca)
assert np.abs(scores_all[fit_idx] - pca["scores_fit"]).max() < 1e-6 * max(1, np.abs(pca["scores_fit"]).max())
evr = pca["evr"]; cum = np.cumsum(evr)
below = np.where(evr < 0.001)[0]
n_rule = int(below[0]) if len(below) else len(evr)
print("PC1..10 variance %:", np.round(100 * evr[:10], 2).tolist())
print(f"0.1% rule keeps {n_rule} PCs ({100 * cum[n_rule - 1]:.1f}% variance)")

o = obs.iloc[fit_idx]
logd = np.log10(o["nCount_RNA"].values)
rows = []
for j in range(10):
    s = pca["scores_fit"][:, j]
    rows.append(dict(PC=j + 1, pct_variance=100 * evr[j], cell_line=C.eta_squared(s, o["cell_line"]),
                     time=C.eta_squared(s, o["time"]), mix=C.eta_squared(s, o["mix"]),
                     phase=C.eta_squared(s, o["phase"]),
                     log_depth=float(np.corrcoef(s, logd)[0, 1] ** 2),
                     pct_mito=float(np.corrcoef(s, o["pct_mito"].values)[0, 1] ** 2),
                     pseudotime=float(np.corrcoef(s, o["pseudotime"].values)[0, 1] ** 2)))
diag = pd.DataFrame(rows); diag.to_csv(out / "pc_drivers.csv", index=False)
print("\nShare of each PC explained by each factor (0-1; fit cells):")
print(diag.round(3).to_string(index=False))
top = {}
for j in range(5):
    v = pca["components"][j]
    top[f"PC{j+1}+"] = [genes[i] for i in np.argsort(v)[::-1][:15]]
    top[f"PC{j+1}-"] = [genes[i] for i in np.argsort(v)[:15]]
pd.DataFrame(top).to_csv(out / "pc_top_genes.csv", index=False)
np.savez_compressed(out / "pca.npz", scores_all=scores_all, fit_idx=fit_idx, mean=pca["mean"],
                    components=pca["components"], evr=evr[:200], genes=np.array(genes),
                    barcodes=np.array(obs.index, dtype=str))
C.save_json(dict(geneset=args.geneset, tag=args.tag, lines=args.lines, n_genes=len(genes),
                 n_fit_cells=int(len(fit_idx)), n_pcs_rule=n_rule, cumvar_at_rule=float(cum[n_rule - 1])),
            out / "pca_info.json")

fig, ax = plt.subplots(figsize=(5, 3.4))
ax.plot(np.arange(1, 51), 100 * evr[:50], "o-", ms=3); ax.axhline(0.1, ls="--", c="0.5", lw=.8)
ax.set_yscale("log"); ax.set_xlabel("PC"); ax.set_ylabel("% variance")
ax.set_title(f"Scree ({args.geneset}, {args.tag}); 0.1% rule -> {n_rule} PCs", fontsize=9)
fig.savefig(out / "scree.png", dpi=200, bbox_inches="tight"); plt.close(fig)

s = pca["scores_fit"]
fig, axes = plt.subplots(1, 4, figsize=(19, 4.3))
for ax, col in zip(axes, ["cell_line", "time", "mix", "phase"]):
    lab = o[col].astype(str).values
    levels = [t for t in C.TIME_ORDER if t in set(lab)] if col == "time" else sorted(set(lab))
    for x in levels:
        m = lab == x
        ax.scatter(s[m, 0], s[m, 1], s=3, alpha=.5, label=x, rasterized=True)
    ax.set_title(col); ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.legend(markerscale=4, fontsize=6)
fig.savefig(out / "pc1_pc2_by_label.png", dpi=170, bbox_inches="tight"); plt.close(fig)
print(f"wrote {out} ({time.time()-t0:.0f}s)")
