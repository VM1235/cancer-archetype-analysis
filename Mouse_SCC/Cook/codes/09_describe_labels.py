#!/usr/bin/env python3
"""Step 9 (Cook): where do the time points sit in a within-line simplex?

For a cleaned per-line fit (07_within_line_clean.py), every cell of that line, including the removal
arm (which was never fitted), is placed in the simplex and called a specialist of its nearest corner or
a generalist by two rules:
  gap     largest barycentric weight minus the second largest >= 0.35
  groves  the 5% of FIT cells closest to each corner (cut-off distances set on the fit cells,
          then applied to the removal-arm guests)
Tables: per time point, % at each corner and % generalist. Figure: fit cells coloured by time,
by specialist/generalist, plus the removal-arm guests.
Time is ground truth here, so this is a clean check of whether late states are corners or interior.

Usage (from Project_1/):
  .venv/bin/python Mouse_SCC/Cook/codes/09_describe_labels.py --line A549 --geneset all --k 3
"""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
import cook_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--line", required=True)
ap.add_argument("--geneset", default="all")
ap.add_argument("--variant", default="clean")
ap.add_argument("--k", type=int, default=3)
args = ap.parse_args()
k = args.k
d = C.RESULTS / "within_line" / f"{args.geneset}_{args.line}_{args.variant}"
out = C.RESULTS / "within_line" / f"describe_{args.geneset}_{args.line}_{args.variant}_k{k}"
out.mkdir(exist_ok=True)
z = np.load(d / "pca.npz"); arcs = np.load(d / f"archetypes_k{k}.npy")
a = C.load_adata(); o = a.obs
gi = [a.var_names.get_loc(g) for g in z["genes"]]
V3 = z["components"][: k - 1]
line_rows = np.where(o["cell_line"].values == args.line)[0]
pts = np.asarray(a.X[line_rows][:, gi] @ V3.T) - z["mean"] @ V3.T
fit_rows = z["rows"]
pos = {r: i for i, r in enumerate(line_rows)}
is_fit = np.isin(line_rows, fit_rows)
assert np.abs(pts[[pos[r] for r in fit_rows]] - z["scores"][:, : k - 1]).max() < 1e-4 * np.abs(z["scores"]).max()
V = [f"V{i + 1}" for i in range(k)]
w = C.barycentric(pts, arcs)
pur, near = C.purity_and_nearest(w)
gap = np.where(pur >= C.PURITY_CUT, np.array(V)[near], "generalist")
dist = np.linalg.norm(pts[:, None, :] - arcs[None, :, : k - 1], axis=2)
cut = np.quantile(dist[is_fit], 0.05, axis=0)
close = dist <= cut
g5 = np.where(close.sum(1) == 0, "generalist", np.array(V)[np.argmin(dist / cut, axis=1)])
t = o["time"].astype(str).values[line_rows]
cells = pd.DataFrame({"time": t, "is_fit": is_fit, "gap": gap, "groves5": g5, "PC1": pts[:, 0],
                      "PC2": pts[:, 1], "pseudotime": o["pseudotime"].values[line_rows],
                      "phase": o["phase"].astype(str).values[line_rows]}, index=o.index[line_rows])
cells.to_csv(out / "cells.csv.gz")
cols = V + ["generalist"]
lines = [f"# {args.line} ({args.geneset}, {args.variant}), k={k}: {int(is_fit.sum())} fit cells, "
         f"{int((~is_fit).sum())} held-out cells (removal arm and cells dropped by cleaning)", ""]
for rule in ("gap", "groves5"):
    tab = pd.crosstab(cells["time"], cells[rule]).reindex(columns=cols, fill_value=0)
    tab = tab.reindex([x for x in C.TIME_ORDER if x in tab.index])
    pct = (100 * tab.div(tab.sum(1), axis=0)).round(1); pct["n"] = tab.sum(1)
    lines += [f"## {rule} rule: % of each time point at each corner / generalist", "", pct.to_string(), ""]
    pct.to_csv(out / f"{rule}_time.csv")
    ph = pd.crosstab(cells["phase"], cells[rule]).reindex(columns=cols, fill_value=0)
    lines += [f"## {rule} rule: cell-cycle phase (%)", "", (100 * ph.div(ph.sum(1), axis=0)).round(1).to_string(), ""]
(out / "summary.md").write_text("\n".join(lines))
print("\n".join(lines))

rng = np.random.default_rng(0)
tcol = {"0d": "#2E5E8C", "8h": "#5C9BC4", "1d": "#B9B9B9", "3d": "#E3A04F", "7d": "#B5473A",
        "8h_rm": "#4C9F70", "1d_rm": "#2F7A52", "3d_rm": "#1B4D33"}
vcol = {"V1": "#B5473A", "V2": "#1F6F8B", "V3": "#D9822B", "V4": "#4C9F70", "generalist": "#D5D8DC"}
tri = np.vstack([arcs[:, :2], arcs[:1, :2]])
fig, axes = plt.subplots(1, 3, figsize=(17, 5.4), constrained_layout=True)
F = np.where(cells["is_fit"].values & ~cells["time"].str.endswith("_rm").values)[0]; F = rng.permutation(F)
for ax, lab in zip(axes[:2], ["time", "gap"]):
    cm = tcol if lab == "time" else vcol
    ax.scatter(cells["PC1"].values[F], cells["PC2"].values[F], s=5, c=[cm[v] for v in cells[lab].values[F]], alpha=.7, rasterized=True)
    for l, c in cm.items():
        if lab == "time" and l.endswith("_rm"): continue
        if l.startswith("V") and int(l[1:]) > k: continue
        ax.scatter([], [], s=25, color=c, label=l)
    ax.set_title("fit cells by time" if lab == "time" else "specialist (gap rule) or generalist")
ax = axes[2]
G = np.where(cells["time"].str.endswith("_rm").values)[0]
ax.scatter(cells["PC1"].values[F], cells["PC2"].values[F], s=4, c="#E3E5E8", rasterized=True)
for tt in ["8h_rm", "1d_rm", "3d_rm"]:
    m = G[cells["time"].values[G] == tt]
    ax.scatter(cells["PC1"].values[m], cells["PC2"].values[m], s=9, c=tcol[tt], alpha=.8, label=f"{tt} (n={len(m)})", rasterized=True)
ax.set_title("held-out removal arm over the fit cells")
for ax in axes:
    ax.plot(tri[:, 0], tri[:, 1], "k-", lw=1.4)
    for i in range(k): ax.annotate(V[i], arcs[i, :2], fontsize=13, fontweight="bold", xytext=(6, 6), textcoords="offset points")
    ax.legend(fontsize=9, frameon=False, markerscale=1.2); ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
fig.suptitle(f"{args.line}, {args.geneset} genes, {args.variant}, k={k}", fontsize=14)
fig.savefig(out / "triangle_by_time.png", dpi=170); plt.close(fig)
