#!/usr/bin/env python3
"""Step 10 (PDAC): where do the EMT states sit in a within-library simplex?

For a cleaned within-library fit (08_within_unit.py --clean), place the unit's cells in the simplex and
call each cell a specialist of its nearest corner or a generalist (interior) by two rules:
  gap     largest barycentric weight minus the second largest >= 0.35 (Pal work)
  groves  the 5% of cells closest to each corner (Groves 2022); every other cell is a generalist
Tables: per EMT state and per clusterU, % at each corner and % generalist, by each rule.
Figure: the triangle with cells coloured by (1) EMT state, (2) specialist corner or generalist,
(3) clusterU. States are KS-quantile bins, so their share at a corner is partly built in.

Usage (from Project_1/):
  .venv/bin/python "Breast Cancer/GSE114687_pseudospace/codes"/10_describe_labels.py --unit M1-SS --k 3
"""
import argparse
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
import scc_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--unit", required=True)
ap.add_argument("--k", type=int, default=3)
ap.add_argument("--geneset", default="all")
args = ap.parse_args()
k = args.k
d = C.RESULTS / "within_unit" / f"{args.geneset}_{args.unit}_clean"
out = C.RESULTS / "within_unit" / f"describe_{args.geneset}_{args.unit}_clean_k{k}"
out.mkdir(exist_ok=True)
z = np.load(d / "pca.npz"); arcs = np.load(d / f"archetypes_k{k}.npy")
a = C.load_adata(); o = a.obs.iloc[z["rows"]]
pts = z["scores"][:, : k - 1]
V = [f"V{i + 1}" for i in range(k)]
w = C.barycentric(pts, arcs)
pur, near = C.purity_and_nearest(w)
gap = np.where(pur >= C.PURITY_CUT, np.array(V)[near], "generalist")
dist = np.linalg.norm(pts[:, None, :] - arcs[None, :, : k - 1], axis=2)
cut = np.quantile(dist, 0.05, axis=0)
close = dist <= cut
g5 = np.where(close.sum(1) == 0, "generalist", np.array(V)[np.argmin(dist / cut, axis=1)])
st = o["cluster_annotation"].astype(str).values; cu = o["clusterU"].astype(str).values
cells = pd.DataFrame({"state": st, "clusterU": cu, "position": o["position"].astype(str).values, "arm": o["arm"].astype(str).values, "gap": gap, "groves5": g5,
                      "PC1": pts[:, 0], "PC2": pts[:, 1] if k > 2 else 0}, index=o.index)
cells.to_csv(out / "cells.csv.gz")
lines = [f"# {args.unit} cleaned, genes {args.geneset}, k={k}: {len(cells)} cells", ""]
cols = V + ["generalist"]
for rule in ("gap", "groves5"):
    for lab, order in (("state", C.STATE_ORDER), ("position", ["inner", "outer"]), ("arm", ["Mock", "TGFB"]), ("clusterU", C.CLUSTERU_ORDER)):
        t = pd.crosstab(cells[lab], cells[rule]).reindex(columns=cols, fill_value=0)
        t = t.reindex([x for x in order if x in t.index])
        pct = (100 * t.div(t.sum(1), axis=0)).round(1); pct["n"] = t.sum(1)
        lines += [f"## {rule} rule: % of each {lab} at each corner / generalist", "", pct.to_string(), ""]
        pct.to_csv(out / f"{rule}_{lab}.csv")
    sp = (cells[rule] != "generalist").mean() * 100
    lines += [f"Specialists by {rule}: {sp:.1f}% of cells", ""]
(out / "summary.md").write_text("\n".join(lines))
print("\n".join(lines))

rng = np.random.default_rng(0); sh = rng.permutation(len(cells))
scol = dict(zip(C.STATE_ORDER, ["#2E5E8C", "#5C9BC4", "#E3A04F", "#B5473A"]))
vcol = {"V1": "#B5473A", "V2": "#1F6F8B", "V3": "#D9822B", "V4": "#4C9F70", "generalist": "#D5D8DC"}
fig, axes = plt.subplots(1, 3, figsize=(17, 5.4), constrained_layout=True)
tri = np.vstack([arcs[:, :2], arcs[:1, :2]])
for ax, lab, cm in zip(axes, ["state", "gap", "clusterU"], [scol, vcol, None]):
    vals = cells[lab].values[sh]
    if cm is None:
        levels = sorted(set(vals)); cmap = plt.get_cmap("tab10"); cm = {l: cmap(i) for i, l in enumerate(levels)}
    ax.scatter(cells["PC1"].values[sh], cells["PC2"].values[sh], s=4, c=[cm[v] for v in vals], alpha=.7, rasterized=True)
    ax.plot(tri[:, 0], tri[:, 1], "k-", lw=1.4)
    for i in range(k): ax.annotate(V[i], arcs[i, :2], fontsize=13, fontweight="bold", xytext=(6, 6), textcoords="offset points")
    for l, c in cm.items():
        if l.startswith("V") and l[1:].isdigit() and int(l[1:]) > k: continue
        ax.scatter([], [], s=25, color=c, label=l)
    ax.legend(fontsize=9, frameon=False, markerscale=1.2)
    ax.set_title({"state": "EMT state (KS quantiles)", "gap": "specialist (gap rule) or generalist", "clusterU": "clusterU (KS-gene k-means)"}[lab])
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
fig.suptitle(f"{args.unit}, cleaned, k={k}", fontsize=14)
fig.savefig(out / "triangle_by_label.png", dpi=170); plt.close(fig)
