#!/usr/bin/env python3
"""Step 10 (Cook): draw the pooled k=4 simplex (KS genes), by cell line and by time.

The pooled fit is significant at every k (p=0.005) but each vertex is one cell line, so it is a
simplex of cell identity, not of EMT. Output: results/mail_figures/cook_pooled_tetrahedron.png
Usage (from Project_1/):  .venv/bin/python Mouse_SCC/Cook/codes/10_pooled_tetrahedron.py
"""
import sys, itertools
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cook_common as C
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13})
a = C.load_adata(); o = a.obs
d = C.RESULTS / "ks_pooled"; z = np.load(d / "pca.npz"); S = z["scores_all"]; arcs = np.load(d / "archetypes_k4.npy")
fwd = ~o.removal_arm.values; idx = np.random.default_rng(0).permutation(np.where(fwd)[0])
lc = {"A549": "#1F6F8B", "DU145": "#D9822B", "MCF7": "#4C9F70", "OVCA420": "#B5473A"}
tc = {"0d": "#2E5E8C", "8h": "#5C9BC4", "1d": "#9E9E9E", "3d": "#E3A04F", "7d": "#B5473A"}
w = C.barycentric(S[:, :3], arcs); pur, near = C.purity_and_nearest(w)
names = {v: pd.Series(o.cell_line.values[fwd & (pur >= .35) & (near == v)]).value_counts().index[0] for v in range(4)}
fig, axes = plt.subplots(1, 2, figsize=(16, 6.6))
for ax, mode in zip(axes, ["line", "time"]):
    cmap = lc if mode == "line" else tc
    lab = o.cell_line.values if mode == "line" else o.time.astype(str).values
    ax.scatter(S[idx, 0], S[idx, 1], s=3, c=[cmap[lab[i]] for i in idx], rasterized=True)
    for p, q in itertools.combinations(range(4), 2): ax.plot(arcs[[p, q], 0], arcs[[p, q], 1], "k-", lw=1.3, zorder=4)
    ax.scatter(arcs[:, 0], arcs[:, 1], s=70, c="k", zorder=6)
    cen = arcs[:, :2].mean(0)
    for v in range(4):
        dv = arcs[v, :2] - cen; dv = dv / np.linalg.norm(dv)
        ax.annotate(f"V{v + 1}  {names[v]}", arcs[v, :2], fontsize=13, fontweight="bold", xytext=tuple(14 * dv),
                    textcoords="offset points", ha="left" if dv[0] > 0.3 else ("right" if dv[0] < -0.3 else "center"),
                    va="bottom" if dv[1] > 0 else "top", bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="0.7"), zorder=7)
    for k_, c_ in cmap.items(): ax.scatter([], [], s=40, color=c_, label=k_)
    ax.legend(loc="lower left", fontsize=11, framealpha=.9)
    xl = ax.get_xlim(); yl = ax.get_ylim(); wd = xl[1] - xl[0]; ht = yl[1] - yl[0]
    ax.set_xlim(xl[0] - .2 * wd, xl[1] + .2 * wd); ax.set_ylim(yl[0] - .12 * ht, yl[1] + .14 * ht)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("Pooled fit, k = 4: coloured by cell line" if mode == "line" else "Same simplex: coloured by time after TGF-β", fontsize=15)
fig.text(0.5, 0.01, "All four lines, KS genes, forward arm (k = 4, p = 0.005). Each vertex is one cell line (100% of its specialists); time varies only a little inside each line's cloud.", ha="center", fontsize=12.5)
fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(C.RESULTS / "mail_figures" / "cook_pooled_tetrahedron.png", dpi=170)
print(names)
