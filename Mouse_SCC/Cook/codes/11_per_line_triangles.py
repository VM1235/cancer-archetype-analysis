#!/usr/bin/env python3
"""Step 11 (Cook): the four per-line k=3 simplices side by side, cells coloured by time.

Cleaned fits (07_within_line_clean.py, variant 'clean', all genes). The triangle is drawn whether or not
it is significant; the shuffle-test p-values (all genes, KS genes) are in each title.
Output: results/mail_figures/cook_per_line_triangles.png
Usage (from Project_1/):  .venv/bin/python Mouse_SCC/Cook/codes/11_per_line_triangles.py
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cook_common as C
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13})
a = C.load_adata(); o = a.obs
tc = {"0d": "#2E5E8C", "8h": "#5C9BC4", "1d": "#9E9E9E", "3d": "#E3A04F", "7d": "#B5473A"}
fig, axes = plt.subplots(2, 2, figsize=(13.5, 11), constrained_layout=True)
for ax, L in zip(axes.ravel(), C.LINES):
    d = C.RESULTS / "within_line" / f"all_{L}_clean"
    z = np.load(d / "pca.npz"); arcs = np.load(d / "archetypes_k3.npy"); S = z["scores"]
    t = o["time"].astype(str).values[z["rows"]]
    pa = json.load(open(d / "fit_k3.json"))["p_value"]
    pk = json.load(open(C.RESULTS / "within_line" / f"ks_{L}_clean" / "fit_k3.json"))["p_value"]
    idx = np.random.default_rng(0).permutation(len(t))
    ax.scatter(S[idx, 0], S[idx, 1], s=5, c=[tc[x] for x in t[idx]], rasterized=True)
    T = np.vstack([arcs[:, :2], arcs[:1, :2]]); ax.plot(T[:, 0], T[:, 1], "k-", lw=1.8, zorder=4)
    ax.scatter(arcs[:, 0], arcs[:, 1], s=60, c="k", zorder=6)
    for v in range(3):
        ax.annotate(f"V{v + 1}", arcs[v, :2], fontsize=14, fontweight="bold", xytext=(8, 8), textcoords="offset points",
                    bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="0.7"), zorder=7)
    sig = "significant" if min(pa, pk) < 0.05 else "not significant"
    ax.set_title(f"{L} (n={len(t)})\nk=3: p={pa:.2f} all genes, p={pk:.2f} KS genes", fontsize=13)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.spines[["top", "right"]].set_visible(False)
for k_, c_ in tc.items(): axes[0, 0].scatter([], [], s=40, color=c_, label=k_)
axes[0, 0].legend(loc="upper left", fontsize=11, framealpha=.9, title="time after TGF-β")
fig.suptitle("Cook: each cell line's own k = 3 simplex (cleaned, all genes), coloured by time", fontsize=15, fontweight="bold")
fig.savefig(C.RESULTS / "mail_figures" / "cook_per_line_triangles.png", dpi=170)
print("ok")
