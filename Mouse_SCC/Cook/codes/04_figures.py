#!/usr/bin/env python3
"""Step 4 (Cook): per-fit PC figures coloured by time, Mix, phase, with the k=3 and elbow-k simplex drawn.

Top row: PC1-PC2 coloured by time (forward arm), Mix, cell-cycle phase.
Bottom row: PC1-PC3, PC2-PC3 by time; removal arms (_rm) as guests over the grey forward cells;
            mean position of each time point along PC1/PC2.
Usage (from Project_1/): .venv/bin/python Mouse_SCC/Cook/codes/04_figures.py --geneset all --tag A549
"""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
import cook_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--geneset", default="all"); ap.add_argument("--tag", required=True)
args = ap.parse_args()
run = C.run_dir(args.geneset, args.tag)
z = np.load(run / "pca.npz")
obs = C.load_adata().obs
assert list(obs.index) == list(z["barcodes"])
S, fi = z["scores_all"], z["fit_idx"]
esv = pd.read_csv(run / "esv_curve.csv"); ke = int(pd.read_csv(run / "fit_summary.csv").iloc[0]["k"]) if False else None
import json
summ = pd.read_csv(run / "fit_summary.csv")
fwd = np.zeros(len(obs), bool); fwd[fi] = True
rm = (obs["time"].str.endswith("_rm")).values
if args.tag != "pooled":
    rm = rm & (obs["cell_line"] == args.tag).values      # guests: removal arms of this line only
tcol = dict(zip(["0d", "8h", "1d", "3d", "7d"], plt.get_cmap("viridis")(np.linspace(0, .95, 5))))
rcol = dict(zip(["8h_rm", "1d_rm", "3d_rm"], ["#fdae61", "#f46d43", "#a50026"]))

def sc(ax, i, j, mask, colours, label_map, s=3):
    for lab in label_map:
        m = mask & (colours == lab)
        if m.any():
            ax.scatter(S[m, i], S[m, j], s=s, alpha=.55, c=[label_map[lab]], label=f"{lab} ({m.sum()})", rasterized=True)
    ax.set_xlabel(f"PC{i+1}"); ax.set_ylabel(f"PC{j+1}"); ax.legend(markerscale=4, fontsize=6)

t = obs["time"].astype(str).values
fig, ax = plt.subplots(2, 3, figsize=(17, 10))
sc(ax[0, 0], 0, 1, fwd, t, tcol); ax[0, 0].set_title("PC1-PC2 by time (fit cells)")
for k in (3, int(summ.loc[(summ.k > 0), "k"].iloc[0]) if False else 3):
    pass
try:
    a3 = np.load(run / "archetypes_k3.npy")
    for i in range(3):
        jx = (i + 1) % 3
        ax[0, 0].plot(a3[[i, jx], 0], a3[[i, jx], 1], "k-", lw=1)
    ax[0, 0].scatter(a3[:, 0], a3[:, 1], c="k", s=40, zorder=5)
except FileNotFoundError:
    pass
mixes = obs["mix"].astype(str).values
sc(ax[0, 1], 0, 1, fwd, mixes, dict(zip(sorted(set(mixes)), plt.get_cmap("tab10")(range(10))))); ax[0, 1].set_title("by Mix (pooling batch)")
ph = obs["phase"].astype(str).values
sc(ax[0, 2], 0, 1, fwd, ph, {"G1": "#1f77b4", "S": "#2ca02c", "G2M": "#d62728"}); ax[0, 2].set_title("by cell-cycle phase")
sc(ax[1, 0], 0, 2, fwd, t, tcol); ax[1, 0].set_title("PC1-PC3 by time")
sc(ax[1, 1], 1, 2, fwd, t, tcol); ax[1, 1].set_title("PC2-PC3 by time")
a = ax[1, 2]
a.scatter(S[fwd, 0], S[fwd, 1], s=2, c="#dddddd", rasterized=True)
sc(a, 0, 1, rm, t, rcol, s=6); a.set_title("removal arms (guests) over fit cells")
for x in ax.ravel(): x.grid(alpha=.2)
fig.suptitle(f"{args.geneset}/{args.tag}", fontsize=11)
fig.tight_layout(); fig.savefig(run / "pc_by_time_mix_phase.png", dpi=150); plt.close(fig)

# mean score per time point (PC1..PC5) incl. removal arms: where does each time sit?
rows = []
for tm in C.TIME_ORDER:
    m = (t == tm) & (obs["cell_line"].isin(C.LINES)).values
    if args.tag != "pooled":
        m &= (obs["cell_line"] == args.tag).values
    if m.sum():
        rows.append(dict(time=tm, n=int(m.sum()), **{f"PC{j+1}": float(S[m, j].mean()) for j in range(5)}))
pd.DataFrame(rows).round(2).to_csv(run / "time_means.csv", index=False)
print(pd.DataFrame(rows).round(2).to_string(index=False))
