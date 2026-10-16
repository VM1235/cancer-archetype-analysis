#!/usr/bin/env python3
"""Step 6 (Cook): specialists vs generalists, Groves 2022 closest-5% rule. Mirrors the mouse script.

Reference cells = the forward-arm fit cells; the 5% cut-off distance per vertex is set on them.
Removal-arm cells (guests) are classified with the same cut-offs. Reported by time point,
cell-cycle phase and Mix (and cell line for the pooled fit).
Usage (from Project_1/): .venv/bin/python Mouse_SCC/Cook/codes/06_groves_specialists.py --geneset all --tag A549 --k 4
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
ap.add_argument("--k", type=int, default=3)
args = ap.parse_args(); k = args.k
run = C.run_dir(args.geneset, args.tag); out = run / f"groves5_k{k}"; out.mkdir(exist_ok=True)
z = np.load(run / "pca.npz"); arcs = np.load(run / f"archetypes_k{k}.npy")
obs = C.load_adata().obs
assert list(obs.index) == list(z["barcodes"])
pts = z["scores_all"][:, : k - 1]; V = [f"V{i+1}" for i in range(k)]
ref = np.zeros(len(obs), bool); ref[z["fit_idx"]] = True
guest = obs["time"].str.endswith("_rm").values
if args.tag != "pooled":
    guest &= (obs["cell_line"] == args.tag).values
d = np.linalg.norm(pts[:, None, :] - arcs[None, :, : k - 1], axis=2)
cut = np.quantile(d[ref], 0.05, axis=0)
cls = np.where((d <= cut).sum(1) == 0, "generalist", np.array(V)[np.argmin(d / cut, axis=1)])
w = C.composition(C.barycentric(pts, arcs)); pur, near = C.purity_and_nearest(w)
cells = pd.DataFrame({"groves5": cls, "weight095": np.where(w.max(1) > .95, np.array(V)[near], "generalist"),
                      "gap035": np.where(pur >= C.PURITY_CUT, np.array(V)[near], "generalist"),
                      "ref": ref, "guest": guest}, index=obs.index)
for c in ["cell_line", "time", "mix", "phase"]:
    cells[c] = obs[c].astype(str).values
cells.to_csv(out / "cells_groves5.csv.gz")
cols = V + ["generalist"]
lines = [f"# Groves closest-5% specialists: {args.geneset}/{args.tag}, k={k}", "",
         f"Reference cells: {int(ref.sum())}; specialists {100*(cells[ref].groves5!='generalist').mean():.1f}%"]
def add(title, df, note=None):
    lines.extend(["", f"## {title}", ""] + ([note, ""] if note else []) + [df.to_string()])
    print(f"\n== {title} ==\n{df.to_string()}")
def tab(sub, key, order, title):
    t = pd.crosstab(sub[key], sub["groves5"]).reindex(columns=cols, fill_value=0)
    t = t.reindex([o for o in order if o in t.index]); t["total"] = t.sum(1)
    add(f"{title}: % of each {key} in each class (rows sum to 100)", (100 * t[cols].div(t["total"], axis=0)).round(1).assign(n=t["total"]))
fwd_order = [x for x in C.TIME_ORDER if not x.endswith("_rm")]
tab(cells[ref], "time", fwd_order, "Fit cells")
if guest.any():
    tab(cells[guest], "time", [x for x in C.TIME_ORDER if x.endswith("_rm")], "Removal-arm guests")
if args.tag == "pooled":
    tab(cells[ref], "cell_line", C.LINES, "Fit cells")
tab(cells[ref], "phase", ["G1", "S", "G2M"], "Fit cells")
tab(cells[ref], "mix", sorted(cells.mix.unique()), "Fit cells")
comp = pd.DataFrame({r: cells[ref][r].value_counts().reindex(cols, fill_value=0) for r in ["groves5", "weight095", "gap035"]})
add("Rule comparison (counts, fit cells)", comp)
sp = pd.DataFrame({r: (100 * (cells[ref][r] != "generalist").groupby(cells[ref]["time"]).mean()).reindex(fwd_order)
                   for r in ["groves5", "weight095", "gap035"]}).round(1)
add("% of each time point called specialist, by rule", sp)
(out / "summary.md").write_text("\n".join(lines) + "\n")

# figure: PC1-PC2 with specialists coloured by vertex and by time
R = cells[ref]; P = pts[ref]; gen = (R.groves5 == "generalist").values
tcol = dict(zip(["0d", "8h", "1d", "3d", "7d"], plt.get_cmap("viridis")(np.linspace(0, .95, 5))))
vcol = dict(zip(V, plt.get_cmap("tab10")(range(k))))
fig, ax = plt.subplots(1, 2, figsize=(13, 5.8))
for a in ax:
    a.scatter(P[gen, 0], P[gen, 1], s=1, c="#d0d0d0", rasterized=True)
    for i in range(k):
        a.scatter(*arcs[i, :2], c="k", s=40, zorder=5); a.annotate(V[i], arcs[i, :2], fontsize=11, fontweight="bold", xytext=(5, 5), textcoords="offset points")
        a.add_patch(plt.Circle(arcs[i, :2], cut[i], fill=False, ls="--", lw=.8, color="0.3"))
    if k == 3:
        for i in range(3):
            j = (i + 1) % 3; a.plot(arcs[[i, j], 0], arcs[[i, j], 1], "k", lw=1)
    a.set_xlabel("PC1"); a.set_ylabel("PC2"); a.set_aspect("equal", "datalim")
for v in V:
    m = (R.groves5 == v).values; ax[0].scatter(P[m, 0], P[m, 1], s=4, c=[vcol[v]], label=f"{v} (n={m.sum()})", rasterized=True)
ax[0].legend(fontsize=8); ax[0].set_title("specialists by vertex")
for tm, c in tcol.items():
    m = (~gen) & (R.time == tm).values
    ax[1].scatter(P[m, 0], P[m, 1], s=4, c=[c], label=f"{tm} (n={m.sum()})", rasterized=True)
ax[1].legend(fontsize=8); ax[1].set_title("specialists by time")
fig.suptitle(f"{args.geneset}/{args.tag} k={k}", fontsize=11); fig.tight_layout()
fig.savefig(out / "specialists.png", dpi=150); plt.close(fig)
print("wrote", out)
