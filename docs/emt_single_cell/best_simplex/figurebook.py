#!/usr/bin/env python3
"""Per-analysis figure book: for each dataset and gene set (A entire / B KS / C Hallmark EMT) one page per 3 units
(pooled first, then every tumour / cell line / library / clone). Per unit: the best simplex (lowest p; ties -> cleaned arm, then smaller k)
drawn on PC1 vs PC2 with cells coloured by EMT state, and by specialist corner / generalist under the gap, closest-5% and weight>0.95 rules.
Run from Project_1/:  .venv/bin/python docs/emt_single_cell/best_simplex/figurebook.py   (reads best_fits.csv; run analyze.py first)
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import analyze as A
import figures as FG
OUT = HERE / "figures" / "book"; OUT.mkdir(parents=True, exist_ok=True)
BG, FGC, MUT = FG.BG, FG.FG, FG.MUT
GS = {"all": "A. Entire gene set", "ks": "B. KS EMT genes", "hallmark": "C. Hallmark EMT genes"}
STATES = {"SCC": ["Epithelial", "Early-Hybrid", "Hybrid", "Late-EMT"]}
SCOL = ["#4299e1", "#ecc94b", "#ed8936", "#e53e3e"]
VCOL = ["#63b3ed", "#f6ad55", "#b794f4", "#68d391", "#fc8181", "#f687b0", "#a0aec0"]
ARM_PRI = {"cleaned": 0, "forward arm": 1, "raw": 1, "all cells": 1, "light QC": 2, "scRNA only": 3, "six clean tumours": 4}
PER = 3
B = pd.read_csv(HERE / "best_fits.csv")

def pick(ds, unit, gs):
    d = B[(B.dataset == ds) & (B.unit == unit) & (B.geneset == gs)].copy()
    if d.empty: return None
    d["pri"] = d.arm.map(lambda a: ARM_PRI.get(a, 9))
    return d.sort_values(["p", "pri", "k"]).iloc[0]

def outline(ax, V, color):
    if len(V) == 3:
        T = np.vstack([V, V[:1]]); ax.plot(T[:, 0], T[:, 1], color=color, lw=1.1)
    else:
        P = V[:, :2]
        try:
            h = ConvexHull(P); T = np.vstack([P[h.vertices], P[h.vertices][:1]]); ax.plot(T[:, 0], T[:, 1], color=color, lw=1.1)
        except Exception: pass
    ax.scatter(V[:, 0], V[:, 1], s=16, color=color, zorder=5)

def draw(ax, pts, V, kind, gen, near, state_idx, stat=None):
    rng = np.random.default_rng(0); n = len(pts); idx = rng.permutation(n)[:3500]
    x, y = pts[idx, 0], pts[idx, 1]
    if kind == "state":
        for j in range(4):
            m = state_idx[idx] == j
            ax.scatter(x[m], y[m], s=2.2, color=SCOL[j], alpha=.7, lw=0)
    else:
        g = gen[idx]
        ax.scatter(x[g], y[g], s=1.6, color="#3a414c", alpha=.7, lw=0)
        for j in range(V.shape[0]):
            m = (~g) & (near[idx] == j)
            ax.scatter(x[m], y[m], s=3.2, color=VCOL[j % 7], alpha=.95, lw=0)
    outline(ax, V[:, :2], "#e8e8e8")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_color("#2a2f37")
    if stat: ax.set_xlabel(stat, fontsize=6.6, color=MUT, labelpad=2)

def page(ds, gs, units, pg, npg, obs):
    fig, axs = plt.subplots(PER, 4, figsize=(13.3, 7.3), gridspec_kw=dict(left=0.01, right=0.99, top=0.90, bottom=0.09, wspace=0.05, hspace=0.30))
    names = ["EMT state (hybrid = middle states)", "gap rule: specialists by corner", "closest-5% rule", "weight > 0.95 rule"]
    for r in range(PER):
        for c in range(4): axs[r, c].set_facecolor(BG); axs[r, c].axis("off") if r >= len(units) else None
    for r, u in enumerate(units):
        row = pick(ds, u, gs)
        if row is None:
            for c in range(4):
                axs[r, c].text(.5, .5, f"{u}: not run" if c == 0 else "", ha="center", va="center", color=MUT, fontsize=9); axs[r, c].set_xticks([]); axs[r, c].set_yticks([])
            continue
        pts, V, o = A.cells_for(row, obs); gen, inside, near = A.rules(pts, V)
        st = o["cluster_annotation"].astype(str).values
        order = STATES.get(ds, ["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"])
        sidx = np.array([order.index(x) if x in order else 0 for x in st])
        hyb = np.isin(st, list(A.HYB))
        ttl = f"{u}  |  {row.arm}  |  k={int(row.k)}  |  p={row.p:.3g}{'' if row.p <= 0.05 else '  (NOT significant)'}  |  {len(o):,} cells"
        for c, (kind, rule) in enumerate((("state", None), ("rule", "gap"), ("rule", "groves5"), ("rule", "w95"))):
            ax = axs[r, c]
            stat = None
            if rule:
                g = gen[rule]; stat = f"{100 * (hyb & g).sum() / max(hyb.sum(), 1):.0f}% of hybrids generalist | {100 * (hyb & g).sum() / max(g.sum(), 1):.0f}% of generalists hybrid"
            else:
                stat = f"hybrids = {100 * hyb.mean():.0f}% of cells; {100 * inside.mean():.0f}% inside simplex"
            draw(ax, pts, V, kind, gen[rule] if rule else None, near, sidx, stat)
            if r == 0 and c > 0: ax.set_title(names[c], fontsize=8.3, color=MUT, pad=14)
            if r == 0 and c == 0: ax.set_title("", pad=14)
        axs[r, 0].text(0.0, 1.015, ttl, transform=axs[r, 0].transAxes, fontsize=8.2, color=FGC if row.p <= 0.05 else "#fc8181", ha="left", va="bottom", fontweight="bold")
    lg = [plt.Line2D([], [], marker="o", ls="", color=SCOL[j], markersize=4, label=n) for j, n in enumerate(STATES.get(ds, ["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"]))]
    lg += [plt.Line2D([], [], marker="o", ls="", color="#3a414c", markersize=4, label="generalist")] + [plt.Line2D([], [], marker="o", ls="", color=VCOL[j], markersize=4, label=f"specialist V{j + 1}") for j in range(6)]
    fig.legend(handles=lg, loc="lower center", ncol=11, frameon=False, fontsize=7, labelcolor=FGC)
    fig.suptitle(f"{FG.NAME[ds]} | {GS[gs]} | page {pg}/{npg}", fontsize=10.5, color=FGC)
    fn = OUT / f"{ds}_{gs}_p{pg}.png"; fig.savefig(fn, dpi=150); plt.close(fig); return fn

def main():
    obs = A.load_obs(); made = []
    for ds in FG.ORDER:
        for gs in ("all", "ks", "hallmark"):
            have = [u for u in FG.ORDER[ds] if pick(ds, u, gs) is not None]
            if not have: continue
            pages = [FG.ORDER[ds][i:i + PER] for i in range(0, len(FG.ORDER[ds]), PER)]
            pages = [[u for u in p if pick(ds, u, gs) is not None or True] for p in pages]
            # drop trailing pages with no fits at all
            pages = [p for p in pages if any(pick(ds, u, gs) is not None for u in p)]
            for i, p in enumerate(pages, 1):
                made.append(page(ds, gs, p, i, len(pages), obs)); print("wrote", made[-1].name, flush=True)
    (OUT / "index.txt").write_text("\n".join(f.name for f in made))

if __name__ == "__main__":
    main()
