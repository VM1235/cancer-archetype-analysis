#!/usr/bin/env python3
"""Figures for the best-simplex deck and report (dark theme). Reads best_fits.csv; example triangles are redrawn from the saved fits.
Run from Project_1/:  .venv/bin/python docs/emt_single_cell/best_simplex/figures.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import analyze as A
FIG = HERE / "figures"; FIG.mkdir(exist_ok=True)
BG, FG, MUT, GRID = "#14171c", "#e8e8e8", "#9aa0a6", "#2a2f37"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "text.color": FG, "axes.labelcolor": FG,
                     "xtick.color": MUT, "ytick.color": MUT, "axes.edgecolor": GRID, "font.size": 10, "axes.titlecolor": FG})
RC = {"gap": "#63b3ed", "groves5": "#f6ad55", "w95": "#b794f4"}
RN = {"gap": "gap rule (purity < 0.35)", "groves5": "closest-5% rule", "w95": "weight > 0.95 rule"}
B = pd.read_csv(HERE / "best_fits.csv")
NAME = A.NAME
ORDER = {"SCC": ["pooled", "SC15", "SC32", "SC46", "SC33", "SC2b1", "SC2b3", "SC35", "SC43", "MU3", "MU4", "MU7", "MU8"],
         "COOK": ["pooled (4 lines)", "A549", "DU145", "MCF7", "OVCA420"],
         "PDAC": ["pooled (in vivo)", "M1", "M1-PT", "M1-SS", "M1-Liver", "M1-Met", "M2", "M2-PT", "M2-Liver", "M2-Met"],
         "T4": ["pooled (6 clones)", "E", "EM1", "EM2", "EM3", "M", "4T1"],
         "GSE": ["pooled (4 samples)", "MCF10A", "MCF10A_Mock", "MCF10A_TGFB", "HuMEC", "HuMEC_Mock", "HuMEC_TGFB"]}

GSL = {"all": "all genes", "ks": "KS genes", "hallmark": "Hallmark"}
GSC = {"all": "#63b3ed", "ks": "#f6ad55", "hallmark": "#68d391"}
def label(r):
    gs = GSL[r.geneset]
    return f"{r.unit}  |  {gs}  |  {r.arm}"

def fig_pheat(ds):
    d = B[B.dataset == ds]
    cols = [("all", None), ("ks", None)]
    units = [u for u in ORDER[ds] if u in set(d.unit)]
    AO = ["light QC", "cleaned", "scRNA only", "six clean tumours"]
    arms = sorted(set(d.arm), key=lambda a: AO.index(a) if a in AO else 9)
    colkeys = [(g, a) for g in ("all", "ks", "hallmark") for a in arms]
    colkeys = [ck for ck in colkeys if len(d[(d.geneset == ck[0]) & (d.arm == ck[1])]) or ck[1] in ("light QC", "cleaned") and not (ds == "SCC" and ck[1] == "cleaned")]
    M = np.full((len(units), len(colkeys)), np.nan); K = np.full(M.shape, "", dtype=object)
    for i, u in enumerate(units):
        for j, (g, a) in enumerate(colkeys):
            r = d[(d.unit == u) & (d.geneset == g) & (d.arm == a)]
            if len(r):
                r = r.iloc[0]; M[i, j] = -np.log10(max(r.p, 1e-3)); K[i, j] = f"k={int(r.k)}\np={r.p:.3g}" + ("\u2020" if r.n_null < 150 else "")
    cmap = LinearSegmentedColormap.from_list("p", ["#222831", "#2b4a6f", "#2f855a", "#9ae6b4"])
    fig, ax = plt.subplots(figsize=(1.25 * len(colkeys) + 2.6, 0.62 * len(units) + 1.9))
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=2.3, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if K[i, j]:
                sig = M[i, j] >= -np.log10(0.05)
                ax.text(j, i, K[i, j], ha="center", va="center", fontsize=8, color="#0b0d10" if M[i, j] > 1.6 else FG, fontweight="bold" if sig else "normal")
            else:
                na = (ds == "T4" and units[i] == "E" and colkeys[j][1] == "cleaned")
                pooled_only = colkeys[j][1] in ("scRNA only", "six clean tumours") and not str(units[i]).startswith("pooled")
                ax.text(j, i, "n/a\n(142 cells)" if na else ("n/a\n(pooled-only arm)" if pooled_only else "not run"), ha="center", va="center", fontsize=7, color=MUT)
    ax.set_xticks(range(len(colkeys))); ax.set_xticklabels([GSL[g] + "\n" + a for g, a in colkeys], fontsize=8.5)
    ax.xaxis.tick_top(); ax.set_yticks(range(len(units))); ax.set_yticklabels(units)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title(f"{NAME[ds]}: most significant simplex per unit (green = lower p; bold = p < 0.05; \u2020 = 100 shuffles)", fontsize=10.5, pad=34)
    fig.tight_layout(); fig.savefig(FIG / f"p_heat_{ds}.png", dpi=170); plt.close(fig)

def fig_rules(ds):
    d = B[(B.dataset == ds)].copy()
    d["lab"] = d.apply(label, axis=1)
    d["sig"] = d.p <= 0.05
    d = d.sort_values(["sig", "p"], ascending=[False, True]).head(26).iloc[::-1]
    n = len(d); y = np.arange(n)
    fig, axs = plt.subplots(1, 2, figsize=(13.2, 0.36 * n + 2.2), sharey=True)
    for ax, kind, ttl in ((axs[0], "hyb_gen", "% of hybrid cells that are generalists"), (axs[1], "gen_hyb", "% of generalist cells that are hybrid")):
        for r in RC:
            ax.scatter(d[f"{kind}_{r}"], y, s=34, color=RC[r], label=RN[r], zorder=3, alpha=0.95)
        if kind == "gen_hyb":
            ax.scatter(d["hybrid_pct"], y, marker="|", s=170, color="#e8e8e8", label="hybrid share of all cells (baseline)", zorder=4)
        ax.set_xlim(0, 101); ax.set_xlabel("%"); ax.set_title(ttl, fontsize=10)
        ax.grid(axis="x", color=GRID); ax.set_axisbelow(True)
        ax.axvline(50, color=GRID, lw=1)
    axs[0].set_yticks(y); axs[0].set_yticklabels([f"{'' if s else '(ns) '}{l}  k={int(k)}  p={p:.3g}" for l, k, p, s in zip(d.lab, d.k, d.p, d.sig)], fontsize=7.8)
    axs[0].legend(loc="lower left", fontsize=7.5, frameon=False, bbox_to_anchor=(0, -0.2 - 0.9 / max(n, 1)), ncol=3)
    axs[1].legend(loc="lower right", fontsize=7.5, frameon=False, bbox_to_anchor=(1, -0.2 - 0.9 / max(n, 1)))
    fig.suptitle(f"{NAME[ds]}: hybrids and generalists, three rules (best simplex of each fit; significant first, (ns) = p > 0.05)", fontsize=11)
    fig.tight_layout(); fig.savefig(FIG / f"rules_{ds}.png", dpi=170, bbox_inches="tight"); plt.close(fig)

def fig_summary():
    s = B.copy(); s["sig"] = s.p <= 0.05
    rows = []
    for ds in ORDER:
        for g in ("all", "ks", "hallmark"):
            for lv in ("pooled", "within"):
                x = s[(s.dataset == ds) & (s.geneset == g) & (s.level == lv)]
                if len(x): rows.append((ds, g, lv, int(x.sig.sum()), len(x)))
    T = pd.DataFrame(rows, columns=["dataset", "geneset", "level", "n_sig", "n"])
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.2), sharey=True)
    for ax, lv in zip(axs, ("pooled", "within")):
        for gi, (g, col) in enumerate((("all", "#63b3ed"), ("ks", "#f6ad55"), ("hallmark", "#68d391"))):
            v = []
            for ds in ORDER:
                x = T[(T.dataset == ds) & (T.geneset == g) & (T.level == lv)]
                v.append(100 * x.n_sig.sum() / x.n.sum() if len(x) else np.nan)
            xs = np.arange(len(ORDER)) + (gi - 1) * 0.27
            ax.bar(xs, v, 0.25, color=col, label=GSL[g])
            for xi, vi, ds in zip(xs, v, ORDER):
                x = T[(T.dataset == ds) & (T.geneset == g) & (T.level == lv)]
                if len(x): ax.text(xi, (vi if np.isfinite(vi) else 0) + 2, f"{int(x.n_sig.sum())}/{int(x.n.sum())}", ha="center", fontsize=8, color=FG)
        ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([NAME[d].replace(" ", "\n", 1) for d in ORDER], fontsize=8.5)
        ax.set_ylim(0, 115); ax.set_title({"pooled": "Pooled fits (all samples together)", "within": "Within single units"}[lv]); ax.grid(axis="y", color=GRID); ax.set_axisbelow(True)
    axs[0].set_ylabel("% of fits with p ≤ 0.05 at the best k"); axs[0].legend(frameon=False, fontsize=9)
    fig.suptitle("How often does a simplex pass the permutation test? (fractions are fits passing / fits run; light QC and cleaned arms counted separately)", fontsize=10.5)
    fig.tight_layout(); fig.savefig(FIG / "summary_pass.png", dpi=170); plt.close(fig)

def fig_cross_rules():
    s = B[B.p <= 0.05].copy()
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4), sharey=True)
    for ax, r in zip(axs, RC):
        pos = 0; ticks = []; labs = []
        for ds in ORDER:
            for g, col, off in (("all", "#63b3ed", -0.25), ("ks", "#f6ad55", 0.0), ("hallmark", "#68d391", 0.25)):
                v = s[(s.dataset == ds) & (s.geneset == g)][f"hyb_gen_{r}"].values
                ax.scatter(np.full(len(v), pos + off) + np.random.default_rng(0).uniform(-.07, .07, len(v)), v, s=14, color=col, alpha=.8, zorder=3)
                if len(v): ax.hlines(np.median(v), pos + off - .1, pos + off + .1, color="white", lw=1.6, zorder=4)
            ticks.append(pos); labs.append(NAME[ds].replace(" ", "\n", 1)); pos += 1
        ax.set_xticks(ticks); ax.set_xticklabels(labs, fontsize=8); ax.set_title(RN[r], color=RC[r], fontsize=10); ax.grid(axis="y", color=GRID); ax.set_ylim(0, 102)
    axs[0].set_ylabel("% of hybrid cells that are generalists")
    axs[0].scatter([], [], color="#63b3ed", s=14, label="all genes"); axs[0].scatter([], [], color="#f6ad55", s=14, label="KS genes"); axs[0].scatter([], [], color="#68d391", s=14, label="Hallmark"); axs[0].legend(frameon=False, fontsize=8, loc="lower left")
    fig.suptitle("Significant simplices only (p ≤ 0.05): hybrids as generalists, one dot per fit, white bar = median", fontsize=10.5)
    fig.tight_layout(); fig.savefig(FIG / "cross_rules.png", dpi=170); plt.close(fig)

def fig_triangles():
    obs = A.load_obs()
    cand = {}
    for ds in ORDER:
        base = B[(B.dataset == ds) & (B.geneset == "all") & (B.p <= 0.05) & (B.k == 3) & (B.n_hybrid >= 200) & (B.hybrid_pct.between(10, 90))]
        d = base[base.level == "within"]
        if not len(d): d = base
        if len(d): cand[ds] = d.sort_values(["p", "n_cells"], ascending=[True, False]).iloc[0]
    fig, axs = plt.subplots(2, len(cand), figsize=(3.7 * len(cand), 7.4))
    for j, (ds, r) in enumerate(cand.items()):
        pts, V, o = A.cells_for(r, obs); gen, inside, near = A.rules(pts, V)
        hyb = o["cluster_annotation"].astype(str).isin(A.HYB).values
        rng = np.random.default_rng(0); idx = rng.permutation(len(o))[:5000]
        for row, (rule, ttl) in enumerate((("gap", "gap rule"), ("groves5", "closest-5% rule"))):
            ax = axs[row, j]; g = gen[rule]
            ax.scatter(pts[idx][~hyb[idx], 0], pts[idx][~hyb[idx], 1], s=3, color="#4a5568", alpha=.6, label="non-hybrid")
            m = hyb[idx] & g[idx]; ax.scatter(pts[idx][m, 0], pts[idx][m, 1], s=4, color="#f6ad55", alpha=.8, label="hybrid, generalist")
            m = hyb[idx] & ~g[idx]; ax.scatter(pts[idx][m, 0], pts[idx][m, 1], s=6, color="#fc8181", alpha=.95, label="hybrid, specialist")
            T = np.vstack([V, V[:1]]); ax.plot(T[:, 0], T[:, 1], color="#63b3ed", lw=1.4)
            ax.scatter(V[:, 0], V[:, 1], s=40, color="#63b3ed", zorder=5)
            ax.set_xticks([]); ax.set_yticks([])
            if row == 0: ax.set_title(f"{NAME[ds]}\n{r.unit} ({r.arm}), k=3, p={r.p:.3g}", fontsize=9)
            ax.set_xlabel(f"{ttl}: {100 * hyb[g].mean() if False else 100 * (hyb & g).sum() / hyb.sum():.0f}% of hybrids generalist", fontsize=8.5)
    axs[0, 0].legend(frameon=False, fontsize=7, loc="lower left", markerscale=3)
    fig.suptitle("Best all-genes triangle in each dataset: where hybrid cells sit (PC1 vs PC2 of the unit)", fontsize=11)
    fig.tight_layout(); fig.savefig(FIG / "triangles.png", dpi=170); plt.close(fig)

def fig_schematic():
    fig, ax = plt.subplots(figsize=(13, 3.8)); ax.axis("off"); ax.set_xlim(0, 13); ax.set_ylim(0, 3.8)
    steps = [("Saved fit", "PCA scores +\nvertices (k corners)\nper unit / gene set"), ("Barycentric\nweights", "each cell = mix of\nk corners (weights\nsum to 1)"),
             ("3 rules", "gap: top−second ≥ .35\ngroves: nearest 5%\nw95: weight > .95"), ("Generalist\nor specialist", "cell in no\nspecialist set =\ngeneralist"),
             ("Hybrid\nlabel", "Early-Hybrid,\nHybrid, Late-Hybrid\n(Cook: KS bins)"), ("Two numbers", "% hybrids that are\ngeneralists;\n% generalists that\nare hybrids")]
    for i, (t, d) in enumerate(steps):
        x = 0.15 + i * 2.15
        ax.add_patch(plt.Rectangle((x, 0.5), 1.8, 2.7, fc="#1d232b", ec="#63b3ed", lw=1.3))
        ax.text(x + .9, 2.78, t, ha="center", va="center", fontsize=10.5, fontweight="bold", color="#e8e8e8")
        ax.text(x + .9, 1.55, d, ha="center", va="center", fontsize=8.6, color="#b8bec6")
        if i < len(steps) - 1: ax.annotate("", xy=(x + 2.1, 1.85), xytext=(x + 1.82, 1.85), arrowprops=dict(arrowstyle="->", color="#f6ad55", lw=1.6))
    fig.savefig(FIG / "schematic.png", dpi=170, bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    fig_schematic(); fig_summary(); fig_cross_rules()
    for ds in ORDER: fig_pheat(ds); fig_rules(ds)
    fig_triangles()
    print("done")
