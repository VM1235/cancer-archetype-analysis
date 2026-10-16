#!/usr/bin/env python3
"""Generate every figure used in Research_Report.md.

Two kinds of figure, kept apart by file name:
  schem_*.png  schematic / teaching illustrations. Toy data or diagrams. NOT results.
  emp_*.png    empirical figures computed from the saved pipeline outputs in
               Mouse_SCC/results/ and Mouse_SCC/data/processed/. Nothing is refit,
               except emp_v1_longgene.png, which computes a simple marker-gene score
               from the expression matrix (a report-time verification analysis).
Project figures produced by the pipeline itself are COPIED (unchanged) to
figures/project/.

Reads only; never writes into Mouse_SCC/results or Mouse_SCC/codes.

Usage (from Project_1/):
  .venv/bin/python Mouse_SCC/research_report/make_report_figures.py
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SCC = HERE.parent
ROOT = SCC.parent
RES = SCC / "results"
FIG = HERE / "figures"
PROJ = FIG / "project"
FIG.mkdir(exist_ok=True)
PROJ.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCC / "codes"))
from src.archetypes import fit_pcha, simplex_volume  # noqa: E402
from scipy.spatial import ConvexHull  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "white", "savefig.bbox": "tight", "savefig.dpi": 170})
STATE_ORDER = ["Epithelial", "Early-Hybrid", "Hybrid", "Late-EMT"]
SCOL = {"Epithelial": "#1f77b4", "Early-Hybrid": "#ff7f0e", "Hybrid": "#2ca02c", "Late-EMT": "#d62728"}
RUNS = {"all_main": "A. all genes, all cells", "all_scrna": "B. all genes, scRNA cells",
        "ks_main": "C. KS genes, all cells"}


def save(fig, name):
    fig.savefig(FIG / name)
    plt.close(fig)
    print("wrote", name)


def box(ax, xy, w, h, text, fc="#eef3fb", ec="#3b6fa0", fs=9):
    ax.add_patch(FancyBboxPatch(xy, w, h, boxstyle="round,pad=0.02", fc=fc, ec=ec, lw=1.2))
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center", fontsize=fs, wrap=True)


def arrow(ax, a, b, color="0.3"):
    ax.annotate("", xy=b, xytext=a, arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2))


# =========================================================== schematics
def schem_pipeline():
    fig, ax = plt.subplots(figsize=(13, 8.6))
    ax.set_xlim(0, 13); ax.set_ylim(0, 8.6); ax.axis("off")
    B = [
        ((0.2, 7.4), 3.0, 0.9, "GEO GSE288137\nSeurat object (.rds)\n23,959 cells x 20,358 genes", "#f6f6f6"),
        ((0.2, 6.0), 3.0, 0.9, "R export (done before this run)\ncounts.mtx, genes.txt,\nbarcodes.txt, cell_metadata.csv", "#f6f6f6"),
        ((0.2, 4.6), 3.0, 0.9, "01_build_anndata.py\nfilter genes (>=1% cells, no mt/ribo)\nlog1p(CP10k) -> 14,502 genes", "#eef3fb"),
        ((4.0, 7.4), 3.6, 0.9, "02_pca.py  (per version A/B/C)\nbalanced fit cells (<=800/tumour)\nexact PCA, project all cells", "#eef3fb"),
        ((4.0, 6.0), 3.6, 0.9, "02b_tumour_within_state.py\nhow much of each PC is tumour\nbeyond EMT state", "#eef3fb"),
        ((4.0, 4.6), 3.6, 0.9, "03_fit_and_test.py (run_fits.sh)\nPCHA k=3..7, t-ratio,\n200-shuffle permutation test, ESV", "#eef3fb"),
        ((8.4, 7.4), 4.3, 0.9, "choose k* (pre-registered rule)\nsmallest k with p<0.01, <= ESV elbow\n-> k*=3 for A and B; none for C", "#fff4e0"),
        ((8.4, 6.0), 4.3, 0.9, "04_describe_simplex.py --k 3\nweights, purity, label tables,\nenrichment, marker genes, figures", "#eef3fb"),
        ((8.4, 4.6), 4.3, 0.9, "05_robustness.py --k 3\nleave-one-tumour-out, scRNA-only,\nphase removed, new subsample", "#eef3fb"),
        ((8.4, 3.2), 4.3, 0.9, "06_groves_specialists.py --k 3\nGroves closest-5% specialists,\nrule comparison, triangle figure", "#eef3fb"),
        ((4.0, 3.2), 3.6, 0.9, "scc_common.py (shared code)\nloading, gene sets, balanced cells,\nPCA, barycentric, fit_simplex", "#eaf6ea"),
        ((0.2, 3.2), 3.0, 0.9, "src/archetypes.py (repo engine)\nfit_pcha, fit_pcha_best,\nsimplex_volume  -> py_pcha", "#eaf6ea"),
    ]
    for xy, w, h, t, c in B:
        box(ax, xy, w, h, t, fc=c)
    arrow(ax, (1.7, 7.4), (1.7, 6.9)); arrow(ax, (1.7, 6.0), (1.7, 5.5))
    arrow(ax, (3.2, 5.05), (4.0, 7.6)); arrow(ax, (5.8, 7.4), (5.8, 6.9)); arrow(ax, (5.8, 6.0), (5.8, 5.5))
    arrow(ax, (7.6, 5.05), (8.4, 7.6)); arrow(ax, (10.55, 7.4), (10.55, 6.9))
    arrow(ax, (10.55, 6.0), (10.55, 5.5)); arrow(ax, (10.55, 4.6), (10.55, 4.1))
    arrow(ax, (3.2, 3.65), (4.0, 3.65))
    ax.text(5.8, 2.95, "imported by 02-06", ha="center", fontsize=8, color="0.35")
    ax.text(0.2, 2.4, "Output folders: Mouse_SCC/results/<geneset>_<tag>/   e.g. all_main (A), all_scrna (B), ks_main (C)",
            fontsize=9)
    ax.text(0.2, 2.0, "Version A = all 14,502 genes, fit on all 12 tumours | B = all genes, fit on the 8 scRNA tumours only |"
            " C = 176 KS EMT genes, all tumours", fontsize=9)
    ax.text(0.2, 1.6, "Every version places all 23,959 cells into the fitted simplex after fitting.", fontsize=9)
    ax.text(0.2, 1.0, "Schematic. Grey = data provenance; blue = pipeline scripts; green = shared code; orange = decision.",
            fontsize=8, color="0.4")
    save(fig, "schem_pipeline.png")


def schem_emt_states():
    fig, ax = plt.subplots(figsize=(13, 4.6))
    ax.set_xlim(0, 13); ax.set_ylim(0, 4.6); ax.axis("off")
    st = [("Epithelial", "Epcam+, Cdh1+, Krt14+,\nEsrp1+, S100a14+;\nno Vim"),
          ("Early hybrid", "Epcam, Cdh1 lost;\nKrt14, Krt5, Itga6 kept;\nVim up"),
          ("Hybrid", "Krt14, Krt5 lost;\nKrt8 induced;\nmore mesenchymal genes"),
          ("Late EMT", "fibroblast-like:\nCol3a1, Pdgfrb,\nAspn, Postn")]
    cols = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
    for i, ((n, g), c) in enumerate(zip(st, cols)):
        x = 0.3 + i * 3.2
        box(ax, (x, 2.3), 2.7, 1.9, f"{n}\n\n{g}", fc=c + "22", ec=c, fs=9)
        if i < 3:
            arrow(ax, (x + 2.75, 3.25), (x + 3.15, 3.25))
    ax.text(0.3, 1.75, "Prrx1, Twist1 rise progressively from early hybrid to late EMT; Zeb1 is expressed in all states.",
            fontsize=9)
    ax.text(0.3, 1.35, "Slingshot (on Harmony embedding, epithelial root): two branches from epithelial: one through "
            "hybrid states to late EMT, one ending at the hybrid EMT state.", fontsize=9)
    ax.text(0.3, 0.95, "Cycling cells: similar proportions in epithelial, early hybrid and hybrid; markedly fewer in late EMT. "
            "No cluster dominated by cell-cycle genes.", fontsize=9)
    ax.text(0.3, 0.35, "Schematic summarising Perez-Gonzalez et al. 2026 (Nat Commun 17:8561), Results text for Fig. 1. "
            "Arrows show the ordering the paper infers, not measured time.", fontsize=8, color="0.4")
    save(fig, "schem_emt_states.png")


def schem_pareto():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    x = np.linspace(0, 1, 50)
    a1, a2 = np.array([0.2, 0.8]), np.array([0.8, 0.25])
    for r in (0.15, 0.3, 0.45):
        ax.add_patch(plt.Circle(a1, r, fill=False, color="#1f77b4", alpha=0.6))
        ax.add_patch(plt.Circle(a2, r, fill=False, color="#d62728", alpha=0.6))
    ax.plot(*zip(a1, a2), "k-", lw=2.5, label="Pareto front (line)")
    ax.scatter(*a1, s=80, c="#1f77b4", zorder=5); ax.scatter(*a2, s=80, c="#d62728", zorder=5)
    ax.text(a1[0] - 0.05, a1[1] + 0.06, "archetype 1\n(best at task 1)", ha="right", fontsize=9)
    ax.text(a2[0] + 0.03, a2[1] - 0.12, "archetype 2\n(best at task 2)", fontsize=9)
    p = np.array([0.55, 0.65])
    ax.scatter(*p, c="0.4", s=40); ax.text(p[0] + 0.02, p[1] + 0.02, "off the front:\nworse at both tasks", fontsize=8)
    ax.set_xlim(0, 1.1); ax.set_ylim(0, 1.1); ax.set_aspect("equal")
    ax.set_xlabel("trait 1 (e.g. gene-expression axis)"); ax.set_ylabel("trait 2")
    ax.set_title("Two tasks: optimal phenotypes lie on a line", fontsize=10)
    ax.legend(loc="lower left", fontsize=8)
    ax = axes[1]
    V = np.array([[0, 0], [1, 0], [0.5, 0.87]])
    ax.add_patch(Polygon(V, fill=False, lw=2))
    rng = np.random.default_rng(1)
    w = rng.dirichlet([1.2, 1.2, 1.2], 600)
    P = w @ V
    pur = np.sort(w, 1)[:, -1] - np.sort(w, 1)[:, -2]
    ax.scatter(P[:, 0], P[:, 1], c=np.where(pur >= 0.35, "#e6550d", "0.7"), s=6)
    for v, n in zip(V, ["task A", "task B", "task C"]):
        ax.text(v[0], v[1] + (0.05 if v[1] > 0 else -0.09), n, ha="center", fontsize=9, fontweight="bold")
    ax.scatter([], [], c="#e6550d", label="specialists (near a corner)")
    ax.scatter([], [], c="0.7", label="generalists (interior mixtures)")
    ax.legend(loc="upper right", fontsize=8); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("Three tasks: phenotypes fill a triangle (simplex)", fontsize=10)
    fig.text(0.5, -0.02, "Schematic (toy data). Shoval et al. 2012; Hart et al. 2015.", ha="center",
             fontsize=8, color="0.4")
    save(fig, "schem_pareto.png")


def schem_pca():
    rng = np.random.default_rng(0)
    X = rng.multivariate_normal([0, 0], [[3, 1.6], [1.6, 1.4]], 400)
    Xc = X - X.mean(0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = S ** 2 / (len(X) - 1)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    ax = axes[0]
    ax.scatter(X[:, 0], X[:, 1], s=6, c="0.6")
    for j, c in zip(range(2), ["#d62728", "#1f77b4"]):
        v = Vt[j] * 2 * np.sqrt(var[j])
        ax.annotate("", xy=v, xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=c, lw=2.5))
        ax.text(*(v * 1.1), f"PC{j + 1}: {100 * var[j] / var.sum():.0f}% of variance", color=c, fontsize=9)
    ax.set_aspect("equal"); ax.set_xlabel("gene 1 (centred)"); ax.set_ylabel("gene 2 (centred)")
    ax.set_title("Cells in gene space; PCs = directions of largest spread", fontsize=10)
    ax = axes[1]
    sc = U * S
    ax.scatter(sc[:, 0], sc[:, 1], s=6, c="0.6"); ax.set_aspect("equal")
    ax.set_xlabel("PC1 score"); ax.set_ylabel("PC2 score")
    ax.set_title("Same cells, rotated: PC scores", fontsize=10)
    fig.text(0.5, -0.03, "Schematic (toy 2-gene data). With 14,502 genes the idea is identical, "
             "only in more dimensions.", ha="center", fontsize=8, color="0.4")
    save(fig, "schem_pca.png")


def bary(points, V):
    a = np.vstack([V.T, np.ones((1, len(V)))])
    b = np.vstack([np.asarray(points).T, np.ones((1, len(points)))])
    return np.linalg.solve(a, b).T


def schem_barycentric():
    V = np.array([[0.0, 0.0], [4.0, 0.0], [1.0, 3.0]])
    pts = np.array([[1.4, 1.0], [3.4, 1.2]])
    W = bary(pts, V)
    fig, ax = plt.subplots(figsize=(7.5, 5.4))
    ax.add_patch(Polygon(V, fill=False, lw=2))
    for v, n in zip(V, ["V1 (0,0)", "V2 (4,0)", "V3 (1,3)"]):
        ax.scatter(*v, c="k", s=50); ax.text(v[0] - 0.1, v[1] + 0.15, n, fontsize=10, fontweight="bold")
    for p, w, c in zip(pts, W, ["#1f77b4", "#d62728"]):
        ax.scatter(*p, c=c, s=60, zorder=5)
        lab = ", ".join(f"{x:+.2f}" for x in w)
        ax.text(p[0] + 0.1, p[1] + 0.1, f"w = ({lab})", color=c, fontsize=10)
    ax.set_aspect("equal"); ax.set_xlim(-0.6, 5.6); ax.set_ylim(-0.6, 3.6)
    ax.set_title("Barycentric weights: a point as a weighted mix of the corners", fontsize=10)
    ax.text(-0.5, -0.5, "Blue: inside, all weights >= 0. Red: outside, one negative weight "
            "(the pipeline clips it to 0 and rescales).", fontsize=8.5)
    save(fig, "schem_barycentric.png")
    return W


def toy_triangle(rng, n=900, noise=0.06):
    V = np.array([[0, 0], [1, 0], [0.5, 0.87]])
    w = rng.dirichlet([0.7, 0.7, 0.7], n)
    return w @ V + rng.normal(0, noise, (n, 2)), V


def schem_pcha_and_tratio():
    rng = np.random.default_rng(3)
    X, Vtrue = toy_triangle(rng)
    np.random.seed(0)
    best = None
    for _ in range(10):
        a, w, e = fit_pcha(X, 3, delta=0.0)
        vol = simplex_volume(a)
        if best is None or vol > best[1]:
            best = (a, vol, e)
    arcs, vol, esv = best
    hull = ConvexHull(X)
    t_real = vol / hull.volume
    Xs = X.copy()
    for j in range(2):
        rng.shuffle(Xs[:, j])
    best_s = None
    for _ in range(10):
        a, w, e = fit_pcha(Xs, 3, delta=0.0)
        v = simplex_volume(a)
        if best_s is None or v > best_s[1]:
            best_s = (a, v)
    hull_s = ConvexHull(Xs)
    t_shuf = best_s[1] / hull_s.volume
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    for ax, D, H, A, t, ttl in [(axes[0], X, hull, arcs, t_real, "toy triangle data"),
                                (axes[1], Xs, hull_s, best_s[0], t_shuf, "same data, each axis shuffled")]:
        ax.scatter(D[:, 0], D[:, 1], s=4, c="0.65")
        hv = D[H.vertices]
        ax.add_patch(Polygon(hv, fill=False, ec="#3b6fa0", lw=1.5, ls="--"))
        ax.add_patch(Polygon(A, fill=True, fc="#e6550d33", ec="#e6550d", lw=2))
        ax.scatter(A[:, 0], A[:, 1], c="#e6550d", s=60, zorder=5)
        ax.set_aspect("equal")
        ax.set_title(f"{ttl}\nt-ratio = area(orange) / area(blue dashed) = {t:.2f}", fontsize=10)
    axes[0].plot(*np.vstack([Vtrue, Vtrue[0]]).T, "k:", lw=1, label="true corners (toy)")
    axes[0].legend(fontsize=8, loc="upper right")
    fig.text(0.5, -0.02, "Schematic computed with the project's own PCHA engine on toy data "
             f"(delta=0, best of 10 starts). ESV of the left fit = {esv:.3f}.", ha="center",
             fontsize=8, color="0.4")
    save(fig, "schem_pcha_tratio.png")
    return dict(t_real=t_real, t_shuf=t_shuf, esv=esv)


def schem_rules():
    rng = np.random.default_rng(3)
    X, V = toy_triangle(rng)
    np.random.seed(0)
    arcs = max((fit_pcha(X, 3, delta=0.0) for _ in range(10)), key=lambda r: simplex_volume(r[0]))[0]
    w = bary(X, arcs)
    c = np.clip(w, 0, None); c = c / c.sum(1, keepdims=True)
    s = np.sort(c, 1); gap = s[:, -1] - s[:, -2]; near = np.argmax(c, 1)
    d = np.linalg.norm(X[:, None] - arcs[None], axis=2)
    cut = np.quantile(d, 0.05, axis=0)
    g5 = (d <= cut).any(1)
    rules = [("purity gap >= 0.35 (Pal work)", gap >= 0.35),
             ("closest 5% per corner (Groves Fig. 3E)", g5),
             ("weight > 0.95 (Groves Methods)", c.max(1) > 0.95)]
    cols = np.array(["#7b3294", "#d7301f", "#1a9850"])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, (n, m) in zip(axes, rules):
        ax.scatter(X[~m, 0], X[~m, 1], s=4, c="0.8")
        ax.scatter(X[m, 0], X[m, 1], s=5, c=cols[near[m]])
        ax.add_patch(Polygon(arcs, fill=False, lw=1.5))
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"{n}\n{100 * m.mean():.0f}% of cells are specialists", fontsize=10)
    fig.text(0.5, -0.02, "Schematic (same toy data as the t-ratio figure). Grey = generalists; colour = "
             "specialist of that corner.", ha="center", fontsize=8, color="0.4")
    save(fig, "schem_rules.png")


def schem_confounding():
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(0, 10); ax.set_ylim(0, 5); ax.axis("off")
    box(ax, (0.3, 3.6), 2.4, 0.9, "Tumour\n(12 mice/tumours)", fc="#f6f6f6")
    box(ax, (0.3, 2.0), 2.4, 0.9, "Assay\nscRNA (cells) vs\nmultiome (nuclei)", fc="#f6f6f6")
    box(ax, (0.3, 0.4), 2.4, 0.9, "EMT state\n(the biology of interest)", fc="#fff4e0")
    box(ax, (4.0, 2.0), 2.5, 0.9, "Expression profile\nof each cell", fc="#eef3fb")
    box(ax, (7.4, 2.0), 2.3, 0.9, "Leading PCs ->\nsimplex corners", fc="#eaf6ea")
    for y in (4.05, 2.45, 0.85):
        arrow(ax, (2.7, y), (4.0, 2.45))
    arrow(ax, (6.5, 2.45), (7.4, 2.45))
    ax.annotate("", xy=(1.5, 2.9), xytext=(1.5, 3.6), arrowprops=dict(arrowstyle="<->", color="#d62728"))
    ax.annotate("", xy=(1.0, 1.3), xytext=(1.0, 2.0), arrowprops=dict(arrowstyle="<->", color="#d62728"))
    ax.annotate("", xy=(2.2, 1.3), xytext=(2.2, 3.6), arrowprops=dict(arrowstyle="<->", color="#d62728",
                                                                      connectionstyle="arc3,rad=-0.3"))
    ax.text(3.0, 4.4, "red arrows = the three factors are correlated in this dataset\n"
            "(each tumour has one assay; tumours differ in state mix)", fontsize=8.5, color="#d62728")
    ax.text(0.3, 0.05, "Schematic. A corner can arise from any of the three inputs; the design cannot fully separate them.",
            fontsize=8, color="0.4")
    save(fig, "schem_confounding.png")


# =========================================================== empirical
def emp_data():
    m = pd.read_csv(SCC / "data" / "processed" / "cell_metadata.csv")
    t = pd.crosstab(m["sample"], m["cluster_annotation"]).reindex(columns=STATE_ORDER)
    order = ["MU3", "MU4", "MU7", "MU8", "SC15", "SC32", "SC33", "SC46", "SC2b1", "SC2b3", "SC35", "SC43"]
    t = t.reindex(order)
    p = 100 * t.div(t.sum(1), axis=0)
    q = m.groupby("sample")[["nCount_RNA", "nFeature_RNA", "pct_counts_mt"]].median().reindex(order)
    assay = m.groupby("sample")["assay"].first().reindex(order)
    fig, axes = plt.subplots(1, 4, figsize=(18, 5), gridspec_kw=dict(width_ratios=[2.2, 1, 1, 1]))
    ax = axes[0]
    left = np.zeros(len(p))
    for s in STATE_ORDER:
        ax.barh(range(len(p)), p[s], left=left, color=SCOL[s], label=s)
        left += p[s].values
    ax.set_yticks(range(len(p)))
    ax.set_yticklabels([f"{i} ({assay[i]}, n={t.loc[i].sum()})" for i in p.index], fontsize=8)
    ax.invert_yaxis(); ax.set_xlabel("% of tumour's cells"); ax.legend(fontsize=8, ncol=4, loc="lower center",
                                                                       bbox_to_anchor=(0.5, 1.0))
    ax.set_title("EMT state mix per tumour", pad=24)
    for ax, col, lab in zip(axes[1:], ["nCount_RNA", "nFeature_RNA", "pct_counts_mt"],
                            ["median UMI counts", "median genes detected", "median % mitochondrial"]):
        ax.barh(range(len(q)), q[col], color=["#8c6bb1" if assay[i] == "multiome" else "#41ab5d" for i in q.index])
        ax.set_yticks(range(len(q))); ax.set_yticklabels(q.index, fontsize=8); ax.invert_yaxis()
        ax.set_title(lab)
    axes[3].text(0.98, -0.12, "purple = multiome, green = scRNA", transform=axes[3].transAxes, ha="right", fontsize=8)
    save(fig, "emp_data_composition.png")
    return t, q


def emp_pc_drivers():
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
    for n, (ax, (run, lab)) in enumerate(zip(axes, RUNS.items())):
        a = pd.read_csv(RES / run / "pc_drivers.csv").set_index("PC")
        b = pd.read_csv(RES / run / "pc_drivers_within_state.csv").set_index("PC")
        M = pd.DataFrame({"state": a["state"], "clusterU": a["clusterU"], "tumour": a["tumour"],
                          "tumour within state": b["tumour_in_state"], "assay": a["assay"],
                          "phase": a["phase"], "log depth": a["log_depth"], "% mito": a["pct_mito"]})
        im = ax.imshow(M.T.values, vmin=0, vmax=1, cmap="magma_r", aspect="auto")
        ax.set_yticks(range(M.shape[1]))
        ax.set_yticklabels(M.columns if n == 0 else [], fontsize=8)
        ax.set_xticks(range(10)); ax.set_xticklabels([f"PC{i}\n{v:.1f}%" for i, v in
                                                      zip(M.index, a["pct_variance"])], fontsize=7)
        for i in range(M.shape[1]):
            for j in range(10):
                v = M.T.values[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if v > 0.5 else "black")
        ax.set_title(lab, fontsize=10)
    fig.colorbar(im, ax=axes, shrink=0.8, label="share of PC variance explained (eta^2 or R^2)")
    save(fig, "emp_pc_drivers.png")


def emp_scree():
    fig, ax = plt.subplots(figsize=(7, 4))
    for run, lab in RUNS.items():
        evr = np.load(RES / run / "pca.npz")["evr"]
        ax.plot(range(1, 31), 100 * evr[:30], "o-", ms=3, label=lab)
    ax.axhline(0.1, ls="--", c="0.5", lw=0.8); ax.set_yscale("log")
    ax.set_xlabel("principal component"); ax.set_ylabel("% of total variance (log scale)")
    ax.legend(fontsize=8); ax.set_title("Scree curves (fit cells)")
    save(fig, "emp_scree.png")


def emp_null_and_esv():
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5))
    for col, (run, lab) in enumerate(RUNS.items()):
        ax = axes[0, col]
        nl = np.load(RES / run / "null_k3.npy")
        f = json.load(open(RES / run / "fit_k3.json"))
        ax.hist(nl, bins=25, color="0.7", edgecolor="white")
        ax.axvline(f["t_ratio"], color="#d62728", lw=2.5, label=f"observed t = {f['t_ratio']:.3f}")
        ax.set_title(f"{lab}\nk = 3: p = {f['p_value']:.3f} (200 shuffles)", fontsize=10)
        ax.set_xlabel("t-ratio of shuffled data"); ax.legend(fontsize=8)
        ax = axes[1, col]
        s = pd.read_csv(RES / run / "fit_summary.csv")
        e = pd.read_csv(RES / run / "esv_curve.csv")
        ax.plot(e["k"], e["esv"], "o-", c="#3b6fa0", label="ESV (first 10 PCs)")
        ax.set_xlabel("k (number of corners)"); ax.set_ylabel("ESV")
        ax2 = ax.twinx()
        ax2.plot(s["k"], -np.log10(s["p_value"]), "s--", c="#d62728", label="-log10 p")
        ax2.axhline(2, color="#d62728", lw=0.7, ls=":")
        ax2.set_ylabel("-log10 p (dotted: p = 0.01)", color="#d62728")
        ax.set_title("ESV curve and permutation p by k", fontsize=10)
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower right")
    fig.tight_layout()
    save(fig, "emp_null_esv.png")


def emp_where():
    c = pd.read_csv(RES / "all_scrna" / "groves5_k3" / "cells_groves5.csv.gz", index_col=0)
    r = c[c["reference_cell"]]
    fig, axes = plt.subplots(1, 3, figsize=(17, 4))
    cl = {"V1": "#7b3294", "V2": "#d7301f", "V3": "#1a9850", "generalist": "#cccccc"}
    for ax, rule, lab in zip(axes, ["gap035", "groves5", "weight095"],
                             ["purity gap >= 0.35", "Groves closest 5%", "weight > 0.95"]):
        t = pd.crosstab(r["cluster_annotation"], r[rule]).reindex(STATE_ORDER)
        t = t.reindex(columns=["V1", "V2", "V3", "generalist"], fill_value=0)
        p = 100 * t.div(t.sum(1), axis=0)
        left = np.zeros(len(p))
        for k in p.columns:
            ax.barh(range(len(p)), p[k], left=left, color=cl[k], label=k)
            for i, (lv, v) in enumerate(zip(left, p[k])):
                if v >= 6:
                    ax.text(lv + v / 2, i, f"{v:.0f}", ha="center", va="center", fontsize=8)
            left += p[k].values
        ax.set_yticks(range(len(p))); ax.set_yticklabels(p.index); ax.invert_yaxis()
        ax.set_xlabel("% of the state's scRNA cells"); ax.set_title(lab)
    axes[0].legend(fontsize=8, ncol=4, loc="lower center", bbox_to_anchor=(1.6, 1.08))
    save(fig, "emp_state_classes.png")


def emp_robustness():
    rows = []
    for run, var in [("all_scrna", "loto"), ("all_scrna", "nophase"), ("all_scrna", "seed"), ("all_main", "loto")]:
        d = pd.read_csv(RES / run / "robustness" / f"k3_{var}.csv")
        for _, x in d.iterrows():
            rows.append(dict(fit="B scRNA" if run == "all_scrna" else "A all cells",
                             variant=x["variant"], agree=x["agreement_pct"], t=x["t_ratio"]))
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(13, 4.3))
    colr = ["#41ab5d" if f.startswith("B") else "#8c6bb1" for f in df["fit"]]
    ax.bar(range(len(df)), df["agree"], color=colr)
    ax.set_xticks(range(len(df))); ax.set_xticklabels([f"{f[0]}: {v}" for f, v in zip(df["fit"], df["variant"])],
                                                     rotation=70, fontsize=7)
    ax.set_ylim(75, 100); ax.set_ylabel("% of all 23,959 cells with the same\nnearest corner as the main fit")
    ax.set_title("Robustness of the k = 3 fits (green: version B refits; purple: version A refits)")
    for i, (a, t) in enumerate(zip(df["agree"], df["t"])):
        ax.text(i, a + 0.3, f"{a:.0f}", ha="center", fontsize=7)
    save(fig, "emp_robustness.png")
    return df


def emp_v1_longgene():
    import anndata as ad
    mk = pd.read_csv(RES / "all_scrna" / "describe_k3" / "vertex_marker_genes.csv")
    out = {}
    adm = ad.read_h5ad(SCC / "data" / "processed" / "expr_log_filtered.h5ad")
    for v in ["V1", "V2", "V3"]:
        g = mk[f"{v}_up"].tolist()[:20]
        idx = [adm.var_names.get_loc(x) for x in g]
        out[v] = np.asarray(adm.X[:, idx].mean(1)).ravel()
    o = adm.obs[["sample", "assay", "cluster_annotation"]].astype(str).copy()
    for v in out:
        o[v] = out[v]
    order = ["MU3", "MU4", "MU7", "MU8", "SC15", "SC46", "SC32", "SC33", "SC2b1", "SC2b3", "SC35", "SC43"]
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.6), gridspec_kw=dict(width_ratios=[2, 1]))
    ax = axes[0]
    data = [o.loc[o["sample"] == s, "V1"].values for s in order]
    bp = ax.boxplot(data, showfliers=False, patch_artist=True)
    for patch, s in zip(bp["boxes"], order):
        patch.set_facecolor("#8c6bb1" if s.startswith("MU") else "#41ab5d")
    ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(order)
    ax.set_ylabel("mean log-expression of\nV1's top-20 marker genes")
    ax.set_title("V1 'long-gene' score per tumour (purple = multiome, green = scRNA)")
    ax = axes[1]
    med = o.groupby(["assay", "cluster_annotation"])["V1"].median().unstack().reindex(columns=STATE_ORDER)
    for i, (assay, row) in enumerate(med.iterrows()):
        ax.bar(np.arange(4) + i * 0.4, row.values, width=0.4,
               color="#8c6bb1" if assay == "multiome" else "#41ab5d", label=assay)
    ax.set_xticks(np.arange(4) + 0.2); ax.set_xticklabels(STATE_ORDER, fontsize=8)
    ax.set_title("median V1 score by state and assay"); ax.legend(fontsize=8)
    fig.text(0.5, -0.04, "Report-time verification analysis (not part of the pipeline): genes = top 20 of "
             "results/all_scrna/describe_k3/vertex_marker_genes.csv, column V1_up.", ha="center", fontsize=8,
             color="0.4")
    save(fig, "emp_v1_longgene.png")
    return o.groupby(["assay", "sample"])["V1"].median().round(3)


def copy_project():
    pairs = {
        "all_main/pc1_pc2_by_label.png": "A_pc1_pc2_by_label.png",
        "all_scrna/pc1_pc2_by_label.png": "B_pc1_pc2_by_label.png",
        "ks_main/pc1_pc2_by_label.png": "C_pc1_pc2_by_label.png",
        "all_main/fit_and_null.png": "A_fit_and_null.png",
        "all_scrna/fit_and_null.png": "B_fit_and_null.png",
        "ks_main/fit_and_null.png": "C_fit_and_null.png",
        "all_main/describe_k3/simplex_PC1_PC2.png": "A_k3_simplex.png",
        "all_main/describe_k4/simplex_PC1_PC2.png": "A_k4_simplex_PC1_PC2.png",
        "all_scrna/describe_k3/simplex_PC1_PC2.png": "B_k3_simplex.png",
        "all_scrna/describe_k4/simplex_PC1_PC2.png": "B_k4_simplex_PC1_PC2.png",
        "all_scrna/describe_k4/simplex_PC1_PC3.png": "B_k4_simplex_PC1_PC3.png",
        "ks_main/describe_k3/simplex_PC1_PC2.png": "C_k3_simplex.png",
        "all_scrna/groves5_k3/triangle_groves5.png": "B_k3_groves5.png",
        "all_main/scree.png": "A_scree.png",
    }
    for src, dst in pairs.items():
        shutil.copy2(RES / src, PROJ / dst)
    print("copied", len(pairs), "project figures")


if __name__ == "__main__":
    schem_pipeline(); schem_emt_states(); schem_pareto(); schem_pca()
    W = schem_barycentric(); print("barycentric toy weights:", np.round(W, 3).tolist())
    print("toy t-ratios:", schem_pcha_and_tratio())
    schem_rules(); schem_confounding()
    emp_data(); emp_pc_drivers(); emp_scree(); emp_null_and_esv(); emp_where()
    print(emp_robustness().round(2).to_string())
    print(emp_v1_longgene().to_string())
    copy_project()
