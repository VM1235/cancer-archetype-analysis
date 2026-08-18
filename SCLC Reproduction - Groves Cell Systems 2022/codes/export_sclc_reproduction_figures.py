#!/usr/bin/env python3
"""Write three official figures matching p2.pdf Figure 1A, 1B, and 1C."""

from pathlib import Path
import sys

SCLC = Path(__file__).resolve().parents[1]
ROOT = SCLC.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.enrichment import SUBTYPES

OUT = SCLC / "figures"
OUT.mkdir(parents=True, exist_ok=True)

PANEL_A = SCLC / "results" / "panel_a"
PANEL_B = SCLC / "results" / "panel_b"
PANEL_C = SCLC / "results" / "panel_c"

ARC_COLORS = {
    0: "#4C78A8",
    1: "#F58518",
    2: "#E45756",
    3: "#72B7B2",
    4: "#54A24B",
}


def style():
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.bbox": "tight",
            "savefig.dpi": 220,
        }
    )


def _polytope(ax, scores, arcs):
    ax.scatter(scores[:, 0], scores[:, 1], s=18, c="#4C78A8", alpha=0.75, label="cell lines")
    ax.scatter(arcs[:, 0], arcs[:, 1], s=80, c="#E45756", zorder=3, label="archetypes")
    for i in range(5):
        for j in range(i + 1, 5):
            ax.plot([arcs[i, 0], arcs[j, 0]], [arcs[i, 1], arcs[j, 1]], color="#E45756", lw=1, zorder=2)
        ax.annotate(f"A{i+1}", (arcs[i, 0], arcs[i, 1]), xytext=(3, 3), textcoords="offset points", fontsize=8)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("k = 5 polytope")
    ax.legend(frameon=False, loc="best", fontsize=8)


def fig_1a():
    esv = pd.read_csv(PANEL_A / "esv_curve.csv")
    tab = pd.read_csv(PANEL_A / "t_ratio_sanity_parti.csv")
    scores = np.load(PANEL_A / "pc_scores_12.npy")
    arcs = np.load(PANEL_A / "archetypes_k5_parti.npy")

    fig = plt.figure(figsize=(11.2, 6.4))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1.15, 1, 1], hspace=0.38, wspace=0.38)
    ax_poly = fig.add_subplot(gs[:, 0])
    ax_esv = fig.add_subplot(gs[0, 1])
    ax_elb = fig.add_subplot(gs[0, 2])
    ax_t = fig.add_subplot(gs[1, 1])
    ax_p = fig.add_subplot(gs[1, 2])

    _polytope(ax_poly, scores, arcs)

    ax_esv.plot(esv["k"], 100 * esv["esv"], marker="o", color="#4C78A8")
    ax_esv.axvline(5, color="0.5", ls="--", lw=1)
    ax_esv.set_xlabel("k")
    ax_esv.set_ylabel("ESV (%)")
    ax_esv.set_title("Explained sample variance")

    ax_elb.plot(esv["k"].iloc[1:], 100 * esv["delta_esv"].iloc[1:], marker="o", color="#F58518")
    ax_elb.axvline(5, color="0.5", ls="--", lw=1)
    ax_elb.set_xlabel("k")
    ax_elb.set_ylabel("ΔESV (pp)")
    ax_elb.set_title("Elbow")

    x = np.arange(len(tab))
    w = 0.35
    ax_t.bar(x - w / 2, tab["paper_t_ratio"], w, color="#9ecae1", label="paper")
    ax_t.bar(x + w / 2, tab["t_ratio"], w, color="#4C78A8", label="ours")
    ax_t.set_xticks(x)
    ax_t.set_xticklabels([f"k={int(k)}" for k in tab["k"]])
    ax_t.set_ylabel("t-ratio")
    ax_t.set_title("t-ratio")
    ax_t.legend(frameon=False, fontsize=8)

    ax_p.bar(x - w / 2, tab["paper_p"], w, color="#fdbb84", label="paper")
    ax_p.bar(x + w / 2, tab["p_value"], w, color="#E45756", label="ours")
    ax_p.axhline(0.05, color="0.3", ls="--", lw=1)
    ax_p.set_xticks(x)
    ax_p.set_xticklabels([f"k={int(k)}" for k in tab["k"]])
    ax_p.set_ylabel("p-value")
    ax_p.set_title("permutation p (k=5 smallest significant)")
    ax_p.legend(frameon=False, fontsize=8)

    fig.suptitle("Figure 1A  —  how many archetypes?", y=0.98, fontsize=13)
    fig.savefig(OUT / "Figure_1A.png")
    plt.close(fig)


def fig_1b():
    table = pd.read_csv(PANEL_B / "enrichment_author.csv")
    fig, axes = plt.subplots(1, 5, figsize=(13.5, 3.4), sharey=True)
    for ax, subtype in zip(axes, SUBTYPES):
        sub = table[table["subtype"] == subtype]
        for arc in sorted(sub["archetype"].unique()):
            arc_tab = sub[sub["archetype"] == arc].sort_values("bin")
            ax.plot(
                arc_tab["bin"],
                arc_tab["fold_enrichment"],
                marker="o",
                color=ARC_COLORS[int(arc)],
                label=f"Arc {int(arc) + 1}",
            )
            sig = arc_tab[arc_tab["sig_peak_at_bin0"]]
            if len(sig):
                ax.scatter(
                    sig["bin"],
                    sig["fold_enrichment"],
                    s=70,
                    facecolors="none",
                    edgecolors=ARC_COLORS[int(arc)],
                    linewidths=1.6,
                    zorder=3,
                )
        ax.axhline(1.0, color="0.7", lw=1, ls="--")
        ax.set_title(f"SCLC-{subtype}")
        ax.set_xlabel("distance bin")
        ax.set_xticks(range(10))
    axes[0].set_ylabel("fold enrichment")
    axes[-1].legend(loc="upper right", fontsize=7, frameon=False)
    fig.suptitle("Figure 1B  —  each known subtype enriches at one archetype (bin 0)", y=1.06, fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_1B.png")
    plt.close(fig)


def fig_1c():
    scores = pd.read_csv(PANEL_C / "combined_pca_scores.csv", index_col=0)
    arcs = pd.read_csv(PANEL_C / "archetypes_in_combined_pca.csv", index_col=0)
    sources = pd.read_csv(PANEL_C / "sample_sources.csv", index_col=0).iloc[:, 0]
    var = pd.read_csv(PANEL_C / "variance_explained.csv")
    null = np.load(PANEL_C / "ev_null_shuffles.npy")
    is_tumor = sources.eq("Tumor").values
    is_cell = ~is_tumor
    a = arcs.values
    ratio5 = float(var.loc[var["n_components"] == 5, "ratio_combined_over_tumor"].iloc[0])

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.2))
    ax = axes[0]
    ax.scatter(scores.values[is_cell, 0], scores.values[is_cell, 1], s=18, c="#9BD770", label="cell lines")
    ax.scatter(scores.values[is_tumor, 0], scores.values[is_tumor, 1], s=18, c="#2E7D32", label="tumors")
    ax.scatter(a[:, 0], a[:, 1], s=80, c="#E45756", zorder=3, label="cell-line archetypes")
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            ax.plot([a[i, 0], a[j, 0]], [a[i, 1], a[j, 1]], color="#E45756", lw=1, zorder=2)
        ax.annotate(f"A{i+1}", (a[i, 0], a[i, 1]), xytext=(3, 3), textcoords="offset points", fontsize=8)
    ax.set_xlabel("PC1 (combined)")
    ax.set_ylabel("PC2 (combined)")
    ax.set_title("Tumors lie inside the cell-line polytope")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    xs = var["n_components"]
    for row in null:
        ax.plot(xs, np.cumsum(row), color="0.75", lw=0.8, alpha=0.85)
    ax.plot(xs, var["cum_tumor_only"], marker="o", color="#4C78A8", label="tumor-only PCA")
    ax.plot(xs, var["cum_combined_on_tumors"], marker="o", color="#F58518", label="combined PCA on tumors")
    ax.axvline(5, color="0.5", ls="--", lw=1)
    ax.set_xlabel("number of components")
    ax.set_ylabel("cumulative % tumor variance")
    ax.set_title(f"5 combined PCs = {100 * ratio5:.0f}% of tumor-only ceiling")
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle("Figure 1C  —  the same space describes human tumors", y=1.02, fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "Figure_1C.png")
    plt.close(fig)


def main():
    style()
    fig_1a()
    fig_1b()
    fig_1c()
    keep = {"Figure_1A.png", "Figure_1B.png", "Figure_1C.png"}
    for path in OUT.glob("*.png"):
        if path.name not in keep:
            path.unlink()
    print("Wrote", sorted(p.name for p in OUT.glob("*.png")))


if __name__ == "__main__":
    main()
