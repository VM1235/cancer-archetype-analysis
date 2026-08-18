#!/usr/bin/env python3
"""Figure 1A in the Groves et al. layout, using our ParTI-matched analysis."""

from pathlib import Path
import sys

SCLC = Path(__file__).resolve().parents[1]
ROOT = SCLC.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.io import load_expression_csv
from src.pca import align_pca_signs, fit_pca, inverse_transform_scores

OUT = SCLC / "figures"
PANEL_A = SCLC / "results" / "panel_a"

# Paper-like vertex colors (A, A2, N, P, Y)
SUBTYPE_COLOR = {
    "A": "#E07A3D",
    "A2": "#3AA6A1",
    "N": "#E2C44F",
    "P": "#3D6FA8",
    "Y": "#D989B5",
}
# Panel B: archetype index -> subtype
ARC_TO_SUBTYPE = {0: "N", 1: "A2", 2: "A", 3: "Y", 4: "P"}


def style():
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.bbox": "tight",
            "savefig.dpi": 240,
        }
    )


def incremental_esv():
    esv = pd.read_csv(PANEL_A / "esv_curve.csv")
    pca_frac = pd.read_csv(PANEL_A / "pca_variance.csv")["cumulative"].iloc[-1]
    # ParTI plots extra variance as a percent of *total gene-space* variance
    gene_esv = 100.0 * esv["esv"].values * pca_frac
    k = esv["k"].values
    delta = np.diff(np.concatenate([[0.0], gene_esv]))
    return k, delta


def fig_1a():
    k_esv, delta = incremental_esv()
    tab_path = PANEL_A / "t_ratio_parti_1000.csv"
    if not tab_path.exists():
        tab_path = PANEL_A / "t_ratio_sanity_parti.csv"
    tab = pd.read_csv(tab_path).sort_values("k")

    expr = load_expression_csv()
    X = expr.T.values
    saved = np.load(PANEL_A / "pc_scores_12.npy")
    pca12, scores12 = fit_pca(X, n_components=12)
    pca12, _, _ = align_pca_signs(scores12, saved, pca12)
    arcs4 = np.load(PANEL_A / "archetypes_k5_parti.npy")
    arcs12 = np.zeros((5, 12))
    arcs12[:, : arcs4.shape[1]] = arcs4
    gene_arcs = inverse_transform_scores(pca12, arcs12)
    pca2, scores2 = fit_pca(X, n_components=2)
    arcs2 = pca2.transform(gene_arcs)

    fig = plt.figure(figsize=(6.6, 7.4))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 1.55], hspace=0.32, wspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(k_esv, delta, "-o", color="#3B6FA0", ms=5, lw=1.4)
    ax.axvline(5, color="0.35", ls="--", lw=1)
    ax.annotate(
        "Suggested number of\narchetypes by elbow",
        xy=(5, float(delta[k_esv == 5][0])),
        xytext=(7.2, max(delta) * 0.72),
        fontsize=7.5,
        color="0.25",
        arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8),
    )
    ax.set_xlim(1.5, 15.5)
    ax.set_ylim(0, max(12, delta.max() * 1.08))
    ax.set_xticks(range(2, 16))
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("% ESV on top of N−1 model")
    ax.set_title("Explained sample variance (ESV)")

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(tab["k"], tab["t_ratio"], "-o", color="#3B6FA0", ms=6, lw=1.4)
    ax.axvline(5, color="0.35", ls="--", lw=1)
    ymax = max(0.62, float(tab["t_ratio"].max()) * 1.15)
    ax.set_ylim(0, ymax)
    ax.set_xlim(2.5, 6.5)
    ax.set_xticks(tab["k"].astype(int))
    for _, row in tab.iterrows():
        p = row["p_value"]
        label = f"p = {p:.2f}" if p >= 0.01 else f"p = {p:.3f}"
        if int(row["k"]) == 5:
            label = f"p = {p:.3f}"
        ax.annotate(
            label,
            (row["k"], row["t_ratio"]),
            textcoords="offset points",
            xytext=(8, 6),
            fontsize=8,
            color="#222",
        )
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("t-ratio")
    ax.set_title("t-ratio of polytopes by number of vertices")

    ax = fig.add_subplot(gs[1, :])
    ax.scatter(scores2[:, 0], scores2[:, 1], s=14, c="#B0B0B0", zorder=1, linewidths=0)
    for i in range(5):
        for j in range(i + 1, 5):
            ax.plot(
                [arcs2[i, 0], arcs2[j, 0]],
                [arcs2[i, 1], arcs2[j, 1]],
                color="#888888",
                lw=0.9,
                zorder=2,
            )
    for i in range(5):
        subtype = ARC_TO_SUBTYPE[i]
        ax.scatter(
            arcs2[i, 0],
            arcs2[i, 1],
            s=90,
            c=SUBTYPE_COLOR[subtype],
            edgecolors="k",
            linewidths=0.4,
            zorder=3,
        )
        ax.annotate(
            f"{i+1}  SCLC-{subtype}",
            (arcs2[i, 0], arcs2[i, 1]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=8,
            color=SUBTYPE_COLOR[subtype],
            fontweight="bold",
        )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Archetype space (120 human SCLC cell lines)")

    fig.suptitle("Figure 1A  —  Archetype analysis on human SCLC cell lines", y=0.995, fontsize=12)
    fig.savefig(OUT / "Figure_1A.png")
    plt.close(fig)
    print("Wrote", OUT / "Figure_1A.png")


if __name__ == "__main__":
    style()
    fig_1a()
