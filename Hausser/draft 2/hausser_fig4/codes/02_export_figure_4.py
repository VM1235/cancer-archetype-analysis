#!/usr/bin/env python3
"""Export Hausser Fig. 4a and 4c style panels.

4a — single cancer cells inside the METABRIC tetrahedron (PC1–PC3)
4c — same cells on METABRIC PC1, PC2, PC50 (alignment of ITH with ITD)

Matches scBC.R plotting: cells only + archetype spheres + grey edges
(no bulk-tumor cloud in panel a).

Usage:
  .venv/bin/python -u "Hausser/draft 2/hausser_fig4/codes/02_export_figure_4.py"
"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

from _paths import ANALYSIS

RESULTS = ANALYSIS / "results" / "projection"
FIGURES = ANALYSIS / "figures"

# scBC.R breast archetype colors
ARCH_COLORS = ["#007eb1", "#019e59", "#111111", "#BBBBBB"]
EDGE_COLOR = "#888888"


def draw_edges(ax, pts, color=EDGE_COLOR, lw=1.2):
    k = pts.shape[0]
    for i in range(k):
        for j in range(i + 1, k):
            ax.plot(
                [pts[i, 0], pts[j, 0]],
                [pts[i, 1], pts[j, 1]],
                [pts[i, 2], pts[j, 2]],
                color=color,
                lw=lw,
                alpha=0.85,
            )


def fit_view(ax, points, pad_frac=0.08):
    spans = np.ptp(points, axis=0)
    spans = np.maximum(spans, 1e-6)
    for i, setter in enumerate([ax.set_xlim, ax.set_ylim, ax.set_zlim]):
        lo, hi = points[:, i].min(), points[:, i].max()
        pad = spans[i] * pad_frac
        setter(lo - pad, hi + pad)
    try:
        ax.set_box_aspect(spans / spans.max())
    except AttributeError:
        pass


def tumor_colors(tumor_ids: pd.Series):
    cats = tumor_ids.astype("category")
    codes = cats.cat.codes.to_numpy()
    n = len(cats.cat.categories)
    # Distinct hues similar to paper (magenta/red/cyan/yellow/green/blue)
    base = plt.cm.tab10(np.linspace(0, 1, max(n, 10)))[:n]
    return codes, ListedColormap(base), list(cats.cat.categories)


def panel_3d(
    ax,
    sc,
    arcs,
    tumor_ids,
    title,
    axis_labels,
    elev=18,
    azim=-60,
    paper_style=True,
):
    """paper_style=True matches scBC.R / Fig 4: no axes, cells + tetrahedron only."""
    codes, cmap, labels = tumor_colors(tumor_ids)
    span = float(np.ptp(arcs))
    cell_size = max(4.0, 18.0 * (span / max(np.ptp(sc), 1e-6)) ** 0)  # keep readable
    ax.scatter(
        sc[:, 0],
        sc[:, 1],
        sc[:, 2],
        s=12 if paper_style else 8,
        c=codes,
        cmap=cmap,
        alpha=0.9,
        depthshade=True,
        linewidths=0,
    )
    draw_edges(ax, arcs, lw=1.6 if paper_style else 1.2)
    ax.scatter(
        arcs[:, 0],
        arcs[:, 1],
        arcs[:, 2],
        s=320 if paper_style else 180,
        c=ARCH_COLORS[: arcs.shape[0]],
        depthshade=False,
        edgecolors="white",
        linewidths=0.9,
        zorder=10,
    )
    # View uses archetypes as bounding box (like scBC.R xlim=range(archProj))
    fit_view(ax, arcs, pad_frac=0.12)
    ax.view_init(elev=elev, azim=azim)
    if paper_style:
        ax.set_axis_off()
    else:
        ax.set_xlabel(axis_labels[0], labelpad=-2, fontsize=8)
        ax.set_ylabel(axis_labels[1], labelpad=-2, fontsize=8)
        ax.set_zlabel(axis_labels[2], labelpad=-2, fontsize=8)
        ax.tick_params(labelsize=6, pad=-2)
    ax.set_title(title, fontsize=10, pad=2)
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=cmap(i),
            markersize=5,
            label=lab,
        )
        for i, lab in enumerate(labels)
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        fontsize=6,
        frameon=False,
        title="tumor",
        title_fontsize=7,
    )


def main():
    needed = [
        RESULTS / "sc_pc123.csv",
        RESULTS / "sc_pc1_pc2_pc50.csv",
        RESULTS / "archetypes_pc123.csv",
        RESULTS / "archetypes_pc1_pc2_pc50.csv",
        RESULTS / "sc_cell_metadata.csv",
    ]
    if not all(p.is_file() for p in needed):
        print("Missing projection outputs; run 01_project_sc_onto_metabric_shape.py first.")
        return 1

    meta = pd.read_csv(RESULTS / "sc_cell_metadata.csv")
    sc_a = pd.read_csv(RESULTS / "sc_pc123.csv", index_col=0).values
    arcs_a = pd.read_csv(RESULTS / "archetypes_pc123.csv", index_col=0).values
    sc_c = pd.read_csv(RESULTS / "sc_pc1_pc2_pc50.csv", index_col=0).values
    arcs_c = pd.read_csv(RESULTS / "archetypes_pc1_pc2_pc50.csv", index_col=0).values
    tumors = meta["tumor_id"]

    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"figure.facecolor": "white", "savefig.dpi": 300})

    # Combined 4a | 4c
    fig = plt.figure(figsize=(14, 6.5))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    panel_3d(
        ax1,
        sc_a,
        arcs_a,
        tumors,
        title="Fig. 4a — single cells on METABRIC tetrahedron (PC1–3)",
        axis_labels=("PC1", "PC2", "PC3"),
        elev=22,
        azim=-55,
    )
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    panel_3d(
        ax2,
        sc_c,
        arcs_c,
        tumors,
        title="Fig. 4c — single cells on METABRIC PC1, PC2, PC50",
        axis_labels=("PC1", "PC2", "PC50"),
        elev=18,
        azim=-70,
    )
    fig.suptitle(
        "Hausser et al. 2019 Fig. 4 reproduction\n"
        "(Karaayvaz scRNA-seq projected into original METABRIC ParTI shape)",
        fontsize=12,
        y=0.98,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIGURES / "Figure_4_hausser_reproduction.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print("Wrote", out)

    # Also save separate panels
    for name, sc, arcs, labels, elev, azim in [
        ("Figure_4a_tetrahedron.png", sc_a, arcs_a, ("PC1", "PC2", "PC3"), 22, -55),
        ("Figure_4c_pc1_pc2_pc50.png", sc_c, arcs_c, ("PC1", "PC2", "PC50"), 18, -70),
    ]:
        fig = plt.figure(figsize=(7, 6.5))
        ax = fig.add_subplot(111, projection="3d")
        panel_3d(ax, sc, arcs, tumors, title=name.replace(".png", ""), axis_labels=labels, elev=elev, azim=azim)
        fig.savefig(FIGURES / name, bbox_inches="tight")
        plt.close(fig)
        print("Wrote", FIGURES / name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
