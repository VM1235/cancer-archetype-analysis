#!/usr/bin/env python3
"""Export a Figure-1d-style panel: one 3D polytope plot per cancer type.

Hausser Fig. 1d caption: tumors plotted in the space spanned by the first
three gene expression PCs. Archetype vertices and polyhedron edges are drawn
in that same PC1–PC3 space (k=3 archetypes lie in the PC1–PC2 plane).

For each type, uses the k in {3,4,5} with the lowest permutation p-value.

Usage:

    .venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/04_export_figure_1d.py"
"""

from __future__ import annotations

from pathlib import Path
import sys

from _paths import ANALYSIS, HERE

sys.path.insert(0, str(HERE))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _registry import CANCER_TYPES, fig1d_types

RESULTS = ANALYSIS / "results"
FIGURES = ANALYSIS / "figures"

ARC_COLOR = "#E45756"
DATA_COLOR = "#4C78A8"
VIEW_ELEV = 22
VIEW_AZIM = -58


def style():
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 9,
            "figure.facecolor": "white",
            "savefig.dpi": 300,
        }
    )


def p_label(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "p n/a"
    if p == 0:
        return "p < 1/n_perm"
    return f"p = {p:.3f}"


def display_k_from_summary(summary: pd.DataFrame) -> int | None:
    valid = summary.dropna(subset=["p_value"])
    if valid.empty:
        return None
    return int(valid.loc[valid["p_value"].idxmin(), "k"])


def archetypes_pc123(arcs: np.ndarray) -> np.ndarray:
    """First three PC coordinates for 3D display; pad with 0 if k=3 (2-D fit)."""
    out = np.zeros((arcs.shape[0], 3), dtype=float)
    ncol = min(3, arcs.shape[1])
    out[:, :ncol] = arcs[:, :ncol]
    return out


def draw_polyhedron_edges_3d(ax, arcs_3d: np.ndarray, color: str, lw: float = 1.0):
    k = arcs_3d.shape[0]
    for i in range(k):
        for j in range(i + 1, k):
            ax.plot(
                [arcs_3d[i, 0], arcs_3d[j, 0]],
                [arcs_3d[i, 1], arcs_3d[j, 1]],
                [arcs_3d[i, 2], arcs_3d[j, 2]],
                color=color,
                lw=lw,
                alpha=0.9,
            )


def fit_3d_view(ax, points: np.ndarray, pad_frac: float = 0.03):
    """Tight per-PC limits so the cloud fills the panel."""
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


def panel(ax, code):
    info = CANCER_TYPES[code]
    panel_dir = RESULTS / code / "panel_a"
    pc_path = panel_dir / "pc_scores.csv"
    summary_path = panel_dir / "t_ratio_summary.csv"

    if not (pc_path.is_file() and summary_path.is_file()):
        ax.text2D(
            0.5,
            0.5,
            f"{info['label']}\n(no results yet)",
            ha="center",
            va="center",
            fontsize=8,
            color="0.5",
            transform=ax.transAxes,
        )
        ax.set_axis_off()
        return

    scores = pd.read_csv(pc_path, index_col=0).values
    summary = pd.read_csv(summary_path)
    display_k = display_k_from_summary(summary)
    pts = scores[:, :3]

    ax.scatter(
        pts[:, 0],
        pts[:, 1],
        pts[:, 2],
        s=3,
        c=DATA_COLOR,
        alpha=0.45,
        depthshade=True,
        linewidths=0,
    )

    if display_k is not None:
        arcs = pd.read_csv(panel_dir / f"archetypes_k{display_k}.csv", index_col=0).values
        arcs_3d = archetypes_pc123(arcs)
        draw_polyhedron_edges_3d(ax, arcs_3d, ARC_COLOR)
        ax.scatter(
            arcs_3d[:, 0],
            arcs_3d[:, 1],
            arcs_3d[:, 2],
            s=40,
            c=ARC_COLOR,
            depthshade=False,
            edgecolors="white",
            linewidths=0.5,
            zorder=10,
        )
        row = summary[summary["k"] == display_k].iloc[0]
        subtitle = f"k={display_k}, t={row['t_ratio']:.2f}, {p_label(row['p_value'])}"
        plot_pts = np.vstack([pts, arcs_3d])
    else:
        subtitle = "no fitted k"
        plot_pts = pts

    fit_3d_view(ax, plot_pts)
    ax.view_init(elev=VIEW_ELEV, azim=VIEW_AZIM)
    ax.set_axis_off()
    ax.set_title(f"{info['label']}\n{subtitle}", fontsize=9, pad=2)


def main():
    style()
    FIGURES.mkdir(parents=True, exist_ok=True)
    codes = fig1d_types()

    ncols = 4
    nrows = 2
    fig = plt.figure(figsize=(22, 11))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.02, wspace=0.02, hspace=0.06)

    for i, code in enumerate(codes):
        ax = fig.add_subplot(nrows, ncols, i + 1, projection="3d")
        panel(ax, code)

    fig.suptitle(
        "Reproduction attempt: Hausser et al. 2019 Fig. 1d "
        "(tumor transcriptomes fall on per-cancer-type polyhedra)",
        y=0.96,
        fontsize=12,
    )
    out_path = FIGURES / "Figure_1D_hausser_reproduction.png"
    fig.savefig(out_path, pad_inches=0.08)
    plt.close(fig)
    print("Wrote", out_path)


if __name__ == "__main__":
    raise SystemExit(main())
