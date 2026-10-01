#!/usr/bin/env python3
"""Export a Figure-1d-style panel: one small polytope plot per cancer type,
tumors in PC1 x PC2 with the archetype vertices for that type's smallest
significant k overlaid, annotated with our t-ratio/p vs. Hausser's reported p.

Run after 03_run_all_cancer_types.py (or after running 02_... for individual
types) -- picks up whichever results/<CODE>/panel_a/ folders exist.

Usage (from repository root):

    .venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/04_export_figure_1d.py"
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
ROOT = ANALYSIS.parent
sys.path.insert(0, str(HERE))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from _registry import CANCER_TYPES, fig1d_types

RESULTS = ANALYSIS / "results"
FIGURES = ANALYSIS / "figures"

ARC_COLOR = "#E45756"
DATA_COLOR = "#4C78A8"


def style():
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.bbox": "tight",
            "savefig.dpi": 240,
        }
    )


def p_label(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "p n/a"
    if p == 0:
        return "p < 1/n_perm"
    return f"p = {p:.3f}"


def panel(ax, code):
    info = CANCER_TYPES[code]
    panel_dir = RESULTS / code / "panel_a"
    best_k_path = panel_dir / "best_k.json"
    pc_path = panel_dir / "pc_scores.csv"
    summary_path = panel_dir / "t_ratio_summary.csv"

    if not (best_k_path.is_file() and pc_path.is_file() and summary_path.is_file()):
        ax.text(0.5, 0.5, f"{info['label']}\n(no results yet)", ha="center", va="center",
                fontsize=8, color="0.5", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        return

    best_k = json.loads(best_k_path.read_text())["best_k"]
    scores = pd.read_csv(pc_path, index_col=0).values
    summary = pd.read_csv(summary_path)

    ax.scatter(scores[:, 0], scores[:, 1], s=8, c=DATA_COLOR, alpha=0.5, linewidths=0)

    if best_k is not None:
        arcs = pd.read_csv(panel_dir / f"archetypes_k{best_k}.csv", index_col=0).values
        ax.scatter(arcs[:, 0], arcs[:, 1], s=60, c=ARC_COLOR, zorder=3, edgecolors="white", linewidths=0.6)
        k = arcs.shape[0]
        for i in range(k):
            for j in range(i + 1, k):
                ax.plot([arcs[i, 0], arcs[j, 0]], [arcs[i, 1], arcs[j, 1]],
                        color=ARC_COLOR, lw=0.9, zorder=2, alpha=0.8)
        row = summary[summary["k"] == best_k].iloc[0]
        subtitle = f"k={best_k}, t={row['t_ratio']:.2f}, {p_label(row['p_value'])}"
    else:
        subtitle = "no significant k found"

    ax.set_title(f"{info['label']}\n{subtitle}\n(Hausser: {info.get('hausser_p', 'n/a')})", fontsize=8)
    ax.set_xlabel("PC1", fontsize=7)
    ax.set_ylabel("PC2", fontsize=7)
    ax.tick_params(labelsize=6)


def main():
    style()
    FIGURES.mkdir(parents=True, exist_ok=True)
    codes = fig1d_types()

    n = len(codes)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.1 * ncols, 3.1 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, code in zip(axes, codes):
        panel(ax, code)
    for ax in axes[len(codes):]:
        ax.axis("off")

    fig.suptitle(
        "Reproduction attempt: Hausser et al. 2019 Fig. 1d "
        "(tumor transcriptomes fall on per-cancer-type polyhedra)",
        y=1.02,
        fontsize=11,
    )
    fig.tight_layout()
    out_path = FIGURES / "Figure_1D_hausser_reproduction.png"
    fig.savefig(out_path)
    plt.close(fig)
    print("Wrote", out_path)


if __name__ == "__main__":
    raise SystemExit(main())
