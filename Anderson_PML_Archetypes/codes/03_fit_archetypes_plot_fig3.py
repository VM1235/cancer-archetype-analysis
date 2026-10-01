"""Fit k=4 PCHA on GSVA scores and export Fig 3A / 3B plots."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA

ROOT = Path("/Users/apple/Desktop/thesis_iisc/Project_1")
sys.path.insert(0, str(ROOT))
from src.archetypes import fit_pcha_best

HERE = ROOT / "Anderson_PML_Archetypes"
PROC = HERE / "data/processed"
OUT = HERE / "results/fig3"

ARCH_NAMES = [
    "A1: Normal-like",
    "A2: Inflammation",
    "A3: Cell adhesion",
    "A4: Proliferation",
]
# Sampled from Anderson et al. Fig. 3 (Mol Cancer Res 2026)
ARCH_COLORS = {
    "A1: Normal-like": "#416c3e",  # forest green
    "A2: Inflammation": "#67b3bb",  # teal
    "A3: Cell adhesion": "#843385",  # purple
    "A4: Proliferation": "#ddb137",  # gold / amber
    "Secondary": "#d0d0d0",
}


def archetype_cmap(hex_color: str) -> LinearSegmentedColormap:
    """Paper Fig 3A: pale → archetype hue (far → close)."""
    return LinearSegmentedColormap.from_list(
        "arch", ["#f5f5f5", "#e8e8e8", hex_color, hex_color], N=256
    )


def assign_octile(distances: np.ndarray, prop: float = 0.125) -> np.ndarray:
    n, k = distances.shape
    primary = np.full(n, -1, dtype=int)
    in_bin = np.zeros((n, k), dtype=bool)
    for a in range(k):
        thr = np.quantile(distances[:, a], prop)
        in_bin[:, a] = distances[:, a] <= thr
    for i in range(n):
        hits = np.where(in_bin[i])[0]
        if hits.size == 1:
            primary[i] = int(hits[0])
    return primary


def map_archetypes_to_paper_labels(arcs_pc: np.ndarray, score_df: pd.DataFrame) -> list[str]:
    """Order fitted vertices to paper names using module biology when possible.

    Fallback: keep fit order as A1..A4.
    """
    # Prefer matching by distance of archetype profile to module axes:
    # Normal-like ~ endothelial proliferation high / cell cycle low
    # Proliferation ~ cell cycle high
    # Inflammation ~ inflammatory / cytokine
    # Cell adhesion ~ ECM
    cols = list(score_df.columns)
    # arcs in PC space only — use mean GSVA of nearest samples later; for now keep order
    return ARCH_NAMES


def main() -> None:
    score_df = pd.read_csv(OUT / "gsva_scores.csv", index_col=0)
    info = pd.read_csv(PROC / "PML_sample_info.csv").set_index("sample_id")
    info = info.loc[score_df.index]

    X = score_df.to_numpy(dtype=float)
    X = X - X.mean(axis=0, keepdims=True)
    pca = PCA(n_components=min(6, X.shape[1], X.shape[0] - 1), random_state=42792)
    pcs = pca.fit_transform(X)
    pd.DataFrame(pcs, index=score_df.index, columns=[f"PC{i+1}" for i in range(pcs.shape[1])]).to_csv(
        OUT / "pc_scores.csv"
    )

    k = 4
    arcs, weights, varexpl, vol, n_ok = fit_pcha_best(
        pcs, k=k, n_init=150, delta=0.0, conv_crit=1e-6, maxiter=500
    )
    print(f"PCHA k={k}: varexpl={varexpl:.3f} vol={vol:.4g} n_ok={n_ok}")
    print("PC variance:", np.round(pca.explained_variance_ratio_[:3], 3))

    data_k = pcs[:, : k - 1]
    arcs_k = arcs[:, : k - 1]
    dists = np.linalg.norm(data_k[:, None, :] - arcs_k[None, :, :], axis=2)
    signed = -dists

    # Label vertices using nearest samples' author labels when available
    labels_fit = ARCH_NAMES.copy()
    if "author_combined" in info.columns:
        author = info["author_combined"].astype(str).str.replace("Cell Adhesion", "Cell adhesion")
        # for each fitted vertex, among closest 12.5%, majority author primary label
        prop = 0.125
        mapped = []
        used = set()
        for a in range(k):
            thr = np.quantile(dists[:, a], prop)
            near = dists[:, a] <= thr
            labs = author[near]
            labs = labs[labs != "Secondary"]
            if len(labs):
                maj = labs.value_counts().index[0]
            else:
                maj = f"A{a+1}"
            mapped.append(maj)
        # resolve collisions by assigning unique paper labels greedily
        paper = ARCH_NAMES.copy()
        final = ["Secondary"] * k
        # score each (vertex, paper_label) by count of author matches in octile
        scores = np.zeros((k, 4))
        for a in range(k):
            thr = np.quantile(dists[:, a], prop)
            near = dists[:, a] <= thr
            labs = author[near]
            for j, name in enumerate(paper):
                scores[a, j] = (labs == name).sum()
        # greedy max matching
        assigned_v = set()
        assigned_l = set()
        pairs = [((i, j), scores[i, j]) for i in range(k) for j in range(4)]
        pairs.sort(key=lambda x: -x[1])
        final = [None] * k
        for (i, j), sc in pairs:
            if i in assigned_v or j in assigned_l:
                continue
            if sc <= 0 and len(assigned_v) < k:
                continue
            final[i] = paper[j]
            assigned_v.add(i)
            assigned_l.add(j)
        for i in range(k):
            if final[i] is None:
                for j, name in enumerate(paper):
                    if j not in assigned_l:
                        final[i] = name
                        assigned_l.add(j)
                        break
        labels_fit = final
        print("Mapped archetype vertices ->", labels_fit)

    # reorder columns by paper order
    order = [labels_fit.index(n) for n in ARCH_NAMES]
    dists = dists[:, order]
    signed = signed[:, order]
    arcs = arcs[order]
    arcs_k = arcs_k[order]

    primary = assign_octile(dists, 0.125)
    labels = np.array(["Secondary"] * len(score_df), dtype=object)
    for i, p in enumerate(primary):
        if p >= 0:
            labels[i] = ARCH_NAMES[p]

    assign = info.copy()
    assign["primary"] = labels
    for i, name in enumerate(ARCH_NAMES):
        assign[f"dist_{name}"] = dists[:, i]
        assign[f"signed_{name}"] = signed[:, i]
    assign["PC1"] = pcs[:, 0]
    assign["PC2"] = pcs[:, 1]
    assign.to_csv(OUT / "sample_assignments.csv")

    if "author_combined" in assign.columns:
        a = assign["author_combined"].astype(str).str.replace("Cell Adhesion", "Cell adhesion")
        b = assign["primary"].astype(str)
        both = (a != "Secondary") & (b != "Secondary")
        if both.any():
            print(f"Primary agreement vs authors: {(a[both]==b[both]).mean():.1%} (n={both.sum()})")
        print("Ours:\n", pd.Series(labels).value_counts())
        print("Authors:\n", a.value_counts())

    arc_pc = arcs[:, :2]
    # shared color limits for 3A (signed -distance; higher/less-negative = closer)
    vmin = float(np.percentile(signed, 2))
    vmax = float(np.percentile(signed, 98))

    # Fig 3A — each panel uses that archetype's monochrome scale (paper style)
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), sharex=True, sharey=True)
    for i, name in enumerate(ARCH_NAMES):
        ax = axes.ravel()[i]
        cmap = archetype_cmap(ARCH_COLORS[name])
        sc = ax.scatter(
            pcs[:, 0],
            pcs[:, 1],
            c=signed[:, i],
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            s=28,
            edgecolors="0.35",
            linewidths=0.25,
        )
        # polytope outline
        order_cycle = list(range(k)) + [0]
        ax.plot(
            arc_pc[order_cycle, 0],
            arc_pc[order_cycle, 1],
            color="0.25",
            lw=1.0,
            zorder=3,
        )
        ax.scatter(arc_pc[:, 0], arc_pc[:, 1], c="0.15", s=40, marker="o", zorder=4)
        for j, lab in enumerate(["1", "2", "3", "4"]):
            ax.text(
                arc_pc[j, 0],
                arc_pc[j, 1],
                lab,
                fontsize=14,
                fontweight="bold",
                color="0.35",
                ha="center",
                va="center",
                zorder=5,
            )
        ax.set_title(name, color=ARCH_COLORS[name], fontweight="bold")
        ax.set_xlabel(f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}% variance)")
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=7)
    fig.suptitle("Fig 3A-style: distance to each PML archetype (rebuild)", y=1.01)
    fig.tight_layout()
    fig.savefig(OUT / "fig3A_distance_panels.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig3A_distance_panels.pdf", bbox_inches="tight")
    plt.close(fig)

    # Fig 3B — categorical paper colors; AAH=circle, AIS/MIA=triangle
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    # draw secondary first (underneath)
    for hist, marker in [("AAH", "o"), ("AIS/MIA", "^")]:
        m = (assign["Histology"].to_numpy() == hist) & (
            assign["primary"].to_numpy() == "Secondary"
        )
        if m.any():
            ax.scatter(
                pcs[m, 0],
                pcs[m, 1],
                c=ARCH_COLORS["Secondary"],
                marker=marker,
                s=28,
                edgecolors="0.55",
                linewidths=0.3,
                alpha=0.85,
                label=f"Secondary | {hist}",
                zorder=1,
            )
    for hist, marker in [("AAH", "o"), ("AIS/MIA", "^")]:
        m_hist = assign["Histology"].to_numpy() == hist
        for lab in ARCH_NAMES:
            mm = m_hist & (assign["primary"].to_numpy() == lab)
            if not mm.any():
                continue
            ax.scatter(
                pcs[mm, 0],
                pcs[mm, 1],
                c=ARCH_COLORS[lab],
                marker=marker,
                s=55,
                edgecolors="0.2",
                linewidths=0.35,
                label=f"{lab} | {hist}",
                zorder=2,
            )
    order_cycle = list(range(k)) + [0]
    ax.plot(
        arc_pc[order_cycle, 0],
        arc_pc[order_cycle, 1],
        color="0.15",
        lw=1.8,
        zorder=3,
    )
    ax.scatter(arc_pc[:, 0], arc_pc[:, 1], c="0.1", s=50, marker="o", zorder=4)
    for j, lab in enumerate(["1", "2", "3", "4"]):
        ax.text(
            arc_pc[j, 0],
            arc_pc[j, 1],
            lab,
            fontsize=16,
            fontweight="bold",
            color="0.2",
            ha="center",
            va="center",
            zorder=5,
        )
    ax.set_xlabel(f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}% variance)")
    ax.set_title("Fig 3B-style: primary archetype (top 12.5%) vs secondary")
    ax.legend(fontsize=7, loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig3B_primary_assignments.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig3B_primary_assignments.pdf", bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(
        arcs, index=ARCH_NAMES, columns=[f"PC{i+1}" for i in range(arcs.shape[1])]
    ).to_csv(OUT / "archetypes_pc.csv")
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
