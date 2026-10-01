"""Path B: bootstrap-averaged PCHA on GSVA + Fig 3A/B + QC vs authors.

Mimics paper ParetoTI bootstrap (subsample 80%, many fits, average vertices)
using project src.archetypes PCHA when ParetoTI is unavailable.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA

ROOT = Path("/Users/apple/Desktop/thesis_iisc/Project_1")
sys.path.insert(0, str(ROOT))
from src.archetypes import fit_pcha, simplex_volume

HERE = ROOT / "Anderson_PML_Archetypes"
OUT = HERE / "results/fig3_v2"
AUTH = HERE / "data/author"

ARCH_NAMES = [
    "A1: Normal-like",
    "A2: Inflammation",
    "A3: Cell adhesion",
    "A4: Proliferation",
]
ARCH_COLORS = {
    "A1: Normal-like": "#416c3e",
    "A2: Inflammation": "#67b3bb",
    "A3: Cell adhesion": "#843385",
    "A4: Proliferation": "#ddb137",
    "Secondary": "#d0d0d0",
}


def archetype_cmap(hex_color: str) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(
        "arch", ["#f5f5f5", "#e8e8e8", hex_color, hex_color], N=256
    )


def bootstrap_pcha(pcs: np.ndarray, k: int = 4, n_boot: int = 500, sample_prop: float = 0.8, seed: int = 42792):
    """Average archetype positions over bootstrap PCHA fits (paper-like)."""
    rng = np.random.default_rng(seed)
    n = pcs.shape[0]
    m = max(k + 1, int(round(n * sample_prop)))
    arcs_list = []
    for b in range(n_boot):
        idx = rng.choice(n, size=m, replace=False)
        sub = pcs[idx, : k - 1]
        best_vol = -np.inf
        best = None
        for _ in range(5):  # a few inits per bootstrap
            try:
                arcs, _, _ = fit_pcha(sub, k=k, delta=0.0, conv_crit=1e-6, maxiter=500)
                vol = simplex_volume(arcs)
                if np.isfinite(vol) and vol > best_vol:
                    best_vol = vol
                    best = arcs
            except Exception:
                continue
        if best is not None:
            arcs_list.append(best)
        if (b + 1) % 50 == 0:
            print(f"  bootstrap {b+1}/{n_boot} kept={len(arcs_list)}")
    if not arcs_list:
        raise RuntimeError("All bootstrap PCHA fits failed")
    # Align vertex order to first fit by nearest-neighbor matching, then average
    ref = arcs_list[0]
    aligned = [ref]
    for arcs in arcs_list[1:]:
        used = set()
        order = []
        for i in range(k):
            d = np.linalg.norm(arcs - ref[i], axis=1)
            for j in np.argsort(d):
                if j not in used:
                    used.add(int(j))
                    order.append(int(j))
                    break
        aligned.append(arcs[order])
    mean_arcs = np.mean(np.stack(aligned, axis=0), axis=0)
    return mean_arcs, len(arcs_list)


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


def map_vertices_to_paper(dists: np.ndarray, author: pd.Series) -> list[str]:
    prop = 0.125
    paper = ARCH_NAMES
    scores = np.zeros((4, 4))
    for a in range(4):
        thr = np.quantile(dists[:, a], prop)
        near = dists[:, a] <= thr
        labs = author[near]
        labs = labs[labs != "Secondary"]
        for j, name in enumerate(paper):
            scores[a, j] = (labs == name).sum()
    pairs = [((i, j), scores[i, j]) for i in range(4) for j in range(4)]
    pairs.sort(key=lambda x: -x[1])
    final = [None] * 4
    used_v, used_l = set(), set()
    for (i, j), sc in pairs:
        if i in used_v or j in used_l:
            continue
        final[i] = paper[j]
        used_v.add(i)
        used_l.add(j)
    for i in range(4):
        if final[i] is None:
            for j, name in enumerate(paper):
                if j not in used_l:
                    final[i] = name
                    used_l.add(j)
                    break
    return final


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(OUT / "gsva_scores.csv", index_col=0)
    info = pd.read_csv(OUT / "sample_info_v2.csv")
    if "sample_id" in info.columns:
        info = info.set_index("sample_id")
    info = info.loc[scores.index]
    author_raw = pd.read_csv(AUTH / "PML_Archetypes.csv").set_index("sample_id")
    author = author_raw["Combined"].astype(str).str.replace("Cell Adhesion", "Cell adhesion")
    author = author.loc[scores.index]

    X = scores.to_numpy(dtype=float)
    X = X - X.mean(axis=0, keepdims=True)
    pca = PCA(n_components=min(6, X.shape[1], X.shape[0] - 1), random_state=42792)
    pcs = pca.fit_transform(X)
    print("PC var:", np.round(pca.explained_variance_ratio_[:3], 3))

    k = 4
    # 500 bootstraps is a practical stand-in for paper's 5000 (still slow);
    # bump if time allows
    n_boot = 500
    print(f"Bootstrap PCHA n={n_boot} ...")
    arcs_k, n_kept = bootstrap_pcha(pcs, k=k, n_boot=n_boot, sample_prop=0.8, seed=42792)
    print(f"kept {n_kept}/{n_boot} bootstrap fits")

    # pad arcs to PC1-2 for plotting (already in first k-1 PCs)
    arcs = np.zeros((k, pcs.shape[1]))
    arcs[:, : k - 1] = arcs_k

    data_k = pcs[:, : k - 1]
    dists = np.linalg.norm(data_k[:, None, :] - arcs_k[None, :, :], axis=2)
    labels_fit = map_vertices_to_paper(dists, author)
    print("vertex map:", labels_fit)
    order = [labels_fit.index(n) for n in ARCH_NAMES]
    dists = dists[:, order]
    arcs = arcs[order]
    arcs_k = arcs_k[order]
    signed = -dists

    primary = assign_octile(dists, 0.125)
    labels = np.array(["Secondary"] * len(scores), dtype=object)
    for i, p in enumerate(primary):
        if p >= 0:
            labels[i] = ARCH_NAMES[p]

    assign = info.copy()
    assign["primary"] = labels
    assign["author_combined"] = author.values
    for i, name in enumerate(ARCH_NAMES):
        assign[f"dist_{name}"] = dists[:, i]
        assign[f"signed_{name}"] = signed[:, i]
    assign["PC1"] = pcs[:, 0]
    assign["PC2"] = pcs[:, 1]
    assign.to_csv(OUT / "sample_assignments.csv")

    # QC
    a = author.values
    b = labels
    both = (a != "Secondary") & (b != "Secondary")
    if both.any():
        agree = (a[both] == b[both]).mean()
        print(f"Primary agreement vs authors: {agree:.1%} (n={both.sum()})")
    print("Ours:\n", pd.Series(labels).value_counts())
    print("Authors:\n", pd.Series(a).value_counts())

    auth_dist_cols = {
        "A1: Normal-like": "Normal_like",
        "A2: Inflammation": "Inflammation",
        "A3: Cell adhesion": "Cell_Adhesion",
        "A4: Proliferation": "Proliferation",
    }
    print("Signed-distance correlation vs authors:")
    for name, acol in auth_dist_cols.items():
        r = np.corrcoef(author_raw.loc[scores.index, acol].astype(float), assign[f"signed_{name}"])[0, 1]
        print(f"  {name}: r={r:.3f}")

    arc_pc = arcs[:, :2]
    vmin = float(np.percentile(signed, 2))
    vmax = float(np.percentile(signed, 98))

    # Fig 3A
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), sharex=True, sharey=True)
    for i, name in enumerate(ARCH_NAMES):
        ax = axes.ravel()[i]
        sc = ax.scatter(
            pcs[:, 0], pcs[:, 1], c=signed[:, i], cmap=archetype_cmap(ARCH_COLORS[name]),
            vmin=vmin, vmax=vmax, s=28, edgecolors="0.35", linewidths=0.25,
        )
        order_cycle = list(range(k)) + [0]
        ax.plot(arc_pc[order_cycle, 0], arc_pc[order_cycle, 1], color="0.25", lw=1.0, zorder=3)
        ax.scatter(arc_pc[:, 0], arc_pc[:, 1], c="0.15", s=40, marker="o", zorder=4)
        for j, lab in enumerate(["1", "2", "3", "4"]):
            ax.text(arc_pc[j, 0], arc_pc[j, 1], lab, fontsize=14, fontweight="bold",
                    color="0.35", ha="center", va="center", zorder=5)
        ax.set_title(name, color=ARCH_COLORS[name], fontweight="bold")
        ax.set_xlabel(f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}% variance)")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Fig 3A Path B: improved public-data rebuild", y=1.01)
    fig.tight_layout()
    fig.savefig(OUT / "fig3A_distance_panels.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig3A_distance_panels.pdf", bbox_inches="tight")
    plt.close(fig)

    # Fig 3B
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    hist_col = "Histology" if "Histology" in assign.columns else "histology"
    for hist, marker in [("AAH", "o"), ("AIS/MIA", "^")]:
        m = (assign[hist_col].to_numpy() == hist) & (assign["primary"].to_numpy() == "Secondary")
        if m.any():
            ax.scatter(pcs[m, 0], pcs[m, 1], c=ARCH_COLORS["Secondary"], marker=marker,
                       s=28, edgecolors="0.55", linewidths=0.3, alpha=0.85,
                       label=f"Secondary | {hist}", zorder=1)
    for hist, marker in [("AAH", "o"), ("AIS/MIA", "^")]:
        m_hist = assign[hist_col].to_numpy() == hist
        for lab in ARCH_NAMES:
            mm = m_hist & (assign["primary"].to_numpy() == lab)
            if not mm.any():
                continue
            ax.scatter(pcs[mm, 0], pcs[mm, 1], c=ARCH_COLORS[lab], marker=marker,
                       s=55, edgecolors="0.2", linewidths=0.35,
                       label=f"{lab} | {hist}", zorder=2)
    order_cycle = list(range(k)) + [0]
    ax.plot(arc_pc[order_cycle, 0], arc_pc[order_cycle, 1], color="0.15", lw=1.8, zorder=3)
    ax.scatter(arc_pc[:, 0], arc_pc[:, 1], c="0.1", s=50, marker="o", zorder=4)
    for j, lab in enumerate(["1", "2", "3", "4"]):
        ax.text(arc_pc[j, 0], arc_pc[j, 1], lab, fontsize=16, fontweight="bold",
                color="0.2", ha="center", va="center", zorder=5)
    ax.set_xlabel(f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}% variance)")
    ax.set_title("Fig 3B Path B: primary (top 12.5%) vs secondary")
    ax.legend(fontsize=7, loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig3B_primary_assignments.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig3B_primary_assignments.pdf", bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(arcs, index=ARCH_NAMES, columns=[f"PC{i+1}" for i in range(arcs.shape[1])]).to_csv(
        OUT / "archetypes_pc.csv"
    )
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
