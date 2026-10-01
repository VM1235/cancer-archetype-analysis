"""Re-fit Path B correctly: PCHA + distances in 9-module GSVA space (paper Step4).

Bug in codes/08: fitted PCHA on first (k-1) PCs and used PC-space distances.
Author XC is 9 modules × 4 archetypes; merge_arch_dist is Euclidean in GSVA space.
PCA is for plotting only.
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
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA

ROOT = Path("/Users/apple/Desktop/thesis_iisc/Project_1")
sys.path.insert(0, str(ROOT))
from src.archetypes import fit_pcha, simplex_volume

HERE = ROOT / "Anderson_PML_Archetypes"
V2 = HERE / "results/fig3_v2"
OUT = HERE / "results/fig3_v2_refit9d"
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


def bootstrap_pcha_gsva(
    X: np.ndarray, k: int = 4, n_boot: int = 500, sample_prop: float = 0.8, seed: int = 42792
):
    """Bootstrap-average archetypes in full module space (n_samples × n_modules)."""
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    m = max(k + 1, int(round(n * sample_prop)))
    arcs_list = []
    for b in range(n_boot):
        idx = rng.choice(n, size=m, replace=False)
        sub = X[idx]
        best_vol = -np.inf
        best = None
        for _ in range(5):
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
            print(f"  bootstrap {b + 1}/{n_boot} kept={len(arcs_list)}")
    if not arcs_list:
        raise RuntimeError("All bootstrap PCHA fits failed")
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
    return np.mean(np.stack(aligned, axis=0), axis=0), len(arcs_list)


def match_to_author_xc(arcs: np.ndarray, xc_paper: np.ndarray) -> list[str]:
    """Map each fitted vertex -> paper name by nearest author XC column.

    xc_paper columns already ordered as ARCH_NAMES.
    """
    k = arcs.shape[0]
    scores = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            scores[i, j] = -np.linalg.norm(arcs[i] - xc_paper[:, j])
    pairs = [((i, j), scores[i, j]) for i in range(k) for j in range(k)]
    pairs.sort(key=lambda x: -x[1])
    final = [None] * k
    used_v, used_l = set(), set()
    for (i, j), _ in pairs:
        if i in used_v or j in used_l:
            continue
        final[i] = ARCH_NAMES[j]
        used_v.add(i)
        used_l.add(j)
    for i in range(k):
        if final[i] is None:
            for j, name in enumerate(ARCH_NAMES):
                if j not in used_l:
                    final[i] = name
                    used_l.add(j)
                    break
    return final


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


def align_pc_signs(pcs: np.ndarray, arc_pc: np.ndarray, target_arc_pc: np.ndarray):
    """Flip PC1/PC2 so our projected arcs best match author arcs in PC space."""
    best = pcs.copy()
    best_arc = arc_pc.copy()
    best_err = np.inf
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            trial_arc = arc_pc.copy()
            trial_arc[:, 0] *= s1
            trial_arc[:, 1] *= s2
            # match vertex order already ARCH_NAMES
            err = np.sum((trial_arc[:, :2] - target_arc_pc[:, :2]) ** 2)
            if err < best_err:
                best_err = err
                best = pcs.copy()
                best[:, 0] *= s1
                best[:, 1] *= s2
                best_arc = trial_arc
    return best, best_arc, best_err


def hull_edges(points_2d: np.ndarray):
    if len(points_2d) < 3:
        return []
    hull = ConvexHull(points_2d)
    return list(hull.simplices)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(V2 / "gsva_scores.csv", index_col=0)
    xc = pd.read_csv(AUTH / "PML_arcA4_XC.csv", index_col=0)
    scores = scores.loc[:, xc.index]
    info = pd.read_csv(V2 / "sample_info_v2.csv").set_index("sample_id").loc[scores.index]
    author_raw = pd.read_csv(AUTH / "PML_Archetypes.csv").set_index("sample_id")
    author = (
        author_raw["Combined"]
        .astype(str)
        .str.replace("Cell Adhesion", "Cell adhesion")
        .loc[scores.index]
    )

    # Author XC in ARCH_NAMES column order (Step4 map on raw ParetoTI cols)
    # Raw XC cols = Arch1..4 → A1, A4, A2, A3
    xc_raw = xc.to_numpy(dtype=float)
    paper_order_from_raw = [
        "A1: Normal-like",
        "A4: Proliferation",
        "A2: Inflammation",
        "A3: Cell adhesion",
    ]
    xc_paper = np.column_stack(
        [xc_raw[:, paper_order_from_raw.index(n)] for n in ARCH_NAMES]
    )  # modules × 4 in ARCH_NAMES order

    X = scores.to_numpy(dtype=float)
    print(f"GSVA matrix: {X.shape[0]} samples × {X.shape[1]} modules")
    print("Bootstrap PCHA in full GSVA space (not PC subspace) ...")
    arcs, n_kept = bootstrap_pcha_gsva(X, k=4, n_boot=500, sample_prop=0.8, seed=42792)
    print(f"kept {n_kept}/500; arcs shape {arcs.shape}")

    labels_fit = match_to_author_xc(arcs, xc_paper)
    print("vertex map (nearest author XC):", labels_fit)
    order = [labels_fit.index(n) for n in ARCH_NAMES]
    arcs = arcs[order]

    # distances in 9D GSVA space (paper-like)
    dists = np.linalg.norm(X[:, None, :] - arcs[None, :, :], axis=2)
    signed = -dists
    primary = assign_octile(dists, 0.125)
    labels = np.array(["Secondary"] * len(scores), dtype=object)
    for i, p in enumerate(primary):
        if p >= 0:
            labels[i] = ARCH_NAMES[p]

    # Compare fitted XC to author XC
    print("Euclidean distance fitted vs author XC (same name):")
    for i, name in enumerate(ARCH_NAMES):
        d = np.linalg.norm(arcs[i] - xc_paper[:, i])
        print(f"  {name}: {d:.3f}")

    # PCA for plotting only
    Xc = X - X.mean(axis=0, keepdims=True)
    pca = PCA(n_components=min(6, X.shape[1], X.shape[0] - 1), random_state=42792)
    pcs = pca.fit_transform(Xc)
    mean = X.mean(axis=0)
    arc_pc = (arcs - mean) @ pca.components_.T
    author_arc_pc = (xc_paper.T - mean) @ pca.components_.T
    print("PC var:", np.round(pca.explained_variance_ratio_[:3], 3))

    pcs, arc_pc, flip_err = align_pc_signs(pcs, arc_pc, author_arc_pc)
    # also flip author target for reporting
    author_arc_pc2 = author_arc_pc.copy()
    # recompute which signs were applied by matching
    print(f"PC sign alignment residual vs author arcs: {flip_err:.3f}")
    print("Our arcs PC1-2 after sign align:")
    for name, row in zip(ARCH_NAMES, arc_pc):
        print(f"  {name}: ({row[0]:.3f}, {row[1]:.3f})")

    assign = info.copy()
    assign["primary"] = labels
    assign["author_combined"] = author.values
    for i, name in enumerate(ARCH_NAMES):
        assign[f"dist_{name}"] = dists[:, i]
        assign[f"signed_{name}"] = signed[:, i]
    assign["PC1"] = pcs[:, 0]
    assign["PC2"] = pcs[:, 1]
    assign.to_csv(OUT / "sample_assignments.csv")
    pd.DataFrame(arcs, index=ARCH_NAMES, columns=list(scores.columns)).to_csv(OUT / "archetypes_gsva.csv")
    pd.DataFrame(arc_pc, index=ARCH_NAMES, columns=[f"PC{i+1}" for i in range(arc_pc.shape[1])]).to_csv(
        OUT / "archetypes_pc.csv"
    )

    a = author.values
    b = labels
    both = (a != "Secondary") & (b != "Secondary")
    if both.any():
        print(f"Primary agreement vs authors: {(a[both] == b[both]).mean():.1%} (n={int(both.sum())})")
    print("Ours:\n", pd.Series(labels).value_counts())
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

    arc2 = arc_pc[:, :2]
    k = 4
    vmin, vmax = float(np.percentile(signed, 2)), float(np.percentile(signed, 98))
    edges = hull_edges(arc2)

    fig, axes = plt.subplots(2, 2, figsize=(10, 9), sharex=True, sharey=True)
    for i, name in enumerate(ARCH_NAMES):
        ax = axes.ravel()[i]
        sc = ax.scatter(
            pcs[:, 0], pcs[:, 1], c=signed[:, i], cmap=archetype_cmap(ARCH_COLORS[name]),
            vmin=vmin, vmax=vmax, s=28, edgecolors="0.35", linewidths=0.25,
        )
        for e0, e1 in edges:
            ax.plot(arc2[[e0, e1], 0], arc2[[e0, e1], 1], color="0.25", lw=1.0, zorder=3)
        ax.scatter(arc2[:, 0], arc2[:, 1], c="0.15", s=40, marker="o", zorder=4)
        for j, lab in enumerate(["1", "2", "3", "4"]):
            ax.text(arc2[j, 0], arc2[j, 1], lab, fontsize=14, fontweight="bold",
                    color="0.35", ha="center", va="center", zorder=5)
        ax.set_title(name, color=ARCH_COLORS[name], fontweight="bold")
        ax.set_xlabel(f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}% variance)")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Fig 3A: re-fit PCHA in 9D GSVA (fixed vs codes/08)", y=1.01)
    fig.tight_layout()
    fig.savefig(OUT / "fig3A_distance_panels.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig3A_distance_panels.pdf", bbox_inches="tight")
    plt.close(fig)

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
    for e0, e1 in edges:
        ax.plot(arc2[[e0, e1], 0], arc2[[e0, e1], 1], color="0.15", lw=1.8, zorder=3)
    ax.scatter(arc2[:, 0], arc2[:, 1], c="0.1", s=50, marker="o", zorder=4)
    for j, lab in enumerate(["1", "2", "3", "4"]):
        ax.text(arc2[j, 0], arc2[j, 1], lab, fontsize=16, fontweight="bold",
                color="0.2", ha="center", va="center", zorder=5)
    ax.set_xlabel(f"PC1 ({100*pca.explained_variance_ratio_[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({100*pca.explained_variance_ratio_[1]:.1f}% variance)")
    ax.set_title("Fig 3B: re-fit in 9D GSVA + PC sign align to author XC")
    ax.legend(fontsize=7, loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "fig3B_primary_assignments.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "fig3B_primary_assignments.pdf", bbox_inches="tight")
    plt.close(fig)

    (OUT / "method.txt").write_text(
        "fix vs codes/08:\n"
        "- PCHA bootstrap on full 9-module GSVA (not first 3 PCs)\n"
        "- distances Euclidean in 9D GSVA\n"
        "- vertex names by nearest author XC column\n"
        "- PCA for plot only; PC1/PC2 signs flipped to match author XC projection\n"
    )
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
