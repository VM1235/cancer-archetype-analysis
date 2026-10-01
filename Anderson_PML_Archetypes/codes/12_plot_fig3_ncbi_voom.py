"""Fig 3A/3B from NCBI-count voom GSVA, distances to author-fixed XC (no re-fit)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from sklearn.decomposition import PCA

HERE = Path("/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes")
IN = HERE / "results/fig3_ncbi_voom"
OUT = HERE / "results/fig3_ncbi_voom"
AUTH = HERE / "data/author"

ARCH_NAMES = [
    "A1: Normal-like",
    "A2: Inflammation",
    "A3: Cell adhesion",
    "A4: Proliferation",
]
XC_TO_PAPER = [
    "A1: Normal-like",
    "A4: Proliferation",
    "A2: Inflammation",
    "A3: Cell adhesion",
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


def plot_panels(pcs, signed, arc_pc, pca, labels, assign, out_dir: Path, tag: str, flip_note: str):
    k = 4
    vmin = float(np.percentile(signed, 2))
    vmax = float(np.percentile(signed, 98))
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), sharex=True, sharey=True)
    for i, name in enumerate(ARCH_NAMES):
        ax = axes.ravel()[i]
        sc = ax.scatter(
            pcs[:, 0],
            pcs[:, 1],
            c=signed[:, i],
            cmap=archetype_cmap(ARCH_COLORS[name]),
            vmin=vmin,
            vmax=vmax,
            s=28,
            edgecolors="0.35",
            linewidths=0.25,
        )
        order_cycle = list(range(k)) + [0]
        ax.plot(arc_pc[order_cycle, 0], arc_pc[order_cycle, 1], color="0.25", lw=1.0, zorder=3)
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
        ax.set_xlabel(f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}% variance)")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(f"Fig 3A: NCBI-count voom GSVA → author XC{flip_note}", y=1.01)
    fig.tight_layout()
    fig.savefig(out_dir / f"fig3A_distance_panels{tag}.png", dpi=200, bbox_inches="tight")
    fig.savefig(out_dir / f"fig3A_distance_panels{tag}.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    hist_col = "Histology" if "Histology" in assign.columns else "histology"
    for hist, marker in [("AAH", "o"), ("AIS/MIA", "^")]:
        m = (assign[hist_col].to_numpy() == hist) & (assign["primary"].to_numpy() == "Secondary")
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
        m_hist = assign[hist_col].to_numpy() == hist
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
    ax.plot(arc_pc[order_cycle, 0], arc_pc[order_cycle, 1], color="0.15", lw=1.8, zorder=3)
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
    ax.set_xlabel(f"PC1 ({100 * pca.explained_variance_ratio_[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({100 * pca.explained_variance_ratio_[1]:.1f}% variance)")
    ax.set_title(f"Fig 3B: NCBI-count voom + author-fixed archetypes{flip_note}")
    ax.legend(fontsize=7, loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / f"fig3B_primary_assignments{tag}.png", dpi=200, bbox_inches="tight")
    fig.savefig(out_dir / f"fig3B_primary_assignments{tag}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(IN / "gsva_scores.csv", index_col=0)
    info = pd.read_csv(IN / "sample_info.csv").set_index("sample_id").loc[scores.index]
    author_raw = pd.read_csv(AUTH / "PML_Archetypes.csv").set_index("sample_id")
    author = (
        author_raw["Combined"]
        .astype(str)
        .str.replace("Cell Adhesion", "Cell adhesion")
        .loc[scores.index]
    )

    xc = pd.read_csv(AUTH / "PML_arcA4_XC.csv", index_col=0)
    missing = [m for m in xc.index if m not in scores.columns]
    if missing:
        raise SystemExit(f"GSVA missing modules for XC: {missing}")
    scores = scores.loc[:, xc.index]
    X = scores.to_numpy(dtype=float)
    XC = xc.to_numpy(dtype=float)

    dists_xc = np.linalg.norm(X[:, None, :] - XC.T[None, :, :], axis=2)
    order = [XC_TO_PAPER.index(n) for n in ARCH_NAMES]
    dists = dists_xc[:, order]
    XC_ord = XC[:, order]
    signed = -dists

    primary = assign_octile(dists, 0.125)
    labels = np.array(["Secondary"] * len(scores), dtype=object)
    for i, p in enumerate(primary):
        if p >= 0:
            labels[i] = ARCH_NAMES[p]

    Xc = X - X.mean(axis=0, keepdims=True)
    pca = PCA(n_components=min(6, X.shape[1], X.shape[0] - 1), random_state=42792)
    pcs = pca.fit_transform(Xc)
    arc_pc_full = (XC_ord.T - X.mean(axis=0)) @ pca.components_.T
    print("PC var:", np.round(pca.explained_variance_ratio_[:3], 3))
    print("Author XC projected PC1-2:")
    for name, row in zip(ARCH_NAMES, arc_pc_full):
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

    a = author.values
    b = labels
    both = (a != "Secondary") & (b != "Secondary")
    if both.any():
        agree = (a[both] == b[both]).mean()
        print(f"Primary agreement vs authors: {agree:.1%} (n={int(both.sum())})")
    print("Ours:\n", pd.Series(labels).value_counts())
    print("Authors:\n", pd.Series(a).value_counts())
    print("Confusion (rows=author, cols=ours) both-primary:")
    print(pd.crosstab(pd.Series(a[both], name="author"), pd.Series(b[both], name="ours")))

    auth_dist_cols = {
        "A1: Normal-like": "Normal_like",
        "A2: Inflammation": "Inflammation",
        "A3: Cell adhesion": "Cell_Adhesion",
        "A4: Proliferation": "Proliferation",
    }
    print("Signed-distance correlation vs authors (9D to fixed XC):")
    for name, acol in auth_dist_cols.items():
        r = np.corrcoef(
            author_raw.loc[scores.index, acol].astype(float), assign[f"signed_{name}"]
        )[0, 1]
        print(f"  {name}: r={r:.3f}")

    log_path = OUT / "plot_log.txt"
    # capture key QC into the log via this print stream; also write a summary file at end

    arc_pc = arc_pc_full[:, :2]
    plot_panels(pcs, signed, arc_pc, pca, labels, assign, OUT, "", "")

    # Cosmetic PC1 flip if A1 is on the left (paper has A1 on the right)
    a1_pc1 = float(arc_pc_full[0, 0])
    if a1_pc1 < 0:
        pcs_f = pcs.copy()
        pcs_f[:, 0] *= -1
        arc_f = arc_pc.copy()
        arc_f[:, 0] *= -1
        assign_f = assign.copy()
        assign_f["PC1"] = pcs_f[:, 0]
        assign_f.to_csv(OUT / "sample_assignments_pc1flip.csv")
        plot_panels(
            pcs_f,
            signed,
            arc_f,
            pca,
            labels,
            assign_f,
            OUT,
            "_pc1flip",
            " (PC1 flipped)",
        )
        print("Wrote PC1-flipped panels (A1 was on the left)")

    pd.DataFrame(
        arc_pc_full, index=ARCH_NAMES, columns=[f"PC{i + 1}" for i in range(arc_pc_full.shape[1])]
    ).to_csv(OUT / "archetypes_pc.csv")
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
