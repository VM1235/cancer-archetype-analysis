#!/usr/bin/env python3
"""GBM Wang2017-restricted figures using k=6 (smallest t-ratio p), not k=5.

Does not overwrite panel_a_genelist/suggested_k.txt (still DimensionFinder k=5)
or Figure_1*_gbm_genelist.png. Writes *_k6.png and panel_{b,c}_genelist_k6/.
k=6 p=0.066 is the smallest p in that table; it is still not p<0.05.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
GBM = HERE.parent
ROOT = GBM.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.utils import shuffle

from src.enrichment import distance_bins, hypergeometric_enrichment
from src.gene_lists import rename_tumor_genes_to_cell_line
from src.io import load_expression_csv
from src.pca import align_pca_signs, cumulative_variance, fit_pca, inverse_transform_scores

K = 6
PANEL_A = GBM / "results" / "panel_a_genelist"
MATRIX = GBM / "data" / "processed" / "input_panelA_wang2017_genelist.csv"
LABELS_SRC = GBM / "results" / "panel_b" / "verhaak_labels_panelA.csv"
OUT_B = GBM / "results" / "panel_b_genelist_k6"
OUT_C = GBM / "results" / "panel_c_genelist_k6"
COMBAT_SRC = GBM / "results" / "panel_c_genelist" / "CCLE_TCGA_COMBAT.csv"
TCGA_EXPR = GBM / "results" / "panel_c_tcga" / "tcga_primary_log_shared_genes.csv"
TCGA_META = GBM / "results" / "panel_c_tcga" / "tcga_primary_metadata.csv"
FIG = GBM / "figures"

KEEP = ("Classical", "Mesenchymal", "Proneural")
N_BINS = 5
FDR = 0.1
N_PCS_VAR = 20
SEED = 0
N_SHUFFLE = 20

SUBTYPE_COLOR = {
    "Classical": "#3B6FA0",
    "Mesenchymal": "#E07A3D",
    "Proneural": "#54A24B",
    "Neural": "#D989B5",
}
ARC_COLORS = ["#4C78A8", "#F58518", "#E45756", "#72B7B2", "#54A24B", "#B279A2"]


def p_label(p):
    if not np.isfinite(p):
        return "p = NA"
    if p == 0:
        return "p < 0.002"
    if p < 0.01:
        return f"p = {p:.3f}"
    return f"p = {p:.2f}"


def barycentric_weights(points, vertices):
    k = vertices.shape[0]
    a = np.vstack([vertices.T, np.ones((1, k))])
    b = np.vstack([points.T, np.ones((1, points.shape[0]))])
    w, *_ = np.linalg.lstsq(a, b, rcond=None)
    return w.T


def cumulative(values):
    return np.cumsum(np.asarray(values, dtype=float))


def tumor_variance_on_pcs(tumor_X, pca, total_var):
    scores = pca.transform(tumor_X)
    return 100.0 * scores.var(axis=0, ddof=1) / total_var


def draw_panel_a():
    expr = load_expression_csv(MATRIX)
    X = expr.T.values
    n_samples = X.shape[0]
    n_pcs = int((PANEL_A / "n_pcs.txt").read_text().strip().splitlines()[0])
    pca, scores = fit_pca(X, n_components=n_pcs)
    saved = np.load(PANEL_A / "pc_scores.npy")
    pca, scores, _ = align_pca_signs(scores, saved, pca)
    esv = pd.read_csv(PANEL_A / "esv_curve.csv")
    tab = pd.read_csv(PANEL_A / "t_ratio_parti_500.csv")
    tot_esv = esv["esv"].values * float(cumulative_variance(pca)[-1])
    gene_esv = 100.0 * tot_esv
    k_esv = esv["k"].values
    delta_esv = np.diff(np.concatenate([[0.0], gene_esv]))

    pca2, scores2 = fit_pca(X, n_components=2)
    arcs = np.load(PANEL_A / f"archetypes_k{K}_parti.npy")
    arcs_full = np.zeros((arcs.shape[0], n_pcs))
    arcs_full[:, : arcs.shape[1]] = arcs
    arcs2 = pca2.transform(inverse_transform_scores(pca, arcs_full))

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
    fig = plt.figure(figsize=(6.6, 7.4))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 1.55], hspace=0.32, wspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(k_esv, delta_esv, "-o", color="#3B6FA0", ms=5, lw=1.4)
    ax.axvline(K, color="0.35", ls="--", lw=1)
    ax.annotate(
        "k=6 (smallest p, still NS)",
        xy=(K, float(delta_esv[k_esv == K][0])),
        xytext=(8.0, max(float(delta_esv.max()) * 0.72, 4)),
        fontsize=7.5,
        color="0.25",
        arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8),
    )
    ax.set_xlim(1.5, 15.5)
    ax.set_ylim(0, max(12, float(delta_esv.max()) * 1.08))
    ax.set_xticks(range(2, 16))
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("% ESV on top of N−1 model")
    ax.set_title("Explained sample variance (ESV)")

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(tab["k"], tab["t_ratio"], "-o", color="#3B6FA0", ms=6, lw=1.4)
    ax.axvline(K, color="0.35", ls="--", lw=1)
    ymax = max(0.62, float(tab["t_ratio"].max()) * 1.15)
    ax.set_ylim(0, ymax)
    ax.set_xlim(2.5, 7.5)
    ax.set_xticks(tab["k"].astype(int))
    for _, row in tab.iterrows():
        ax.annotate(
            p_label(row["p_value"]),
            (row["k"], row["t_ratio"]),
            textcoords="offset points",
            xytext=(8, 6),
            fontsize=8,
            color="#222",
        )
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("t-ratio")
    ax.set_title("t-ratio (500 shuffles, PCHA numIter=50)")

    ax = fig.add_subplot(gs[1, :])
    ax.scatter(scores2[:, 0], scores2[:, 1], s=18, c="#B0B0B0", zorder=1, linewidths=0)
    for i in range(arcs2.shape[0]):
        for j in range(i + 1, arcs2.shape[0]):
            ax.plot(
                [arcs2[i, 0], arcs2[j, 0]],
                [arcs2[i, 1], arcs2[j, 1]],
                color="#888888",
                lw=0.9,
                zorder=2,
            )
    ax.scatter(
        arcs2[:, 0], arcs2[:, 1], s=90, c="#3B6FA0", edgecolors="k", linewidths=0.4, zorder=3
    )
    for i in range(arcs2.shape[0]):
        ax.annotate(
            str(i + 1),
            (arcs2[i, 0], arcs2[i, 1]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=8,
            fontweight="bold",
        )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"Archetype space ({n_samples} lines, Wang 2017 genes, k={K})")

    fig.suptitle(
        f"Figure 1A  —  GBM gene-list, k={K} (smallest p=0.066, not <0.05)",
        y=0.995,
        fontsize=11,
    )
    path = FIG / "Figure_1A_gbm_genelist_k6.png"
    fig.savefig(path)
    plt.close(fig)
    print("Wrote", path)


def run_panel_b():
    OUT_B.mkdir(parents=True, exist_ok=True)
    expr = load_expression_csv(MATRIX)
    sample_ids = list(expr.columns.astype(str))
    scores = np.load(PANEL_A / "pc_scores.npy")
    arcs = np.load(PANEL_A / f"archetypes_k{K}_parti.npy")
    n_vol = arcs.shape[1]
    labels = pd.read_csv(LABELS_SRC)
    labels["cell_line"] = labels["cell_line"].astype(str)
    labels = labels.set_index("cell_line").reindex(sample_ids)
    keep_mask = labels["verhaak_subtype"].isin(KEEP)
    subtypes = labels.loc[keep_mask, "verhaak_subtype"]
    scores_u = scores[keep_mask.values][:, :n_vol]
    bin_ids, distances = distance_bins(scores_u, arcs, n_bins=N_BINS)
    n_arcs = arcs.shape[0]
    pd.DataFrame(
        distances, index=subtypes.index, columns=[f"arc{i+1}" for i in range(n_arcs)]
    ).to_csv(OUT_B / "distances_verhaak.csv")
    pd.DataFrame(
        bin_ids, index=subtypes.index, columns=[f"arc{i+1}" for i in range(n_arcs)]
    ).to_csv(OUT_B / "bins_verhaak.csv")
    table = hypergeometric_enrichment(
        subtypes, bin_ids, fdr=FDR, subtype_levels=KEEP
    )
    table.to_csv(OUT_B / "enrichment_verhaak.csv", index=False)

    fig, axes = plt.subplots(1, len(KEEP), figsize=(10.2, 3.4), sharey=True)
    for ax, subtype in zip(np.atleast_1d(axes), KEEP):
        sub = table[table["subtype"] == subtype]
        for arc in sorted(sub["archetype"].unique()):
            arc_tab = sub[sub["archetype"] == arc].sort_values("bin")
            color = ARC_COLORS[int(arc) % len(ARC_COLORS)]
            ax.plot(
                arc_tab["bin"],
                arc_tab["fold_enrichment"],
                marker="o",
                color=color,
                label=f"Arc {int(arc) + 1}",
            )
            sig = arc_tab[arc_tab["sig_peak_at_bin0"]]
            if len(sig):
                ax.scatter(
                    sig["bin"],
                    sig["fold_enrichment"],
                    s=70,
                    facecolors="none",
                    edgecolors=color,
                    linewidths=1.6,
                    zorder=3,
                )
        ax.axhline(1.0, color="0.7", lw=1, ls="--")
        ax.set_title(subtype)
        ax.set_xlabel("distance bin")
        ax.set_xticks(range(N_BINS))
    axes[0].set_ylabel("fold enrichment")
    axes[-1].legend(loc="upper right", fontsize=6, frameon=False)
    fig.suptitle(
        f"Figure 1B  —  GBM gene-list k={K} marker enrichment",
        y=1.06,
        fontsize=13,
    )
    fig.tight_layout()
    path = FIG / "Figure_1B_gbm_genelist_k6.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(OUT_B / "enrichment_verhaak.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    hits = table[table["sig_peak_at_bin0"]]
    print("Wrote", path)
    print("bin-0 peaks:")
    if hits.empty:
        print("  none")
    else:
        print(hits[["subtype", "archetype", "fold_enrichment", "q_value"]].to_string(index=False))
    subtype_to_arc = {}
    for subtype in KEEP:
        sub_hits = hits[hits["subtype"] == subtype]
        if sub_hits.empty:
            print(f"  {subtype}: no significant bin-0 peak")
        else:
            arcs_hit = sorted(int(a) + 1 for a in sub_hits["archetype"].unique())
            subtype_to_arc[subtype] = arcs_hit
            print(f"  {subtype}: peaks at archetype(s) {arcs_hit}")
    empty = [i + 1 for i in range(n_arcs) if i + 1 not in {a for v in subtype_to_arc.values() for a in v}]
    print("  unmatched arcs:", empty)


def run_panel_c():
    OUT_C.mkdir(parents=True, exist_ok=True)
    n_pcs = int((PANEL_A / "n_pcs.txt").read_text().strip().splitlines()[0])
    expr_cl = load_expression_csv(MATRIX)
    saved = np.load(PANEL_A / "pc_scores.npy")
    arcs = np.load(PANEL_A / f"archetypes_k{K}_parti.npy")
    n_vol = arcs.shape[1]
    k = arcs.shape[0]
    pca_cl, scores_cl = fit_pca(expr_cl.T.values, n_components=n_pcs)
    pca_cl, scores_cl, _ = align_pca_signs(scores_cl, saved, pca_cl)
    arcs_full = np.zeros((k, n_pcs))
    arcs_full[:, :n_vol] = arcs
    gene_arcs = pd.DataFrame(
        inverse_transform_scores(pca_cl, arcs_full).T,
        index=expr_cl.index.astype(str),
        columns=[f"arc{i+1}" for i in range(k)],
    )
    tumors = pd.read_csv(TCGA_EXPR, index_col=0)
    tumors.index = tumors.index.astype(str)
    tumors.columns = tumors.columns.astype(str)
    rename, _ = rename_tumor_genes_to_cell_line(tumors.index, expr_cl.index)
    if rename:
        tumors = tumors.rename(index=rename)
    meta = pd.read_csv(TCGA_META, index_col=0)
    meta.index = meta.index.astype(str)
    tumors = tumors.loc[:, [c for c in tumors.columns if c in meta.index]]
    meta = meta.reindex(tumors.columns)
    shared = expr_cl.index.astype(str).intersection(tumors.index)
    cl_s = expr_cl.loc[shared]
    cl_s.columns = cl_s.columns.astype(str)
    tu_s = tumors.loc[shared]
    if not COMBAT_SRC.is_file():
        raise FileNotFoundError(COMBAT_SRC)
    combined = pd.read_csv(COMBAT_SRC, index_col=0)
    combined.columns = combined.columns.astype(str)
    combined = combined.loc[shared]
    gene_arcs = gene_arcs.loc[shared]
    is_cell = np.array([True] * cl_s.shape[1] + [False] * tu_s.shape[1])
    is_tumor = ~is_cell
    n_comp = min(N_PCS_VAR, combined.shape[1] - 1, combined.shape[0])
    pca_all, scores_all = fit_pca(combined.T.values, n_components=n_comp)
    arc_scores = pca_all.transform(gene_arcs.T.values)
    tumor_X = combined.loc[:, is_tumor].T.values
    tumor_total_var = float(tumor_X.var(axis=0, ddof=1).sum())
    pca_tumor, _ = fit_pca(tumor_X, n_components=n_comp)
    ev_tumor_only = 100.0 * pca_tumor.explained_variance_ratio_
    ev_combined = tumor_variance_on_pcs(tumor_X, pca_all, tumor_total_var)
    t3 = scores_all[is_tumor][:, :n_vol]
    a3 = arc_scores[:, :n_vol]
    dist = np.sqrt(((t3[:, None, :] - a3[None, :, :]) ** 2).sum(axis=2))
    nearest = dist.argmin(axis=1) + 1
    weights = barycentric_weights(t3, a3)
    inside = (weights >= -1e-6).all(axis=1)
    print(f"Tumors inside {n_vol}-D k={k} simplex: {int(inside.sum())}/{inside.size}")

    SUB_COL = "GeneExp_Subtype"
    tumor_ids = combined.columns[is_tumor]
    table = meta.copy()
    if SUB_COL in table.columns:
        table[SUB_COL] = table[SUB_COL].astype(str)
    score_cols = [f"PC{i+1}" for i in range(n_comp)]
    for i, c in enumerate(score_cols[:8]):
        table[c] = scores_all[is_tumor, i]
    dist_df = pd.DataFrame(dist, index=tumor_ids, columns=[f"dist_arc{i+1}" for i in range(k)])
    dist_df["nearest_archetype"] = nearest
    dist_df["inside_simplex"] = inside
    table = table.join(dist_df)
    table.to_csv(OUT_C / "tcga_sample_subtype_archetype.csv")
    labeled = table[table[SUB_COL].isin(SUBTYPE_COLOR)] if SUB_COL in table.columns else table.iloc[0:0]
    if len(labeled):
        pd.crosstab(labeled[SUB_COL], labeled["nearest_archetype"], margins=True).to_csv(
            OUT_C / "subtype_by_nearest_archetype_counts.csv"
        )
        pd.crosstab(labeled[SUB_COL], labeled["inside_simplex"], margins=True).to_csv(
            OUT_C / "subtype_inside_simplex_counts.csv"
        )

    cum_t = cumulative(ev_tumor_only)
    cum_c = cumulative(ev_combined)
    rng = np.random.default_rng(SEED)
    rows = []
    for i in range(N_SHUFFLE):
        ran = shuffle(combined, random_state=int(rng.integers(0, 2**31 - 1)))
        ran.index = combined.index
        pca_rand, _ = fit_pca(ran.T.values, n_components=n_comp)
        rows.append(tumor_variance_on_pcs(tumor_X, pca_rand, tumor_total_var))
    ev_null = np.vstack(rows)

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.3))
    ax = axes[0]
    ax.scatter(
        scores_all[is_cell, 0],
        scores_all[is_cell, 1],
        s=22,
        c="#B8B8B8",
        zorder=1,
        label=f"cell lines (n={int(is_cell.sum())})",
    )
    unlabeled = table[~table.index.isin(labeled.index)] if len(labeled) else table
    if len(unlabeled):
        ax.scatter(
            unlabeled["PC1"], unlabeled["PC2"], s=10, c="#DDDDDD", alpha=0.55, zorder=2,
            label=f"tumor, subtype missing (n={len(unlabeled)})",
        )
    for gname, color in SUBTYPE_COLOR.items():
        sub = labeled[labeled[SUB_COL] == gname] if len(labeled) else labeled
        if len(sub):
            ax.scatter(sub["PC1"], sub["PC2"], s=12, c=color, alpha=0.8, zorder=3, label=f"{gname} (n={len(sub)})")
    for i in range(k):
        for j in range(i + 1, k):
            ax.plot(
                [arc_scores[i, 0], arc_scores[j, 0]],
                [arc_scores[i, 1], arc_scores[j, 1]],
                color="#555555",
                lw=0.9,
                zorder=4,
            )
        ax.scatter(
            arc_scores[i, 0],
            arc_scores[i, 1],
            s=90,
            c=ARC_COLORS[i % len(ARC_COLORS)],
            edgecolors="k",
            linewidths=0.4,
            zorder=5,
        )
        ax.annotate(
            f"Arc {i+1}",
            (arc_scores[i, 0], arc_scores[i, 1]),
            textcoords="offset points",
            xytext=(5, 5),
            fontsize=8,
            fontweight="bold",
        )
    ax.set_xlabel("PC1 (combined, ComBat)")
    ax.set_ylabel("PC2 (combined, ComBat)")
    ax.set_title("TCGA tumors in cell-line archetype space")
    ax.legend(frameon=False, fontsize=7, loc="best")

    ax = axes[1]
    xs = np.arange(1, n_comp + 1)
    for row in ev_null:
        ax.plot(xs, cumulative(row), color="0.75", lw=0.8, alpha=0.8)
    ax.plot(xs, cum_t, marker="o", color="#4C78A8", label="tumor-only PCA (ceiling)")
    ax.plot(xs, cum_c, marker="o", color="#F58518", label="combined PCA on tumors")
    ax.axvline(k, color="0.5", ls="--", lw=1)
    ax.set_xlabel("number of components")
    ax.set_ylabel("cumulative % tumor variance")
    ratio_cum = cum_c / np.maximum(cum_t, 1e-12)
    pct = 100.0 * float(ratio_cum[k - 1])
    ax.set_title(f"{k} combined PCs = {pct:.0f}% of tumor-only ceiling")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle(f"Figure 1C  —  TCGA-GBM, Wang2017-restricted k={K}", y=1.02)
    fig.tight_layout()
    path = FIG / "Figure_1C_gbm_genelist_k6.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print("Wrote", path)
    print(f"inside simplex={int(inside.sum())}/{inside.size}")


def main():
    if not (PANEL_A / f"archetypes_k{K}_parti.npy").is_file():
        print("Missing", PANEL_A / f"archetypes_k{K}_parti.npy")
        return 1
    FIG.mkdir(parents=True, exist_ok=True)
    draw_panel_a()
    run_panel_b()
    run_panel_c()
    print("k=5 genelist figures left in place (Figure_1*_gbm_genelist.png).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
