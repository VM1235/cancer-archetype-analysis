#!/usr/bin/env python3
"""Panel C: project cell-line archetypes into combined cell-line + tumor PCA."""

from pathlib import Path
import sys

SCLC = Path(__file__).resolve().parents[1]
ROOT = SCLC.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.utils import shuffle

from src.io import load_combined_combat, load_expression_csv, sample_source_from_name
from src.pca import align_pca_signs, fit_pca, inverse_transform_scores

PANEL_A = SCLC / "results" / "panel_a"
OUT = SCLC / "results" / "panel_c"
OUT.mkdir(parents=True, exist_ok=True)

N_COMPONENTS = 20
N_SHUFFLE = 20
SEED = 0


def cumulative(values):
    return np.cumsum(np.asarray(values, dtype=float))


def tumor_variance_on_pcs(tumor_X, pca, total_var):
    scores = pca.transform(tumor_X)
    return 100.0 * scores.var(axis=0, ddof=1) / total_var


def shuffled_tumor_ev(combined, tumor_X, n_components, n_shuffles, total_var, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_shuffles):
        # Match the Groves notebook: permute gene rows, keep sample columns.
        ran = shuffle(combined, random_state=int(rng.integers(0, 2**31 - 1)))
        ran.index = combined.index
        pca_rand, _ = fit_pca(ran.T.values, n_components=n_components)
        rows.append(tumor_variance_on_pcs(tumor_X, pca_rand, total_var))
        print(f"    shuffle {i + 1}/{n_shuffles}", flush=True)
    return np.vstack(rows)


def main():
    expr_cl = load_expression_csv()
    combined = load_combined_combat()
    saved_scores = np.load(PANEL_A / "pc_scores_12.npy")
    archetypes = np.load(PANEL_A / "archetypes_k5.npy")

    pca12, scores12 = fit_pca(expr_cl.T.values, n_components=12)
    pca12, scores12, signs = align_pca_signs(scores12, saved_scores, pca12)
    max_diff = np.max(np.abs(scores12 - saved_scores))
    print(f"Aligned 12-PC scores; max abs diff vs saved = {max_diff:.4g}")
    np.save(PANEL_A / "pca12_components.npy", pca12.components_)
    np.save(PANEL_A / "pca12_mean.npy", pca12.mean_)

    gene_arcs = pd.DataFrame(
        inverse_transform_scores(pca12, archetypes).T,
        index=expr_cl.index,
        columns=[f"arc{i + 1}" for i in range(5)],
    )

    shared = gene_arcs.index.intersection(combined.index)
    print(f"Combined matrix: {combined.shape[0]} genes × {combined.shape[1]} samples")
    print(f"Shared genes with cell-line archetypes: {len(shared)}")
    combined = combined.loc[shared]
    gene_arcs = gene_arcs.loc[shared]

    sources = pd.Series(
        [sample_source_from_name(c) for c in combined.columns],
        index=combined.columns,
        name="source",
    )
    is_tumor = sources.eq("Tumor")
    is_cell = ~is_tumor
    print(sources.value_counts().to_string())

    X_all = combined.T.values
    pca_all, scores_all = fit_pca(X_all, n_components=N_COMPONENTS)
    arc_scores = pca_all.transform(gene_arcs.T.values)

    tumor = combined.loc[:, is_tumor]
    tumor_X = tumor.T.values
    tumor_total_var = float(tumor_X.var(axis=0, ddof=1).sum())
    pca_tumor, _ = fit_pca(tumor_X, n_components=N_COMPONENTS)

    ev_tumor_only = 100.0 * pca_tumor.explained_variance_ratio_
    ev_combined = tumor_variance_on_pcs(tumor_X, pca_all, tumor_total_var)
    print("Shuffled null PCA fits ...", flush=True)
    ev_null = shuffled_tumor_ev(
        combined, tumor_X, N_COMPONENTS, N_SHUFFLE, tumor_total_var, SEED
    )

    cum_tumor = cumulative(ev_tumor_only)
    cum_combined = cumulative(ev_combined)
    ratio = cum_combined / cum_tumor
    print("Cumulative combined / tumor-only at k=1..8:")
    for k, r in enumerate(ratio[:8], start=1):
        print(f"  {k}: {100 * r:.1f}%")

    pd.DataFrame(
        {
            "n_components": np.arange(1, N_COMPONENTS + 1),
            "ev_tumor_only": ev_tumor_only,
            "ev_combined_on_tumors": ev_combined,
            "cum_tumor_only": cum_tumor,
            "cum_combined_on_tumors": cum_combined,
            "ratio_combined_over_tumor": ratio,
        }
    ).to_csv(OUT / "variance_explained.csv", index=False)
    np.save(OUT / "ev_null_shuffles.npy", ev_null)
    pd.DataFrame(scores_all, index=combined.columns, columns=[f"PC{i+1}" for i in range(N_COMPONENTS)]).to_csv(
        OUT / "combined_pca_scores.csv"
    )
    pd.DataFrame(arc_scores, index=gene_arcs.columns, columns=[f"PC{i+1}" for i in range(N_COMPONENTS)]).to_csv(
        OUT / "archetypes_in_combined_pca.csv"
    )
    gene_arcs.to_csv(OUT / "archetypes_gene_space_shared.csv")
    sources.to_csv(OUT / "sample_sources.csv")

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(
        scores_all[is_cell.values, 0],
        scores_all[is_cell.values, 1],
        s=22,
        c="#9BD770",
        alpha=0.85,
        label="cell lines",
    )
    ax.scatter(
        scores_all[is_tumor.values, 0],
        scores_all[is_tumor.values, 1],
        s=22,
        c="#2E7D32",
        alpha=0.9,
        label="tumors",
    )
    ax.scatter(arc_scores[:, 0], arc_scores[:, 1], s=90, c="#E45756", zorder=3, label="archetypes")
    for i in range(5):
        for j in range(i + 1, 5):
            ax.plot(
                [arc_scores[i, 0], arc_scores[j, 0]],
                [arc_scores[i, 1], arc_scores[j, 1]],
                color="#E45756",
                lw=1,
                zorder=2,
            )
        ax.annotate(f"A{i+1}", (arc_scores[i, 0], arc_scores[i, 1]), fontsize=8, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("PC1 (combined)")
    ax.set_ylabel("PC2 (combined)")
    ax.set_title("Panel C: tumors in cell-line archetype space")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "combined_pca_scatter.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    xs = np.arange(1, N_COMPONENTS + 1)
    for row in ev_null:
        ax.plot(xs, cumulative(row), color="0.75", lw=0.8, alpha=0.8)
    ax.plot(xs, cum_tumor, marker="o", color="#4C78A8", label="tumor-only PCA (ceiling)")
    ax.plot(xs, cum_combined, marker="o", color="#F58518", label="combined PCA on tumors")
    ax.axvline(5, color="0.5", ls="--", lw=1)
    ax.set_xlabel("number of components")
    ax.set_ylabel("cumulative % tumor variance")
    ax.set_title(f"5 combined PCs capture {100 * ratio[4]:.0f}% of tumor-only variance")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "variance_explained.png", dpi=150)
    plt.close(fig)

    print(f"\nWrote outputs to {OUT}")
    print(f"5-component ratio vs tumor-only ceiling: {100 * ratio[4]:.1f}% (paper ~80%)")


if __name__ == "__main__":
    main()
