#!/usr/bin/env python3
"""Panel B: subtype clustering vs distance-to-archetype enrichment."""

from pathlib import Path
import sys

SCLC = Path(__file__).resolve().parents[1]
ROOT = SCLC.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score

from src.io import load_author_subtypes, load_expression_csv
from src.enrichment import (
    SUBTYPES,
    distance_bins,
    hypergeometric_enrichment,
    name_clusters_by_overlap,
    name_clusters_by_tfs,
    spearman_average_clusters,
)

PANEL_A = SCLC / "results" / "panel_a"
OUT = SCLC / "results" / "panel_b"
OUT.mkdir(parents=True, exist_ok=True)

ARC_COLORS = {
    0: "#4C78A8",
    1: "#F58518",
    2: "#E45756",
    3: "#72B7B2",
    4: "#54A24B",
}


def plot_enrichment(table, title, path):
    fig, axes = plt.subplots(1, 5, figsize=(14, 3.2), sharey=True)
    for ax, subtype in zip(axes, SUBTYPES):
        sub = table[table["subtype"] == subtype]
        for arc in sorted(sub["archetype"].unique()):
            arc_tab = sub[sub["archetype"] == arc].sort_values("bin")
            ax.plot(
                arc_tab["bin"],
                arc_tab["fold_enrichment"],
                marker="o",
                color=ARC_COLORS[int(arc)],
                label=f"Arc {int(arc) + 1}",
            )
            sig = arc_tab[arc_tab["sig_peak_at_bin0"]]
            if len(sig):
                ax.scatter(
                    sig["bin"],
                    sig["fold_enrichment"],
                    s=60,
                    facecolors="none",
                    edgecolors=ARC_COLORS[int(arc)],
                    linewidths=1.5,
                    zorder=3,
                )
        ax.axhline(1.0, color="0.7", lw=1, ls="--")
        ax.set_title(f"SCLC-{subtype}")
        ax.set_xlabel("distance bin")
        ax.set_xticks(range(10))
    axes[0].set_ylabel("fold enrichment")
    axes[-1].legend(loc="upper right", fontsize=8, frameon=False)
    fig.suptitle(title, y=1.04)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run_one(subtypes, scores, archetypes, tag):
    bin_ids, distances = distance_bins(scores, archetypes, n_bins=10)
    pd.DataFrame(distances, index=subtypes.index, columns=[f"arc{i+1}" for i in range(5)]).to_csv(
        OUT / f"distances_{tag}.csv"
    )
    pd.DataFrame(bin_ids, index=subtypes.index, columns=[f"arc{i+1}" for i in range(5)]).to_csv(
        OUT / f"bins_{tag}.csv"
    )
    table = hypergeometric_enrichment(subtypes, bin_ids, fdr=0.1)
    table.to_csv(OUT / f"enrichment_{tag}.csv", index=False)
    plot_enrichment(
        table,
        title=f"Panel B ({tag} labels)",
        path=OUT / f"enrichment_{tag}.png",
    )
    hits = table[table["sig_peak_at_bin0"]][
        ["subtype", "archetype", "fold_enrichment", "p_value", "q_value"]
    ]
    print(f"\nSignificant bin-0 peaks ({tag}):")
    if hits.empty:
        print("  none")
    else:
        print(hits.to_string(index=False))
    return table


def main():
    expr = load_expression_csv()
    scores = np.load(PANEL_A / "pc_scores_12.npy")
    archetypes = np.load(PANEL_A / "archetypes_k5.npy")
    sample_names = list(expr.columns)
    if scores.shape[0] != len(sample_names):
        raise ValueError("PCA scores and expression samples are misaligned")

    author = load_author_subtypes().reindex(sample_names)
    if author.isna().any():
        missing = author[author.isna()].index.tolist()
        raise ValueError(f"Author labels missing for: {missing[:5]}")

    clusters = spearman_average_clusters(expr, n_clusters=5)
    tf_named, tf_map, tf_means = name_clusters_by_tfs(expr, clusters)
    overlap_named, overlap_map = name_clusters_by_overlap(clusters, author)

    labels = pd.DataFrame(
        {
            "cluster": clusters,
            "subtype_tf": tf_named,
            "subtype_overlap": overlap_named,
            "author_NEW_10_2020": author,
        }
    )
    labels.to_csv(OUT / "subtype_labels.csv")
    tf_means.to_csv(OUT / "cluster_tf_means.csv")

    ari_tf = adjusted_rand_score(author, tf_named)
    ari_overlap = adjusted_rand_score(author, overlap_named)
    print("Cluster → TF map:", tf_map)
    print("Cluster → overlap map:", overlap_map)
    print(f"ARI vs author labels (TF naming): {ari_tf:.3f}")
    print(f"ARI vs author labels (overlap naming): {ari_overlap:.3f}")
    print("Our TF-named counts:\n", tf_named.value_counts().to_string())
    print("Author counts:\n", author.value_counts().to_string())

    run_one(tf_named, scores, archetypes, tag="tf")
    run_one(author, scores, archetypes, tag="author")


if __name__ == "__main__":
    main()
