#!/usr/bin/env python3
"""GBM Panel B on Wang2017-restricted archetypes. Writes panel_b_genelist/."""

from collections import defaultdict
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
GBM = HERE.parent
ROOT = GBM.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.enrichment import distance_bins, hypergeometric_enrichment
from src.io import load_expression_csv

MATRIX = GBM / "data" / "processed" / "input_panelA_wang2017_genelist.csv"
SCORES = GBM / "results" / "panel_a_genelist" / "pc_scores.npy"
LABELS_SRC = GBM / "results" / "panel_b" / "verhaak_labels_panelA.csv"
OUT = GBM / "results" / "panel_b_genelist"
FIG = GBM / "figures"
PANEL_A = GBM / "results" / "panel_a_genelist"

DROP = ()
KEEP = ("Classical", "Mesenchymal", "Proneural")
N_BINS = 5
FDR = 0.1

ARC_COLORS = {
    0: "#4C78A8",
    1: "#F58518",
    2: "#E45756",
    3: "#72B7B2",
    4: "#54A24B",
}


def plot_enrichment(table, path, n_bins):
    fig, axes = plt.subplots(1, len(KEEP), figsize=(10.2, 3.4), sharey=True)
    for ax, subtype in zip(np.atleast_1d(axes), KEEP):
        sub = table[table["subtype"] == subtype]
        for arc in sorted(sub["archetype"].unique()):
            arc_tab = sub[sub["archetype"] == arc].sort_values("bin")
            color = ARC_COLORS.get(int(arc), "#333333")
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
        ax.set_xticks(range(n_bins))
    axes[0].set_ylabel("fold enrichment")
    axes[-1].legend(loc="upper right", fontsize=7, frameon=False)
    fig.suptitle(
        "Figure 1B  —  Verhaak-style markers on Wang2017-restricted archetypes",
        y=1.06,
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    suggested = (PANEL_A / "suggested_k.txt").read_text().splitlines()
    k_star = int(suggested[0])
    print(f"Suggested k from Panel A genelist: {k_star} ({suggested[1] if len(suggested) > 1 else ''})")
    arcs_path = PANEL_A / f"archetypes_k{k_star}_parti.npy"

    expr = load_expression_csv(MATRIX)
    sample_ids = list(expr.columns.astype(str))
    scores = np.load(SCORES)
    arcs = np.load(arcs_path)
    if scores.shape[0] != len(sample_ids):
        raise ValueError(
            f"PCA scores n={scores.shape[0]} vs matrix samples n={len(sample_ids)}"
        )
    n_vol = arcs.shape[1]
    print(f"Using scores {SCORES.name} {scores.shape} -> first {n_vol} PCs")
    print(f"Using archetypes {arcs_path.name} {arcs.shape}")

    labels = pd.read_csv(LABELS_SRC)
    labels["cell_line"] = labels["cell_line"].astype(str)
    labels = labels.set_index("cell_line").reindex(sample_ids)
    labels.to_csv(OUT / "verhaak_labels_all.csv")

    dropped = labels[labels["verhaak_subtype"].isin(DROP)] if DROP else labels.iloc[0:0]
    if DROP:
        dropped.to_csv(OUT / "verhaak_dropped.csv")
        print(f"Dropped {list(DROP)}: {len(dropped)} lines")
    else:
        print("No subtypes dropped from enrichment")

    keep_mask = labels["verhaak_subtype"].isin(KEEP)
    n_keep = int(keep_mask.sum())
    print(f"Kept for enrichment: {n_keep}/{len(sample_ids)}")
    labels.loc[keep_mask].to_csv(OUT / "verhaak_labels_panelB_input.csv")

    subtypes = labels.loc[keep_mask, "verhaak_subtype"]
    scores_u = scores[keep_mask.values][:, :n_vol]
    print("Verhaak counts after drop:\n", subtypes.value_counts().to_string())

    n_bins = N_BINS
    if n_keep < n_bins * 2:
        n_bins = max(3, n_keep // 8)
        print(f"Reduced n_bins to {n_bins} because n={n_keep}")

    bin_ids, distances = distance_bins(scores_u, arcs, n_bins=n_bins)
    print("\nBin sizes (rows=archetypes, cols=bins):")
    for j in range(arcs.shape[0]):
        sizes = [(bin_ids[:, j] == b).sum() for b in range(n_bins)]
        print(f"  Arc {j+1}: {sizes}  min={min(sizes)} max={max(sizes)}")
        if min(sizes) <= 1:
            raise ValueError(f"Degenerate bin for archetype {j+1}: {sizes}")

    n_arcs = arcs.shape[0]
    pd.DataFrame(
        distances, index=subtypes.index, columns=[f"arc{i+1}" for i in range(n_arcs)]
    ).to_csv(OUT / "distances_verhaak.csv")
    pd.DataFrame(
        bin_ids, index=subtypes.index, columns=[f"arc{i+1}" for i in range(n_arcs)]
    ).to_csv(OUT / "bins_verhaak.csv")

    table = hypergeometric_enrichment(
        subtypes, bin_ids, fdr=FDR, subtype_levels=KEEP
    )
    keep_cols = [
        "archetype",
        "subtype",
        "bin",
        "fold_enrichment",
        "p_value",
        "q_value",
        "sig_peak_at_bin0",
        "k_in_bin",
        "bin_size",
        "K_subtype",
        "significant",
        "peak_bin",
    ]
    table = table[keep_cols]
    table.to_csv(OUT / "enrichment_verhaak.csv", index=False)
    plot_enrichment(table, FIG / "Figure_1B_gbm_genelist.png", n_bins)
    plot_enrichment(table, OUT / "enrichment_verhaak.png", n_bins)

    hits = table[table["sig_peak_at_bin0"]]
    print("\nSignificant bin-0 peaks (q < 0.1 and peak at bin 0):")
    if hits.empty:
        print("  none")
    else:
        print(
            hits[
                ["subtype", "archetype", "fold_enrichment", "p_value", "q_value"]
            ].to_string(index=False)
        )

    print("\n=== GBM Panel B summary ===")
    subtype_to_arc = {}
    for subtype in KEEP:
        sub_hits = hits[hits["subtype"] == subtype]
        if sub_hits.empty:
            subtype_to_arc[subtype] = None
            print(f"  {subtype}: no significant bin-0 peak")
        else:
            arcs_hit = sorted(int(a) + 1 for a in sub_hits["archetype"].unique())
            subtype_to_arc[subtype] = arcs_hit
            print(f"  {subtype}: peaks at archetype(s) {arcs_hit}")

    arc_to_sub = defaultdict(list)
    for subtype, arcs_hit in subtype_to_arc.items():
        if arcs_hit:
            for a in arcs_hit:
                arc_to_sub[a].append(subtype)
    empty_arcs = [i + 1 for i in range(n_arcs) if (i + 1) not in arc_to_sub]
    if empty_arcs:
        print(f"  Archetypes with no significant subtype match: {empty_arcs}")
    print(
        f"  Enrichment n={n_keep} (dropped {len(dropped)}). "
        f"Panel A polytope used all {len(sample_ids)} lines."
    )
    print("Wrote", OUT / "enrichment_verhaak.csv")
    print("Wrote", FIG / "Figure_1B_gbm_genelist.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
