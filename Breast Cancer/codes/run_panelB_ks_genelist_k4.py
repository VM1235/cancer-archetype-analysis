#!/usr/bin/env python3
"""Breast Panel B on KS genelist-restricted archetypes, FORCED to k=4.

Uses the extended-k Panel A run (which force-fit enough PCs to also test
k=4,5,6,7) since the default panel_a_ks_genelist/ only fit 2 PCs (enough
for k=3 alone). Writes panel_b_ks_genelist_k4/, does not touch panel_b_ks_genelist/.
"""

from collections import defaultdict
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.enrichment import distance_bins, hypergeometric_enrichment
from src.io import load_expression_csv

MATRIX = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
SCORES = BREAST / "results" / "panel_a_ks_genelist_extendedk" / "pc_scores.npy"
LABELS_SRC = BREAST / "results" / "panel_b" / "pam50_labels_panelA.csv"
OUT = BREAST / "results" / "panel_b_ks_genelist_k4"
FIG = BREAST / "figures"
PANEL_A = BREAST / "results" / "panel_a_ks_genelist_extendedk"

DROP = ("LumA", "Normal")
KEEP = ("LumB", "Her2", "Basal")
N_BINS = 5
FDR = 0.1

ARC_COLORS = {
    0: "#4C78A8",
    1: "#F58518",
    2: "#E45756",
    3: "#72B7B2",
}


def plot_enrichment(table, path, n_bins, n_arcs):
    fig, axes = plt.subplots(1, len(KEEP), figsize=(10.2, 3.4), sharey=True)
    for ax, subtype in zip(axes, KEEP):
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
                    s=70,
                    facecolors="none",
                    edgecolors=ARC_COLORS[int(arc)],
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
        "Figure 1B  —  PAM50 enrichment (KS gene-list archetypes; LumA/Normal removed)",
        y=1.06,
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    k_star = 4  # forced: exploring the k=4 solution specifically, not the
                # smallest-significant-k default (which is k=3)
    print(f"Forcing k={k_star} (Panel A also found k=3 and k=5 significant; see extendedk sweep)")
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
    labels.to_csv(OUT / "pam50_labels_all63.csv")

    dropped = labels[labels["pam50_subtype"].isin(DROP)]
    dropped.to_csv(OUT / "pam50_dropped_lumA_normal.csv")
    print(f"Dropped LumA/Normal: {len(dropped)} lines")
    print(dropped[["pam50_subtype", "confidence_score"]].to_string())

    keep_mask = labels["pam50_subtype"].isin(KEEP)
    n_keep = int(keep_mask.sum())
    print(f"Kept for enrichment: {n_keep}/{len(sample_ids)}")
    labels.loc[keep_mask].to_csv(OUT / "pam50_labels_panelB_input.csv")

    subtypes = labels.loc[keep_mask, "pam50_subtype"]
    scores_u = scores[keep_mask.values][:, :n_vol]
    print("PAM50 counts after drop:\n", subtypes.value_counts().to_string())

    bin_ids, distances = distance_bins(scores_u, arcs, n_bins=N_BINS)
    print("\nBin sizes (rows=archetypes, cols=bins 0-4):")
    for j in range(arcs.shape[0]):
        sizes = [(bin_ids[:, j] == b).sum() for b in range(N_BINS)]
        print(f"  Arc {j+1}: {sizes}  min={min(sizes)} max={max(sizes)}")
        if min(sizes) <= 1:
            raise ValueError(f"Degenerate bin for archetype {j+1}: {sizes}")

    n_arcs = arcs.shape[0]
    pd.DataFrame(
        distances, index=subtypes.index, columns=[f"arc{i+1}" for i in range(n_arcs)]
    ).to_csv(OUT / "distances_pam50.csv")
    pd.DataFrame(
        bin_ids, index=subtypes.index, columns=[f"arc{i+1}" for i in range(n_arcs)]
    ).to_csv(OUT / "bins_pam50.csv")

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
    table.to_csv(OUT / "enrichment_pam50.csv", index=False)
    plot_enrichment(table, FIG / "Figure_1B_breast_ks_genelist_k4.png", N_BINS, n_arcs)
    plot_enrichment(table, OUT / "enrichment_pam50.png", N_BINS, n_arcs)

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

    print("\n=== Run_2 Panel B summary ===")
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
    collapsed = {a: s for a, s in arc_to_sub.items() if len(s) > 1}
    if collapsed:
        print("  Subtypes sharing an archetype:")
        for a, s in sorted(collapsed.items()):
            print(f"    Arc {a}: {s}")
    else:
        print("  No two subtypes collapsed onto the same archetype.")
    empty_arcs = [i + 1 for i in range(n_arcs) if (i + 1) not in arc_to_sub]
    if empty_arcs:
        print(f"  Archetypes with no significant subtype match: {empty_arcs}")
    else:
        print("  Every archetype has at least one significant subtype match.")
    print(
        f"  Enrichment n={n_keep} (dropped {len(dropped)} LumA/Normal from 63). "
        "Panel A polytope still used all 63 lines."
    )
    print("Wrote", OUT / "enrichment_pam50.csv")
    print("Wrote", FIG / "Figure_1B_breast_ks_genelist_k4.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
