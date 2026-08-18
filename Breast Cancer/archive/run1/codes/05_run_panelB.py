#!/usr/bin/env python3
"""Panel B: PAM50 enrichment vs distance to k=4 breast archetypes.

Uses the same 8-PC scores as Panel A, but distances are in the (k-1)=3-D
space the k=4 simplex was actually fit in (do not pad PCs 4–8 with zeros).
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parents[2]
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))
RUN1 = HERE.parent

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.enrichment import distance_bins, hypergeometric_enrichment
from src.io import load_expression_csv

MATRIX = BREAST / "data" / "processed" / "input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv"
SCORES = RUN1 / "results" / "panel_a" / "pc_scores.npy"
ARCS = RUN1 / "results" / "panel_a" / "archetypes_k4_parti.npy"
LABELS = BREAST / "results" / "panel_b" / "pam50_labels_panelA.csv"
OUT = RUN1 / "results" / "panel_b"
FIG = RUN1 / "figures"

PAM50 = ("LumA", "LumB", "Her2", "Basal", "Normal")
N_BINS = 5
FDR = 0.1

ARC_COLORS = {
    0: "#4C78A8",
    1: "#F58518",
    2: "#E45756",
    3: "#72B7B2",
}


def plot_enrichment(table, path, n_bins):
    fig, axes = plt.subplots(1, 5, figsize=(13.5, 3.4), sharey=True)
    for ax, subtype in zip(axes, PAM50):
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
        "Figure 1B  —  PAM50 subtype enrichment at k=4 archetypes (bin 0 = closest)",
        y=1.06,
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    expr = load_expression_csv(MATRIX)
    sample_ids = list(expr.columns.astype(str))
    scores = np.load(SCORES)
    arcs = np.load(ARCS)
    if scores.shape[0] != len(sample_ids):
        raise ValueError(
            f"PCA scores n={scores.shape[0]} vs matrix samples n={len(sample_ids)}"
        )
    if arcs.shape != (4, 3):
        raise ValueError(
            f"Expected breast k=4 archetypes shape (4, 3); got {arcs.shape}. "
            "Refusing to use a leftover SCLC file."
        )
    if scores.shape[1] < 3:
        raise ValueError(f"PCA scores have {scores.shape[1]} PCs; need >= 3")

    # Fit space for k=4 is the first 3 PCs.
    scores3 = scores[:, :3]
    print(f"Using scores {SCORES.name} {scores.shape} -> first 3 PCs {scores3.shape}")
    print(f"Using archetypes {ARCS.name} {arcs.shape}")

    labels = pd.read_csv(LABELS)
    labels["cell_line"] = labels["cell_line"].astype(str)
    labels = labels.set_index("cell_line").reindex(sample_ids)
    usable = labels["pam50_subtype"].isin(PAM50)
    n_usable = int(usable.sum())
    print(f"Usable PAM50 labels: {n_usable}/{len(sample_ids)}")
    if n_usable < int(0.8 * len(sample_ids)):
        print("STOP: usable PAM50 coverage < 80%.")
        return 2

    subtypes = labels.loc[usable, "pam50_subtype"]
    scores_u = scores3[usable.values]
    print("PAM50 counts:\n", subtypes.value_counts().to_string())

    bin_ids, distances = distance_bins(scores_u, arcs, n_bins=N_BINS)
    print("\nBin sizes (rows=archetypes 1-4, cols=bins 0-4):")
    for j in range(arcs.shape[0]):
        sizes = [(bin_ids[:, j] == b).sum() for b in range(N_BINS)]
        print(f"  Arc {j+1}: {sizes}  min={min(sizes)} max={max(sizes)}")
        if min(sizes) <= 1:
            raise ValueError(f"Degenerate bin for archetype {j+1}: {sizes}")

    pd.DataFrame(
        distances, index=subtypes.index, columns=[f"arc{i+1}" for i in range(4)]
    ).to_csv(OUT / "distances_pam50.csv")
    pd.DataFrame(
        bin_ids, index=subtypes.index, columns=[f"arc{i+1}" for i in range(4)]
    ).to_csv(OUT / "bins_pam50.csv")

    table = hypergeometric_enrichment(
        subtypes, bin_ids, fdr=FDR, subtype_levels=PAM50
    )
    keep = [
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
    table = table[keep]
    table.to_csv(OUT / "enrichment_pam50.csv", index=False)
    plot_enrichment(table, FIG / "Figure_1B_breast.png", N_BINS)
    plot_enrichment(table, OUT / "enrichment_pam50.png", N_BINS)

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

    print("\n=== Panel B summary ===")
    subtype_to_arc = {}
    for subtype in PAM50:
        sub_hits = hits[hits["subtype"] == subtype]
        if sub_hits.empty:
            subtype_to_arc[subtype] = None
            print(f"  {subtype}: no significant bin-0 peak")
        else:
            arcs_hit = sorted(int(a) + 1 for a in sub_hits["archetype"].unique())
            subtype_to_arc[subtype] = arcs_hit
            print(f"  {subtype}: peaks at archetype(s) {arcs_hit}")

    from collections import defaultdict

    arc_to_sub = defaultdict(list)
    for subtype, arcs_hit in subtype_to_arc.items():
        if arcs_hit:
            for a in arcs_hit:
                arc_to_sub[a].append(subtype)
    collapsed = {a: s for a, s in arc_to_sub.items() if len(s) > 1}
    if collapsed:
        print("  Subtypes sharing an archetype (expected with 5 vs 4):")
        for a, s in sorted(collapsed.items()):
            print(f"    Arc {a}: {s}")
    else:
        print("  No two subtypes collapsed onto the same archetype.")

    empty_arcs = [i + 1 for i in range(4) if (i + 1) not in arc_to_sub]
    if empty_arcs:
        print(f"  Archetypes with no significant subtype match: {empty_arcs}")
    else:
        print("  Every archetype has at least one significant subtype match.")
    print(f"  Final match coverage: {n_usable}/{len(sample_ids)} lines with a usable PAM50 label")
    print("Wrote", OUT / "enrichment_pam50.csv")
    print("Wrote", FIG / "Figure_1B_breast.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
