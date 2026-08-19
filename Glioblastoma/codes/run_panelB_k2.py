#!/usr/bin/env python3
"""Panel B for hypothesis-driven k=2: Proneural vs Mesenchymal only.

Polytope = all 54 lines (k=2 PCHA). Enrichment drops Classical.
Distances are Euclidean in the 1-D fit space. Equal-count bins (3–5).
"""

from collections import defaultdict
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
GBM = HERE.parent
ROOT = GBM.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.enrichment import distance_bins, hypergeometric_enrichment
from src.io import load_expression_csv
from k2_utils import SUBTYPE_COLOR, corr_table, style_mpl

MATRIX = GBM / "data" / "processed" / "input_panelA_glioblastoma_ccle_logtpm_filtered.csv"
SCORES = GBM / "results" / "panel_a" / "pc_scores.npy"
ARCS = GBM / "results" / "panel_a_k2" / "archetypes_k2_parti.npy"
WEIGHTS = GBM / "results" / "panel_a_k2" / "S_k2_parti.npy"
POLES = GBM / "results" / "panel_a_k2" / "pole_assignment.txt"
LABELS_SRC = GBM / "results" / "panel_b" / "verhaak_labels_panelA.csv"
SIG = GBM / "results" / "panel_a_k2" / "signature_scores_cell_lines.csv"
OUT = GBM / "results" / "panel_b_k2"
FIG = GBM / "figures"

DROP = ("Classical",)
KEEP = ("Proneural", "Mesenchymal")
N_BINS_TRY = (5, 4, 3)
FDR = 0.1

ARC_COLORS = {0: "#4C78A8", 1: "#E07A3D"}


def read_poles():
    mes_idx, pn_idx = 1, 0
    for line in POLES.read_text().splitlines():
        if line.startswith("mes_idx="):
            mes_idx = int(line.split("=", 1)[1])
        if line.startswith("pn_idx="):
            pn_idx = int(line.split("=", 1)[1])
    return mes_idx, pn_idx


def plot_enrichment(table, path, n_bins, mes_idx, pn_idx):
    style_mpl(plt)
    fig, axes = plt.subplots(1, len(KEEP), figsize=(8.4, 3.6), sharey=True)
    names = {pn_idx: "PN pole", mes_idx: "MES pole"}
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
                label=f"Arc {int(arc) + 1} ({names.get(int(arc), '')})",
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
        ax.set_xlabel("distance bin (0 = closest)")
        ax.set_xticks(range(n_bins))
    axes[0].set_ylabel("fold enrichment")
    axes[-1].legend(loc="upper right", fontsize=7, frameon=False)
    fig.suptitle(
        "Figure 1B (k=2)  —  PN vs MES enrichment on the 1-D axis (not Groves k)",
        y=1.06,
        fontsize=12,
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
    weights = np.load(WEIGHTS)
    n_vol = arcs.shape[1]
    if n_vol != 1:
        raise ValueError(f"Expected 1-D archetypes for k=2, got {arcs.shape}")
    print(f"1-D fit space: scores {scores.shape} -> first {n_vol} PC; archetypes {arcs.shape}")
    print("Enrichment: Proneural vs Mesenchymal only; Classical dropped from the test.")
    print("Polytope still uses all 54 lines.")

    mes_idx, pn_idx = read_poles()
    print(f"Poles from Panel A k=2: Arc {mes_idx + 1}=MES, Arc {pn_idx + 1}=PN")

    labels = pd.read_csv(LABELS_SRC)
    labels["cell_line"] = labels["cell_line"].astype(str)
    labels = labels.set_index("cell_line").reindex(sample_ids)
    labels.to_csv(OUT / "verhaak_labels_all.csv")

    dropped = labels[labels["verhaak_subtype"].isin(DROP)]
    dropped.to_csv(OUT / "verhaak_dropped.csv")
    print(f"Dropped {list(DROP)}: {len(dropped)} lines")

    keep_mask = labels["verhaak_subtype"].isin(KEEP)
    n_keep = int(keep_mask.sum())
    labels.loc[keep_mask].to_csv(OUT / "verhaak_labels_panelB_input.csv")
    subtypes = labels.loc[keep_mask, "verhaak_subtype"]
    scores_u = scores[keep_mask.values][:, :n_vol]
    print(f"Kept for enrichment: {n_keep}/{len(sample_ids)}")
    print("Counts:\n", subtypes.value_counts().to_string())

    chosen = None
    bin_ids = distances = None
    for n_bins in N_BINS_TRY:
        if n_keep < n_bins * 2:
            print(f"Skip n_bins={n_bins}: n={n_keep} too small")
            continue
        bid, dist = distance_bins(scores_u, arcs, n_bins=n_bins)
        sizes_ok = True
        print(f"\nTrying n_bins={n_bins} (rows=archetypes, cols=bins):")
        for j in range(arcs.shape[0]):
            sizes = [(bid[:, j] == b).sum() for b in range(n_bins)]
            print(f"  Arc {j + 1}: {sizes}  min={min(sizes)} max={max(sizes)}")
            if min(sizes) <= 1:
                sizes_ok = False
        if sizes_ok:
            chosen = n_bins
            bin_ids, distances = bid, dist
            break
        print(f"  degenerate bins at n_bins={n_bins}")
    if chosen is None:
        raise ValueError("Could not find non-degenerate equal-count bins in 3–5")
    n_bins = chosen
    print(f"Using n_bins={n_bins}")

    n_arcs = arcs.shape[0]
    pd.DataFrame(
        distances, index=subtypes.index, columns=[f"arc{i + 1}" for i in range(n_arcs)]
    ).to_csv(OUT / "distances_verhaak.csv")
    pd.DataFrame(
        bin_ids, index=subtypes.index, columns=[f"arc{i + 1}" for i in range(n_arcs)]
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
    plot_enrichment(table, FIG / "Figure_1B_gbm_k2.png", n_bins, mes_idx, pn_idx)
    plot_enrichment(table, OUT / "enrichment_verhaak.png", n_bins, mes_idx, pn_idx)

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

    print("\n=== 1-to-1 PN vs MES vs vertices? ===")
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

    expected = {"Proneural": pn_idx + 1, "Mesenchymal": mes_idx + 1}
    one_to_one = (
        subtype_to_arc.get("Proneural") == [expected["Proneural"]]
        and subtype_to_arc.get("Mesenchymal") == [expected["Mesenchymal"]]
    )
    print(f"  Signature-based poles: PN=Arc {expected['Proneural']}, MES=Arc {expected['Mesenchymal']}")
    print(f"  Groves-style 1-to-1 (sig bin-0 peaks match those poles): {one_to_one}")

    arc_to_sub = defaultdict(list)
    for subtype, arcs_hit in subtype_to_arc.items():
        if arcs_hit:
            for a in arcs_hit:
                arc_to_sub[a].append(subtype)
    empty_arcs = [i + 1 for i in range(n_arcs) if (i + 1) not in arc_to_sub]
    if empty_arcs:
        print(f"  Archetypes with no significant subtype match: {empty_arcs}")

    if SIG.exists() and WEIGHTS.exists():
        sig = pd.read_csv(SIG, index_col=0)
        sig.index = sig.index.astype(str)
        mix = pd.DataFrame(
            {
                "w_MES": weights[:, mes_idx],
                "w_PN": weights[:, pn_idx],
            },
            index=sample_ids,
        )
        cor = pd.concat(
            [
                corr_table(mix["w_MES"], sig, "w_MES"),
                corr_table(mix["w_PN"], sig, "w_PN"),
            ]
        )
        cor.to_csv(OUT / "signature_vs_mixture_cell_lines.csv", index=False)
        print("\nCell-line signature vs PCHA mixture (Pearson r):")
        print(cor[["weight", "signature", "pearson_r", "pearson_p"]].to_string(index=False))

    print(
        f"  Enrichment n={n_keep} (dropped {len(dropped)} Classical). "
        f"k=2 polytope used all {len(sample_ids)} lines."
    )
    print("Wrote", OUT / "enrichment_verhaak.csv")
    print("Wrote", FIG / "Figure_1B_gbm_k2.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
