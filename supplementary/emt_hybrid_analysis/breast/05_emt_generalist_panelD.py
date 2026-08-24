#!/usr/bin/env python3
"""Panel D: are hybrid E/M breast cell lines archetype "generalists"?

Combines two things already in this repo/paper set:

1. George, Jolly, Xu, Somarelli & Levine (Cancer Res 2017) - a
   2-predictor (CLDN7, VIM/CDH1) ordinal logistic model that scores any
   sample's EMT status as a continuous mu in [0, 2] and bins it into
   E / hybrid E/M / M.
2. This repo's own Panel A PCHA fit on the 63 DepMap invasive breast
   cell lines - each line's convex-combination weight over the k=4
   archetypes (task specialists), already saved as S_k4_parti.npy.

Step (2)'s weights give a standard "generalist score" (normalized
Shannon entropy of the weight row: 0 = pure specialist at one
archetype, 1 = uniform generalist across all archetypes - the
specialist/generalist readout used in Hausser 2019 / Groves 2022).

The question: do hybrid E/M lines (mu in the paper's [0.5, 1.5] band)
have a significantly different generalist score than E lines and than
M lines? We test this with two independent two-sample
Kolmogorov-Smirnov tests:

    KS-epi = KS_2samp(generalist_score[E],  generalist_score[E/M])
    KS-mes = KS_2samp(generalist_score[M],  generalist_score[E/M])

If both are significant and the hybrid group's median is higher, that
supports "hybrid E/M cells are archetype generalists, not specialists
at either the epithelial or mesenchymal task" - directly testable with
what is already in this repo, no refit needed.

Does NOT touch Panel A/B/C outputs. Writes only to results/panel_d_emt/
and figures/Figure_1D_breast_emt_generalist.png.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from src.io import load_expression_csv
from src.emt_score import score_cell_lines, sanity_check_against_labels
from src.enrichment import (
    distance_bins,
    archetype_weight_entropy,
    ks_epi_mes,
)

MATRIX = BREAST / "data" / "processed" / "input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv"
PC_SCORES = BREAST / "results" / "panel_a" / "pc_scores.npy"
SUGGESTED_K = BREAST / "results" / "panel_a" / "suggested_k.txt"
PAM50_LABELS = BREAST / "results" / "panel_b" / "pam50_labels_panelA.csv"
PANEL_B_ENRICHMENT = BREAST / "results" / "panel_b" / "enrichment_pam50.csv"
OUT = BREAST / "results" / "panel_d_emt"
FIG = BREAST / "figures"

N_BINS = 5  # match Panel B, for the archetype-identity cross-reference only
CATEGORY_COLORS = {"E": "#4C78A8", "E/M": "#F58518", "M": "#E45756"}

# Independent biological prior used ONLY to pick the model's sign
# convention (see src/emt_score.py docstring): Basal/TNBC breast lines
# are well established to be more mesenchymal/EMT-active than
# Luminal-B or Her2 lines. This is not circular with the archetype
# analysis below - PAM50 calls do not depend on CLDN7/CDH1/VIM at all.
SIGN_CHECK_MESENCHYMAL_LABEL = "Basal"
SIGN_CHECK_EPITHELIAL_LABELS = ["LumB", "Her2"]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    # ---- 1. EMT score (George et al. 2017) for all 63 cell lines ----
    expr = load_expression_csv(MATRIX)
    sample_ids = list(expr.columns.astype(str))
    scored = score_cell_lines(expr)
    scored = scored.reindex(sample_ids)

    labels = pd.read_csv(PAM50_LABELS)
    labels["cell_line"] = labels["cell_line"].astype(str)
    labels = labels.set_index("cell_line")["pam50_subtype"].reindex(sample_ids)

    check = sanity_check_against_labels(
        scored, labels,
        mesenchymal_label=SIGN_CHECK_MESENCHYMAL_LABEL,
        epithelial_labels=SIGN_CHECK_EPITHELIAL_LABELS,
    )
    print("=== Sign-convention sanity check (George et al. Eq. 2) ===")
    for sign in ("paper", "flipped"):
        r = check[sign]
        print(
            f"  {sign:8s}: mean mu[{SIGN_CHECK_MESENCHYMAL_LABEL}]="
            f"{r['mesenchymal_label_mean_mu']:.3f}  "
            f"mean mu[{'+'.join(SIGN_CHECK_EPITHELIAL_LABELS)}]="
            f"{r['epithelial_labels_mean_mu']:.3f}  "
            f"diff={r['difference_mes_minus_epi']:+.3f}  "
            f"agrees_with_biology={r['agrees_with_biology']}"
        )
    chosen = check["chosen_sign"]
    if chosen is None:
        raise RuntimeError(
            "Neither sign convention gives Basal a higher mean mu than "
            "LumB/Her2. Do not trust mu until this is resolved by hand "
            "against the original paper's fitted model or NCI-60 data."
        )
    print(f"  -> using sign='{chosen}' for all downstream results.\n")

    mu = scored[f"mu_{chosen}"].rename("mu")
    category = scored[f"category_{chosen}"].rename("category")
    tidy_emt = pd.concat(
        [scored[["CLDN7", "CDH1", "VIM", "X1", "X2"]], mu, category, labels], axis=1
    )
    tidy_emt.index.name = "cell_line"
    tidy_emt.to_csv(OUT / "emt_scores.csv")
    print("EMT category counts (chosen sign):")
    print(category.value_counts().to_string())
    print()

    # ---- 2. Archetype generalist score (existing Panel A weights) ----
    k_star = int(SUGGESTED_K.read_text().splitlines()[0])
    weights = np.load(BREAST / "results" / "panel_a" / f"S_k{k_star}_parti.npy")
    arcs = np.load(BREAST / "results" / "panel_a" / f"archetypes_k{k_star}_parti.npy")
    pc_scores = np.load(PC_SCORES)
    if weights.shape[0] != len(sample_ids):
        raise ValueError(
            f"S weights n={weights.shape[0]} vs matrix samples n={len(sample_ids)}"
        )
    generalist = pd.Series(
        archetype_weight_entropy(weights), index=sample_ids, name="generalist_score"
    )
    generalist_df = pd.DataFrame(
        weights, index=sample_ids, columns=[f"weight_arc{i+1}" for i in range(k_star)]
    )
    generalist_df["generalist_score"] = generalist
    generalist_df.to_csv(OUT / "archetype_weights_and_generalist_score.csv")
    print(f"Generalist score computed from S_k{k_star}_parti.npy "
          f"(0=specialist at one archetype, 1=uniform generalist).\n")

    # ---- 3. KS-epi / KS-mes ----
    ks_table = ks_epi_mes(generalist, category, hybrid_label="E/M", epi_label="E", mes_label="M")
    ks_table.to_csv(OUT / "ks_epi_mes.csv")
    print("=== KS-epi / KS-mes: is E/M's generalist-score distribution ===")
    print("=== different from E's (KS-epi) and from M's (KS-mes)?      ===")
    print(ks_table.to_string())
    print()
    for name in ("KS-epi", "KS-mes"):
        row = ks_table.loc[name]
        verdict = "SIGNIFICANT" if row["p_value"] < 0.05 else "not significant"
        print(
            f"  {name}: D={row['ks_statistic']:.3f}, p={row['p_value']:.4f} "
            f"({verdict}); median generalist score "
            f"{row['group_a']}={row['median_a']:.3f} vs "
            f"{row['group_b']}={row['median_b']:.3f} "
            f"({row['direction']})"
        )
    print()

    # ---- 4. Continuous robustness check: mu vs generalist score ----
    common = mu.index.intersection(generalist.index)
    rho, p_rho = spearmanr(mu.loc[common], generalist.loc[common])
    print(
        f"Robustness check (no thresholding): Spearman rho(mu, "
        f"generalist_score) = {rho:.3f}, p = {p_rho:.4f}"
    )
    print(
        "  If hybrid E/M cells are generalists, mu (peaking at 1 = "
        "maximal hybrid) and generalist_score should be positively, "
        "non-linearly related (both extremes of mu -> low generalist "
        "score); a plain linear Spearman rho is a coarse check only, "
        "see the figure for the actual (inverted-U) shape.\n"
    )

    # ---- 5. Which archetype looks epithelial vs mesenchymal? ----
    # Cross-reference with Panel B's own PAM50 x archetype enrichment,
    # purely descriptive (does not feed back into steps 1-4 above).
    if PANEL_B_ENRICHMENT.exists():
        panelb = pd.read_csv(PANEL_B_ENRICHMENT)
        hits = panelb[panelb["sig_peak_at_bin0"]]
        print("=== Cross-reference: Panel B PAM50-vs-archetype peaks ===")
        if hits.empty:
            print("  No significant bin-0 PAM50 peaks in Panel B; skipping.")
        else:
            for _, row in hits.iterrows():
                print(
                    f"  {row['subtype']} peaks at archetype "
                    f"{int(row['archetype']) + 1} (fold={row['fold_enrichment']:.2f}, "
                    f"q={row['q_value']:.3f})"
                )
            print(
                "  Basal is the field's standard proxy for high EMT/"
                "mesenchymal activity in breast cancer; LumB/Her2 for "
                "epithelial. Compare the archetype indices above to "
                "which archetype each hybrid-scoring cell line sits "
                "closest to, in distances_pam50.csv (Panel B) or by "
                "rerunning distance_bins() on this script's pc_scores."
            )
        print()

    # ---- 6. Figure ----
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4))
    ax = axes[0]
    for cat, color in CATEGORY_COLORS.items():
        mask = category == cat
        ax.scatter(
            mu[mask], generalist.reindex(mu.index)[mask],
            color=color, label=cat, s=45, edgecolor="white", linewidth=0.5,
        )
    ax.set_xlabel("EMT score $\\mu$ (George et al. 2017)")
    ax.set_ylabel("Generalist score (archetype-weight entropy)")
    ax.axvline(0.5, color="0.7", lw=1, ls="--")
    ax.axvline(1.5, color="0.7", lw=1, ls="--")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(f"$\\rho$={rho:.2f}, p={p_rho:.3f}")

    ax2 = axes[1]
    order = ["E", "E/M", "M"]
    data = [generalist.reindex(mu.index)[category == c].dropna().values for c in order]
    try:
        bp = ax2.boxplot(data, tick_labels=order, patch_artist=True)
    except TypeError:  # matplotlib < 3.9 used `labels`
        bp = ax2.boxplot(data, labels=order, patch_artist=True)
    for patch, c in zip(bp["boxes"], order):
        patch.set_facecolor(CATEGORY_COLORS[c])
        patch.set_alpha(0.6)
    ax2.set_ylabel("Generalist score")
    ax2.set_title("KS-epi p={:.3f}, KS-mes p={:.3f}".format(
        ks_table.loc["KS-epi", "p_value"], ks_table.loc["KS-mes", "p_value"]
    ))
    fig.suptitle(
        "Panel D — EMT score vs archetype-generalist score, "
        f"breast Panel A (k={k_star})",
        y=1.03,
    )
    fig.tight_layout()
    fig.savefig(FIG / "Figure_1D_breast_emt_generalist.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / "emt_generalist.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("Wrote", OUT / "emt_scores.csv")
    print("Wrote", OUT / "archetype_weights_and_generalist_score.csv")
    print("Wrote", OUT / "ks_epi_mes.csv")
    print("Wrote", FIG / "Figure_1D_breast_emt_generalist.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
