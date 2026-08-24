#!/usr/bin/env python3
"""Panel D (GBM analogue): are Proneural/Mesenchymal "hybrid" glioblastoma
cell lines archetype "generalists" rather than specialists?

This mirrors Breast Cancer/codes/05_emt_generalist_panelD.py in structure
and statistics, but uses a GBM-appropriate continuous phenotype axis
instead of the George et al. 2017 CLDN7/CDH1/VIM epithelial-EMT score.

WHY NOT REUSE THE BREAST EMT SCORE LITERALLY
----------------------------------------------
George et al.'s model needs CDH1 (E-cadherin, an epithelial adherens-
junction gene) as a live, varying predictor. Glioblastoma is not an
epithelial tumor (it arises from glial/neural lineage, no CDH1+
epithelium), so CDH1 sits near the assay floor in essentially every
GBM cell line while VIM sits near the ceiling in almost all of them
regardless of Proneural/Mesenchymal identity. Run literally on this
repo's GBM matrix (see docs/panel_d_pnmes_gbm.md for the numbers),
the model pins 53/54 lines at mu=2.00 (maximal "mesenchymal") with
essentially zero variance - i.e. it is not measuring GBM biology, it
is measuring the fact that glial cells lack E-cadherin.

GBM's own literature (and this repo's README, "Hypothesis-driven k=2
(PN-MES axis)") already identifies the correct antagonistic pair for
glioblastoma: Verhaak's Proneural vs Mesenchymal transcriptional
programs (Wang et al. 2017 marker sets), not an epithelial-cadherin
axis. So Panel D here uses the *existing* Panel B marker z-scores
(results/panel_b/verhaak_labels_panelA.csv, already computed by
02_assign_verhaak.py - not recomputed) as the continuous phenotype
axis, and asks the same generalist/specialist question as breast
Panel D, but framed as PN vs MES vs Hybrid.

Inputs (nothing recomputed):
- results/panel_b/verhaak_labels_panelA.csv  (score_Mesenchymal,
  score_Proneural per cell line - independent gene-marker z-scores,
  already computed, does not touch PCHA)
- results/panel_a/S_k{k}_parti.npy           (PCHA mixture weights,
  same k as Panel B uses, i.e. suggested_k.txt)

Outputs: results/panel_d_pnmes/*.csv, figures/Figure_1D_gbm_pnmes_generalist.png
Does not touch Panel A/B/C. Safe to rerun.
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
from scipy.stats import spearmanr

from src.enrichment import archetype_weight_entropy, ks_epi_mes

VERHAAK_LABELS = GBM / "results" / "panel_b" / "verhaak_labels_panelA.csv"
SUGGESTED_K = GBM / "results" / "panel_a" / "suggested_k.txt"
PANEL_A = GBM / "results" / "panel_a"
PANEL_B_ENRICHMENT = GBM / "results" / "panel_b" / "enrichment_verhaak.csv"
OUT = GBM / "results" / "panel_d_pnmes"
FIG = GBM / "figures"

CATEGORY_COLORS = {"PN": "#4C78A8", "Hybrid": "#F58518", "MES": "#E45756"}


def classify_tertiles(score):
    """PN / Hybrid / MES via data-driven tertiles of the continuous
    Mesenchymal-minus-Proneural marker z-score.

    Unlike George et al.'s mu (which has paper-defined 0.5/1.5
    thresholds from an NCI-60-trained ordinal model), there is no
    published absolute threshold for this marker-based score, so we
    use tertiles: bottom third = PN-leaning, middle third = Hybrid,
    top third = MES-leaning. This is a transparent, data-driven
    choice, not a literature-derived one - stated explicitly here so
    it is not mistaken for a validated clinical cutoff.
    """
    q1, q2 = score.quantile([1 / 3, 2 / 3])
    cat = pd.Series(
        np.where(score <= q1, "PN", np.where(score >= q2, "MES", "Hybrid")),
        index=score.index,
        name="category",
    )
    return cat, float(q1), float(q2)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    # ---- 1. PN-MES continuous score (already-computed marker z-scores) ----
    labels = pd.read_csv(VERHAAK_LABELS)
    labels["cell_line"] = labels["cell_line"].astype(str)
    labels = labels.set_index("cell_line")
    sample_ids = list(labels.index)

    pnmes_score = (labels["score_Mesenchymal"] - labels["score_Proneural"]).rename("pnmes_score")
    category, q1, q2 = classify_tertiles(pnmes_score)
    print("=== PN-MES continuous score (score_Mesenchymal - score_Proneural) ===")
    print(f"  tertile cut points: q1={q1:.3f}, q2={q2:.3f}")
    print("Category counts (tertile split):")
    print(category.value_counts().to_string())
    print()

    tidy = pd.concat(
        [labels[["verhaak_subtype", "score_Classical", "score_Mesenchymal", "score_Proneural"]],
         pnmes_score, category], axis=1
    )
    tidy.index.name = "cell_line"
    tidy.to_csv(OUT / "pnmes_scores.csv")

    # Cross-tab against the existing argmax Verhaak call, purely descriptive
    print("Tertile category vs argmax verhaak_subtype (existing Panel B label):")
    print(pd.crosstab(category, labels["verhaak_subtype"]).to_string())
    print()

    # ---- 2. Archetype generalist score (existing Panel A weights, same k as Panel B) ----
    k_star = int(SUGGESTED_K.read_text().splitlines()[0])
    weights = np.load(PANEL_A / f"S_k{k_star}_parti.npy")
    if weights.shape[0] != len(sample_ids):
        raise ValueError(f"S weights n={weights.shape[0]} vs labels n={len(sample_ids)}")
    generalist = pd.Series(
        archetype_weight_entropy(weights), index=sample_ids, name="generalist_score"
    )
    generalist_df = pd.DataFrame(
        weights, index=sample_ids, columns=[f"weight_arc{i+1}" for i in range(k_star)]
    )
    generalist_df["generalist_score"] = generalist
    generalist_df.to_csv(OUT / "archetype_weights_and_generalist_score.csv")
    print(f"Generalist score computed from S_k{k_star}_parti.npy "
          f"(0=specialist at one archetype, 1=uniform generalist).")
    print(f"NOTE: k={k_star} is Panel A's *elbow-selected* k, not a statistically "
          f"significant one - t_ratio_parti_500.csv shows p>0.05 for every k tested "
          f"(k=3..7; best p=0.068 at k=7). Panel B's own subtype-vs-archetype "
          f"enrichment also found zero significant peaks at this k. Treat the "
          f"'specialist/generalist' framing below as exploratory, more so than "
          f"in the breast analysis (which was at least borderline, p=0.082).\n")

    # ---- 3. KS-PN / KS-MES (reuses the same generic function as breast) ----
    ks_table = ks_epi_mes(generalist, category, hybrid_label="Hybrid", epi_label="PN", mes_label="MES")
    ks_table = ks_table.rename(index={"KS-epi": "KS-PN", "KS-mes": "KS-MES"})
    ks_table.to_csv(OUT / "ks_pn_mes.csv")
    print("=== KS-PN / KS-MES: is Hybrid's generalist-score distribution ===")
    print("=== different from PN's (KS-PN) and from MES's (KS-MES)?      ===")
    print(ks_table.to_string())
    print()
    for name in ("KS-PN", "KS-MES"):
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

    # ---- 4. Continuous robustness check ----
    common = pnmes_score.index.intersection(generalist.index)
    rho, p_rho = spearmanr(pnmes_score.loc[common], generalist.loc[common])
    print(
        f"Robustness check (no thresholding): Spearman rho(pnmes_score, "
        f"generalist_score) = {rho:.3f}, p = {p_rho:.4f}\n"
    )

    # ---- 5. Cross-reference Panel B enrichment (descriptive only) ----
    if PANEL_B_ENRICHMENT.exists():
        panelb = pd.read_csv(PANEL_B_ENRICHMENT)
        hits = panelb[panelb["sig_peak_at_bin0"]] if "sig_peak_at_bin0" in panelb.columns else panelb.iloc[0:0]
        print("=== Cross-reference: Panel B Verhaak-vs-archetype peaks ===")
        if hits.empty:
            print("  No significant bin-0 Verhaak peaks in Panel B at this k - "
                  "no archetype can be confidently called 'the PN vertex' or "
                  "'the MES vertex' the way breast Panel B could call an "
                  "archetype Basal-leaning.")
        print()

    # ---- 6. Figure ----
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4))
    ax = axes[0]
    for cat, color in CATEGORY_COLORS.items():
        mask = category == cat
        ax.scatter(
            pnmes_score[mask], generalist.reindex(pnmes_score.index)[mask],
            color=color, label=cat, s=45, edgecolor="white", linewidth=0.5,
        )
    ax.set_xlabel("PN-MES marker score (Mesenchymal z - Proneural z)")
    ax.set_ylabel("Generalist score (archetype-weight entropy)")
    ax.axvline(q1, color="0.7", lw=1, ls="--")
    ax.axvline(q2, color="0.7", lw=1, ls="--")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title(f"$\\rho$={rho:.2f}, p={p_rho:.3f}")

    ax2 = axes[1]
    order = ["PN", "Hybrid", "MES"]
    data = [generalist.reindex(pnmes_score.index)[category == c].dropna().values for c in order]
    try:
        bp = ax2.boxplot(data, tick_labels=order, patch_artist=True)
    except TypeError:
        bp = ax2.boxplot(data, labels=order, patch_artist=True)
    for patch, c in zip(bp["boxes"], order):
        patch.set_facecolor(CATEGORY_COLORS[c])
        patch.set_alpha(0.6)
    ax2.set_ylabel("Generalist score")
    ax2.set_title("KS-PN p={:.3f}, KS-MES p={:.3f}".format(
        ks_table.loc["KS-PN", "p_value"], ks_table.loc["KS-MES", "p_value"]
    ))
    fig.suptitle(
        f"Panel D — PN-MES score vs archetype-generalist score, GBM Panel A (k={k_star})",
        y=1.03,
    )
    fig.tight_layout()
    fig.savefig(FIG / "Figure_1D_gbm_pnmes_generalist.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / "pnmes_generalist.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("Wrote", OUT / "pnmes_scores.csv")
    print("Wrote", OUT / "archetype_weights_and_generalist_score.csv")
    print("Wrote", OUT / "ks_pn_mes.csv")
    print("Wrote", FIG / "Figure_1D_gbm_pnmes_generalist.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
