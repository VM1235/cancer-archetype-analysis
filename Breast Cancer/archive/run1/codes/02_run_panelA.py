#!/usr/bin/env python3
"""Panel A for invasive breast carcinoma CCLE lines (ParTI-matched PCHA).

Same call as the SCLC reproduction:
  - PCA on samples × genes; keep enough PCs to reach ~50% variance
    (Groves used 12 PCs for SCLC because that was ~50%)
  - ESV vs k in that PC space (delta=0)
  - t-ratio: PCHA on first (k-1) PCs, delta=0, 150 inits, max volume
  - permutation: 100 shuffles × 15 inits (same sanity scale as SCLC)

Outputs go in Breast Cancer/results/panel_a/
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
from scipy.spatial import ConvexHull, QhullError

from src.io import load_expression_csv
from src.pca import fit_pca, cumulative_variance, inverse_transform_scores
from src.archetypes import (
    _shuffle_columns,
    esv_curve,
    fit_pcha_best,
    simplex_volume,
)

MATRIX = BREAST / "data" / "processed" / "input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv"
OUT = RUN1 / "results" / "panel_a"
FIG = RUN1 / "figures"

TARGET_CUMVAR = 0.50
K_ESV = range(2, 16)
K_PERM = (3, 4, 5, 6)
N_INIT_OBS = 150
N_INIT_NULL = 15
N_PERM = 100
DELTA = 0.0
SEED = 0


def choose_n_pcs(cumvar, target=TARGET_CUMVAR):
    """Smallest n with cumulative variance >= target; else all computed PCs."""
    hit = np.where(cumvar >= target)[0]
    if len(hit):
        return int(hit[0] + 1)
    return int(len(cumvar))


def hull_t_ratio(scores, archetypes):
    k = archetypes.shape[0]
    data = np.asarray(scores, dtype=float)[:, : k - 1]
    return float(simplex_volume(archetypes) / ConvexHull(data).volume)


def style():
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


def suggested_k(tab, k_esv, delta_esv):
    sig = tab.loc[tab["p_value"] < 0.05].sort_values("k")
    if len(sig):
        return int(sig["k"].iloc[0]), "smallest k with p < 0.05"
    k_esv = np.asarray(k_esv)
    delta_esv = np.asarray(delta_esv, dtype=float)
    if len(k_esv) >= 3:
        # largest drop in extra ESV after k=3
        drops = -np.diff(delta_esv)
        idx = int(np.argmax(drops)) + 1
        return int(k_esv[idx]), "ESV elbow (no k with p<0.05)"
    return 5, "default"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    if not MATRIX.exists():
        print("Missing", MATRIX)
        print("Run 01_build_input_panelA.py first.")
        return 1

    expr = load_expression_csv(MATRIX)
    X = expr.T.values
    n_samples, n_genes = X.shape
    print(f"Loaded {n_genes} genes × {n_samples} cell lines")

    n_fit = min(n_samples - 1, n_genes, 40)
    pca_all, scores_all = fit_pca(X, n_components=n_fit)
    cum_all = cumulative_variance(pca_all)
    n_pcs = choose_n_pcs(cum_all)
    pca, scores = fit_pca(X, n_components=n_pcs)
    cumvar = cumulative_variance(pca)
    print(
        f"PCA: {n_pcs} PCs explain {100 * cumvar[-1]:.1f}% variance "
        f"(target ~{100 * TARGET_CUMVAR:.0f}%; SCLC used 12 PCs for ~50%)"
    )
    np.save(OUT / "pc_scores.npy", scores)
    pd.DataFrame(
        {
            "pc": np.arange(1, n_pcs + 1),
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "cumulative": cumvar,
        }
    ).to_csv(OUT / "pca_variance.csv", index=False)
    (OUT / "n_pcs.txt").write_text(str(n_pcs) + "\n")

    print(f"ESV curve k=2..15 in {n_pcs}-D (delta={DELTA}) ...")
    fits = esv_curve(scores, K_ESV, delta=DELTA, seed=SEED)
    esv = pd.DataFrame([{"k": f["k"], "esv": f["esv"]} for f in fits])
    esv["delta_esv"] = esv["esv"].diff()
    esv.to_csv(OUT / "esv_curve.csv", index=False)
    print(esv.to_string(index=False))

    gene_esv = 100.0 * esv["esv"].values * float(cumvar[-1])
    k_esv = esv["k"].values
    delta_esv = np.diff(np.concatenate([[0.0], gene_esv]))

    rows = []
    for k in K_PERM:
        print(
            f"\n=== k={k}: {N_INIT_OBS} observed inits (delta={DELTA}, {k-1} PCs) ===",
            flush=True,
        )
        archetypes, weights, varexpl, vol, n_ok = fit_pcha_best(
            scores, k, n_init=N_INIT_OBS, delta=DELTA
        )
        observed = hull_t_ratio(scores, archetypes)
        np.save(OUT / f"archetypes_k{k}_parti.npy", archetypes)
        np.save(OUT / f"S_k{k}_parti.npy", weights)
        print(
            f"  inits_ok={n_ok}/{N_INIT_OBS}  ESV={varexpl:.3f}  t={observed:.4f}",
            flush=True,
        )

        print(f"  permutation: {N_PERM} shuffles x {N_INIT_NULL} inits ...", flush=True)
        rng = np.random.default_rng(SEED + k)
        null = []
        n_fail = 0
        for i in range(N_PERM):
            shuffled = _shuffle_columns(scores[:, : k - 1], rng)
            try:
                arch, _, _, _, _ = fit_pcha_best(
                    shuffled, k, n_init=N_INIT_NULL, delta=DELTA
                )
                null.append(hull_t_ratio(shuffled, arch))
            except (QhullError, ValueError, RuntimeError):
                n_fail += 1
            if (i + 1) % 20 == 0:
                print(f"    {i + 1}/{N_PERM}", flush=True)
        null = np.asarray(null, dtype=float)
        null = null[np.isfinite(null)]
        p_value = float(np.mean(null >= observed)) if len(null) else np.nan
        np.save(OUT / f"null_t_ratios_k{k}_parti_n{N_PERM}.npy", null)
        rows.append(
            {
                "k": k,
                "t_ratio": observed,
                "p_value": p_value,
                "n_perm": N_PERM,
                "n_init_obs": N_INIT_OBS,
                "n_init_null": N_INIT_NULL,
                "n_success": int(len(null)),
                "n_fail": n_fail,
                "n_pcs_fit": n_pcs,
            }
        )
        print(
            f"  k={k}: t={observed:.4f}, p={p_value:.3f} from {len(null)}/{N_PERM}",
            flush=True,
        )

    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "t_ratio_sanity_parti.csv", index=False)
    print("\nBreast Panel A t-ratio table:")
    print(tab.to_string(index=False))

    k_star, why = suggested_k(tab, k_esv, delta_esv)
    print(f"\nSuggested k={k_star} ({why})")
    (OUT / "suggested_k.txt").write_text(f"{k_star}\n{why}\n")

    pca2, scores2 = fit_pca(X, n_components=2)
    arc_path = OUT / f"archetypes_k{k_star}_parti.npy"
    if not arc_path.exists():
        archetypes, _, _, _, _ = fit_pcha_best(scores, k_star, n_init=N_INIT_OBS, delta=DELTA)
        np.save(arc_path, archetypes)
    arcs = np.load(arc_path)
    arcs_full = np.zeros((arcs.shape[0], n_pcs))
    arcs_full[:, : arcs.shape[1]] = arcs
    gene_arcs = inverse_transform_scores(pca, arcs_full)
    arcs2 = pca2.transform(gene_arcs)

    style()
    fig = plt.figure(figsize=(6.6, 7.4))
    from matplotlib.gridspec import GridSpec

    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 1.55], hspace=0.32, wspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(k_esv, delta_esv, "-o", color="#3B6FA0", ms=5, lw=1.4)
    ax.axvline(k_star, color="0.35", ls="--", lw=1)
    ax.set_xlim(1.5, 15.5)
    ax.set_ylim(0, max(12, float(delta_esv.max()) * 1.08))
    ax.set_xticks(range(2, 16))
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("% ESV on top of N−1 model")
    ax.set_title("Explained sample variance (ESV)")

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(tab["k"], tab["t_ratio"], "-o", color="#3B6FA0", ms=6, lw=1.4)
    ax.axvline(k_star, color="0.35", ls="--", lw=1)
    ymax = max(0.62, float(tab["t_ratio"].max()) * 1.15)
    ax.set_ylim(0, ymax)
    ax.set_xlim(2.5, 6.5)
    ax.set_xticks(tab["k"].astype(int))
    for _, row in tab.iterrows():
        p = row["p_value"]
        label = f"p = {p:.2f}" if p >= 0.01 else f"p = {p:.3f}"
        ax.annotate(
            label,
            (row["k"], row["t_ratio"]),
            textcoords="offset points",
            xytext=(8, 6),
            fontsize=8,
            color="#222",
        )
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("t-ratio")
    ax.set_title("t-ratio of polytopes by number of vertices")

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
        arcs2[:, 0],
        arcs2[:, 1],
        s=90,
        c="#3B6FA0",
        edgecolors="k",
        linewidths=0.4,
        zorder=3,
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
    ax.set_title(f"Archetype space ({n_samples} invasive breast carcinoma cell lines)")

    fig.suptitle(
        "Figure 1A  —  Archetype analysis on invasive breast carcinoma CCLE lines",
        y=0.995,
        fontsize=11,
    )
    fig.savefig(FIG / "Figure_1A_breast.png")
    plt.close(fig)
    print("Wrote", FIG / "Figure_1A_breast.png")
    print(f"Wrote outputs to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
