#!/usr/bin/env python3
"""Fit KS Panel A-style archetypes DIRECTLY on GSE161529 (reverse direction).

Every other script in this project that calls fit_pcha_best/esv_curve does so
on the 63 DepMap cell lines only -- GSE176078/GSE173634/GSE161529 have only
ever been PROJECTED onto that already-fit polytope (see
run_panelA_ks_project_gse*_sc*.py, which just load archetypes_k{k}_parti.npy
from panel_a_ks_genelist/ or _extendedk/ and never refit). This script runs
the fitting machinery itself, for the first time, on single-cell data.

Deliberately holds back the 500-shuffle permutation/t-ratio significance test
(see run_panelA_ks_genelist_extendedk.py for that block): this project has
never run PCHA at ~85k-sample scale before (vs. 63 lines everywhere else), so
runtime here is genuinely unknown. This script gets the ESV elbow curve and
point-estimate archetype coordinates first; decide on the permutation test
after seeing how long this takes.

Reuses the MAGIC-imputed KS-gene cache already built by
1run_panelA_ks_project_gse161529_sc.py -- no need to recompute MAGIC.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/run_panelA_ks_genelist_gse161529.py"
"""

from __future__ import annotations

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.pca import fit_pca, cumulative_variance
from src.archetypes import esv_curve, fit_pcha_best

# Output of 1run_panelA_ks_project_gse161529_sc.py (RUN_MAGIC=True) -- the
# already-QC'd, already-MAGIC-imputed 206 KS genes x 84,602 epithelial cells
# (ER-0001 already excluded there).
MAGIC_CACHE = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse161529_magic" / "ks_genes_magic_imputed.npz"
)
OUT = BREAST / "results" / "panel_a_ks_genelist_gse161529"
FIG = BREAST / "figures"

TARGET_CUMVAR = 0.50
K_FIT = (3, 4, 5, 6, 7)
NUM_ITER = 5  # ParTI.m algNum==5, same protocol as the cell-line runs
N_INIT_OBS = 3 * NUM_ITER
DELTA = 0.0
SEED = 0
# Set an integer (e.g. 10000) to stratified*-subsample cells before fitting
# if full-scale turns out to be too slow. *not actually stratified yet --
# plain random subsample; upgrade this if we end up needing it.
MAX_CELLS = None


def choose_n_pcs(cumvar, target=TARGET_CUMVAR):
    hit = np.where(cumvar >= target)[0]
    if len(hit):
        return int(hit[0] + 1)
    return int(len(cumvar))


def dimension_finder(esv):
    """ParTI DimensionFinder.m: farthest point from the chord of the ESV curve."""
    esv = np.asarray(esv, dtype=float)
    n = len(esv)
    slope = (esv[-1] - esv[0]) / (n - 1)
    intercept = esv[0] - slope * 1.0
    di = np.empty(n)
    for i in range(1, n + 1):
        s = -1.0 / slope
        inter2 = esv[i - 1] - s * i
        x = (inter2 - intercept) / (slope - s)
        y = s * x + inter2
        di[i - 1] = np.hypot(x - i, y - esv[i - 1])
    return int(np.argmax(di) + 2)


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


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    if not MAGIC_CACHE.is_file():
        print(f"Missing {MAGIC_CACHE}")
        print("Run 1run_panelA_ks_project_gse161529_sc.py (RUN_MAGIC=True) first.")
        return 1

    print(f"Loading cached MAGIC-imputed KS genes: {MAGIC_CACHE}")
    z = np.load(MAGIC_CACHE, allow_pickle=True)
    genes = list(z["genes"])
    barcodes = list(z["barcodes"])
    expr = pd.DataFrame(z["imputed"], index=genes, columns=barcodes)  # genes x cells
    print(f"Loaded {expr.shape[0]} genes x {expr.shape[1]} cells")

    rng = np.random.default_rng(SEED)
    if MAX_CELLS is not None and expr.shape[1] > MAX_CELLS:
        keep = np.sort(rng.choice(expr.shape[1], size=MAX_CELLS, replace=False))
        expr = expr.iloc[:, keep]
        print(f"Subsampled to {expr.shape[1]} cells (MAX_CELLS={MAX_CELLS})")

    barcodes_used = expr.columns.tolist()
    X = expr.T.values  # cells (samples) x genes -- orientation fit_pca expects
    n_samples, n_genes = X.shape
    print(f"GSE161529 KS-gene matrix: {n_genes} genes x {n_samples} cells")
    print(
        f"ParTI PCHA (POINT ESTIMATE ONLY -- no permutation test): "
        f"numIter={NUM_ITER}, n_init={N_INIT_OBS}, delta={DELTA}"
    )
    print("First time this fitting machinery has run on single-cell-scale data in this project.")

    n_fit = min(n_samples - 1, n_genes, 40)
    pca_all, _ = fit_pca(X, n_components=n_fit)
    cum_all = cumulative_variance(pca_all)
    n_pcs = max(choose_n_pcs(cum_all), max(K_FIT) - 1)  # force enough PCs for k up to 7
    pca, scores = fit_pca(X, n_components=n_pcs)
    cumvar = cumulative_variance(pca)
    print(
        f"dim={n_pcs} PCs explain {100 * cumvar[-1]:.1f}% variance "
        f"(cell-line Panel A needed only dim=2; forced here to cover k up to {max(K_FIT)})"
    )

    pd.DataFrame(
        {
            "pc": np.arange(1, n_fit + 1),
            "explained_variance_ratio": pca_all.explained_variance_ratio_,
            "cumulative": cum_all,
        }
    ).to_csv(OUT / "pca_variance_full.csv", index=False)
    np.save(OUT / "pc_scores.npy", scores)
    pd.Series(barcodes_used).to_csv(OUT / "pc_scores_barcodes.csv", index=False, header=["barcode"])
    (OUT / "n_pcs.txt").write_text(str(n_pcs) + "\n")

    print(f"ESV curve k=2..{n_pcs + 1} in {n_pcs}-D (delta={DELTA}) ...", flush=True)
    k_esv_vals = list(range(2, n_pcs + 2))
    fits = esv_curve(scores, k_esv_vals, delta=DELTA, seed=SEED)
    esv = pd.DataFrame([{"k": f["k"], "esv": f["esv"]} for f in fits])
    esv["delta_esv"] = esv["esv"].diff()
    esv.to_csv(OUT / "esv_curve.csv", index=False)
    print(esv.to_string(index=False))

    tot_esv = esv["esv"].values * float(cumvar[-1])
    k_from_finder = dimension_finder(tot_esv)
    print(f"ParTI DimensionFinder suggests k={k_from_finder} (elbow only -- NOT significance-tested yet)")
    (OUT / "k_from_finder.txt").write_text(str(k_from_finder) + "\n")

    gene_esv = 100.0 * tot_esv
    k_esv = esv["k"].values
    delta_esv = np.diff(np.concatenate([[0.0], gene_esv]))

    print("\nFitting point-estimate archetypes at each k (no permutation test):")
    for k in K_FIT:
        if k - 1 > n_pcs:
            print(f"  skip k={k}: needs {k - 1} PCs, only have {n_pcs}")
            continue
        print(f"  k={k}: {N_INIT_OBS} inits, {k - 1} PCs ...", flush=True)
        archetypes, weights, varexpl, vol, n_ok = fit_pcha_best(
            scores, k, n_init=N_INIT_OBS, delta=DELTA
        )
        np.save(OUT / f"archetypes_k{k}_parti.npy", archetypes)
        np.save(OUT / f"S_k{k}_parti.npy", weights)
        print(f"    inits_ok={n_ok}/{N_INIT_OBS}  ESV={varexpl:.3f}  vol={vol:.4g}", flush=True)

    print(
        f"\nDone -- point estimates + elbow only, no p-values yet. "
        f"Elbow suggests k={k_from_finder}; decide on the permutation test next."
    )

    k_plot = min(k_from_finder, max(K_FIT))
    arcs_path = OUT / f"archetypes_k{k_plot}_parti.npy"
    style()
    fig = plt.figure(figsize=(11, 4.6))
    gs = GridSpec(1, 2, figure=fig, wspace=0.32)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(k_esv, delta_esv, "-o", color="#3B6FA0", ms=5, lw=1.4)
    ax.axvline(k_from_finder, color="0.35", ls="--", lw=1)
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("% ESV on top of N-1 model")
    ax.set_title("GSE161529: explained sample variance (ESV)")

    ax = fig.add_subplot(gs[0, 1])
    take = rng.choice(scores.shape[0], size=min(20000, scores.shape[0]), replace=False)
    ax.scatter(scores[take, 0], scores[take, 1], s=3, c="#B0B0B0", alpha=0.3, linewidths=0, zorder=1)
    if arcs_path.is_file():
        arcs2 = np.load(arcs_path)[:, :2]
        pts = list(arcs2) + [arcs2[0]]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#333333", lw=1.3, zorder=3)
        ax.scatter(arcs2[:, 0], arcs2[:, 1], s=90, c="#F58518", edgecolors="k", linewidths=0.5, zorder=4)
        for i in range(arcs2.shape[0]):
            ax.annotate(
                str(i + 1), (arcs2[i, 0], arcs2[i, 1]), textcoords="offset points",
                xytext=(6, 6), fontsize=8, fontweight="bold",
            )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"GSE161529-fit archetypes, k={k_plot} (point estimate, provisional)")

    fig.suptitle(
        "GSE161529 self-fit archetypes -- ESV elbow + point estimate only (no significance test yet)",
        fontsize=11,
        y=1.03,
    )
    fig.savefig(FIG / "Figure_1A_gse161529_provisional.png")
    plt.close(fig)
    print("Wrote", FIG / "Figure_1A_gse161529_provisional.png")
    print("Wrote outputs to", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
