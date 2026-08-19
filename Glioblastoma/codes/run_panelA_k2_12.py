#!/usr/bin/env python3
"""Hypothesis-driven k=2 PCHA in the full 12-PC space (not k−1=1).

The ParTI (k−1)-PC rule is for a k-simplex that lives in (k−1)-D. At k=2
that forces a 1-D fit, where delta=0 PCHA sits on the hull endpoints and
t-ratio is identically 1. Here we relax the rule: fit two archetypes in the
same 12-PC matrix used for GBM k=3–7, keep the max-length pair, and define

    t = ||a1 − a2|| / extent of the data along that axis

i.e. 1-D ParTI volumes in the affine span of the pair after a 12-D fit.
Shuffles still permute each PC independently (Groves). 150/50 inits, 500
shuffles, delta=0.

Writes results/panel_a_k2_12/ and figures/Figure_1A_gbm_k2_12.png.
Does not overwrite panel_a/, panel_a_k2/, or Figure_1A_gbm.png.
"""

from pathlib import Path
import argparse
import os
import sys
import time

HERE = Path(__file__).resolve().parent
GBM = HERE.parent
ROOT = GBM.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from src.io import load_expression_csv
from src.pca import align_pca_signs, cumulative_variance, fit_pca, inverse_transform_scores
from src.archetypes import _shuffle_columns, fit_pcha

MATRIX = GBM / "data" / "processed" / "input_panelA_glioblastoma_ccle_logtpm_filtered.csv"
SCORES = GBM / "results" / "panel_a" / "pc_scores.npy"
N_PCS_PATH = GBM / "results" / "panel_a" / "n_pcs.txt"
ESV_PATH = GBM / "results" / "panel_a" / "esv_curve.csv"
OUT = GBM / "results" / "panel_a_k2_12"
FIG = GBM / "figures"

K = 2
NUM_ITER = 50
N_INIT_OBS = 3 * NUM_ITER
N_INIT_NULL = NUM_ITER
N_PERM = 500
N_CHECKPOINT = 25
DELTA = 0.0
SEED = 0


def pair_length(archetypes):
    a = np.asarray(archetypes, dtype=float)
    return float(np.linalg.norm(a[0] - a[1]))


def axis_extent(scores, archetypes):
    """Length of the data cloud along the line through the two archetypes."""
    a = np.asarray(archetypes, dtype=float)
    x = np.asarray(scores, dtype=float)
    direction = a[1] - a[0]
    norm = float(np.linalg.norm(direction))
    if norm <= 0:
        return np.nan
    u = direction / norm
    proj = x @ u
    return float(np.max(proj) - np.min(proj))


def t_ratio_pair(scores, archetypes):
    extent = axis_extent(scores, archetypes)
    if not np.isfinite(extent) or extent <= 0:
        return np.nan
    return pair_length(archetypes) / extent


def fit_pcha_best_full(pc_scores, k, n_init, delta=DELTA, progress_every=0):
    """PCHA on all provided PCs (do not slice to k−1). Keep max pair length."""
    scores = np.asarray(pc_scores, dtype=float)
    best = None
    best_vol = -np.inf
    n_ok = 0
    t0 = time.time()
    for i in range(int(n_init)):
        try:
            archetypes, weights, varexpl = fit_pcha(scores, k, delta=delta)
            vol = pair_length(archetypes)
            n_ok += 1
            if np.isfinite(vol) and vol > best_vol:
                best_vol = vol
                best = (archetypes, weights, varexpl, vol)
        except RuntimeError:
            continue
        if progress_every and (i + 1) % progress_every == 0:
            print(
                f"    observed inits {i + 1}/{n_init}  n_ok={n_ok}  "
                f"{time.time() - t0:.1f}s",
                flush=True,
            )
    if best is None:
        raise RuntimeError(f"PCHA failed all {n_init} inits for k={k}")
    return best[0], best[1], best[2], best[3], n_ok


def p_label(p):
    if not np.isfinite(p):
        return "p = NA"
    if p == 0:
        return f"p < {1.0 / N_PERM:.3f}"
    if p < 0.01:
        return f"p = {p:.3f}"
    return f"p = {p:.2f}"


def permute_t_ratio(scores, archetypes, n_perm):
    observed = t_ratio_pair(scores, archetypes)
    null_path = OUT / f"null_t_ratios_k2_12pc_n{n_perm}.npy"
    if null_path.exists():
        null = [float(x) for x in np.load(null_path).ravel()]
        print(f"  resume {null_path.name}: {len(null)}/{n_perm}", flush=True)
    else:
        null = []
    n_fail = sum(1 for x in null if not np.isfinite(x))
    t0 = time.time()
    for i in range(len(null), n_perm):
        shuffled = _shuffle_columns(scores, np.random.default_rng([SEED, 12, int(i)]))
        try:
            arch, _, _, _, _ = fit_pcha_best_full(
                shuffled, K, n_init=N_INIT_NULL, delta=DELTA
            )
            value = t_ratio_pair(shuffled, arch)
            if not np.isfinite(value):
                n_fail += 1
                value = np.nan
        except (ValueError, RuntimeError):
            n_fail += 1
            value = np.nan
        null.append(value)
        if (i + 1) % N_CHECKPOINT == 0 or (i + 1) == n_perm:
            np.save(null_path, np.asarray(null, dtype=float))
            finite = np.asarray([x for x in null if np.isfinite(x)], dtype=float)
            running_p = float(np.mean(finite >= observed)) if len(finite) else np.nan
            print(
                f"    k=2 (12-PC)  {i + 1}/{n_perm}  n_ok={len(finite)}  "
                f"n_fail={n_fail}  running p={running_p:.4f}  {time.time() - t0:.0f}s",
                flush=True,
            )
    null_arr = np.asarray(null[:n_perm], dtype=float)
    np.save(null_path, null_arr)
    finite = null_arr[np.isfinite(null_arr)]
    p_value = float(np.mean(finite >= observed)) if len(finite) else np.nan
    return observed, p_value, int(len(finite)), int(n_fail), null_arr


def draw_figure(expr, pca12, n_pcs, archetypes, n_samples, t_row, esv, null_arr):
    X = expr.T.values
    tot_esv = esv["esv"].values * float(cumulative_variance(pca12)[-1])
    gene_esv = 100.0 * tot_esv
    k_esv = esv["k"].values
    delta_esv = np.diff(np.concatenate([[0.0], gene_esv]))

    pca2, scores2 = fit_pca(X, n_components=2)
    gene_arcs = inverse_transform_scores(pca12, archetypes)
    arcs2 = pca2.transform(gene_arcs)

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
    fig = plt.figure(figsize=(6.6, 7.4))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 1.55], hspace=0.32, wspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(k_esv, delta_esv, "-o", color="#3B6FA0", ms=5, lw=1.4)
    ax.axvline(K, color="0.35", ls="--", lw=1)
    ax.set_xlim(1.5, 15.5)
    ax.set_ylim(0, max(12, float(delta_esv.max()) * 1.08))
    ax.set_xticks(range(2, 16))
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("% ESV on top of N−1 model")
    ax.set_title("Explained sample variance (ESV)")

    ax = fig.add_subplot(gs[0, 1])
    finite = null_arr[np.isfinite(null_arr)]
    observed = float(t_row["t_ratio"])
    ax.hist(finite, bins=25, color="#C5D4E8", edgecolor="#3B6FA0", linewidth=0.6)
    ax.axvline(observed, color="#C44E52", lw=1.6, label=f"observed t={observed:.3f}")
    ax.set_xlabel("t-ratio (null shuffles)")
    ax.set_ylabel("count")
    ax.set_title(f"k=2 in 12-PC space  {p_label(float(t_row['p_value']))}")
    ax.legend(frameon=False, fontsize=7)

    ax = fig.add_subplot(gs[1, :])
    ax.scatter(scores2[:, 0], scores2[:, 1], s=18, c="#B0B0B0", zorder=1, linewidths=0)
    ax.plot(arcs2[:, 0], arcs2[:, 1], color="#888888", lw=0.9, zorder=2)
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
    ax.set_title(
        f"k=2 pair in 12-PC space, shown in PC1–PC2  "
        f"({n_samples} GBM cell lines)"
    )

    fig.suptitle(
        "Figure 1A (k=2, 12 PCs)  —  pair fit in full PCA space "
        "(not k−1=1; 500 shuffles, numIter=50)",
        y=0.995,
        fontsize=10,
    )
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "Figure_1A_gbm_k2_12.png")
    plt.close(fig)
    print("Wrote", FIG / "Figure_1A_gbm_k2_12.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-perm", type=int, default=N_PERM)
    parser.add_argument("--force-fit", action="store_true")
    args = parser.parse_args()

    for key in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[key] = "1"

    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    scores = np.load(SCORES)
    n_pcs = int(N_PCS_PATH.read_text().strip().splitlines()[0])
    expr = load_expression_csv(MATRIX)
    sample_ids = list(expr.columns.astype(str))
    if scores.shape != (len(sample_ids), n_pcs):
        raise ValueError(f"pc_scores {scores.shape} vs {len(sample_ids)} × {n_pcs}")

    pca12, rebuilt = fit_pca(expr.T.values, n_components=n_pcs)
    pca12, rebuilt, _ = align_pca_signs(rebuilt, scores, pca12)
    print(f"Reused {SCORES}; max |score diff| = {np.max(np.abs(rebuilt - scores)):.4g}")
    print(
        f"k=2 PCHA in full {n_pcs}-PC space (not first 1 PC). "
        f"numIter={NUM_ITER}, obs inits={N_INIT_OBS}, null inits={N_INIT_NULL}, "
        f"shuffles={int(args.n_perm)}, delta={DELTA}"
    )
    print(
        "t = pair length / data extent along that axis. "
        "This relaxes ParTI k−1 because k=2 in 1-D is degenerate."
    )

    (OUT / "METHOD.txt").write_text(
        "Hypothesis-driven k=2 fit in the full 12-PC cell-line space.\n"
        "Deviation from ParTI: do not restrict to (k-1)=1 PC. That convention\n"
        "is for k>=3 simplices; at k=2 it forces a 1-D hull-endpoint fit with\n"
        "t-ratio identically 1. Here PCHA uses all 12 PCs; t-ratio is the\n"
        "1-D volume ratio in the affine span of the two archetypes after that\n"
        "12-D fit. Inits/shuffles match GBM Panel A (150 / 50 x 500, delta=0).\n"
        "Not Groves-selected k. Does not overwrite panel_a or panel_a_k2.\n"
    )

    arch_path = OUT / "archetypes_k2_12pc.npy"
    s_path = OUT / "S_k2_12pc.npy"
    row_path = OUT / "t_ratio_k2_12pc.csv"
    n_perm = int(args.n_perm)

    if arch_path.exists() and s_path.exists() and not args.force_fit:
        archetypes = np.load(arch_path)
        weights = np.load(s_path)
        print(f"Reused {arch_path.name} {archetypes.shape}")
        if archetypes.shape != (K, n_pcs):
            archetypes = None
    else:
        archetypes = None

    if archetypes is None:
        print(f"Observed fit: {N_INIT_OBS} inits on {n_pcs} PCs, max pair length ...", flush=True)
        archetypes, weights, varexpl, vol, n_ok = fit_pcha_best_full(
            scores, K, n_init=N_INIT_OBS, delta=DELTA, progress_every=25
        )
        np.save(arch_path, archetypes)
        np.save(s_path, weights)
        print(
            f"  inits_ok={n_ok}/{N_INIT_OBS}  ESV={varexpl:.3f}  "
            f"length={vol:.4g}  t={t_ratio_pair(scores, archetypes):.4f}",
            flush=True,
        )

    observed = t_ratio_pair(scores, archetypes)
    print(f"  observed t={observed:.4f}", flush=True)

    null_path = OUT / f"null_t_ratios_k2_12pc_n{n_perm}.npy"
    need_perm = True
    t_row = None
    if row_path.exists() and null_path.exists() and np.load(null_path).size >= n_perm:
        prev = pd.read_csv(row_path)
        print(f"Already complete: {row_path.name}", flush=True)
        need_perm = False
        t_row = prev.iloc[0].to_dict()
        null_arr = np.load(null_path)

    if need_perm:
        print(f"Permutation: {n_perm} shuffles × {N_INIT_NULL} inits", flush=True)
        observed, p_value, n_ok, n_fail, null_arr = permute_t_ratio(
            scores, archetypes, n_perm
        )
        t_row = {
            "k": K,
            "n_pcs_fit": n_pcs,
            "t_ratio": observed,
            "p_value": p_value,
            "n_perm": n_perm,
            "n_init_obs": N_INIT_OBS,
            "n_init_null": N_INIT_NULL,
            "numIter": NUM_ITER,
            "n_success": n_ok,
            "n_fail": n_fail,
            "pair_length": pair_length(archetypes),
            "axis_extent": axis_extent(scores, archetypes),
            "note": "k=2 fit in 12 PCs; t = length / axis extent; not Groves k",
        }
        pd.DataFrame([t_row]).to_csv(row_path, index=False)
        print(
            f"  k=2 (12-PC) done: t={observed:.4f}, p={p_value:.4f} "
            f"from {n_ok}/{n_perm} (n_fail={n_fail})",
            flush=True,
        )

    esv = pd.read_csv(ESV_PATH)
    draw_figure(expr, pca12, n_pcs, archetypes, len(sample_ids), t_row, esv, null_arr)
    print("Wrote outputs to", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
