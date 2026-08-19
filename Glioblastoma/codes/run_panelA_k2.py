#!/usr/bin/env python3
"""Hypothesis-driven k=2 Panel A: t-ratio with Groves/ParTI PCHA settings.

Same protocol as Glioblastoma/codes/run_panelA.py (k=3–7):
  numIter=50, 150 observed inits (max volume), 50 inits per shuffle,
  500 shuffles, delta=0, fit in (k−1)=1 PC of the saved 12-PC scores.

py_pcha can hang on 1-D line search; this script uses a bounded copy of the
same PCHA updates so the permutation can finish. t-ratio in 1-D is segment
length / hull length (max−min). This is NOT Groves-selected k.

Writes only results/panel_a_k2/ and figures/Figure_1A_gbm_k2.png.
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
from src.archetypes import _shuffle_columns, simplex_volume
from k2_utils import MARKERS, assign_poles, signature_scores, style_mpl
from pcha_bounded import PCHA as PCHA_BOUNDED

MATRIX = GBM / "data" / "processed" / "input_panelA_glioblastoma_ccle_logtpm_filtered.csv"
SCORES = GBM / "results" / "panel_a" / "pc_scores.npy"
N_PCS_PATH = GBM / "results" / "panel_a" / "n_pcs.txt"
ESV_PATH = GBM / "results" / "panel_a" / "esv_curve.csv"
T_OTHER = GBM / "results" / "panel_a" / "t_ratio_parti_500.csv"
LABELS = GBM / "results" / "panel_b" / "verhaak_labels_panelA.csv"
OUT = GBM / "results" / "panel_a_k2"
FIG = GBM / "figures"

K = 2
NUM_ITER = 50
N_INIT_OBS = 3 * NUM_ITER
N_INIT_NULL = NUM_ITER
N_PERM = 500
N_CHECKPOINT = 25
DELTA = 0.0
SEED = 0


def fit_pcha_once(pc_scores, k, delta=DELTA):
    scores = np.asarray(pc_scores, dtype=float)[:, : int(k) - 1]
    X = scores.T
    XC, S, C, SSE, varexpl = PCHA_BOUNDED(X, noc=int(k), delta=delta)
    archetypes = np.asarray(XC).T
    weights = np.asarray(S).T
    return archetypes, weights, float(varexpl)


def fit_pcha_best(pc_scores, k, n_init, delta=DELTA, progress_every=0):
    best = None
    best_vol = -np.inf
    n_ok = 0
    t0 = time.time()
    for i in range(int(n_init)):
        try:
            archetypes, weights, varexpl = fit_pcha_once(pc_scores, k, delta=delta)
            vol = simplex_volume(archetypes)
            n_ok += 1
            if np.isfinite(vol) and vol > best_vol:
                best_vol = vol
                best = (archetypes, weights, varexpl, vol)
        except Exception:
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


def hull_t_ratio(scores, archetypes):
    """1-D hull length is max−min (scipy ConvexHull needs ≥2-D)."""
    k = archetypes.shape[0]
    data = np.asarray(scores, dtype=float)[:, : k - 1]
    if k == 2:
        hull = float(np.max(data) - np.min(data))
    else:
        from scipy.spatial import ConvexHull

        hull = float(ConvexHull(data).volume)
    if hull <= 0:
        return np.nan
    return float(simplex_volume(archetypes) / hull)


def p_label(p):
    if not np.isfinite(p):
        return "p = NA"
    if p == 0:
        return f"p < {1.0 / N_PERM:.3f}"
    if p < 0.01:
        return f"p = {p:.3f}"
    return f"p = {p:.2f}"


def permute_t_ratio(scores, archetypes, n_perm):
    observed = hull_t_ratio(scores, archetypes)
    null_path = OUT / f"null_t_ratios_k2_parti_n{n_perm}.npy"
    if null_path.exists():
        null = [float(x) for x in np.load(null_path).ravel()]
        print(f"  resume {null_path.name}: {len(null)}/{n_perm}", flush=True)
    else:
        null = []
    n_fail = sum(1 for x in null if not np.isfinite(x))
    start = len(null)
    t0 = time.time()
    for i in range(start, n_perm):
        shuffled = _shuffle_columns(
            scores[:, : K - 1], np.random.default_rng([SEED, int(K), int(i)])
        )
        try:
            arch, _, _, _, _ = fit_pcha_best(shuffled, K, n_init=N_INIT_NULL, delta=DELTA)
            value = hull_t_ratio(shuffled, arch)
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
                f"    k=2  {i + 1}/{n_perm}  n_ok={len(finite)}  n_fail={n_fail}  "
                f"running p={running_p:.4f}  {time.time() - t0:.0f}s",
                flush=True,
            )
    null_arr = np.asarray(null[:n_perm], dtype=float)
    np.save(null_path, null_arr)
    finite = null_arr[np.isfinite(null_arr)]
    p_value = float(np.mean(finite >= observed)) if len(finite) else np.nan
    return observed, p_value, int(len(finite)), int(n_fail)


def draw_figure(expr, pca12, n_pcs, archetypes, n_samples, t_row, esv, tab_all):
    X = expr.T.values
    tot_esv = esv["esv"].values * float(cumulative_variance(pca12)[-1])
    gene_esv = 100.0 * tot_esv
    k_esv = esv["k"].values
    delta_esv = np.diff(np.concatenate([[0.0], gene_esv]))

    pca2, scores2 = fit_pca(X, n_components=2)
    arcs_full = np.zeros((K, n_pcs))
    arcs_full[:, : archetypes.shape[1]] = archetypes
    gene_arcs = inverse_transform_scores(pca12, arcs_full)
    arcs2 = pca2.transform(gene_arcs)

    style_mpl(plt)
    fig = plt.figure(figsize=(6.6, 7.4))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 1.55], hspace=0.32, wspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(k_esv, delta_esv, "-o", color="#3B6FA0", ms=5, lw=1.4)
    ax.axvline(K, color="0.35", ls="--", lw=1)
    ax.annotate(
        "Hypothesis-driven k=2\n(not Groves-selected k)",
        xy=(K, float(delta_esv[k_esv == K][0])),
        xytext=(5.5, max(float(delta_esv.max()) * 0.72, 4)),
        fontsize=7.5,
        color="0.25",
        arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8),
    )
    ax.set_xlim(1.5, 15.5)
    ax.set_ylim(0, max(12, float(delta_esv.max()) * 1.08))
    ax.set_xticks(range(2, 16))
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("% ESV on top of N−1 model")
    ax.set_title("Explained sample variance (ESV)")

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(tab_all["k"], tab_all["t_ratio"], "-o", color="#3B6FA0", ms=6, lw=1.4)
    ax.axvline(K, color="0.35", ls="--", lw=1)
    ymax = max(1.05, float(tab_all["t_ratio"].max()) * 1.15)
    ax.set_ylim(0, ymax)
    ax.set_xlim(1.5, 7.5)
    ax.set_xticks(tab_all["k"].astype(int))
    for _, row in tab_all.iterrows():
        ax.annotate(
            p_label(row["p_value"]),
            (row["k"], row["t_ratio"]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=7.5,
            color="#222",
        )
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("t-ratio")
    ax.set_title("t-ratio (500 shuffles, PCHA numIter=50)")

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
    t = float(t_row["t_ratio"])
    p = float(t_row["p_value"])
    ax.set_title(
        f"Archetype space ({n_samples} glioblastoma cell lines)  "
        f"k=2  t={t:.3f}  {p_label(p)}"
    )

    fig.suptitle(
        "Figure 1A (k=2)  —  GBM CCLE, hypothesis-driven PN–MES axis "
        "(ParTI PCHA, 500 shuffles, numIter=50)",
        y=0.995,
        fontsize=10,
    )
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "Figure_1A_gbm_k2.png")
    plt.close(fig)
    print("Wrote", FIG / "Figure_1A_gbm_k2.png")


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
    print(f"Reused {SCORES}; max |score diff| vs rebuilt PCA = {np.max(np.abs(rebuilt - scores)):.4g}")
    print(
        f"k=2 ParTI PCHA on first 1 PC: numIter={NUM_ITER}, observed inits={N_INIT_OBS}, "
        f"null inits={N_INIT_NULL}, shuffles={int(args.n_perm)}, delta={DELTA}"
    )
    print("Not Groves smallest-significant-k. t-ratio is computed first.")

    arch_path = OUT / "archetypes_k2_parti.npy"
    s_path = OUT / "S_k2_parti.npy"
    row_path = OUT / "t_ratio_k2.csv"
    n_perm = int(args.n_perm)

    if arch_path.exists() and s_path.exists() and not args.force_fit:
        archetypes = np.load(arch_path)
        weights = np.load(s_path)
        print(f"Reused {arch_path.name} {archetypes.shape}")
        if archetypes.shape != (2, 1):
            print("Saved archetypes are not a 1-D PCHA fit; refitting.")
            archetypes = None
    else:
        archetypes = None

    if archetypes is None:
        print(f"Observed fit: {N_INIT_OBS} inits, keep max volume ...", flush=True)
        archetypes, weights, varexpl, vol, n_ok = fit_pcha_best(
            scores, K, n_init=N_INIT_OBS, delta=DELTA, progress_every=25
        )
        np.save(arch_path, archetypes)
        np.save(s_path, weights)
        print(
            f"  inits_ok={n_ok}/{N_INIT_OBS}  ESV={varexpl:.3f}  vol={vol:.4g}  "
            f"arcs={archetypes.ravel()}",
            flush=True,
        )

    observed = hull_t_ratio(scores, archetypes)
    print(f"  observed t={observed:.4f}", flush=True)

    null_path = OUT / f"null_t_ratios_k2_parti_n{n_perm}.npy"
    need_perm = True
    if row_path.exists() and null_path.exists() and np.load(null_path).size >= n_perm:
        prev = pd.read_csv(row_path)
        if "n_init_obs" in prev.columns and str(prev["n_init_obs"].iloc[0]) == str(N_INIT_OBS):
            print(f"Already complete: {row_path.name}", flush=True)
            need_perm = False
            t_row = prev.iloc[0].to_dict()

    if need_perm:
        print(
            f"Permutation: {n_perm} shuffles × {N_INIT_NULL} inits (same as k=3–7)",
            flush=True,
        )
        observed, p_value, n_ok, n_fail = permute_t_ratio(scores, archetypes, n_perm)
        t_row = {
            "k": K,
            "t_ratio": observed,
            "p_value": p_value,
            "n_perm": n_perm,
            "n_init_obs": N_INIT_OBS,
            "n_init_null": N_INIT_NULL,
            "numIter": NUM_ITER,
            "n_success": n_ok,
            "n_fail": n_fail,
            "n_pcs_fit": 1,
            "note": "hypothesis-driven k=2; ParTI PCHA; 1-D hull = max-min",
        }
        pd.DataFrame([t_row]).to_csv(row_path, index=False)
        print(
            f"  k=2 done: t={observed:.4f}, p={p_value:.4f} from {n_ok}/{n_perm} "
            f"(n_fail={n_fail})",
            flush=True,
        )
        print("  Do not treat this as Groves k selection.")

    sig, present = signature_scores(expr)
    sig.to_csv(OUT / "signature_scores_cell_lines.csv")
    print("Marker genes present on Panel A matrix:")
    for subtype, genes in MARKERS.items():
        print(f"  {subtype}: {len(present[subtype])}/{len(genes)} {present[subtype]}")
    poles = assign_poles(weights, sig)
    print(
        f"Pole assignment from corr(S, MES−PN): "
        f"Arc {poles['mes_idx'] + 1}=MES, Arc {poles['pn_idx'] + 1}=PN"
    )
    (OUT / "pole_assignment.txt").write_text(
        f"hypothesis_driven_k=2\n"
        f"mes_idx={poles['mes_idx']}\n"
        f"pn_idx={poles['pn_idx']}\n"
        f"corr_S_with_MES_minus_PN={poles['corr_S_with_MES_minus_PN']}\n"
        f"t_ratio={t_row['t_ratio']}\n"
        f"p_value={t_row['p_value']}\n"
        f"note=not Groves-selected k\n"
    )
    mix = pd.DataFrame(
        {
            "cell_line": sample_ids,
            "S_arc1": weights[:, 0],
            "S_arc2": weights[:, 1],
            "w_MES": weights[:, poles["mes_idx"]],
            "w_PN": weights[:, poles["pn_idx"]],
        }
    )
    mix = mix.merge(sig.reset_index().rename(columns={"sample": "cell_line"}), on="cell_line")
    mix.to_csv(OUT / "mixture_weights_with_signatures.csv", index=False)

    esv = pd.read_csv(ESV_PATH)
    other = pd.read_csv(T_OTHER)
    tab_all = pd.concat(
        [pd.DataFrame([{"k": t_row["k"], "t_ratio": t_row["t_ratio"], "p_value": t_row["p_value"]}]),
         other[["k", "t_ratio", "p_value"]]],
        ignore_index=True,
    ).sort_values("k")
    labels = pd.read_csv(LABELS)
    labels["cell_line"] = labels["cell_line"].astype(str)
    labels = labels.set_index("cell_line")
    _ = labels
    draw_figure(expr, pca12, n_pcs, archetypes, len(sample_ids), t_row, esv, tab_all)
    print("Wrote outputs to", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
