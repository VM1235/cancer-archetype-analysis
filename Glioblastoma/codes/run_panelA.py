#!/usr/bin/env python3
"""GBM Panel A: Groves/ParTI PCHA on DepMap glioblastoma cell lines.

Calling convention matches Groves / ParTI.m algNum=5 (same as SCLC):
  numIter=50
  observed: 3*numIter = 150 inits, keep max volume
  each shuffle: numIter = 50 inits, keep max t-ratio
  maxRuns: 500 (user request; Groves SCLC used 1000)
  delta=0, fit in (k-1) PCs

dim: smallest n_pcs with >=50% gene-space variance on THIS GBM matrix.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import os
import sys

HERE = Path(__file__).resolve().parent
GBM = HERE.parent
ROOT = GBM.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.spatial import ConvexHull, QhullError

from src.io import load_expression_csv
from src.pca import fit_pca, cumulative_variance, inverse_transform_scores
from src.archetypes import (
    _shuffle_columns,
    esv_curve,
    fit_pcha_best,
    simplex_volume,
)

MATRIX = GBM / "data" / "processed" / "input_panelA_glioblastoma_ccle_logtpm_filtered.csv"
OUT = GBM / "results" / "panel_a"
FIG = GBM / "figures"

TARGET_CUMVAR = 0.50
K_PERM = (3, 4, 5, 6, 7)
NUM_ITER = 50
N_INIT_OBS = 3 * NUM_ITER
N_INIT_NULL = NUM_ITER
N_PERM = 500
N_CHECKPOINT = 25
DELTA = 0.0
SEED = 0


def choose_n_pcs(cumvar, target=TARGET_CUMVAR):
    hit = np.where(cumvar >= target)[0]
    if len(hit):
        return int(hit[0] + 1)
    return int(len(cumvar))


def dimension_finder(esv):
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


def p_label(p):
    if not np.isfinite(p):
        return "p = NA"
    if p == 0:
        return f"p < {1.0 / N_PERM:.3f}"
    if p < 0.01:
        return f"p = {p:.3f}"
    return f"p = {p:.2f}"


def permute_one_k(k):
    for key in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[key] = "1"

    scores = np.load(OUT / "pc_scores.npy")
    n_pcs = int((OUT / "n_pcs.txt").read_text().strip().splitlines()[0])
    if k - 1 > n_pcs:
        print(f"skip k={k}: needs {k-1} PCs, only have {n_pcs}", flush=True)
        return None

    null_path = OUT / f"null_t_ratios_k{k}_parti_n{N_PERM}.npy"
    print(
        f"\n=== k={k}: {N_INIT_OBS} observed inits, "
        f"{N_PERM} shuffles x {N_INIT_NULL} inits, delta={DELTA} ===",
        flush=True,
    )
    arch_path = OUT / f"archetypes_k{k}_parti.npy"
    if arch_path.exists():
        archetypes = np.load(arch_path)
        print(f"  reused {arch_path.name}", flush=True)
    else:
        archetypes, weights, varexpl, vol, n_ok = fit_pcha_best(
            scores, k, n_init=N_INIT_OBS, delta=DELTA
        )
        np.save(arch_path, archetypes)
        np.save(OUT / f"S_k{k}_parti.npy", weights)
        print(
            f"  inits_ok={n_ok}/{N_INIT_OBS}  ESV={varexpl:.3f}  vol={vol:.4g}",
            flush=True,
        )
    observed = hull_t_ratio(scores, archetypes)
    print(f"  observed t={observed:.4f}", flush=True)

    if null_path.exists():
        null = [float(x) for x in np.load(null_path).ravel()]
        print(f"  resume {null_path.name}: {len(null)}/{N_PERM}", flush=True)
    else:
        null = []
    n_fail = sum(1 for x in null if not np.isfinite(x))
    start = len(null)
    for i in range(start, N_PERM):
        shuffled = _shuffle_columns(
            scores[:, : k - 1], np.random.default_rng([SEED, int(k), int(i)])
        )
        try:
            arch, _, _, _, _ = fit_pcha_best(
                shuffled, k, n_init=N_INIT_NULL, delta=DELTA
            )
            value = hull_t_ratio(shuffled, arch)
            if not np.isfinite(value):
                n_fail += 1
                value = np.nan
        except (QhullError, ValueError, RuntimeError):
            n_fail += 1
            value = np.nan
        null.append(value)
        if (i + 1) % N_CHECKPOINT == 0 or (i + 1) == N_PERM:
            np.save(null_path, np.asarray(null, dtype=float))
            finite = np.asarray([x for x in null if np.isfinite(x)], dtype=float)
            running_p = float(np.mean(finite >= observed)) if len(finite) else np.nan
            print(
                f"    {i + 1}/{N_PERM}  n_ok={len(finite)}  n_fail={n_fail}  "
                f"running p={running_p:.4f}",
                flush=True,
            )

    null_arr = np.asarray(null[:N_PERM], dtype=float)
    np.save(null_path, null_arr)
    finite = null_arr[np.isfinite(null_arr)]
    p_value = float(np.mean(finite >= observed)) if len(finite) else np.nan
    row = {
        "k": k,
        "t_ratio": observed,
        "p_value": p_value,
        "n_perm": N_PERM,
        "n_init_obs": N_INIT_OBS,
        "n_init_null": N_INIT_NULL,
        "numIter": NUM_ITER,
        "n_success": int(len(finite)),
        "n_fail": int(n_fail),
        "n_pcs_fit": n_pcs,
    }
    pd.DataFrame([row]).to_csv(OUT / f"t_ratio_parti_500_k{k}.csv", index=False)
    print(
        f"  k={k} done: t={observed:.4f}, p={p_value:.4f} from {len(finite)}/{N_PERM}",
        flush=True,
    )
    return row


def merge_rows():
    rows = []
    for k in K_PERM:
        path = OUT / f"t_ratio_parti_500_k{k}.csv"
        if path.exists():
            rows.extend(pd.read_csv(path).to_dict("records"))
    if not rows:
        return None
    tab = pd.DataFrame(rows).sort_values("k")
    tab.to_csv(OUT / "t_ratio_parti_500.csv", index=False)
    return tab


def draw_figure(expr, pca, scores, n_pcs, esv, tab, k_star, n_samples):
    X = expr.T.values
    tot_esv = esv["esv"].values * float(cumulative_variance(pca)[-1])
    gene_esv = 100.0 * tot_esv
    k_esv = esv["k"].values
    delta_esv = np.diff(np.concatenate([[0.0], gene_esv]))

    pca2, scores2 = fit_pca(X, n_components=2)
    arcs = np.load(OUT / f"archetypes_k{k_star}_parti.npy")
    arcs_full = np.zeros((arcs.shape[0], n_pcs))
    arcs_full[:, : arcs.shape[1]] = arcs
    gene_arcs = inverse_transform_scores(pca, arcs_full)
    arcs2 = pca2.transform(gene_arcs)

    style()
    fig = plt.figure(figsize=(6.6, 7.4))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 1.55], hspace=0.32, wspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(k_esv, delta_esv, "-o", color="#3B6FA0", ms=5, lw=1.4)
    ax.axvline(k_star, color="0.35", ls="--", lw=1)
    ax.annotate(
        "Suggested number of\narchetypes by elbow",
        xy=(k_star, float(delta_esv[k_esv == k_star][0])),
        xytext=(7.2, max(float(delta_esv.max()) * 0.72, 4)),
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
    ax.plot(tab["k"], tab["t_ratio"], "-o", color="#3B6FA0", ms=6, lw=1.4)
    ax.axvline(k_star, color="0.35", ls="--", lw=1)
    ymax = max(0.62, float(tab["t_ratio"].max()) * 1.15)
    ax.set_ylim(0, ymax)
    ax.set_xlim(2.5, 7.5)
    ax.set_xticks(tab["k"].astype(int))
    for _, row in tab.iterrows():
        ax.annotate(
            p_label(row["p_value"]),
            (row["k"], row["t_ratio"]),
            textcoords="offset points",
            xytext=(8, 6),
            fontsize=8,
            color="#222",
        )
    ax.set_xlabel("Number of archetypes (N)")
    ax.set_ylabel("t-ratio")
    ax.set_title("t-ratio (500 shuffles, PCHA numIter=50)")

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
    ax.set_title(f"Archetype space ({n_samples} glioblastoma cell lines)")

    fig.suptitle(
        "Figure 1A  —  GBM CCLE (ParTI PCHA, 500 shuffles, numIter=50)",
        y=0.995,
        fontsize=11,
    )
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "Figure_1A_gbm.png")
    plt.close(fig)
    print("Wrote", FIG / "Figure_1A_gbm.png")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    expr = load_expression_csv(MATRIX)
    X = expr.T.values
    n_samples, n_genes = X.shape
    print(f"GBM matrix: {n_genes} genes × {n_samples} cell lines")
    print(
        f"ParTI PCHA: numIter={NUM_ITER}, observed inits={N_INIT_OBS}, "
        f"null inits={N_INIT_NULL}, shuffles={N_PERM}, delta={DELTA}"
    )

    n_fit = min(n_samples - 1, n_genes, 40)
    pca_all, _ = fit_pca(X, n_components=n_fit)
    cum_all = cumulative_variance(pca_all)
    n_pcs = choose_n_pcs(cum_all)
    pca, scores = fit_pca(X, n_components=n_pcs)
    cumvar = cumulative_variance(pca)
    print(
        f"dim={n_pcs} PCs explain {100 * cumvar[-1]:.1f}% variance "
        f"(rule = first n with ≥50%)"
    )
    pd.DataFrame(
        {
            "pc": np.arange(1, n_fit + 1),
            "explained_variance_ratio": pca_all.explained_variance_ratio_,
            "cumulative": cum_all,
        }
    ).to_csv(OUT / "pca_variance_full.csv", index=False)
    np.save(OUT / "pc_scores.npy", scores)
    pd.DataFrame(
        {
            "pc": np.arange(1, n_pcs + 1),
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "cumulative": cumvar,
        }
    ).to_csv(OUT / "pca_variance.csv", index=False)
    (OUT / "n_pcs.txt").write_text(str(n_pcs) + "\n")

    print(f"ESV curve k=2..{n_pcs + 1} in {n_pcs}-D (delta={DELTA}) ...")
    k_esv_vals = list(range(2, n_pcs + 2))
    fits = esv_curve(scores, k_esv_vals, delta=DELTA, seed=SEED)
    esv = pd.DataFrame([{"k": f["k"], "esv": f["esv"]} for f in fits])
    esv["delta_esv"] = esv["esv"].diff()
    esv.to_csv(OUT / "esv_curve.csv", index=False)
    print(esv.to_string(index=False))
    tot_esv = esv["esv"].values * float(cumvar[-1])
    k_from_finder = dimension_finder(tot_esv)
    print(f"ParTI DimensionFinder suggests k={k_from_finder}")

    pending = []
    for k in K_PERM:
        row_path = OUT / f"t_ratio_parti_500_k{k}.csv"
        null_path = OUT / f"null_t_ratios_k{k}_parti_n{N_PERM}.npy"
        if row_path.exists() and null_path.exists() and np.load(null_path).size >= N_PERM:
            print(f"Already complete: k={k}", flush=True)
        else:
            pending.append(k)

    if pending:
        n_workers = min(len(pending), os.cpu_count() or 1)
        print(f"Permutation workers={n_workers} for k={pending}", flush=True)
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futs = {pool.submit(permute_one_k, k): k for k in pending}
            for fut in as_completed(futs):
                fut.result()
                merge_rows()

    tab = merge_rows()
    if tab is None or tab.empty:
        print("No t-ratio table; permutation did not finish")
        return 1
    print("\nGBM Panel A t-ratio (ParTI PCHA, 500 shuffles, numIter=50):")
    print(tab.to_string(index=False))

    sig = tab.loc[tab["p_value"] < 0.05].sort_values("k")
    if len(sig):
        k_star = int(sig["k"].iloc[0])
        why = "smallest k with p < 0.05"
    else:
        k_star = k_from_finder
        why = f"DimensionFinder elbow (no k with p<0.05); finder={k_from_finder}"
    print(f"\nSuggested k={k_star} ({why})")
    (OUT / "suggested_k.txt").write_text(f"{k_star}\n{why}\n")

    draw_figure(expr, pca, scores, n_pcs, esv, tab, k_star, n_samples)
    print(f"Wrote outputs to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
