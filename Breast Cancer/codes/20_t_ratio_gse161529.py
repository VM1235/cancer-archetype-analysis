#!/usr/bin/env python3
"""Sanity-check t-ratio permutation on GSE161529 self-fit archetypes.

Reuses frozen Pal PC scores and the 15-init point-estimate vertices from
17.9_run_panelA_ks_genelist_gse161529.py. Does NOT refit the observed
simplex. Null: shuffle each of the first (k-1) PCs independently, refit
PCHA with N_INIT_NULL inits, keep max volume (same as DepMap Panel A).

This is a 50-shuffle check at k=3 and k=4 only — not the paper protocol
(500 shuffles). MAGIC-imputed cells make a small p easier than on bulk
DepMap lines; interpret accordingly.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/20_t_ratio_gse161529.py"
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
from scipy.spatial import ConvexHull, QhullError

from src.archetypes import _shuffle_columns, fit_pcha_best, simplex_volume

PAL_FIT = BREAST / "results" / "panel_a_ks_genelist_gse161529"
FIG = BREAST / "figures"

K_PERM = (3, 4)
NUM_ITER = 5
N_INIT_OBS = 3 * NUM_ITER  # recorded only; observed vertices already saved
N_INIT_NULL = NUM_ITER
N_PERM = 50
N_CHECKPOINT = 5
DELTA = 0.0
SEED = 0


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


def main():
    scores_path = PAL_FIT / "pc_scores.npy"
    if not scores_path.is_file():
        print(f"Missing {scores_path}")
        return 1
    scores = np.load(scores_path)
    n_pcs = int((PAL_FIT / "n_pcs.txt").read_text().strip().splitlines()[0])
    print(
        f"Pal scores {scores.shape[0]} cells x {scores.shape[1]} PCs "
        f"(n_pcs.txt={n_pcs})"
    )
    print(
        f"Permutation: {N_PERM} shuffles x {N_INIT_NULL} null inits, "
        f"delta={DELTA}, k={list(K_PERM)}"
    )
    print("Observed vertices reused from 17.9 (not refit).")

    rows = []
    for k in K_PERM:
        arcs_path = PAL_FIT / f"archetypes_k{k}_parti.npy"
        if not arcs_path.is_file():
            print(f"Missing {arcs_path}")
            return 1
        archetypes = np.load(arcs_path)
        observed = hull_t_ratio(scores, archetypes)
        print(
            f"\n=== k={k}: observed t={observed:.4f}  "
            f"vertices {archetypes.shape} ===",
            flush=True,
        )

        null_path = PAL_FIT / f"null_t_ratios_k{k}_parti_n{N_PERM}.npy"
        null: list[float] = []
        n_fail = 0
        if null_path.is_file():
            prev = np.load(null_path)
            null = [float(x) for x in prev]
            n_fail = int(np.sum(~np.isfinite(prev)))
            print(f"  resuming {len(null)}/{N_PERM} shuffles from {null_path.name}", flush=True)

        rng = np.random.default_rng(SEED + 1000 + k)
        # Replay RNG to the resume point so later shuffles stay reproducible.
        for _ in range(len(null)):
            _shuffle_columns(scores[:, : k - 1], rng)

        print(f"  permutation: {N_PERM} shuffles x {N_INIT_NULL} inits ...", flush=True)
        start = len(null)
        for i in range(start, N_PERM):
            shuffled = _shuffle_columns(scores[:, : k - 1], rng)
            try:
                arch, _, _, _, _ = fit_pcha_best(
                    shuffled, k, n_init=N_INIT_NULL, delta=DELTA
                )
                value = hull_t_ratio(shuffled, arch)
                if not np.isfinite(value):
                    n_fail += 1
                    value = np.nan
            except (QhullError, ValueError, RuntimeError) as err:
                n_fail += 1
                value = np.nan
                print(f"    shuffle {i + 1} failed: {err}", flush=True)
            null.append(value)
            np.save(null_path, np.asarray(null, dtype=float))
            finite = np.asarray([x for x in null if np.isfinite(x)], dtype=float)
            running_p = float(np.mean(finite >= observed)) if len(finite) else np.nan
            if (i + 1) % N_CHECKPOINT == 0 or (i + 1) == N_PERM or (i + 1) <= start + 1:
                print(
                    f"    {i + 1}/{N_PERM}  n_ok={len(finite)}  n_fail={n_fail}  "
                    f"running p={running_p:.4f}",
                    flush=True,
                )

        null_arr = np.asarray(null, dtype=float)
        finite = null_arr[np.isfinite(null_arr)]
        p_value = float(np.mean(finite >= observed)) if len(finite) else np.nan
        rows.append(
            {
                "k": k,
                "t_ratio": observed,
                "p_value": p_value,
                "n_perm": N_PERM,
                "n_init_obs": N_INIT_OBS,
                "n_init_null": N_INIT_NULL,
                "numIter": NUM_ITER,
                "n_success": int(len(finite)),
                "n_fail": n_fail,
                "n_pcs_fit": n_pcs,
                "n_cells": int(scores.shape[0]),
                "magic": True,
            }
        )
        print(
            f"  k={k}: t={observed:.4f}, p={p_value:.4f} from {len(finite)}/{N_PERM}",
            flush=True,
        )

    tab = pd.DataFrame(rows)
    out_csv = PAL_FIT / "t_ratio_parti_50.csv"
    tab.to_csv(out_csv, index=False)
    print("\nGSE161529 t-ratio (50 shuffles, MAGIC Pal cells):")
    print(tab.to_string(index=False))
    print("Wrote", out_csv)

    style()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, k in zip(axes, K_PERM):
        observed = float(tab.loc[tab["k"] == k, "t_ratio"].iloc[0])
        p_value = float(tab.loc[tab["k"] == k, "p_value"].iloc[0])
        null = np.load(PAL_FIT / f"null_t_ratios_k{k}_parti_n{N_PERM}.npy")
        finite = null[np.isfinite(null)]
        ax.hist(finite, bins=15, color="#9ecae1", edgecolor="white")
        ax.axvline(observed, color="#E45756", lw=1.6, label=f"observed t={observed:.3f}")
        ax.set_xlabel("t-ratio (null)")
        ax.set_ylabel("shuffles")
        ax.set_title(f"k={k}  p={p_value:.3f}  (n={len(finite)})")
        ax.legend(fontsize=7)
    fig.suptitle("GSE161529 Pal self-fit t-ratio, 50 shuffles (not 500)", y=1.03)
    fig.savefig(FIG / "Figure_1A_gse161529_t_ratio_50.png")
    plt.close(fig)
    print("Wrote", FIG / "Figure_1A_gse161529_t_ratio_50.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
