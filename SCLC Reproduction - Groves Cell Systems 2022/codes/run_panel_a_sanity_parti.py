#!/usr/bin/env python3
"""Sanity-check Panel A using ParTI's PCHA protocol.

Differences vs our first run:
  - fit in (k-1) PCs, not 12 PCs
  - delta=0, not 0.1
  - many random inits; keep the maximum-volume simplex (ParTI numIter)

Observed fit: 150 inits (3 * ParTI numIter=50).
Permutation: 100 shuffles x 15 inits (sanity, not the full 1000 x 50).
"""

from pathlib import Path
import sys

SCLC = Path(__file__).resolve().parents[1]
ROOT = SCLC.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from src.archetypes import (
    _shuffle_columns,
    fit_pcha_best,
    simplex_volume,
    t_ratio,
)
from scipy.spatial import ConvexHull, QhullError

OUT = SCLC / "results" / "panel_a"
PAPER = {
    4: {"p": 0.059, "t_ratio": 0.247},
    5: {"p": 0.034, "t_ratio": 0.107},
    6: {"p": 0.016, "t_ratio": 0.043},
}
K_PERM = (4, 5, 6)
N_INIT_OBS = 150
N_INIT_NULL = 15
N_PERM = 100
DELTA = 0.0
SEED = 0


def hull_t_ratio(scores, archetypes):
    k = archetypes.shape[0]
    data = scores[:, : k - 1]
    return float(simplex_volume(archetypes) / ConvexHull(data).volume)


def main():
    scores12 = np.load(OUT / "pc_scores_12.npy")
    rng = np.random.default_rng(SEED)
    rows = []
    for k in K_PERM:
        print(f"\n=== k={k}: {N_INIT_OBS} observed inits (delta={DELTA}, {k-1} PCs) ===", flush=True)
        archetypes, weights, varexpl, vol, n_ok = fit_pcha_best(
            scores12, k, n_init=N_INIT_OBS, delta=DELTA
        )
        observed = hull_t_ratio(scores12, archetypes)
        np.save(OUT / f"archetypes_k{k}_parti.npy", archetypes)
        np.save(OUT / f"S_k{k}_parti.npy", weights)
        paper = PAPER[k]
        print(
            f"  inits_ok={n_ok}/{N_INIT_OBS}  ESV={varexpl:.3f}  "
            f"t={observed:.4f} (paper {paper['t_ratio']:.3f})",
            flush=True,
        )

        print(f"  permutation: {N_PERM} shuffles x {N_INIT_NULL} inits ...", flush=True)
        null = []
        n_fail = 0
        for i in range(N_PERM):
            shuffled = _shuffle_columns(scores12[:, : k - 1], rng)
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
                "paper_p": paper["p"],
                "paper_t_ratio": paper["t_ratio"],
            }
        )
        print(
            f"  k={k}: t={observed:.4f} (paper {paper['t_ratio']:.3f}), "
            f"p={p_value:.3f} from {len(null)}/{N_PERM} (paper p={paper['p']:.3f})",
            flush=True,
        )

    table = pd.DataFrame(rows)
    table.to_csv(OUT / "t_ratio_sanity_parti.csv", index=False)
    print("\nParTI-matched sanity table:")
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
