#!/usr/bin/env python3
"""Full ParTI Panel A t-ratio permutation test.

Matches reference/ParTI/{ParTI.m, findMinSimplex.m, CalculateSimplexTratiosPCHA.m}:
  - PCHA on first (k-1) PCs, delta=0
  - Observed simplex: 3*numIter = 150 random inits, keep MAX volume
    (reuse results/panel_a/archetypes_k{k}_parti.npy if present)
  - Each shuffle: independently permute each PC across samples
  - Per shuffle: numIter=50 PCHA inits, keep MAX t-ratio
  - maxRuns=1000 shuffles
  - p = fraction of finite null t-ratios >= observed t-ratio
  - Volume: ParTI det simplex / ConvexHull of data in (k-1)-D

Nulls and the summary table are written as soon as each k finishes, and
null arrays are checkpointed every N_CHECKPOINT shuffles so a crash does
not lose a completed k or a long partial run.
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import os
import shutil
import sys

SCLC = Path(__file__).resolve().parents[1]
ROOT = SCLC.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, QhullError

from src.archetypes import _shuffle_columns, fit_pcha_best, simplex_volume

OUT = SCLC / "results" / "panel_a"
OFFICIAL = SCLC / "results" / "t_ratio_official.csv"
PAPER = {
    3: {"p": 0.508, "t_ratio": 0.52},
    4: {"p": 0.059, "t_ratio": 0.247},
    5: {"p": 0.034, "t_ratio": 0.107},
    6: {"p": 0.016, "t_ratio": 0.043},
}
K_PERM = (3, 4, 5, 6)
N_INIT_OBS = 150
N_INIT_NULL = 50
N_PERM = 1000
N_CHECKPOINT = 50
DELTA = 0.0
SEED = 0


def hull_t_ratio(scores, archetypes):
    k = archetypes.shape[0]
    data = np.asarray(scores, dtype=float)[:, : k - 1]
    return float(simplex_volume(archetypes) / ConvexHull(data).volume)


def observed_fit(scores12, k):
    path = OUT / f"archetypes_k{k}_parti.npy"
    weights_path = OUT / f"S_k{k}_parti.npy"
    if path.exists():
        archetypes = np.load(path)
        print(f"  reused saved archetypes {path.name} shape={archetypes.shape}", flush=True)
        return archetypes
    print(f"  missing {path.name}; fitting {N_INIT_OBS} observed inits ...", flush=True)
    archetypes, weights, varexpl, vol, n_ok = fit_pcha_best(
        scores12, k, n_init=N_INIT_OBS, delta=DELTA
    )
    np.save(path, archetypes)
    np.save(weights_path, weights)
    print(f"  saved observed fit inits_ok={n_ok}/{N_INIT_OBS} ESV={varexpl:.3f} vol={vol:.4g}", flush=True)
    return archetypes


def load_checkpoint(path):
    if not path.exists():
        return []
    arr = np.load(path)
    values = [float(x) for x in np.asarray(arr, dtype=float).ravel()]
    print(f"  resume {path.name}: {len(values)}/{N_PERM} shuffles already stored", flush=True)
    return values


def write_table(rows):
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "t_ratio_parti_1000.csv", index=False)
    OFFICIAL.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(OUT / "t_ratio_parti_1000.csv", OFFICIAL)
    return table


def permute_one_k(scores12, k):
    paper = PAPER[k]
    null_path = OUT / f"null_t_ratios_k{k}_parti_n{N_PERM}.npy"
    print(
        f"\n=== k={k}: delta={DELTA}, {k-1} PCs, "
        f"{N_PERM} shuffles x {N_INIT_NULL} inits ===",
        flush=True,
    )
    archetypes = observed_fit(scores12, k)
    observed = hull_t_ratio(scores12, archetypes)
    print(
        f"  observed t={observed:.4f} (paper {paper['t_ratio']:.3f})",
        flush=True,
    )

    null = load_checkpoint(null_path)
    n_fail = sum(1 for x in null if not np.isfinite(x))
    start = len(null)
    if start >= N_PERM:
        print(f"  already complete ({start} stored); skipping shuffles", flush=True)
    else:
        print(f"  permutation: starting at shuffle {start + 1}/{N_PERM} ...", flush=True)

    for i in range(start, N_PERM):
        shuffled = _shuffle_columns(
            scores12[:, : k - 1], np.random.default_rng([SEED, int(k), int(i)])
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
        "n_success": int(len(finite)),
        "n_fail": int(n_fail),
        "paper_p": paper["p"],
        "paper_t_ratio": paper["t_ratio"],
    }
    print(
        f"  k={k} done: t={observed:.4f} (paper {paper['t_ratio']:.3f}), "
        f"p={p_value:.4f} from {len(finite)}/{N_PERM} (paper p={paper['p']:.3f})",
        flush=True,
    )
    return row


def k_is_complete(k):
    null_path = OUT / f"null_t_ratios_k{k}_parti_n{N_PERM}.npy"
    row_path = OUT / f"t_ratio_parti_1000_k{k}.csv"
    if not (null_path.exists() and row_path.exists()):
        return False
    return np.load(null_path).size >= N_PERM


def merge_available_rows():
    rows = []
    for k in K_PERM:
        path = OUT / f"t_ratio_parti_1000_k{k}.csv"
        if path.exists():
            rows.extend(pd.read_csv(path).to_dict("records"))
    if not rows:
        return None
    rows.sort(key=lambda r: int(r["k"]))
    return write_table(rows)


def run_one_k(k):
    for key in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[key] = "1"
    scores12 = np.load(OUT / "pc_scores_12.npy")
    row = permute_one_k(scores12, k)
    pd.DataFrame([row]).to_csv(OUT / f"t_ratio_parti_1000_k{k}.csv", index=False)
    return row


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pending = [k for k in K_PERM if not k_is_complete(k)]
    already = [k for k in K_PERM if k_is_complete(k)]
    if already:
        print(f"Already complete: {already}", flush=True)
        merge_available_rows()
    if not pending:
        table = merge_available_rows()
        print("\nParTI-matched 1000-shuffle table:")
        print(table.to_string(index=False), flush=True)
        print("Copied to", OFFICIAL, flush=True)
        return

    print(f"Running k={pending} in parallel (up to {len(pending)} workers)", flush=True)
    with ProcessPoolExecutor(max_workers=len(pending)) as pool:
        futures = {pool.submit(run_one_k, k): k for k in pending}
        for fut in as_completed(futures):
            k = futures[fut]
            fut.result()
            table = merge_available_rows()
            print(f"\nPersisted after k={k}:", flush=True)
            print(table.to_string(index=False), flush=True)

    table = merge_available_rows()
    print("\nParTI-matched 1000-shuffle table:")
    print(table.to_string(index=False), flush=True)
    print("Copied to", OFFICIAL, flush=True)


if __name__ == "__main__":
    main()
