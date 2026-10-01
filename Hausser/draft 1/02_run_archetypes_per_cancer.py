#!/usr/bin/env python3
"""Fit archetypes directly on one cancer type's tumor transcriptomes.

This is the actual Hausser Fig. 1d step: unlike the Groves-style pipeline
elsewhere in this repo, there is no cell-line fit + tumor projection here.
PCA and PCHA are both run straight on the tumor x gene matrix.

Follows the same ParTI / Hausser convention already documented in
src/archetypes.py and used for SCLC/Breast/GBM: fit PCHA in the first (k-1)
PCs, delta=0, many random restarts, keep the maximum-volume simplex; assess
significance with the t-ratio permutation test (shuffle each PC
independently, refit, compare observed vs. null t-ratio).

ASSUMPTIONS TO VERIFY against Hausser Supplementary Methods (not in the main
text PDF bundled in papers/ -- check the SI on the journal site):
  - Gene filtering before PCA: this script defaults to the top
    `--n-genes` (default 5000) most-variable genes. The paper doesn't state
    an exact number in the main text; adjust --n-genes if the SI gives one.
  - Number of PCA components used for fitting: this script computes enough
    PCs to cover k up to 5 (i.e. 4), plus a couple extra for the variance
    curve. Fig. 1d visualizes PC1-3 regardless of which k is used for
    fitting.
  - Number of permutation shuffles: defaults to 100 for a fast pass; bump to
    1000 with --n-perm 1000 to match Hausser's SCLC-equivalent convention
    (see Groves 2022 / this repo's SCLC t-ratio test) once you're ready for
    the final run.

Usage (from repository root):

    .venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/02_run_archetypes_per_cancer.py" --cancer-type THCA
    .venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/02_run_archetypes_per_cancer.py" --cancer-type BRCA_METABRIC --n-perm 1000
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parent
ROOT = ANALYSIS.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

import numpy as np
import pandas as pd

from src.pca import fit_pca
from src.archetypes import fit_pcha_best, t_ratio, permutation_t_ratio
from _registry import CANCER_TYPES

RESULTS = ANALYSIS / "results"

K_VALUES = (3, 4, 5)  # Hausser found 3-5 archetype polyhedra for these cancer types.
# k=2 is deliberately excluded: py_pcha fits it as a 1-D problem (a single PC),
# which is documented elsewhere in this repo to hang
# (see Glioblastoma/codes/run_panelA_k2.py -- "py_pcha is not used for k=2,
# it can hang in 1-D"). t-ratio is also undefined at k=2 (a line has no
# interior), so nothing is lost by skipping it here.
N_PCS_FOR_FIT = max(K_VALUES) - 1 + 2  # enough PCs for k=5 (needs 4), plus 2 extra headroom


def top_variance_genes(expr_genes_by_samples, n_genes):
    """expr: genes x samples. Returns the n_genes highest-variance rows."""
    var = expr_genes_by_samples.var(axis=1, skipna=True)
    keep = var.sort_values(ascending=False).index[: int(n_genes)]
    return expr_genes_by_samples.loc[keep]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancer-type", required=True, choices=list(CANCER_TYPES))
    ap.add_argument("--n-genes", type=int, default=5000, help="top-variance genes kept before PCA")
    ap.add_argument("--n-init", type=int, default=50, help="PCHA random restarts per k (ParTI convention)")
    ap.add_argument("--n-perm", type=int, default=100, help="permutation shuffles for the t-ratio test")
    ap.add_argument("--delta", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    code = args.cancer_type
    info = CANCER_TYPES[code]
    in_dir = RESULTS / code
    out_dir = RESULTS / code / "panel_a"
    out_dir.mkdir(parents=True, exist_ok=True)

    expr_path = in_dir / "tumor_expr_primary.csv"
    if not expr_path.is_file():
        print(f"Missing {expr_path}")
        print("Run 01_prepare_tcga_tumor.py (or 01b_prepare_metabric_full.py for breast) first.")
        return 1

    print(f"=== {code} ({info['label']}) ===")
    expr = pd.read_csv(expr_path, index_col=0)
    print(f"Loaded: {expr.shape[0]} genes x {expr.shape[1]} tumors")

    if expr.shape[1] < 20:
        print("Too few samples to fit anything meaningful -- check the prep step.")
        return 1

    expr_filt = top_variance_genes(expr, args.n_genes)
    print(f"Kept top {expr_filt.shape[0]} variance genes")

    # samples x genes for PCA
    X = expr_filt.T.dropna(axis=1, how="any")
    dropped = expr_filt.shape[0] - X.shape[1]
    if dropped:
        print(f"Dropped {dropped} genes with missing values before PCA")

    n_pcs = min(N_PCS_FOR_FIT, X.shape[0] - 1, X.shape[1])
    pca, scores = fit_pca(X.values, n_components=n_pcs, random_state=args.seed)
    print(f"PCA: {n_pcs} components, cumulative variance explained: "
          f"{np.cumsum(pca.explained_variance_ratio_).round(3).tolist()}")

    pd.DataFrame(
        scores, index=X.index, columns=[f"PC{i+1}" for i in range(n_pcs)]
    ).to_csv(out_dir / "pc_scores.csv")

    summary_rows = []
    for k in K_VALUES:
        n_dims_needed = k - 1
        if n_dims_needed > n_pcs:
            print(f"k={k}: need {n_dims_needed} PCs, only have {n_pcs} -- skipping")
            continue
        if k >= X.shape[0]:
            print(f"k={k}: not enough samples ({X.shape[0]}) -- skipping")
            continue

        print(f"\n--- k={k} ---")
        archetypes, weights, varexpl, vol, n_ok = fit_pcha_best(
            scores, k, n_init=args.n_init, delta=args.delta
        )
        print(f"  PCHA: {n_ok}/{args.n_init} inits converged, varexpl={varexpl:.3f}")

        obs_t_ratio = t_ratio(scores[:, : k - 1], archetypes)

        perm = permutation_t_ratio(
            scores[:, : k - 1],
            k,
            n_perm=args.n_perm,
            delta=args.delta,
            seed=args.seed,
            observed_archetypes=archetypes,
            verbose=True,
        )
        p_value = perm["p_value"]

        pd.DataFrame(
            archetypes, columns=[f"PC{i+1}" for i in range(archetypes.shape[1])]
        ).to_csv(out_dir / f"archetypes_k{k}.csv", index_label="archetype")
        pd.DataFrame(
            weights, index=X.index, columns=[f"archetype_{i+1}" for i in range(k)]
        ).to_csv(out_dir / f"weights_k{k}.csv")

        summary_rows.append(
            {
                "cancer_type": code,
                "label": info["label"],
                "k": k,
                "n_samples": X.shape[0],
                "n_genes": X.shape[1],
                "varexpl": varexpl,
                "t_ratio": obs_t_ratio,
                "p_value": p_value,
                "n_perm": args.n_perm,
                "hausser_reported_p": info.get("hausser_p", ""),
            }
        )
        print(f"  t-ratio={obs_t_ratio if not np.isnan(obs_t_ratio) else 'n/a'}, "
              f"p={p_value if not np.isnan(p_value) else 'n/a'} "
              f"(Hausser reports {info.get('hausser_p', 'n/a')})")

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "t_ratio_summary.csv", index=False)

    sig = summary[summary["p_value"] < 0.05]
    best_k = int(sig["k"].min()) if len(sig) else None
    (out_dir / "best_k.json").write_text(json.dumps({"cancer_type": code, "best_k": best_k}))
    print(f"\nSmallest significant k (p<0.05): {best_k}")
    print("Wrote", out_dir / "t_ratio_summary.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
