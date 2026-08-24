#!/usr/bin/env python3
"""Run 02_run_archetypes_per_cancer.py for every cancer type in Fig. 1d and
aggregate the results into one summary table (Hausser-style Table 1).

Assumes 01_prepare_tcga_tumor.py (per TCGA type) and
01b_prepare_metabric_full.py (breast) have already been run -- this script
does not download or prepare data itself, it only fits archetypes and
collects results for whichever cancer types already have
results/<CODE>/tumor_expr_primary.csv on disk. Missing ones are skipped
with a warning so you can run this incrementally as downloads complete.

Usage (from repository root):

    .venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/03_run_all_cancer_types.py"
    .venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/03_run_all_cancer_types.py" --n-perm 1000  # final run
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
import sys

from _paths import ANALYSIS, HERE

sys.path.insert(0, str(HERE))

import pandas as pd

from _registry import CANCER_TYPES, fig1d_types

RESULTS = ANALYSIS / "results"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--n-genes",
        type=int,
        default=None,
        help="Optional top-variance gene filter; default None = all genes, matching the paper's Methods.",
    )
    ap.add_argument("--n-init", type=int, default=50)
    ap.add_argument(
        "--n-perm",
        type=int,
        default=100,
        help="Fast-pass default. Use --n-perm 1000 to match the paper's Methods for final numbers.",
    )
    args = ap.parse_args()

    if args.n_perm < 1000:
        print(
            f"NOTE: --n-perm={args.n_perm} is a fast dev pass. The paper's Methods "
            "use 1000 shuffles for the t-ratio permutation test -- re-run with "
            "--n-perm 1000 before treating p-values as final.\n"
        )

    all_summaries = []
    skipped = []

    for code in fig1d_types():
        expr_path = RESULTS / code / "tumor_expr_primary.csv"
        if not expr_path.is_file():
            print(f"[skip] {code}: {expr_path} not found -- run the prep step first")
            skipped.append(code)
            continue

        print(f"\n{'=' * 60}\nRunning {code} ({CANCER_TYPES[code]['label']})\n{'=' * 60}")
        cmd = [
            sys.executable,
            str(HERE / "02_run_archetypes_per_cancer.py"),
            "--cancer-type", code,
            "--n-init", str(args.n_init),
            "--n-perm", str(args.n_perm),
        ]
        if args.n_genes is not None:
            cmd += ["--n-genes", str(args.n_genes)]
        ret = subprocess.run(cmd)
        if ret.returncode != 0:
            print(f"[FAILED] {code} exited with code {ret.returncode}")
            continue

        summary_path = RESULTS / code / "panel_a" / "t_ratio_summary.csv"
        if summary_path.is_file():
            all_summaries.append(pd.read_csv(summary_path))

    if skipped:
        print(f"\nSkipped (no prepared data yet): {skipped}")

    if not all_summaries:
        print("\nNo results to aggregate -- nothing ran successfully.")
        return 1

    combined = pd.concat(all_summaries, ignore_index=True)
    out_path = ANALYSIS / "results" / "fig1d_summary_table.csv"
    combined.to_csv(out_path, index=False)
    print(f"\nWrote combined summary: {out_path}")
    print(combined.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
