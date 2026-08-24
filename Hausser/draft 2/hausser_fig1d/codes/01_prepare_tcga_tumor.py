#!/usr/bin/env python3
"""Prepare one TCGA cancer type's primary-tumor expression matrix.

Unlike the Breast/GBM Panel C prep scripts (which project tumors onto a
*cell-line* archetype space), Hausser Fig. 1d fits archetypes directly on the
tumor transcriptomes -- no cell lines involved. This script just loads Xena
HiSeqV2 (already log2 RSEM), keeps primary tumors, and writes a clean
genes x samples CSV per cancer type.

Usage (from repository root, once data is downloaded -- see docs/README.md):

    .venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/01_prepare_tcga_tumor.py" --cancer-type THCA
    .venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/01_prepare_tcga_tumor.py" --cancer-type BLCA
    ... (repeat for LIHC, COAD, LGG, LUAD, HNSC)

Expects, per cancer type code, files placed under:

    data/<CODE>/HiSeqV2                          (Xena "gene RSEM log2" matrix)
    data/<CODE>/*clinicalMatrix*                  (Xena phenotype/clinical matrix)

matching the layout already used in Glioblastoma/data and Breast Cancer/data/tumors.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from _paths import ANALYSIS, HERE, ROOT

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
import pandas as pd
from _registry import CANCER_TYPES

DATA = ANALYSIS / "data"
RESULTS = ANALYSIS / "results"


def tcga_sample_code(sample_id):
    """Barcode suffix (e.g. '01' = primary solid tumor, '11' = normal)."""
    parts = str(sample_id).split("-")
    if len(parts) >= 4:
        return parts[3][:2]
    return ""


def find_files(cancer_dir):
    expr, clin = None, None
    for p in cancer_dir.rglob("*"):
        if not p.is_file():
            continue
        if p.name == "HiSeqV2" or p.name.endswith("HiSeqV2"):
            expr = p
        if "clinicalMatrix" in p.name or "clinical_matrix" in p.name.lower():
            clin = p
    return expr, clin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--cancer-type",
        required=True,
        choices=[c for c in CANCER_TYPES if CANCER_TYPES[c]["source"] == "TCGA"],
    )
    ap.add_argument(
        "--min-samples",
        type=int,
        default=250,
        help="Hausser's inclusion threshold (>=250 primary tumors); warns if not met.",
    )
    args = ap.parse_args()

    code = args.cancer_type
    info = CANCER_TYPES[code]
    cancer_dir = DATA / code
    out_dir = RESULTS / code
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== {code} ({info['label']}) ===")
    print("Expected Xena cohort:", info.get("xena_cohort"))

    expr_path, clin_path = find_files(cancer_dir)
    if expr_path is None or clin_path is None:
        print(f"Could not find HiSeqV2 / clinicalMatrix under {cancer_dir}")
        print(f"Place the Xena download there first -- see docs/README.md.")
        return 1

    print("Expression:", expr_path)
    print("Clinical:  ", clin_path)

    clin = pd.read_csv(clin_path, sep="\t")
    if "sampleID" not in clin.columns:
        raise ValueError(f"{clin_path} has no sampleID column")
    clin["sampleID"] = clin["sampleID"].astype(str)
    clin = clin.set_index("sampleID")

    expr = pd.read_csv(expr_path, sep="\t", index_col=0)
    expr.index = expr.index.astype(str)
    expr.columns = expr.columns.astype(str)
    n_downloaded = expr.shape[1]
    print(f"Xena HiSeqV2: {expr.shape[0]} genes x {n_downloaded} samples")
    print(f"  value range: {float(expr.min().min()):.3f} .. {float(expr.max().max()):.3f}")
    print("  (Xena HiSeqV2 is already log2; not logged again)")

    matched = [s for s in expr.columns if s in clin.index]
    print(f"Expression IDs in clinical table: {len(matched)}/{n_downloaded}")

    codes = pd.Series([tcga_sample_code(s) for s in matched], index=matched)
    print("Sample-type barcode suffix among matched:")
    print(codes.value_counts().to_string())

    # Methods: "we focused our analyses on primary tumors (field 'sample_type'
    # set to 'Primary Tumor' in the TCGA clinical annotation)". Prefer that
    # literal field when Xena's clinicalMatrix has it; fall back to the
    # barcode-01 heuristic (+ sample_type_id==1 if present) only if it doesn't.
    if "sample_type" in clin.columns:
        is_primary = clin.loc[matched, "sample_type"].astype(str).str.strip().eq("Primary Tumor")
        filter_desc = "sample_type == 'Primary Tumor' (matches paper's Methods verbatim)"
    else:
        is_primary = codes.eq("01")
        if "sample_type_id" in clin.columns:
            is_primary = is_primary & clin.loc[matched, "sample_type_id"].eq(1)
        filter_desc = (
            "barcode suffix '01'"
            + (" & sample_type_id=1" if "sample_type_id" in clin.columns else "")
            + " (fallback: no 'sample_type' column in this clinicalMatrix)"
        )
    primary_ids = [s for s in matched if bool(is_primary.loc[s])]
    print(f"Primary tumors ({filter_desc}): {len(primary_ids)}")

    if len(primary_ids) < args.min_samples:
        print(
            f"WARNING: {len(primary_ids)} primary tumors < Hausser's {args.min_samples} "
            "sample threshold. Check the download / consider excluding this type."
        )

    expr_out = expr.loc[:, primary_ids]
    expr_out = expr_out.loc[~expr_out.index.duplicated(keep="first")]
    n_nan = int(expr_out.isna().sum().sum())
    print(f"Output matrix: {expr_out.shape[0]} genes x {expr_out.shape[1]} primary tumors, {n_nan} NaN")

    expr_out.to_csv(out_dir / "tumor_expr_primary.csv")
    clin.loc[primary_ids].to_csv(out_dir / "tumor_clinical_primary.csv")

    report = [
        f"cancer_type={code}",
        f"expression_file={expr_path}",
        f"clinical_file={clin_path}",
        f"n_downloaded={n_downloaded}",
        f"n_genes={expr_out.shape[0]}",
        f"n_primary_tumors={expr_out.shape[1]}",
        f"n_nan={n_nan}",
        "no extra log transform (Xena HiSeqV2 already log2)",
        "",
    ]
    (out_dir / "prepare_report.txt").write_text("\n".join(report))
    print("Wrote", out_dir / "tumor_expr_primary.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
