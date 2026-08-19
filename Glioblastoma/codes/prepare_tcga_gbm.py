#!/usr/bin/env python3
"""Prepare TCGA-GBM tumors for projection into the cell-line Panel A space.

Xena HiSeqV2 is already log2. DepMap Panel A is already log2(TPM+1).
No extra log transform.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
GBM = HERE.parent
ROOT = GBM.parent
sys.path.insert(0, str(ROOT))

import pandas as pd

from src.io import load_expression_csv

DATA = GBM / "data"
CL_MATRIX = GBM / "data" / "processed" / "input_panelA_glioblastoma_ccle_logtpm_filtered.csv"
OUT = GBM / "results" / "panel_c_tcga"

CLIN_KEEP = [
    "sampleID",
    "sample_type",
    "sample_type_id",
    "histological_type",
    "GeneExp_Subtype",
    "G_CIMP_STATUS",
    "gender",
    "vital_status",
]


def find_tcga_files():
    expr = None
    clin = None
    for p in DATA.rglob("*"):
        if not p.is_file():
            continue
        if p.name == "HiSeqV2" or p.name.endswith("HiSeqV2"):
            expr = p
        if "clinicalMatrix" in p.name:
            clin = p
    return expr, clin


def tcga_sample_code(sample_id):
    parts = str(sample_id).split("-")
    if len(parts) >= 4:
        return parts[3][:2]
    return ""


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    expr_path, clin_path = find_tcga_files()
    if expr_path is None or clin_path is None:
        print("Could not find HiSeqV2 / clinicalMatrix under", DATA)
        return 1
    print("Expression:", expr_path)
    print("Clinical:  ", clin_path)

    clin = pd.read_csv(clin_path, sep="\t")
    print("\n=== Clinical columns used ===")
    for c in CLIN_KEEP:
        present = c in clin.columns
        print(f"  {c}: {'OK' if present else 'MISSING'}")
        if present:
            print(clin[c].value_counts(dropna=False).head(8).to_string())
            print()

    tcga = pd.read_csv(expr_path, sep="\t", index_col=0)
    tcga.index = tcga.index.astype(str)
    tcga.columns = tcga.columns.astype(str)
    n_downloaded = tcga.shape[1]
    print(f"TCGA HiSeqV2: {tcga.shape[0]} genes × {n_downloaded} samples")
    print(f"  value range: {float(tcga.min().min()):.3f} .. {float(tcga.max().max()):.3f}")
    print("  (Xena HiSeqV2 is already log2-transformed; not logged again)")

    cl = load_expression_csv(CL_MATRIX)
    cl_genes = cl.index.astype(str)
    shared = cl_genes.intersection(tcga.index)
    print(f"\nPanel A cell-line genes: {len(cl_genes)}")
    print(f"Shared gene symbols: {len(shared)}")

    meta = clin.copy()
    if "sampleID" not in meta.columns:
        raise ValueError("clinical matrix has no sampleID")
    meta["sampleID"] = meta["sampleID"].astype(str)
    meta = meta.set_index("sampleID")
    matched = [s for s in tcga.columns if s in meta.index]
    print(f"Expression IDs in clinical table: {len(matched)}/{n_downloaded}")

    codes = pd.Series([tcga_sample_code(s) for s in matched], index=matched)
    print("TCGA sample-type barcode suffix among matched:")
    print(codes.value_counts().to_string())

    is_primary = codes.eq("01")
    if "sample_type_id" in meta.columns:
        is_primary = is_primary & meta.loc[matched, "sample_type_id"].eq(1)
    primary_ids = [s for s in matched if bool(is_primary.loc[s])]
    print(f"Primary tumors (barcode 01 & sample_type_id=1): {len(primary_ids)}")

    keep_cols = [c for c in CLIN_KEEP if c in clin.columns and c != "sampleID"]
    meta_out = meta.loc[primary_ids, keep_cols].copy()
    meta_out.index.name = "sampleID"
    n_sub = int(meta_out["GeneExp_Subtype"].notna().sum()) if "GeneExp_Subtype" in meta_out else 0
    print(f"Primary with GeneExp_Subtype: {n_sub}")

    expr_out = tcga.loc[shared, primary_ids]
    expr_out = expr_out.loc[~expr_out.index.duplicated(keep="first")]
    n_nan = int(expr_out.isna().sum().sum())
    print(f"Processed tumor matrix: {expr_out.shape[0]} genes × {expr_out.shape[1]} primary tumors")
    print(f"NaN entries: {n_nan}")

    expr_out.to_csv(OUT / "tcga_primary_log_shared_genes.csv")
    meta_out.to_csv(OUT / "tcga_primary_metadata.csv")
    report = [
        f"expression_file={expr_path}",
        f"clinical_file={clin_path}",
        f"n_tcga_downloaded={n_downloaded}",
        f"n_genes_tcga={tcga.shape[0]}",
        f"n_genes_cell_line={len(cl_genes)}",
        f"n_genes_shared={len(shared)}",
        f"n_expr_in_clinical={len(matched)}",
        f"n_primary_tumors={len(primary_ids)}",
        f"n_primary_with_GeneExp_Subtype={n_sub}",
        f"n_nan={n_nan}",
        "no extra log transform (HiSeqV2 already log2)",
        "",
    ]
    (OUT / "prepare_report.txt").write_text("\n".join(report))
    print("Wrote", OUT / "tcga_primary_log_shared_genes.csv")
    print("Wrote", OUT / "tcga_primary_metadata.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
