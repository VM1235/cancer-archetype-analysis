#!/usr/bin/env python3
"""Prepare the full-gene METABRIC breast tumor matrix for Fig. 1d.

Different from Breast Cancer/codes/prepare_metabric_brca.py:
  - that script restricts METABRIC to the 206-gene KS signature and projects
    onto a *cell-line* archetype space (the Groves-style pipeline).
  - this script keeps the full gene set and does NOT touch cell lines at all,
    because Hausser fits archetypes directly on the 1970 METABRIC tumor
    transcriptomes (per the Fig. 1d caption: "TCGA, breast cancer from
    Metabric").

Reuses the same already-downloaded cBioPortal export
(Breast Cancer/brca_metabric/), so if you've already followed
Breast Cancer/brca_metabric/DOWNLOAD.md, no new download is needed here.

Usage (from repository root):

    .venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/01b_prepare_metabric_full.py"
"""

from __future__ import annotations

from pathlib import Path
import sys

from _paths import ANALYSIS, ROOT

sys.path.insert(0, str(ROOT))

import pandas as pd

BREAST = ROOT / "Breast Cancer"
EXPR_PATH = BREAST / "brca_metabric" / "data_mrna_illumina_microarray.txt"
CLIN_SAMPLE = BREAST / "brca_metabric" / "data_clinical_sample.txt"
CLIN_PATIENT = BREAST / "brca_metabric" / "data_clinical_patient.txt"

OUT = ANALYSIS / "results" / "BRCA_METABRIC"

CLIN_KEEP = [
    "PATIENT_ID",
    "SAMPLE_ID",
    "CANCER_TYPE",
    "CANCER_TYPE_DETAILED",
    "SAMPLE_TYPE",
    "ER_STATUS",
    "HER2_STATUS",
    "PR_STATUS",
    "CLAUDIN_SUBTYPE",
]


def read_cbioportal_table(path):
    """Skip leading '#' comment lines; return a DataFrame."""
    skip = 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                skip += 1
            else:
                break
    return pd.read_csv(path, sep="\t", skiprows=skip)


def main():
    if not EXPR_PATH.is_file():
        print(f"Missing {EXPR_PATH}")
        print("Follow Breast Cancer/brca_metabric/DOWNLOAD.md first (cBioPortal export).")
        return 1

    OUT.mkdir(parents=True, exist_ok=True)
    print("Expression:", EXPR_PATH)
    print("Clinical (sample):", CLIN_SAMPLE)
    print("Clinical (patient):", CLIN_PATIENT)

    patient = read_cbioportal_table(CLIN_PATIENT)
    n_pat, n_unique = len(patient), patient["PATIENT_ID"].nunique()
    print(f"\nPATIENT_ID check: {n_unique} unique / {n_pat} rows")
    if n_unique != n_pat:
        dups = patient["PATIENT_ID"].value_counts()
        raise ValueError(f"Duplicate PATIENT_IDs: {dups[dups > 1].head().to_dict()}")

    sample_clin = read_cbioportal_table(CLIN_SAMPLE)
    clin = sample_clin.merge(
        patient[["PATIENT_ID", "CLAUDIN_SUBTYPE"]],
        on="PATIENT_ID",
        how="left",
        validate="one_to_one",
    )
    print(f"Joined clinical: {clin.shape[0]} sample rows x {clin.shape[1]} cols")

    expr = pd.read_csv(EXPR_PATH, sep="\t")
    expr["Hugo_Symbol"] = expr["Hugo_Symbol"].astype(str)
    expr = expr.loc[~expr["Hugo_Symbol"].duplicated(keep="first")]
    expr = expr.set_index("Hugo_Symbol").drop(columns=["Entrez_Gene_Id"], errors="ignore")
    expr.columns = expr.columns.astype(str)
    n_expr_samples = expr.shape[1]
    print(f"\nMETABRIC expression: {expr.shape[0]} genes x {n_expr_samples} samples")
    print(
        f"  value range: {float(expr.min().min()):.3f} .. {float(expr.max().max()):.3f}"
    )
    print("  (Illumina HT-12 log2 intensity; not logged again)")

    clin["SAMPLE_ID"] = clin["SAMPLE_ID"].astype(str)
    clin = clin.set_index("SAMPLE_ID")
    matched = [s for s in expr.columns if s in clin.index]
    print(f"Expression IDs in clinical table: {len(matched)}/{n_expr_samples}")

    # METABRIC is essentially all primary tumors, but filter on SAMPLE_TYPE
    # if present, matching the "primary tumors only" convention used
    # elsewhere in this repo.
    if "SAMPLE_TYPE" in clin.columns:
        is_primary = clin.loc[matched, "SAMPLE_TYPE"].astype(str).str.lower().eq("primary")
        primary_ids = [s for s in matched if bool(is_primary.get(s, True))]
    else:
        primary_ids = matched
    print(f"Primary tumors: {len(primary_ids)}")

    expr_out = expr.loc[:, primary_ids]
    n_nan = int(expr_out.isna().sum().sum())
    print(f"Output matrix: {expr_out.shape[0]} genes x {expr_out.shape[1]} tumors, {n_nan} NaN")

    meta_out = clin.loc[primary_ids, [c for c in CLIN_KEEP if c in clin.columns and c != "SAMPLE_ID"]].copy()
    meta_out.index.name = "SAMPLE_ID"

    expr_out.to_csv(OUT / "tumor_expr_primary.csv")
    meta_out.to_csv(OUT / "tumor_clinical_primary.csv")

    report = [
        f"expression_file={EXPR_PATH}",
        f"n_genes={expr_out.shape[0]}",
        f"n_primary_tumors={expr_out.shape[1]}",
        f"n_nan={n_nan}",
        "full gene set, no KS restriction, no cell-line involvement",
        "no extra log transform (Illumina HT-12 already log2)",
        "",
    ]
    (OUT / "prepare_report.txt").write_text("\n".join(report))
    print("Wrote", OUT / "tumor_expr_primary.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
