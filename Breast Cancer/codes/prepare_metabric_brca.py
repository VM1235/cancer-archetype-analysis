#!/usr/bin/env python3
"""Prepare METABRIC-BRCA tumors for projection into the KS-restricted Panel A space.

Does not fit archetypes. Restricts METABRIC Illumina HT-12 log2 expression to the
Tan 2014 KS cell-line gene list, reporting symbol matches against Hugo_Symbol.
Writes genes×samples table plus clinical metadata (including CLAUDIN_SUBTYPE from
the patient file). No extra log transform (cBioPortal values are already log2).
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

import pandas as pd

from src.gene_lists import resolve_symbols, subset_to_genes
from src.io import load_expression_csv

EXPR_PATH = BREAST / "brca_metabric" / "data_mrna_illumina_microarray.txt"
CLIN_SAMPLE = BREAST / "brca_metabric" / "data_clinical_sample.txt"
CLIN_PATIENT = BREAST / "brca_metabric" / "data_clinical_patient.txt"
CL_MATRIX = BREAST / "data" / "processed" / "input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv"
KS_GENELIST = BREAST / "data" / "processed" / "ks_cellline_signature_tan2014.csv"
CL_KS_MATRIX = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
OUT = BREAST / "results" / "panel_c_metabric_ks_genelist"

CLIN_KEEP = [
    "PATIENT_ID",
    "SAMPLE_ID",
    "CANCER_TYPE",
    "CANCER_TYPE_DETAILED",
    "ER_STATUS",
    "HER2_STATUS",
    "PR_STATUS",
    "SAMPLE_TYPE",
    "CLAUDIN_SUBTYPE",
]


def read_cbioportal_table(path):
    """Skip # comment lines; return a DataFrame."""
    skip = 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                skip += 1
            else:
                break
    return pd.read_csv(path, sep="\t", skiprows=skip)


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    print("Expression:", EXPR_PATH)
    print("Clinical (sample):", CLIN_SAMPLE)
    print("Clinical (patient):", CLIN_PATIENT)

    # --- Patient ID uniqueness (METABRIC is 1:1 patient:sample) ---
    patient = read_cbioportal_table(CLIN_PATIENT)
    n_pat = len(patient)
    n_unique = patient["PATIENT_ID"].nunique()
    print(f"\nPATIENT_ID check: {n_unique} unique / {n_pat} rows")
    if n_unique != n_pat:
        dups = patient["PATIENT_ID"].value_counts()
        dups = dups[dups > 1]
        raise ValueError(f"Duplicate PATIENT_IDs found: {dups.head().to_dict()}")

    sample_clin = read_cbioportal_table(CLIN_SAMPLE)
    clin = sample_clin.merge(
        patient[["PATIENT_ID", "CLAUDIN_SUBTYPE"]],
        on="PATIENT_ID",
        how="left",
        validate="one_to_one",
    )
    print(f"Joined clinical: {clin.shape[0]} sample rows × {clin.shape[1]} cols")

    for c in ["ER_STATUS", "HER2_STATUS", "PR_STATUS", "CLAUDIN_SUBTYPE"]:
        if c in clin.columns:
            print(f"\n{c}:")
            print(clin[c].value_counts(dropna=False).head(10).to_string())

    # --- Expression matrix (genes × samples) ---
    meta_raw = pd.read_csv(EXPR_PATH, sep="\t")
    meta_raw["Hugo_Symbol"] = meta_raw["Hugo_Symbol"].astype(str)
    meta_raw = meta_raw.loc[~meta_raw["Hugo_Symbol"].duplicated(keep="first")]
    meta_raw = meta_raw.set_index("Hugo_Symbol")
    meta_raw = meta_raw.drop(columns=["Entrez_Gene_Id"], errors="ignore")
    meta_raw.columns = meta_raw.columns.astype(str)
    meta_raw.index = meta_raw.index.astype(str)

    n_expr_samples = meta_raw.shape[1]
    print(f"\nMETABRIC expression: {meta_raw.shape[0]} genes × {n_expr_samples} samples")
    print(f"  sample examples: {list(meta_raw.columns[:3])}")
    print(
        f"  value range: {float(meta_raw.min().min()):.3f} .. "
        f"{float(meta_raw.max().max()):.3f}"
    )
    print("  (Illumina HT-12 log2 intensity; not logged again)")

    # --- KS genelist matching against METABRIC Hugo_Symbol ---
    genelist = pd.read_csv(KS_GENELIST)
    requested = [str(x) for x in genelist["gene_symbol"].tolist()]
    status, keep_metabric = resolve_symbols(requested, meta_raw.index)
    status.to_csv(OUT / "ks_genelist_metabric_gene_status.csv", index=False)

    missing = status.loc[status["status"] == "missing", "requested_symbol"].tolist()
    n_exact = int((status["status"] == "present_exact").sum())
    n_alias = int((status["status"] == "present_alias").sum())
    n_matched = n_exact + n_alias
    print(f"\n=== KS genelist matching (218 requested → METABRIC Hugo_Symbol) ===")
    print(f"  present_exact: {n_exact}")
    print(f"  present_alias: {n_alias}")
    print(f"  matched total: {n_matched}/{len(requested)}")
    print(f"  missing: {len(missing)}")
    if missing:
        print("  missing genes:")
        for g in missing:
            print(f"    - {g}")

    # Intersect with cell-line KS matrix (Panel A restricted genes)
    cl_ks = load_expression_csv(CL_KS_MATRIX)
    cl_ks_genes = cl_ks.index.astype(str)
    keep_final = [g for g in keep_metabric if g in cl_ks_genes]
    absent_from_cl = sorted(set(keep_metabric) - set(cl_ks_genes))
    if absent_from_cl:
        print(f"  matched in METABRIC but absent from cell-line KS matrix: {absent_from_cl}")

    cl_full = load_expression_csv(CL_MATRIX)
    shared_full = cl_full.index.astype(str).intersection(meta_raw.index)
    print(f"\nFull CCLE reference genes: {len(cl_full.index)}")
    print(f"Shared Hugo_Symbol (full CCLE ∩ METABRIC): {len(shared_full)}")
    print(f"Cell-line KS matrix genes: {len(cl_ks_genes)}")
    print(f"Final KS genes for output (METABRIC ∩ cell-line KS): {len(keep_final)}")

    # --- Align samples to clinical ---
    clin["SAMPLE_ID"] = clin["SAMPLE_ID"].astype(str)
    clin = clin.set_index("SAMPLE_ID")
    matched_samples = [s for s in meta_raw.columns if s in clin.index]
    print(f"\nExpression IDs in clinical table: {len(matched_samples)}/{n_expr_samples}")
    if len(matched_samples) != n_expr_samples:
        expr_only = set(meta_raw.columns) - set(clin.index)
        clin_only = set(clin.index) - set(meta_raw.columns)
        print(f"  expression-only samples: {len(expr_only)}")
        print(f"  clinical-only samples: {len(clin_only)}")

    expr_out = subset_to_genes(meta_raw, keep_final)
    expr_out = expr_out.loc[:, matched_samples]
    meta_out = clin.loc[matched_samples, [c for c in CLIN_KEEP if c in clin.columns and c != "SAMPLE_ID"]].copy()
    meta_out.index.name = "SAMPLE_ID"

    n_nan = int(expr_out.isna().sum().sum())
    print(f"\nProcessed tumor matrix: {expr_out.shape[0]} genes × {expr_out.shape[1]} tumors")
    print(f"NaN entries: {n_nan}")
    n_claudin = int(meta_out["CLAUDIN_SUBTYPE"].replace("", pd.NA).notna().sum())
    n_er = int(meta_out["ER_STATUS"].replace("NA", pd.NA).notna().sum())
    print(f"Tumors with CLAUDIN_SUBTYPE: {n_claudin}")
    print(f"Tumors with ER_STATUS: {n_er}")

    expr_out.to_csv(OUT / "metabric_log_shared_genes.csv")
    meta_out.to_csv(OUT / "metabric_clinical_metadata.csv")

    report = [
        f"expression_file={EXPR_PATH}",
        f"clinical_sample={CLIN_SAMPLE}",
        f"clinical_patient={CLIN_PATIENT}",
        f"n_patient_rows={n_pat}",
        f"n_unique_patient_id={n_unique}",
        f"n_genes_metabric={meta_raw.shape[0]}",
        f"n_samples_metabric={n_expr_samples}",
        f"n_ks_requested={len(requested)}",
        f"n_ks_present_exact={n_exact}",
        f"n_ks_present_alias={n_alias}",
        f"n_ks_matched_metabric={n_matched}",
        f"n_ks_missing={len(missing)}",
        f"missing_genes={','.join(missing) if missing else ''}",
        f"n_genes_full_ccle_shared={len(shared_full)}",
        f"n_genes_cell_line_ks={len(cl_ks_genes)}",
        f"n_genes_output={expr_out.shape[0]}",
        f"n_tumors_output={expr_out.shape[1]}",
        f"n_nan={n_nan}",
        "no extra log transform (Illumina HT-12 already log2)",
        "",
    ]
    (OUT / "prepare_report.txt").write_text("\n".join(report))
    print("\nWrote", OUT / "metabric_log_shared_genes.csv")
    print("Wrote", OUT / "metabric_clinical_metadata.csv")
    print("Wrote", OUT / "prepare_report.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
