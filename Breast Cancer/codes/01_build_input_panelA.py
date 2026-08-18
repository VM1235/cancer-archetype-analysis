#!/usr/bin/env python3
"""Build Panel A input from DepMap Model.csv + log2(TPM+1) expression.

Output is genes × cell lines, matching Groves' SCLC matrix layout.
Keeps Invasive Breast Carcinoma *cell lines* that have RNA-seq.
Drops genes with log2(TPM+1) < 1 in every line (TPM < 1 everywhere).
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

import pandas as pd

from src.preprocess import drop_low_genes, expression_to_genes_by_samples, pick_breast_models

RAW = BREAST / "data" / "raw"
PROCESSED = BREAST / "data" / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)
OUT_MATRIX = PROCESSED / "input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv"
OUT_MODELS = PROCESSED / "input_panelA_models_used.csv"
OUT_REPORT = PROCESSED / "input_panelA_build_report.txt"


def find_file(folder, names):
    for name in names:
        hit = folder / name
        if hit.is_file():
            return hit
    return None


def main():
    expr_path = find_file(
        RAW,
        [
            "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        ],
    )
    model_path = find_file(RAW, ["Model.csv"])
    if expr_path is None or model_path is None:
        print("Missing files in", RAW)
        return 1

    models = pd.read_csv(model_path)
    id_col = next(
        c for c in ("ModelID", "DepMap_ID", "depmap_id") if c in models.columns
    )
    mask, rule = pick_breast_models(models, prefer_invasive=True)
    if "ModelType" in models.columns:
        mask = mask & models["ModelType"].astype(str).str.contains("cell line", case=False, na=False)
        rule = rule + "; ModelType=Cell Line"

    breast_models = models.loc[mask].copy()
    breast_ids = set(breast_models[id_col].astype(str))
    print(f"Metadata filter ({rule}): {len(breast_ids)} models")
    print(f"Expression: {expr_path.name}")
    print(f"Models:     {model_path.name}")

    expr = pd.read_csv(expr_path, index_col=0)
    expr = expression_to_genes_by_samples(expr)
    keep_ids = [c for c in expr.columns if str(c) in breast_ids]
    missing = sorted(breast_ids - set(map(str, expr.columns)))
    print(f"Of those, {len(keep_ids)} have RNA-seq")
    if missing:
        print(f"  dropped {len(missing)} metadata models with no expression")

    subset = expr[keep_ids]
    n_before = subset.shape[0]
    subset = drop_low_genes(subset, max_log_threshold=1.0)
    subset = subset.loc[~subset.index.duplicated(keep="first")]
    print(f"Genes: {n_before} → {subset.shape[0]} after dropping all-low genes")
    print(f"Panel A matrix: {subset.shape[0]} genes × {subset.shape[1]} cell lines")

    subset.to_csv(OUT_MATRIX)
    used = breast_models[breast_models[id_col].astype(str).isin(keep_ids)].copy()
    used.to_csv(OUT_MODELS, index=False)
    OUT_REPORT.write_text(
        "\n".join(
            [
                f"expression_file: {expr_path}",
                f"model_file: {model_path}",
                f"id_column: {id_col}",
                f"filter: {rule}",
                f"n_models_in_metadata: {len(breast_ids)}",
                f"n_with_rnaseq: {len(keep_ids)}",
                f"n_genes_before_filter: {n_before}",
                f"n_genes_after_filter: {subset.shape[0]}",
                f"output_matrix: {OUT_MATRIX.name}",
                f"output_models: {OUT_MODELS.name}",
                "",
            ]
        )
    )
    print(f"Wrote {OUT_MATRIX}")
    print(f"Wrote {OUT_MODELS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
