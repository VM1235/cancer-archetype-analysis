#!/usr/bin/env python3
"""Build the breast Panel A matrix from DepMap expression + Model.csv."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import pandas as pd

from src.preprocess import drop_low_genes, expression_to_genes_by_samples, pick_breast_models

RAW = ROOT / "Breast Cancer"
BREAST = ROOT / "data" / "breast"
OUT_MATRIX = BREAST / "breast_ccle_logtpm_filtered.csv"
OUT_MODELS = BREAST / "breast_models_used.csv"
OUT_REPORT = BREAST / "build_report.txt"


def find_file(patterns):
    for folder in (RAW, BREAST):
        if not folder.exists():
            continue
        for pat in patterns:
            hits = sorted(
                p
                for p in folder.glob(pat)
                if p.is_file() and ".download" not in p.parts
            )
            if hits:
                return hits[0]
    return None


def main():
    BREAST.mkdir(parents=True, exist_ok=True)
    expr_path = find_file(
        [
            "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
            "*Expression*TPMLogp1*.csv",
            "*OmicsExpression*.csv",
        ]
    )
    model_path = find_file(["Model.csv", "*Model*.csv", "sample_info.csv"])
    if expr_path is None or model_path is None:
        raw_names = [p.name for p in RAW.iterdir()] if RAW.exists() else []
        print("Waiting on the expression file in Breast Cancer/.")
        print("Have Model.csv already." if model_path else "Need Model.csv too.")
        print("Need the finished file (not the .download folder):")
        print("  OmicsExpressionProteinCodingGenesTPMLogp1.csv")
        print(f"Currently in Breast Cancer/: {raw_names}")
        return 1

    print(f"Expression: {expr_path.name}")
    print(f"Models:     {model_path.name}")
    models = pd.read_csv(model_path)
    id_col = next(
        c for c in ("ModelID", "DepMap_ID", "depmap_id", "CCLE_Name") if c in models.columns
    )
    mask, rule = pick_breast_models(models)
    breast_models = models.loc[mask].copy()
    breast_ids = set(breast_models[id_col].astype(str))
    print(f"Metadata filter ({rule}): {len(breast_ids)} models")

    expr = pd.read_csv(expr_path, index_col=0)
    expr = expression_to_genes_by_samples(expr)
    keep_ids = [c for c in expr.columns if str(c) in breast_ids]
    missing = sorted(breast_ids - set(map(str, expr.columns)))
    print(f"Of those, {len(keep_ids)} have RNA-seq in the expression file")
    if missing:
        print(f"  {len(missing)} metadata models have no expression (dropped)")

    subset = expr[keep_ids]
    n_before = subset.shape[0]
    subset = drop_low_genes(subset, max_log_threshold=1.0)
    subset = subset.loc[~subset.index.duplicated(keep="first")]
    print(f"Genes: {n_before} → {subset.shape[0]} after dropping all-low genes")
    print(f"Panel A matrix: {subset.shape[0]} genes × {subset.shape[1]} cell lines")

    subset.to_csv(OUT_MATRIX)
    used = breast_models[breast_models[id_col].astype(str).isin(keep_ids)]
    used.to_csv(OUT_MODELS, index=False)
    OUT_REPORT.write_text(
        "\n".join(
            [
                f"expression_file: {expr_path.name}",
                f"model_file: {model_path.name}",
                f"id_column: {id_col}",
                f"filter: {rule}",
                f"n_models_in_metadata: {len(breast_ids)}",
                f"n_with_rnaseq: {len(keep_ids)}",
                f"n_genes_before_filter: {n_before}",
                f"n_genes_after_filter: {subset.shape[0]}",
                f"output: {OUT_MATRIX.name}",
                "",
            ]
        )
    )
    print(f"Wrote {OUT_MATRIX}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
