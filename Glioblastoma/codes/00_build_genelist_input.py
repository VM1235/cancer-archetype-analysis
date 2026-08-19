#!/usr/bin/env python3
"""Subset the official GBM Panel A matrix to Wang 2017 subtype signatures.

Source: Wang et al., Cancer Cell 32:42–56 (2017), Table S1 sheet
'Subtype Signatures' (mmc2.xlsx). 50 upregulated genes each for
Mesenchymal, Proneural, Classical. Neural is not a tumor-intrinsic class
in that paper. Corrected Table S1 headers (MES/PN not swapped).
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
GBM = HERE.parent
ROOT = GBM.parent
sys.path.insert(0, str(ROOT))

import pandas as pd

from src.gene_lists import resolve_symbols, subset_to_genes
from src.io import load_expression_csv

SRC_MATRIX = GBM / "data" / "processed" / "input_panelA_glioblastoma_ccle_logtpm_filtered.csv"
SIG = GBM / "data" / "processed" / "wang2017_tableS1_signatures.csv"
OUT_MATRIX = GBM / "data" / "processed" / "input_panelA_wang2017_genelist.csv"
OUT_STATUS = GBM / "data" / "processed" / "input_panelA_wang2017_genelist_gene_status.csv"
OUT_REPORT = GBM / "data" / "processed" / "input_panelA_wang2017_genelist_build_report.txt"


def main():
    sig = pd.read_csv(SIG)
    counts = sig.groupby("subtype")["gene_symbol_paper"].nunique()
    requested = list(dict.fromkeys(sig["gene_symbol_paper"].astype(str)))
    expr = load_expression_csv(SRC_MATRIX)
    status, keep = resolve_symbols(requested, expr.index)
    status = status.merge(
        sig.rename(columns={"gene_symbol_paper": "requested_symbol"})[
            ["requested_symbol", "subtype", "rank"]
        ],
        on="requested_symbol",
        how="left",
    )
    status.to_csv(OUT_STATUS, index=False)
    sub = subset_to_genes(expr, keep)
    sub.to_csv(OUT_MATRIX)
    missing = status.loc[status["status"] == "missing"]
    alias = status.loc[status["status"] == "present_alias"]
    lines = [
        "source=Wang et al. Cancer Cell 32:42-56 (2017) Table S1 Subtype Signatures",
        "file=Glioblastoma/data/processed/Wang2017_CancerCell_TableS1_mmc2.xlsx",
        "doi=10.1016/j.ccell.2017.06.003",
        "note=50 upregulated genes per CL/MES/PN; Neural omitted (not tumor-intrinsic)",
        f"n_per_subtype={counts.to_dict()}",
        f"n_requested_unique={len(requested)}",
        f"n_present_exact={int((status.status=='present_exact').sum())}",
        f"n_present_alias={int((status.status=='present_alias').sum())}",
        f"n_missing={len(missing)}",
        "missing_genes=" + ",".join(missing["requested_symbol"].astype(str)),
        "missing_by_subtype="
        + str(missing.groupby("subtype")["requested_symbol"].apply(list).to_dict())
        if len(missing)
        else "missing_by_subtype={}",
        "alias_maps="
        + ";".join(f"{r.requested_symbol}->{r.matched_matrix_symbol}" for r in alias.itertuples()),
        f"src_matrix={SRC_MATRIX.name} {expr.shape[0]}x{expr.shape[1]}",
        f"restricted_matrix={OUT_MATRIX.name} {sub.shape[0]}x{sub.shape[1]}",
        "no re-log; no low-expression refilter",
    ]
    OUT_REPORT.write_text("\n".join(str(x) for x in lines) + "\n")
    print("\n".join(str(x) for x in lines))
    print("Wrote", OUT_MATRIX)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
