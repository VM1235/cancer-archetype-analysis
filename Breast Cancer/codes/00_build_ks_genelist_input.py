#!/usr/bin/env python3
"""Subset the official breast Panel A matrix to the Tan et al 2014 KS
(Kolmogorov-Smirnov) generic cell-line EMT signature (Table S1B).

Source: Tan TZ, Miow QH, Miki Y, Noda T, Mori S, Huang RY-J, Thiery JP.
Epithelial-mesenchymal transition spectrum quantification and its efficacy
in deciphering survival and drug responses of cancer patients. EMBO Mol Med
6(10):1279-1293, 2014. PMID 25214461.

This uses the GENERIC cell-line EMT signature (Table S1B: 170 Epi + 48 Mes
= 218 genes), not a breast-specific list -- the breast/ovarian-specific KS
signatures in that paper were derived in an earlier paper (Akalay et al,
Cancer Res 2013) and are not included in this supplementary file. The
generic cell-line signature is the correct match for our DepMap breast
cell-line Panel A matrix (samples are cell lines, not tumours); Table S1A
(315-gene TUMOUR signature) is the correct match for TCGA-BRCA Panel C
instead -- see 00_build_ks_genelist_input_panelC.py.

Does not re-filter or re-log. Writes a new matrix; does not touch panel_a/.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

from src.gene_lists import resolve_symbols, subset_to_genes
from src.io import load_expression_csv

SRC_MATRIX = BREAST / "data" / "processed" / "input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv"
GENE_LIST = BREAST / "data" / "processed" / "ks_cellline_signature_tan2014.csv"
OUT_MATRIX = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
OUT_STATUS = BREAST / "data" / "processed" / "input_panelA_ks_genelist_gene_status.csv"
OUT_REPORT = BREAST / "data" / "processed" / "input_panelA_ks_genelist_build_report.txt"
OUT_CATEGORY = BREAST / "data" / "processed" / "input_panelA_ks_genelist_categories.csv"


def main():
    if not GENE_LIST.is_file():
        print("Missing", GENE_LIST)
        return 1
    import pandas as pd

    genelist = pd.read_csv(GENE_LIST)
    requested = [str(x) for x in genelist["gene_symbol"].tolist()]
    cat_map = dict(zip(genelist["gene_symbol"], genelist["category"]))

    expr = load_expression_csv(SRC_MATRIX)
    status, keep = resolve_symbols(requested, expr.index)
    status.to_csv(OUT_STATUS, index=False)
    sub = subset_to_genes(expr, keep)
    sub.to_csv(OUT_MATRIX)

    # category file for downstream Panel B enrichment / KS scoring, keyed by
    # the matrix symbol actually used (handles alias renames)
    cat_rows = []
    for r in status.itertuples():
        if r.status in ("present_exact", "present_alias"):
            cat_rows.append({
                "gene_symbol": r.matched_matrix_symbol,
                "category": cat_map.get(r.requested_symbol, ""),
            })
    pd.DataFrame(cat_rows).drop_duplicates("gene_symbol").to_csv(OUT_CATEGORY, index=False)

    missing = status.loc[status["status"] == "missing", "requested_symbol"].tolist()
    alias = status.loc[status["status"] == "present_alias"]
    n_epi = sum(1 for c in cat_map.values() if c == "Epi")
    n_mes = sum(1 for c in cat_map.values() if c == "Mes")
    lines = [
        "source=Tan et al 2014 EMBO Mol Med, Table S1B (generic CELL LINE EMT signature)",
        "citation=Tan TZ, Miow QH, Miki Y, Noda T, Mori S, Huang RY-J, Thiery JP. "
        "EMBO Mol Med 6(10):1279-1293, 2014. PMID 25214461.",
        "note=generic signature, not breast-specific (breast-specific KS list is in "
        "Akalay et al 2013 Cancer Res, not in this supplementary file)",
        f"gene_list_file={GENE_LIST}",
        f"n_requested={len(requested)} (Epi={n_epi}, Mes={n_mes})",
        f"n_present_exact={int((status.status=='present_exact').sum())}",
        f"n_present_alias={int((status.status=='present_alias').sum())}",
        f"n_missing={len(missing)}",
        f"missing_genes={','.join(missing) if missing else ''}",
        "alias_maps=" + ";".join(
            f"{r.requested_symbol}->{r.matched_matrix_symbol}" for r in alias.itertuples()
        ),
        f"src_matrix={SRC_MATRIX.name} {expr.shape[0]}x{expr.shape[1]}",
        f"restricted_matrix={OUT_MATRIX.name} {sub.shape[0]}x{sub.shape[1]}",
        "no re-log; no low-expression refilter",
    ]
    OUT_REPORT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("Wrote", OUT_MATRIX)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
