#!/usr/bin/env python3
"""Subset the official breast Panel A matrix to genefu PAM50 genes.

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
GENE_MAP = BREAST / "data" / "processed" / "pam50_genefu_centroids_map.csv"
OUT_MATRIX = BREAST / "data" / "processed" / "input_panelA_pam50_genelist.csv"
OUT_STATUS = BREAST / "data" / "processed" / "input_panelA_pam50_genelist_gene_status.csv"
OUT_REPORT = BREAST / "data" / "processed" / "input_panelA_pam50_genelist_build_report.txt"


def main():
    if not GENE_MAP.is_file():
        print("Missing", GENE_MAP)
        print("Extract with genefu: data(pam50); write.csv(pam50$centroids.map, ...)")
        return 1
    mapping = __import__("pandas").read_csv(GENE_MAP)
    # genefu pam50$centroids.map: probe is the gene symbol used by the centroids.
    if "probe" not in mapping.columns:
        raise ValueError(f"Unexpected columns in {GENE_MAP}: {list(mapping.columns)}")
    requested = [str(x) for x in mapping["probe"].tolist()]
    expr = load_expression_csv(SRC_MATRIX)
    status, keep = resolve_symbols(requested, expr.index)
    status.to_csv(OUT_STATUS, index=False)
    sub = subset_to_genes(expr, keep)
    sub.to_csv(OUT_MATRIX)
    missing = status.loc[status["status"] == "missing", "requested_symbol"].tolist()
    alias = status.loc[status["status"] == "present_alias"]
    lines = [
        "source=genefu::pam50$centroids.map (same 50 genes as pam50.robust)",
        "citation=Parker et al. J Clin Oncol 2009; genefu pam50 centroids",
        f"genefu_map_file={GENE_MAP}",
        f"n_requested={len(requested)}",
        f"n_unique_requested={len(set(requested))}",
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
