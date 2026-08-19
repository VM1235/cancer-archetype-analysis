"""Curated subtype gene lists and matrix subsetting for the genelist-restricted run.

Do not use these helpers to overwrite official full-transcriptome Panel A/B/C.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# HGNC / common RNA-seq renaming since Parker 2009 and Wang 2017 Table S1.
# Match is applied only after an exact symbol miss, and is logged.
SYMBOL_ALIASES = {
    "CDCA1": "NUF2",
    "KNTC2": "NDC80",
    "ORC6L": "ORC6",
    "HN1": "JPT1",
    "HRASLS": "PLAAT1",
    "PAK7": "PAK5",
    "ZNF643": "ZSCAN31",
    "C14orf159": "DGLUCY",
    "LHFP": "LHFPL6",
    "SEPT11": "SEPTIN11",
    "KIAA0494": "EFCAB14",
    "ACPP": "ACP3",
}


def resolve_symbols(requested, matrix_index):
    """Map a requested gene list onto matrix row names.

    Returns a status DataFrame and the ordered unique matrix names to keep.
    Missing genes are kept in the table; they are never silently omitted from
    the report.
    """
    idx = set(map(str, matrix_index))
    rows = []
    keep = []
    seen = set()
    for gene in requested:
        gene = str(gene).strip()
        if not gene:
            continue
        status = "missing"
        matched = ""
        how = ""
        if gene in idx:
            status, matched, how = "present_exact", gene, "exact"
        elif SYMBOL_ALIASES.get(gene) in idx:
            status = "present_alias"
            matched = SYMBOL_ALIASES[gene]
            how = f"alias:{gene}->{matched}"
        rows.append(
            {
                "requested_symbol": gene,
                "status": status,
                "matched_matrix_symbol": matched,
                "match_how": how,
            }
        )
        if matched and matched not in seen:
            keep.append(matched)
            seen.add(matched)
    return pd.DataFrame(rows), keep


def subset_to_genes(expr, keep_symbols):
    """Subset genes (rows). Sample columns unchanged. No extra filter/log."""
    keep = [g for g in keep_symbols if g in expr.index]
    return expr.loc[keep]


def rename_tumor_genes_to_cell_line(tumor_index, cell_line_index):
    """Rename tumor symbols that match a cell-line gene via HGNC alias.

    Returns (rename_map old_tumor_name -> cell_line_name, unmatched_cell_line_genes).
    """
    cl = set(map(str, cell_line_index))
    tu = set(map(str, tumor_index))
    rename = {}
    for old, new in SYMBOL_ALIASES.items():
        if new in cl and old in tu and new not in tu:
            rename[old] = new
        if old in cl and new in tu and old not in tu:
            rename[new] = old
    matched = (cl & tu) | set(rename.values())
    unmatched = sorted(cl - matched)
    return rename, unmatched