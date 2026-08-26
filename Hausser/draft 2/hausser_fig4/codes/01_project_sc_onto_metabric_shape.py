#!/usr/bin/env python3
"""Project Karaayvaz scRNA-seq onto METABRIC bulk PCA / ParTI archetypes.

Faithful Python port of the Fig. 4 core in:
  Hausser_Original/.../Karaayvaz2018/scBC.R  (lines ~498–707, 729–767)

Paper Methods (Hausser et al. 2019):
  - After Seurat-style filtering: ~650 cancer cells from 6 tumors
  - Project on first 3 METABRIC PCs using ~1964 genes expressed in ≥50% of
    cancer cells and present in METABRIC
  - Fig 4a: single cells inside the METABRIC tetrahedron
  - Fig 4c: same cells on METABRIC PC1, PC2, PC50

Critical detail from scBC.R (why our first attempt looked wrong):
  Archetypes are re-projected with the SAME common-gene loadings used for
  single cells. Using full-gene archetypes + subset loadings collapses the
  sc cloud to the origin.

Inputs (original Hausser release, not our draft-2 METABRIC prep):
  Karaayvaz2018/brca_metabric/expMatrix.csv
  Karaayvaz2018/brca_metabric/geneListExp.list
  Karaayvaz2018/brca_metabric/arcsOrig_genes.csv
  Karaayvaz2018/GSE118389_counts_rsem.txt

Usage (repo root):
  .venv/bin/python -u "Hausser/draft 2/hausser_fig4/codes/01_project_sc_onto_metabric_shape.py"
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from _paths import ANALYSIS, METABRIC_ORIG, SC_COUNTS, SC_HQ_HEADERS

RESULTS = ANALYSIS / "results" / "projection"
LN2 = float(np.log(2))

# Marker thresholds matching scBC.R cancer-cell selection spirit
# (hematopoietic / CAF / endothelial exclusion → cancer).
MARKER_EXCLUDE = {
    "hematopoietic": ("PTPRC", "CD68", "LYZ", "CD3D"),
    "CAF": ("PDGFRA", "ZEB1", "ACTA2", "FAP"),
    "endothelial": ("PECAM1", "VWF"),
}
MARKER_EPITHELIAL = ("KRT8", "KRT18", "KRT19", "EPCAM")

# Paper Fig. 4a archetype colors (scBC.R arcCols / tCols for breast tetrahedron)
ARCH_COLORS = {
    0: "#007eb1",  # cell division (blue)
    1: "#019e59",  # biomass & energy (green)
    2: "#111111",  # invasion (black)
    3: "#BBBBBB",  # HER2 / remaining (grey) — paper panel uses grey not pink
}


def load_metabric():
    """samples × genes matrix + gene-space archetypes from original ParTI run."""
    gene_names = pd.read_csv(
        METABRIC_ORIG / "geneListExp.list", header=None
    )[0].astype(str).tolist()
    filt_path = METABRIC_ORIG / "geneNamesAfterExprFiltering.list"
    if filt_path.is_file():
        gene_filt = pd.read_csv(filt_path, header=None)[0].astype(str).tolist()
    else:
        gene_filt = gene_names

    exp = pd.read_csv(METABRIC_ORIG / "expMatrix.csv", header=None)
    if exp.shape[1] != len(gene_names):
        raise ValueError(
            f"expMatrix cols ({exp.shape[1]}) != geneListExp ({len(gene_names)})"
        )
    exp.columns = gene_names
    # Keep filtered genes (original run used these for ParTI)
    keep = [g for g in gene_filt if g in exp.columns]
    exp = exp.loc[:, keep]

    arcs = pd.read_csv(METABRIC_ORIG / "arcsOrig_genes.csv", header=None)
    if arcs.shape[1] == len(gene_names):
        arcs.columns = gene_names
        arcs = arcs.loc[:, keep]
    elif arcs.shape[1] == len(keep):
        arcs.columns = keep
    else:
        raise ValueError(
            f"arcsOrig cols ({arcs.shape[1]}) match neither full nor filtered gene list"
        )

    return exp, arcs


def fit_bulk_pca(exp_samples_by_genes: pd.DataFrame, n_components: int):
    """Mean-center genes, no SD scaling — matches ade4::dudi.pca(..., scale=F)."""
    X = exp_samples_by_genes.values.astype(float)
    gene_mean = X.mean(axis=0)
    X0 = X - gene_mean
    pca = PCA(n_components=n_components, svd_solver="full", random_state=0)
    scores = pca.fit_transform(X0)
    # loadings: genes × PCs (ade4 c1)
    loadings = pca.components_.T
    return gene_mean, loadings, scores, pca


def log_normalize_counts(counts: pd.DataFrame, scale_factor: float = 1e5) -> pd.DataFrame:
    """Seurat LogNormalize: ln(1 + counts / library * scale_factor)."""
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log1p(counts.div(lib, axis=1) * scale_factor)


def load_high_quality_cell_ids() -> list[str] | None:
    """Cell IDs from GSE118389_norm_data_headers.txt (Karaayvaz high-quality set)."""
    if not SC_HQ_HEADERS.is_file():
        return None
    # File can be huge; only need the header line.
    with open(SC_HQ_HEADERS, encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n")
    cols = [c.strip().strip('"') for c in header.split("\t")]
    return [c for c in cols[1:] if c]


def prepare_cancer_cells(min_features: int = 1500, max_features: int = 8000):
    """Seurat-like QC + marker filter → cancer cells (paper: ~650)."""
    print("Loading single-cell counts...")
    sc = pd.read_csv(SC_COUNTS, sep="\t")
    sc = sc.rename(columns={sc.columns[0]: "gene"}).set_index("gene")
    sc.index = sc.index.astype(str)
    sc = sc.loc[~sc.index.duplicated(keep="first")]
    print(f"  raw: {sc.shape[0]} genes × {sc.shape[1]} cells")

    hq = load_high_quality_cell_ids()
    if hq:
        keep_cells = [c for c in hq if c in sc.columns]
        print(f"  high-quality header overlap: {len(keep_cells)}/{len(hq)}")
        if keep_cells:
            sc = sc.loc[:, keep_cells]

    # genes in ≥5% cells; cells with 1500–8000 detected genes
    n_cells = sc.shape[1]
    gene_ok = (sc > 0).sum(axis=1) >= max(1, int(round(0.05 * n_cells)))
    sc = sc.loc[gene_ok]
    n_feat = (sc > 0).sum(axis=0)
    cell_ok = (n_feat >= min_features) & (n_feat <= max_features)
    sc = sc.loc[:, cell_ok]
    print(f"  after QC: {sc.shape[0]} genes × {sc.shape[1]} cells")

    log_sc = log_normalize_counts(sc)

    def mean_markers(names):
        present = [g for g in names if g in log_sc.index]
        if not present:
            return pd.Series(0.0, index=log_sc.columns)
        return log_sc.loc[present].mean(axis=0)

    hema = mean_markers(MARKER_EXCLUDE["hematopoietic"])
    caf = mean_markers(MARKER_EXCLUDE["CAF"])
    endo = mean_markers(MARKER_EXCLUDE["endothelial"])
    epi = mean_markers(MARKER_EPITHELIAL)

    # Match scBC.R cluster filters: exclude hematopoietic / CAF / endothelial.
    # Epithelial enrichment is soft (Methods emphasize PTPRC/CAF/PECAM1 exclusion;
    # paper ends at ~650 cancer cells from 6 tumors).
    is_cancer = (hema < 1.0) & (caf < 1.0) & (endo < 1.0)
    if (epi > 0).any():
        # drop clear non-epithelial leftovers if epithelial markers exist
        is_cancer = is_cancer & (epi >= epi.median() * 0.25)
    log_cancer = log_sc.loc[:, is_cancer]
    print(f"  cancer cells kept: {log_cancer.shape[1]} (paper ~650)")
    return log_cancer


def project_items(items_by_genes: np.ndarray, loadings_genes_by_pcs: np.ndarray) -> np.ndarray:
    """items × genes  @  genes × PCs  →  items × PCs."""
    return items_by_genes @ loadings_genes_by_pcs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-gene-cell-frac", type=float, default=0.5)
    ap.add_argument("--n-pcs-bulk", type=int, default=50, help="PCs for 4c (needs ≥50)")
    args = ap.parse_args()

    if not (METABRIC_ORIG / "expMatrix.csv").is_file():
        print(f"Missing original METABRIC at {METABRIC_ORIG}")
        return 1
    if not SC_COUNTS.is_file():
        print(f"Missing {SC_COUNTS}")
        return 1

    RESULTS.mkdir(parents=True, exist_ok=True)

    print("=== METABRIC bulk (original ParTI release) ===")
    exp, arcs = load_metabric()
    print(f"Bulk tumors: {exp.shape[0]} × {exp.shape[1]} genes")
    print(f"Archetypes (gene space): {arcs.shape[0]} × {arcs.shape[1]}")

    gene_mean, loadings, bulk_scores, pca = fit_bulk_pca(exp, n_components=args.n_pcs_bulk)
    gene_index = list(exp.columns)
    print(
        "PCA variance (PC1,2,50 %):",
        (100 * pca.explained_variance_ratio_[[0, 1, min(49, args.n_pcs_bulk - 1)]]).round(2),
    )

    print("\n=== Single-cell cancer cells ===")
    log_cancer = prepare_cancer_cells()
    frac = (log_cancer > 0).mean(axis=1)
    sc_genes = frac[frac >= args.min_gene_cell_frac].index.tolist()
    print(f"Genes in ≥{args.min_gene_cell_frac:.0%} of cancer cells: {len(sc_genes)}")

    common = [g for g in gene_index if g in set(sc_genes)]
    print(f"Common METABRIC ∩ sc genes: {len(common)} (paper ~1964)")
    if len(common) < 500:
        print("Too few common genes.")
        return 1

    g_idx = [gene_index.index(g) for g in common]
    L3 = loadings[np.ix_(g_idx, [0, 1, 2])]
    pcs50 = [0, 1, min(49, args.n_pcs_bulk - 1)]
    L50 = loadings[np.ix_(g_idx, pcs50)]

    # sc: genes × cells → center each gene across cells → cells × genes
    sc_mat = log_cancer.loc[common].values.astype(float)  # genes × cells
    sc_centered = sc_mat - sc_mat.mean(axis=1, keepdims=True)
    # Seurat uses ln; METABRIC is log2 → divide by ln(2) (scBC.R)
    sc_log2 = (sc_centered / LN2).T  # cells × genes

    sc_pc123 = project_items(sc_log2, L3)
    sc_pc1250 = project_items(sc_log2, L50)

    # Archetypes: (arc − bulk gene mean) on common genes, same loadings
    arcs_common = arcs.loc[:, common].values.astype(float)  # k × genes
    mean_common = gene_mean[g_idx]
    arcs0 = arcs_common - mean_common
    arch_pc123 = project_items(arcs0, L3)
    arch_pc1250 = project_items(arcs0, L50)

    # Bulk tumors on same common-gene subspace (optional reference; not drawn in 4a)
    bulk_common = exp.loc[:, common].values.astype(float) - mean_common
    bulk_pc123_common = project_items(bulk_common, L3)

    tumor_ids = pd.Index(log_cancer.columns).to_series().str.replace(r"_.*$", "", regex=True).str.upper()
    meta = pd.DataFrame({"cell_id": log_cancer.columns, "tumor_id": tumor_ids.values})

    pd.DataFrame(sc_pc123, index=log_cancer.columns, columns=["PC1", "PC2", "PC3"]).to_csv(
        RESULTS / "sc_pc123.csv"
    )
    pd.DataFrame(sc_pc1250, index=log_cancer.columns, columns=["PC1", "PC2", "PC50"]).to_csv(
        RESULTS / "sc_pc1_pc2_pc50.csv"
    )
    pd.DataFrame(arch_pc123, columns=["PC1", "PC2", "PC3"]).to_csv(
        RESULTS / "archetypes_pc123.csv", index_label="archetype"
    )
    pd.DataFrame(arch_pc1250, columns=["PC1", "PC2", "PC50"]).to_csv(
        RESULTS / "archetypes_pc1_pc2_pc50.csv", index_label="archetype"
    )
    pd.DataFrame(bulk_scores[:, :3], columns=["PC1", "PC2", "PC3"]).to_csv(
        RESULTS / "bulk_pc123_fullgenes.csv", index_label="tumor"
    )
    pd.DataFrame(bulk_pc123_common, columns=["PC1", "PC2", "PC3"]).to_csv(
        RESULTS / "bulk_pc123_commongenes.csv", index_label="tumor"
    )
    meta.to_csv(RESULTS / "sc_cell_metadata.csv", index=False)

    # Quick sanity: sc should span a meaningful fraction of archetype range
    sc_span = np.ptp(sc_pc123, axis=0)
    arch_span = np.ptp(arch_pc123, axis=0)
    ratio = sc_span / np.maximum(arch_span, 1e-9)
    print(f"sc / archetype axis spans (PC1–3): {ratio.round(3).tolist()}")
    if float(ratio.mean()) < 0.05:
        print("WARNING: sc cloud still tiny vs tetrahedron — check projection.")

    report = {
        "n_bulk_tumors": int(exp.shape[0]),
        "n_bulk_genes": int(exp.shape[1]),
        "n_archetypes": int(arcs.shape[0]),
        "n_sc_cancer_cells": int(log_cancer.shape[1]),
        "n_common_genes": int(len(common)),
        "paper_n_cells": 650,
        "paper_n_common_genes": 1964,
        "sc_span_over_arch_span_pc123": ratio.tolist(),
        "projection": (
            "center each sc gene across cells; divide by ln(2); "
            "multiply by METABRIC PCA loadings on common genes; "
            "re-project gene-space archetypes on the same common genes"
        ),
        "sources": {
            "metabric_exp": str(METABRIC_ORIG / "expMatrix.csv"),
            "arcs": str(METABRIC_ORIG / "arcsOrig_genes.csv"),
            "sc": str(SC_COUNTS),
        },
    }
    (RESULTS / "projection_report.json").write_text(json.dumps(report, indent=2))
    print("Wrote", RESULTS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
