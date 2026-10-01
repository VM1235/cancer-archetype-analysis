#!/usr/bin/env python3
"""Project GSE161529 breast tumor scRNA-seq (Pal et al. 2021) into KS Panel A k=3 space.

Hausser Fig. 4-style: archetypes fitted on DepMap breast lines (Panel A);
cancer-epithelial single cells are projected (not refit). Same downstream
logic as the GSE176078 / GSE173634 scripts -- only the loading stage
differs, since GSE161529's raw GEO deposit is 69 per-sample 10x triplets
plus a separate companion cell-type annotation deposit (Figshare), not one
pre-merged matrix.

This script consumes the OUTPUT of the two preprocessing scripts, run first:
  00_extract_metadata_gse161529.R      -> data/gse161529/metadata_gse161529.csv
  00_build_matrix_gse161529.py         -> data/gse161529/processed/*.npz|csv
All four must exist in GSE161529_DIR before running this.

Known data issue -- ER-0001 excluded:
  GSM4909296 ("ER-MH0001" on disk) does not correspond to the same library
  as the Figshare SeuratObject_ERTotalTum.rds "ER_0001" cells: their raw
  10x barcodes overlap at almost exactly the rate you'd expect from two
  *unrelated* libraries drawn from the same whitelist by chance (~30
  observed vs. ~27.6 expected), whereas every other one of the 27 samples
  showed a clean 100% subset match. Provenance couldn't be reconciled from
  what's on GEO, so this sample is dropped rather than silently mismatched.

No cell-state-subset figure (unlike GSE176078's "_subsets" plot): the
Figshare Tum objects only carry a within-subtype numeric Seurat cluster id,
not a named cell-state label the way GSE176078 has celltype_minor. Written
to the output CSV as `subtype_cluster` for future use, but not plotted here
-- same reasoning GSE173634 already followed (cell lines have no
comparable substates either).

MAGIC (Sahoo et al. 2024 protocol): lib-size normalize + sqrt on the FULL
gene x cell matrix, then magic(..., genes=<KS symbols>, solver="approximate").
Set RUN_MAGIC = False to reproduce the pre-MAGIC (log1p CP10k) pipeline.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/1run_panelA_ks_project_gse161529_sc.py"
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

from src.io import load_expression_csv
from src.pca import align_pca_signs, fit_pca, inverse_transform_scores
from magic_utils import run_magic_on_ks_genes

CL_MATRIX = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
PANEL_A = BREAST / "results" / "panel_a_ks_genelist"
FIG = BREAST / "figures"

# Pal et al. 2021 raw GEO + Figshare live under data/gse161529/ (not Sahoo/).
GSE161529_DIR = BREAST / "data" / "gse161529"
COUNTS_PATH = GSE161529_DIR / "processed" / "gse161529_combined_counts.npz"
GENES_PATH = GSE161529_DIR / "processed" / "gse161529_genes.csv"
BARCODES_PATH = GSE161529_DIR / "processed" / "gse161529_barcodes.csv"
META_PATH = GSE161529_DIR / "metadata_gse161529.csv"

MIN_UMI = 500
# See "Known data issue" above -- unreconciled barcode provenance, not a
# biological exclusion.
EXCLUDE_SAMPLES = {"ER-0001"}

# --- MAGIC (Sahoo / Rmagic defaults) ---
RUN_MAGIC = True
MAGIC_KNN = 5
MAGIC_T = "auto"
MAGIC_N_PCA = 100
MAGIC_SOLVER = "approximate"

SUFFIX = "_magic" if RUN_MAGIC else ""
OUT = BREAST / "results" / f"panel_a_ks_genelist_sc_gse161529{SUFFIX}"


def log_normalize_with_lib(counts: np.ndarray, lib: np.ndarray, scale: float = 1e4) -> np.ndarray:
    lib = np.asarray(lib, dtype=float).copy()
    lib[lib <= 0] = np.nan
    return np.log1p(counts / lib * scale)


def draw_triangle(ax, arcs, color="#333333", lw=1.5):
    pts = list(arcs) + [arcs[0]]
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.plot(xs, ys, color=color, lw=lw, zorder=5)
    ax.scatter(
        arcs[:, 0], arcs[:, 1], s=120, c=["#4C78A8", "#F58518", "#E45756"],
        edgecolors="white", linewidths=0.8, zorder=6,
    )


def barycentric(points, vertices):
    a = np.vstack([vertices.T, np.ones((1, vertices.shape[0]))])
    b = np.vstack([points.T, np.ones((1, points.shape[0]))])
    w, *_ = np.linalg.lstsq(a, b, rcond=None)
    return w.T


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    for p in (COUNTS_PATH, GENES_PATH, BARCODES_PATH, META_PATH):
        if not p.is_file():
            print(f"Missing {p}")
            print("Run 00_extract_metadata_gse161529.R and 00_build_matrix_gse161529.py first.")
            return 1
    if not PANEL_A.joinpath("archetypes_k3_parti.npy").is_file():
        print("Missing Panel A k=3 archetypes; run run_panelA_ks_genelist.py first.")
        return 1

    print("=== Panel A cell-line space (k=3) ===")
    expr_cl = load_expression_csv(CL_MATRIX)
    n_pcs = int((PANEL_A / "n_pcs.txt").read_text().strip().splitlines()[0])
    saved = np.load(PANEL_A / "pc_scores.npy")
    arcs = np.load(PANEL_A / "archetypes_k3_parti.npy")  # (3, 2)
    k = arcs.shape[0]
    print(f"Cell lines: {expr_cl.shape[1]} x {expr_cl.shape[0]} KS genes; n_pcs={n_pcs}")
    print(f"Archetypes k={k} in {arcs.shape[1]}-D; Panel A t-ratio p=0.022")
    print(f"MAGIC imputation: {RUN_MAGIC}")

    pca_cl, scores_cl = fit_pca(expr_cl.T.values, n_components=n_pcs)
    pca_cl, scores_cl, _ = align_pca_signs(scores_cl, saved, pca_cl)
    print(f"Rebuilt PCA; max |score diff| vs saved = {np.max(np.abs(scores_cl - saved)):.4g}")

    print("\n=== Load combined matrix + map KS symbols <-> GSE161529 genes ===")
    mat = sp.load_npz(COUNTS_PATH).tocsr()  # genes x cells (all 27 samples, QC-passed)
    genes = pd.read_csv(GENES_PATH)["gene"].astype(str).tolist()
    bc_df = pd.read_csv(BARCODES_PATH)
    barcodes = bc_df["barcode"].astype(str).tolist()
    sample_of = bc_df["sample"].astype(str).to_numpy()
    assert mat.shape == (len(genes), len(barcodes)), "combined matrix / genes / barcodes size mismatch"

    gene_to_idx = {g: i for i, g in enumerate(genes)}
    uniq_idx, uniq_sym = [], []
    for sym in expr_cl.index.astype(str):
        i = gene_to_idx.get(sym)
        if i is not None:
            uniq_idx.append(i)
            uniq_sym.append(sym)
    print(f"KS genes mappable in GSE161529: {len(uniq_sym)}/{expr_cl.shape[0]}")
    missing = [g for g in expr_cl.index.astype(str) if g not in set(uniq_sym)]
    if missing:
        print(f"  missing ({len(missing)}): {missing}")
    if len(uniq_sym) < 50:
        print("Too few overlapping genes.")
        return 1

    print("\n=== Cell-type annotation (Figshare Tum objects) ===")
    meta = pd.read_csv(META_PATH)
    n_excluded = int(meta["sample_name"].isin(EXCLUDE_SAMPLES).sum())
    if n_excluded:
        print(f"Excluding {n_excluded} cells from {sorted(EXCLUDE_SAMPLES)} (unreconciled provenance)")
    meta = meta[~meta["sample_name"].isin(EXCLUDE_SAMPLES)]
    meta_idx = meta.set_index("full_barcode")

    barcodes_arr = np.asarray(barcodes)
    is_epithelial = pd.Index(barcodes_arr).isin(meta_idx.index)
    also_excluded = np.isin(sample_of, list(EXCLUDE_SAMPLES))
    is_epithelial = is_epithelial & ~also_excluded

    cache = OUT / "ks_counts_cache.npz"
    if cache.is_file():
        print(f"Loading cached KS counts: {cache}")
        z = np.load(cache, allow_pickle=True)
        counts_ks_all = z["counts_ks_all"]
        total_umi_all = z["total_umi_all"]
        uniq_sym = list(z["uniq_sym"])
    else:
        counts_ks_all = np.asarray(mat[uniq_idx, :].toarray(), dtype=np.float32)
        total_umi_all = np.asarray(mat.sum(axis=0), dtype=np.float64).ravel()
        np.savez_compressed(
            cache, counts_ks_all=counts_ks_all, total_umi_all=total_umi_all,
            uniq_sym=np.array(uniq_sym),
        )
        print("Cached KS counts to", cache)

    n_ks_genes_all = (counts_ks_all > 0).sum(axis=0)
    keep_cells = is_epithelial & (total_umi_all >= MIN_UMI) & (n_ks_genes_all >= 10)
    print(
        f"Epithelial after QC (UMI>={MIN_UMI}, KS genes>=10): "
        f"{int(keep_cells.sum())}/{int(is_epithelial.sum())} epithelial; "
        f"{int(is_epithelial.sum())}/{len(is_epithelial)} total raw cells"
    )
    barcodes_k = barcodes_arr[keep_cells]
    meta_k = meta_idx.reindex(barcodes_k)
    if meta_k["subtype"].isna().any():
        n_miss = int(meta_k["subtype"].isna().sum())
        print(f"Warning: {n_miss} barcodes missing metadata after reindex")

    shared = [g for g in expr_cl.index.astype(str) if g in set(uniq_sym)]
    cl_genes = expr_cl.loc[shared]

    if RUN_MAGIC:
        magic_cache = OUT / "ks_genes_magic_imputed.npz"
        if magic_cache.is_file():
            print(f"Loading cached MAGIC-imputed KS genes: {magic_cache}")
            z = np.load(magic_cache, allow_pickle=True)
            sc_genes = pd.DataFrame(z["imputed"], index=list(z["genes"]), columns=list(z["barcodes"]))
            sc_genes = sc_genes.loc[shared]
            sc_genes.columns = barcodes_k
        else:
            print("\n=== MAGIC imputation (Sahoo protocol: lib-size + sqrt, full genes) ===")
            keep_cell_idx = np.where(keep_cells)[0]
            full_mat = mat.tocsc()[:, keep_cell_idx].T.tocsr()  # cells x genes
            imputed = run_magic_on_ks_genes(
                full_mat,
                all_gene_names=genes,
                ks_gene_names=shared,
                n_pca=MAGIC_N_PCA,
                knn=MAGIC_KNN,
                t=MAGIC_T,
                solver=MAGIC_SOLVER,
            )
            imputed.index = barcodes_k
            np.savez_compressed(
                magic_cache, imputed=imputed.T.values, genes=np.array(shared),
                barcodes=np.array(barcodes_k),
            )
            print("Cached MAGIC output to", magic_cache)
            sc_genes = imputed.T.loc[shared]
    else:
        log_sc = log_normalize_with_lib(counts_ks_all[:, keep_cells], total_umi_all[keep_cells])
        sc_genes = pd.DataFrame(log_sc, index=uniq_sym, columns=barcodes_k).loc[shared]

    print("\n=== Match sc gene mean/SD to Panel A cell lines (batch correction) ===")
    cl_mean = cl_genes.mean(axis=1).values.copy()
    cl_sd = cl_genes.std(axis=1, ddof=1).values.copy()
    cl_sd[cl_sd < 1e-8] = 1.0
    sc_mean = sc_genes.mean(axis=1).values.copy()
    sc_sd = sc_genes.std(axis=1, ddof=1).values.copy()
    sc_sd[sc_sd < 1e-8] = 1.0
    sc_matched = ((sc_genes.values - sc_mean[:, None]) / sc_sd[:, None]) * cl_sd[:, None] + cl_mean[:, None]
    combined = np.hstack([cl_genes.values, sc_matched])
    is_cell = np.array([True] * cl_genes.shape[1] + [False] * sc_genes.shape[1])

    n_comp = min(n_pcs, combined.shape[1] - 1, combined.shape[0])
    pca_all, scores_all = fit_pca(combined.T, n_components=n_comp)

    scores_cl_share = scores_all[is_cell]
    if scores_cl_share.shape[1] == saved.shape[1]:
        for j in range(scores_cl_share.shape[1]):
            if np.corrcoef(scores_cl_share[:, j], saved[:, j])[0, 1] < 0:
                scores_all[:, j] *= -1
                pca_all.components_[j] *= -1
                scores_cl_share = scores_all[is_cell]

    arcs_full = np.zeros((k, n_pcs))
    arcs_full[:, : arcs.shape[1]] = arcs
    gene_arcs = inverse_transform_scores(pca_cl, arcs_full)
    gene_arcs_df = pd.DataFrame(gene_arcs.T, index=expr_cl.index.astype(str))
    arcs_pc = pca_all.transform(gene_arcs_df.loc[shared].T.values)
    sc_pc = scores_all[~is_cell]
    print(f"Projected {sc_pc.shape[0]} epithelial cells into {sc_pc.shape[1]}-D batch-matched PC space")

    w = barycentric(sc_pc[:, :2], arcs_pc[:, :2])
    inside = (w >= -1e-6).all(axis=1)
    print(f"Cells inside triangle: {int(inside.sum())}/{len(inside)} ({100 * inside.mean():.1f}%)")

    subtype = meta_k["subtype"].astype(str).values
    patient = meta_k["sample_name"].astype(str).values
    subtype_cluster = (meta_k["subtype"].astype(str) + "-" + meta_k["cluster"].astype(str)).values

    out_df = pd.DataFrame(
        {
            "barcode": barcodes_k,
            "patient": patient,
            "subtype": subtype,
            "subtype_cluster": subtype_cluster,  # numeric within-subtype cluster, not a named cell state
            "PC1": sc_pc[:, 0],
            "PC2": sc_pc[:, 1] if sc_pc.shape[1] > 1 else 0.0,
            "inside_simplex": inside,
            "w_arc1": w[:, 0],
            "w_arc2": w[:, 1],
            "w_arc3": w[:, 2],
            "nearest_archetype": w.argmax(axis=1) + 1,
        }
    )
    out_df.to_csv(OUT / "sc_cells_in_panelA_space.csv", index=False)
    pd.DataFrame(
        scores_cl_share, index=expr_cl.columns.astype(str), columns=["PC1", "PC2"]
    ).to_csv(OUT / "cellline_pc_scores_shared_genes.csv")
    pd.DataFrame(
        arcs_pc, columns=[f"PC{i + 1}" for i in range(arcs_pc.shape[1])],
        index=[f"arc{i + 1}" for i in range(k)],
    ).to_csv(OUT / "archetypes_in_shared_pc.csv")
    pd.Series(shared).to_csv(OUT / "shared_ks_genes.csv", index=False, header=["gene"])

    report = {
        "dataset": "GSE161529",
        "cell_filter": "Epithelial (Figshare Tum objects: TNBCTum/HER2Tum/ERTotalTum)",
        "excluded_samples": sorted(EXCLUDE_SAMPLES),
        "magic_imputation": RUN_MAGIC,
        "magic_protocol": "libsize_median+sqrt+approximate" if RUN_MAGIC else None,
        "magic_knn": MAGIC_KNN if RUN_MAGIC else None,
        "magic_n_pca": MAGIC_N_PCA if RUN_MAGIC else None,
        "n_shared_genes": len(shared),
        "n_cells_projected": int(sc_pc.shape[0]),
        "frac_inside_triangle": float(inside.mean()),
        "panel_a_p": 0.022,
        "k": 3,
    }
    pd.Series(report).to_json(OUT / "projection_report.json")
    print("Inside by clinical subtype:")
    for st, g in out_df.groupby("subtype"):
        print(f"  {st}: {g.inside_simplex.mean():.1%} (n={len(g)})")

    subtype_colors = {"ER+": "#4C78A8", "HER2+": "#E45756", "TNBC": "#54A24B"}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    rng = np.random.default_rng(0)

    ax = axes[0]
    ax.scatter(
        scores_cl_share[:, 0], scores_cl_share[:, 1], s=28, c="#888888", alpha=0.9,
        label="DepMap cell lines (Panel A)", zorder=3, edgecolors="white", linewidths=0.3,
    )
    take = rng.choice(sc_pc.shape[0], size=min(8000, sc_pc.shape[0]), replace=False)
    ax.scatter(
        sc_pc[take, 0], sc_pc[take, 1], s=2, c="#4C78A8", alpha=0.15,
        label="GSE161529 epithelial cells", zorder=2,
    )
    draw_triangle(ax, arcs_pc[:, :2])
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(
        f"KS Panel A k=3 triangle\n"
        f"{int(inside.sum())}/{len(inside)} cells inside ({100 * inside.mean():.1f}%)"
    )
    ax.legend(loc="best", fontsize=7, markerscale=1.2)
    ax.set_aspect("equal", adjustable="datalim")

    ax = axes[1]
    for st, col in subtype_colors.items():
        idx = np.where(subtype == st)[0]
        if len(idx) == 0:
            continue
        idx_plot = idx if len(idx) <= 4000 else rng.choice(idx, 4000, replace=False)
        ax.scatter(sc_pc[idx_plot, 0], sc_pc[idx_plot, 1], s=3, c=col, alpha=0.35,
                   label=f"{st} (n={len(idx)})")
    draw_triangle(ax, arcs_pc[:, :2])
    ax.scatter(
        scores_cl_share[:, 0], scores_cl_share[:, 1], s=20, c="black", marker="x",
        linewidths=0.8, label="Panel A lines", zorder=4,
    )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Epithelial cells by clinical subtype (GSE161529)")
    ax.legend(loc="best", fontsize=6, markerscale=2)
    ax.set_aspect("equal", adjustable="datalim")

    magic_note = " (MAGIC: libsize+sqrt)" if RUN_MAGIC else " (no MAGIC)"
    fig.suptitle(
        f"Hausser Fig. 4-style: GSE161529 epithelial sc in KS Panel A k=3{magic_note}\n"
        "(archetypes on DepMap breast lines; cells projected, not refit; ER-0001 excluded)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out_fig = FIG / f"Figure_4_ks_panelA_gse161529_sc{SUFFIX}.png"
    fig.savefig(out_fig, dpi=300)
    plt.close(fig)
    print("Wrote", out_fig)
    print("Wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# Install MAGIC once in the project venv:
#   pip install --no-deps magic-impute graphtools scprep pygsp tasklogger future decorator
#   pip install deprecated
