#!/usr/bin/env python3
"""Project Wu et al. GSE176078 BRCA tumor scRNA-seq into KS Panel A k=3 space.

Hausser Fig. 4–style: archetypes fitted on DepMap breast lines (Panel A);
cancer-epithelial single cells are projected (not refit).

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/run_panelA_ks_project_gse176078_sc.py"
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))  # for magic_utils.py placed alongside this script

from src.io import load_expression_csv
from src.pca import align_pca_signs, fit_pca, inverse_transform_scores
from magic_utils import stream_full_sparse_mtx, run_magic_on_ks_genes

USE_MAGIC = True  # set False to fall back to the original (no-impute) path
MAGIC_KNN = 5
MAGIC_T = "auto"
MAGIC_N_PCA = 100

CL_MATRIX = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
PANEL_A = BREAST / "results" / "panel_a_ks_genelist"
OUT = BREAST / "results" / "panel_a_ks_genelist_sc_gse176078"
FIG = BREAST / "figures"

SC_DIR = ROOT / "Sahoo" / "GSE176078_sc" / "Wu_etal_2021_BRCA_scRNASeq"
MTX_PATH = SC_DIR / "count_matrix_sparse.mtx"
GENES_PATH = SC_DIR / "count_matrix_genes.tsv"
BARCODES_PATH = SC_DIR / "count_matrix_barcodes.tsv"
META_PATH = SC_DIR / "metadata.csv"

N_GENES_FILE = 29733
N_CELLS_FILE = 100064
NNZ = 177994136
MIN_UMI = 500
# Wu atlas: project malignant epithelium into cancer archetypes
KEEP_MAJOR = {"Cancer Epithelial"}


def load_ks_counts(path: Path, keep_gene_idx: list[int], n_genes_keep: int):
    """Stream MatrixMarket; return dense KS genes × cells and total UMI per cell."""
    row_of = {g + 1: i for i, g in enumerate(keep_gene_idx)}  # MTX 1-based
    data = np.zeros((n_genes_keep, N_CELLS_FILE), dtype=np.float32)
    total_umi = np.zeros(N_CELLS_FILE, dtype=np.float64)
    print(f"Streaming UMI matrix for {n_genes_keep} KS genes (+ library sizes) …")
    with open(path) as f:
        assert f.readline().startswith("%%MatrixMarket")
        dims = f.readline().split()
        assert int(dims[0]) == N_GENES_FILE and int(dims[1]) == N_CELLS_FILE
        n_hit = 0
        for k, line in enumerate(f):
            if k >= NNZ:
                break
            a, b, c = line.split()
            gi, ci, v = int(a), int(b), float(c)
            total_umi[ci - 1] += v
            r = row_of.get(gi)
            if r is not None:
                data[r, ci - 1] = v
                n_hit += 1
            if (k + 1) % 20_000_000 == 0:
                print(f"  … {k + 1:,}/{NNZ:,} entries")
    print(f"  nonzero KS hits stored: {n_hit:,}")
    return data, total_umi


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
        arcs[:, 0],
        arcs[:, 1],
        s=120,
        c=["#4C78A8", "#F58518", "#E45756"],
        edgecolors="white",
        linewidths=0.8,
        zorder=6,
    )


def barycentric(points, vertices):
    a = np.vstack([vertices.T, np.ones((1, vertices.shape[0]))])
    b = np.vstack([points.T, np.ones((1, points.shape[0]))])
    w, *_ = np.linalg.lstsq(a, b, rcond=None)
    return w.T


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    for p in (MTX_PATH, GENES_PATH, BARCODES_PATH, META_PATH):
        if not p.is_file():
            print(f"Missing {p}")
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
    print(f"Cell lines: {expr_cl.shape[1]} × {expr_cl.shape[0]} KS genes; n_pcs={n_pcs}")
    print(f"Archetypes k={k} in {arcs.shape[1]}-D; Panel A t-ratio p=0.022")

    pca_cl, scores_cl = fit_pca(expr_cl.T.values, n_components=n_pcs)
    pca_cl, scores_cl, _ = align_pca_signs(scores_cl, saved, pca_cl)
    print(f"Rebuilt PCA; max |score diff| vs saved = {np.max(np.abs(scores_cl - saved)):.4g}")

    print("\n=== Map KS symbols ↔ GSE176078 genes ===")
    genes = pd.read_csv(GENES_PATH, header=None)[0].astype(str).tolist()
    barcodes = pd.read_csv(BARCODES_PATH, header=None)[0].astype(str).tolist()
    assert len(genes) == N_GENES_FILE and len(barcodes) == N_CELLS_FILE

    gene_to_idx = {g: i for i, g in enumerate(genes)}
    uniq_idx, uniq_sym = [], []
    for sym in expr_cl.index.astype(str):
        i = gene_to_idx.get(sym)
        if i is not None:
            uniq_idx.append(i)
            uniq_sym.append(sym)
    print(f"KS genes mappable in Wu atlas: {len(uniq_sym)}/{expr_cl.shape[0]}")
    missing = [g for g in expr_cl.index.astype(str) if g not in set(uniq_sym)]
    if missing:
        print(f"  missing ({len(missing)}): {missing}")
    if len(uniq_sym) < 50:
        print("Too few overlapping genes.")
        return 1

    meta = pd.read_csv(META_PATH)
    # barcode column may be Unnamed: 0
    if "Unnamed: 0" in meta.columns:
        meta = meta.rename(columns={"Unnamed: 0": "barcode"})
    meta = meta.set_index("barcode").reindex(barcodes)
    if meta["celltype_major"].isna().any():
        n_miss = int(meta["celltype_major"].isna().sum())
        print(f"Warning: {n_miss} barcodes missing metadata after reindex")

    # --- Pass 1: KS-gene counts + library size, same as the original script.
    # We still need this pass regardless of MAGIC: total_umi/n_ks_genes drive
    # the QC mask, and computing them requires touching every nonzero entry
    # anyway, so there's no cheaper way to learn which cells to keep.
    cache = OUT / "ks_counts_cache.npz"
    if cache.is_file():
        print(f"Loading cached KS counts: {cache}")
        z = np.load(cache, allow_pickle=True)
        counts = z["counts"]
        total_umi = z["total_umi"]
        uniq_sym = list(z["uniq_sym"])
        barcodes = list(z["barcodes"])
    else:
        counts, total_umi = load_ks_counts(MTX_PATH, uniq_idx, len(uniq_sym))
        np.savez_compressed(
            cache,
            counts=counts,
            total_umi=total_umi,
            uniq_sym=np.array(uniq_sym),
            barcodes=np.array(barcodes),
        )
        print("Cached KS counts to", cache)

    n_ks_genes = (counts > 0).sum(axis=0)
    is_cancer = meta["celltype_major"].isin(KEEP_MAJOR).fillna(False).to_numpy()
    keep_cells = is_cancer & (total_umi >= MIN_UMI) & (n_ks_genes >= 10)
    print(
        f"Cancer epithelial after QC (UMI>={MIN_UMI}, KS genes>=10): "
        f"{int(keep_cells.sum())}/{int(is_cancer.sum())} cancer; "
        f"{int(is_cancer.sum())}/{len(is_cancer)} total atlas"
    )

    barcodes_k = np.asarray(barcodes)[keep_cells]
    meta_k = meta.loc[barcodes_k]
    shared = [g for g in expr_cl.index.astype(str) if g in set(uniq_sym)]
    cl_genes = expr_cl.loc[shared]

    if not USE_MAGIC:
        log_sc = log_normalize_with_lib(counts[:, keep_cells], total_umi[keep_cells])
        sc_genes = pd.DataFrame(log_sc, index=uniq_sym, columns=barcodes_k).loc[shared]
    else:
        # --- Pass 2 (MAGIC path): stream the FULL gene x cell matrix, but
        # only for the cells that already passed QC above -- this is the
        # gene-axis restriction we deliberately avoided during Pass 1.
        magic_cache = OUT / "ks_genes_magic_imputed.parquet"
        if magic_cache.is_file():
            print(f"Loading cached MAGIC-imputed KS genes: {magic_cache}")
            sc_genes = pd.read_parquet(magic_cache).T  # stored cells x genes
            sc_genes = sc_genes.loc[shared]
        else:
            keep_cell_idx = np.where(keep_cells)[0]
            full_mat = stream_full_sparse_mtx(
                MTX_PATH, N_GENES_FILE, N_CELLS_FILE, NNZ, keep_cell_idx
            )  # cells x genes, only QC-passing cells
            imputed = run_magic_on_ks_genes(
                full_mat,
                all_gene_symbols=genes,
                ks_gene_symbols=shared,
                n_pca=MAGIC_N_PCA,
                knn=MAGIC_KNN,
                t=MAGIC_T,
            )  # cells x KS-genes, index 0..n_keep-1 (order matches barcodes_k)
            imputed.index = barcodes_k
            imputed.T.to_parquet(magic_cache)  # store genes x cells like the raw path
            sc_genes = imputed.T.loc[shared]
        sc_genes.columns = barcodes_k

    print("\n=== Match sc gene mean/SD to Panel A cell lines (batch correction) ===")
    cl_mean = cl_genes.mean(axis=1).values
    cl_sd = cl_genes.std(axis=1, ddof=1).values
    cl_sd[cl_sd < 1e-8] = 1.0
    sc_mean = sc_genes.mean(axis=1).values
    sc_sd = sc_genes.std(axis=1, ddof=1).values
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
    print(f"Projected {sc_pc.shape[0]} cancer cells into {sc_pc.shape[1]}-D batch-matched PC space")

    w = barycentric(sc_pc[:, :2], arcs_pc[:, :2])
    inside = (w >= -1e-6).all(axis=1)
    print(f"Cells inside triangle: {int(inside.sum())}/{len(inside)} ({100 * inside.mean():.1f}%)")

    subtype = meta_k["subtype"].astype(str).values
    subset = meta_k["celltype_minor"].astype(str).values

    out_df = pd.DataFrame(
        {
            "barcode": barcodes_k,
            "patient": meta_k["orig.ident"].astype(str).values,
            "subtype": subtype,
            "celltype_minor": subset,
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
        arcs_pc,
        columns=[f"PC{i + 1}" for i in range(arcs_pc.shape[1])],
        index=[f"arc{i + 1}" for i in range(k)],
    ).to_csv(OUT / "archetypes_in_shared_pc.csv")
    pd.Series(shared).to_csv(OUT / "shared_ks_genes.csv", index=False, header=["gene"])

    report = {
        "dataset": "GSE176078",
        "cell_filter": "Cancer Epithelial",
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
    subset_colors = {
        "Cancer LumA SC": "#4C78A8",
        "Cancer LumB SC": "#F58518",
        "Cancer Her2 SC": "#E45756",
        "Cancer Basal SC": "#B279A2",
        "Cancer Cycling": "#72B7B2",
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    rng = np.random.default_rng(0)

    ax = axes[0]
    ax.scatter(
        scores_cl_share[:, 0],
        scores_cl_share[:, 1],
        s=28,
        c="#888888",
        alpha=0.9,
        label="DepMap cell lines (Panel A)",
        zorder=3,
        edgecolors="white",
        linewidths=0.3,
    )
    take = rng.choice(sc_pc.shape[0], size=min(8000, sc_pc.shape[0]), replace=False)
    ax.scatter(
        sc_pc[take, 0],
        sc_pc[take, 1],
        s=2,
        c="#4C78A8",
        alpha=0.15,
        label="GSE176078 cancer epithelial",
        zorder=2,
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
        ax.scatter(
            sc_pc[idx_plot, 0],
            sc_pc[idx_plot, 1],
            s=3,
            c=col,
            alpha=0.35,
            label=f"{st} (n={len(idx)})",
        )
    # also mark other rare labels if any
    other = [s for s in np.unique(subtype) if s not in subtype_colors]
    for st in other:
        idx = np.where(subtype == st)[0]
        ax.scatter(sc_pc[idx, 0], sc_pc[idx, 1], s=3, c="#CCCCCC", alpha=0.3, label=f"{st} (n={len(idx)})")
    draw_triangle(ax, arcs_pc[:, :2])
    ax.scatter(
        scores_cl_share[:, 0],
        scores_cl_share[:, 1],
        s=20,
        c="black",
        marker="x",
        linewidths=0.8,
        label="Panel A lines",
        zorder=4,
    )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Cancer epithelial cells by clinical subtype")
    ax.legend(loc="best", fontsize=6, markerscale=2)
    ax.set_aspect("equal", adjustable="datalim")

    fig.suptitle(
        "Hausser Fig. 4–style: GSE176078 cancer epithelial sc in KS Panel A k=3\n"
        "(archetypes on DepMap breast lines; cells projected, not refit)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out_fig = FIG / "Figure_4_ks_panelA_gse176078_sc.png"
    fig.savefig(out_fig, dpi=300)
    plt.close(fig)
    print("Wrote", out_fig)

    # Optional second figure: cancer cell-state subsets
    fig2, ax = plt.subplots(figsize=(7, 6))
    plotted = set()
    for st, col in subset_colors.items():
        idx = np.where(subset == st)[0]
        if len(idx) == 0:
            continue
        idx_plot = idx if len(idx) <= 3000 else rng.choice(idx, 3000, replace=False)
        ax.scatter(
            sc_pc[idx_plot, 0],
            sc_pc[idx_plot, 1],
            s=4,
            c=col,
            alpha=0.4,
            label=f"{st} (n={len(idx)})",
        )
        plotted.add(st)
    rest = [s for s in np.unique(subset) if s not in plotted]
    for st in rest:
        idx = np.where(subset == st)[0]
        idx_plot = idx if len(idx) <= 1500 else rng.choice(idx, 1500, replace=False)
        ax.scatter(
            sc_pc[idx_plot, 0],
            sc_pc[idx_plot, 1],
            s=3,
            c="#BBBBBB",
            alpha=0.25,
            label=f"{st} (n={len(idx)})",
        )
    draw_triangle(ax, arcs_pc[:, :2])
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(
        f"GSE176078 cancer epithelial by cell-state subset\n"
        f"{int(inside.sum())}/{len(inside)} inside ({100 * inside.mean():.1f}%)"
    )
    ax.legend(loc="best", fontsize=6, markerscale=2)
    ax.set_aspect("equal", adjustable="datalim")
    fig2.tight_layout()
    out_fig2 = FIG / "Figure_4_ks_panelA_gse176078_sc_subsets.png"
    fig2.savefig(out_fig2, dpi=300)
    plt.close(fig2)
    print("Wrote", out_fig2)
    print("Wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
