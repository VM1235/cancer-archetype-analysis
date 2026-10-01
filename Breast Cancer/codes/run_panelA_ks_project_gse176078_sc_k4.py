#!/usr/bin/env python3
"""Project Wu et al. GSE176078 cancer epithelial sc into KS Panel A k=4 tetrahedron.

Reuses KS UMI cache from run_panelA_ks_project_gse176078_sc.py.
Archetypes from panel_a_ks_genelist_extendedk (k=4, p=0.018); not refit.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/run_panelA_ks_project_gse176078_sc_k4.py"
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

from src.io import load_expression_csv
from src.pca import align_pca_signs, fit_pca, inverse_transform_scores

CL_MATRIX = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
PANEL_A = BREAST / "results" / "panel_a_ks_genelist_extendedk"
CACHE = BREAST / "results" / "panel_a_ks_genelist_sc_gse176078" / "ks_counts_cache.npz"
OUT = BREAST / "results" / "panel_a_ks_genelist_sc_gse176078_k4"
FIG = BREAST / "figures"

META_PATH = ROOT / "Sahoo" / "GSE176078_sc" / "Wu_etal_2021_BRCA_scRNASeq" / "metadata.csv"

MIN_UMI = 500
K = 4
KEEP_MAJOR = {"Cancer Epithelial"}
ARC_COLORS = ["#4C78A8", "#F58518", "#E45756", "#72B7B2"]
SUBTYPE_COLORS = {"ER+": "#4C78A8", "HER2+": "#E45756", "TNBC": "#54A24B"}


def log_normalize_with_lib(counts, lib, scale=1e4):
    lib = np.asarray(lib, dtype=float).copy()
    lib[lib <= 0] = np.nan
    return np.log1p(counts / lib * scale)


def barycentric(points, vertices):
    a = np.vstack([vertices.T, np.ones((1, vertices.shape[0]))])
    b = np.vstack([points.T, np.ones((1, points.shape[0]))])
    w, *_ = np.linalg.lstsq(a, b, rcond=None)
    return w.T


def draw_tetrahedron(ax, arcs, colors):
    k = arcs.shape[0]
    for i in range(k):
        for j in range(i + 1, k):
            ax.plot(
                [arcs[i, 0], arcs[j, 0]],
                [arcs[i, 1], arcs[j, 1]],
                [arcs[i, 2], arcs[j, 2]],
                color="#666666",
                lw=1.2,
                alpha=0.9,
            )
    ax.scatter(
        arcs[:, 0],
        arcs[:, 1],
        arcs[:, 2],
        s=160,
        c=colors[:k],
        edgecolors="white",
        linewidths=0.8,
        depthshade=False,
        zorder=10,
    )


def fit_3d_view(ax, pts, pad=0.08):
    spans = np.ptp(pts, axis=0)
    spans = np.maximum(spans, 1e-6)
    for i, setter in enumerate([ax.set_xlim, ax.set_ylim, ax.set_zlim]):
        lo, hi = pts[:, i].min(), pts[:, i].max()
        p = spans[i] * pad
        setter(lo - p, hi + p)
    try:
        ax.set_box_aspect(spans / spans.max())
    except AttributeError:
        pass


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    arcs_path = PANEL_A / f"archetypes_k{K}_parti.npy"
    if not arcs_path.is_file():
        print(f"Missing {arcs_path}")
        return 1
    if not CACHE.is_file():
        print(f"Missing cache {CACHE}; run k=3 GSE176078 script first.")
        return 1

    print("=== Panel A extended k=4 tetrahedron ===")
    expr_cl = load_expression_csv(CL_MATRIX)
    n_pcs = int((PANEL_A / "n_pcs.txt").read_text().strip().splitlines()[0])
    saved = np.load(PANEL_A / "pc_scores.npy")
    arcs = np.load(arcs_path)  # (4, 3)
    n_vol = arcs.shape[1]
    print(f"Cell lines: {expr_cl.shape[1]} × {expr_cl.shape[0]} genes; n_pcs={n_pcs}")
    print(f"Archetypes k={K} in {n_vol}-D; Panel A extended t-ratio p=0.018")

    pca_cl, scores_cl = fit_pca(expr_cl.T.values, n_components=n_pcs)
    pca_cl, scores_cl, _ = align_pca_signs(scores_cl, saved, pca_cl)
    print(f"Rebuilt PCA; max |score diff| vs saved = {np.max(np.abs(scores_cl - saved)):.4g}")

    print(f"Loading KS count cache: {CACHE}")
    z = np.load(CACHE, allow_pickle=True)
    counts = z["counts"]
    total_umi = z["total_umi"]
    uniq_sym = list(z["uniq_sym"])
    barcodes = list(z["barcodes"])

    meta = pd.read_csv(META_PATH)
    if "Unnamed: 0" in meta.columns:
        meta = meta.rename(columns={"Unnamed: 0": "barcode"})
    meta = meta.set_index("barcode").reindex(barcodes)

    n_ks = (counts > 0).sum(axis=0)
    is_cancer = meta["celltype_major"].isin(KEEP_MAJOR).fillna(False).to_numpy()
    keep = is_cancer & (total_umi >= MIN_UMI) & (n_ks >= 10)
    print(
        f"Cancer epithelial after QC: {int(keep.sum())}/{int(is_cancer.sum())} cancer; "
        f"{int(is_cancer.sum())}/{len(is_cancer)} atlas"
    )

    log_sc = log_normalize_with_lib(counts[:, keep], total_umi[keep])
    barcodes_k = np.asarray(barcodes)[keep]
    meta_k = meta.loc[barcodes_k]

    shared = [g for g in expr_cl.index.astype(str) if g in set(uniq_sym)]
    sc_genes = pd.DataFrame(log_sc, index=uniq_sym, columns=barcodes_k).loc[shared]
    cl_genes = expr_cl.loc[shared]
    print(f"Shared KS genes: {len(shared)}")

    print("=== Match sc gene mean/SD to cell lines + combined PCA ===")
    cl_mean = cl_genes.mean(axis=1).values
    cl_sd = cl_genes.std(axis=1, ddof=1).values
    cl_sd[cl_sd < 1e-8] = 1.0
    sc_mean = sc_genes.mean(axis=1).values
    sc_sd = sc_genes.std(axis=1, ddof=1).values
    sc_sd[sc_sd < 1e-8] = 1.0
    sc_matched = ((sc_genes.values - sc_mean[:, None]) / sc_sd[:, None]) * cl_sd[:, None] + cl_mean[:, None]
    combined = np.hstack([cl_genes.values, sc_matched])
    is_line = np.array([True] * cl_genes.shape[1] + [False] * sc_genes.shape[1])

    n_comp = min(max(n_pcs, n_vol), combined.shape[1] - 1, combined.shape[0])
    pca_all, scores_all = fit_pca(combined.T, n_components=n_comp)

    scores_lines = scores_all[is_line]
    n_align = min(scores_lines.shape[1], saved.shape[1])
    for j in range(n_align):
        if np.corrcoef(scores_lines[:, j], saved[:, j])[0, 1] < 0:
            scores_all[:, j] *= -1
            pca_all.components_[j] *= -1
            scores_lines = scores_all[is_line]

    arcs_full = np.zeros((K, n_pcs))
    arcs_full[:, :n_vol] = arcs
    gene_arcs = inverse_transform_scores(pca_cl, arcs_full)
    gene_arcs_df = pd.DataFrame(gene_arcs.T, index=expr_cl.index.astype(str))
    arcs_pc = pca_all.transform(gene_arcs_df.loc[shared].T.values)
    sc_pc = scores_all[~is_line]
    print(f"Projected {sc_pc.shape[0]} cells; archetypes in {arcs_pc.shape[1]}-D")

    w = barycentric(sc_pc[:, :n_vol], arcs_pc[:, :n_vol])
    inside = (w >= -1e-6).all(axis=1)
    print(f"Cells inside tetrahedron: {int(inside.sum())}/{len(inside)} ({100 * inside.mean():.1f}%)")

    subtype = meta_k["subtype"].astype(str).values
    out_df = pd.DataFrame(
        {
            "barcode": barcodes_k,
            "patient": meta_k["orig.ident"].astype(str).values,
            "subtype": subtype,
            "celltype_minor": meta_k["celltype_minor"].astype(str).values,
            "PC1": sc_pc[:, 0],
            "PC2": sc_pc[:, 1],
            "PC3": sc_pc[:, 2],
            "inside_simplex": inside,
            "nearest_archetype": w.argmax(axis=1) + 1,
        }
    )
    for i in range(K):
        out_df[f"w_arc{i + 1}"] = w[:, i]
    out_df.to_csv(OUT / "sc_cells_in_panelA_k4_space.csv", index=False)
    pd.DataFrame(
        scores_lines[:, :3],
        index=expr_cl.columns.astype(str),
        columns=["PC1", "PC2", "PC3"],
    ).to_csv(OUT / "cellline_pc_scores.csv")
    pd.DataFrame(
        arcs_pc[:, :3],
        index=[f"arc{i + 1}" for i in range(K)],
        columns=["PC1", "PC2", "PC3"],
    ).to_csv(OUT / "archetypes_pc123.csv")
    pd.Series(
        {
            "dataset": "GSE176078",
            "cell_filter": "Cancer Epithelial",
            "k": K,
            "n_shared_genes": len(shared),
            "n_cells": int(sc_pc.shape[0]),
            "frac_inside": float(inside.mean()),
            "panel_a_p": 0.018,
            "n_pcs_panel_a": n_pcs,
        }
    ).to_json(OUT / "projection_report.json")

    print("Inside by clinical subtype:")
    for st, g in out_df.groupby("subtype"):
        print(f"  {st}: {g.inside_simplex.mean():.1%} (n={len(g)})")

    rng = np.random.default_rng(0)
    take = rng.choice(sc_pc.shape[0], size=min(10000, sc_pc.shape[0]), replace=False)

    fig = plt.figure(figsize=(13, 6))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax1.scatter(
        scores_lines[:, 0],
        scores_lines[:, 1],
        scores_lines[:, 2],
        s=25,
        c="#888888",
        alpha=0.95,
        depthshade=True,
        label="DepMap lines",
        edgecolors="white",
        linewidths=0.2,
    )
    ax1.scatter(
        sc_pc[take, 0],
        sc_pc[take, 1],
        sc_pc[take, 2],
        s=3,
        c="#4C78A8",
        alpha=0.12,
        depthshade=True,
        label="GSE176078 cancer epithelial",
        linewidths=0,
    )
    draw_tetrahedron(ax1, arcs_pc[:, :3], ARC_COLORS)
    fit_3d_view(ax1, np.vstack([scores_lines[:, :3], arcs_pc[:, :3], sc_pc[take, :3]]))
    ax1.view_init(elev=18, azim=-60)
    ax1.set_xlabel("PC1")
    ax1.set_ylabel("PC2")
    ax1.set_zlabel("PC3")
    ax1.set_title(
        f"KS Panel A k=4 tetrahedron\n"
        f"{int(inside.sum())}/{len(inside)} cells inside ({100 * inside.mean():.1f}%)"
    )
    ax1.legend(loc="upper left", fontsize=7)

    ax2 = fig.add_subplot(1, 2, 2, projection="3d")
    for st, col in SUBTYPE_COLORS.items():
        idx = np.where(subtype == st)[0]
        if len(idx) == 0:
            continue
        idx_plot = idx if len(idx) <= 2500 else rng.choice(idx, 2500, replace=False)
        ax2.scatter(
            sc_pc[idx_plot, 0],
            sc_pc[idx_plot, 1],
            sc_pc[idx_plot, 2],
            s=3,
            c=col,
            alpha=0.3,
            depthshade=True,
            linewidths=0,
            label=f"{st} (n={len(idx)})",
        )
    other = [s for s in np.unique(subtype) if s not in SUBTYPE_COLORS]
    for st in other:
        idx = np.where(subtype == st)[0]
        ax2.scatter(
            sc_pc[idx, 0],
            sc_pc[idx, 1],
            sc_pc[idx, 2],
            s=3,
            c="#CCCCCC",
            alpha=0.25,
            depthshade=True,
            linewidths=0,
            label=f"{st} (n={len(idx)})",
        )
    draw_tetrahedron(ax2, arcs_pc[:, :3], ARC_COLORS)
    ax2.scatter(
        scores_lines[:, 0],
        scores_lines[:, 1],
        scores_lines[:, 2],
        s=18,
        c="black",
        marker="x",
        linewidths=0.7,
        label="Panel A lines",
        depthshade=False,
    )
    fit_3d_view(ax2, np.vstack([scores_lines[:, :3], arcs_pc[:, :3], sc_pc[take, :3]]))
    ax2.view_init(elev=18, azim=-60)
    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")
    ax2.set_zlabel("PC3")
    ax2.set_title("Cancer epithelial by clinical subtype")
    ax2.legend(loc="upper left", fontsize=6)

    fig.suptitle(
        "Hausser Fig. 4–style: GSE176078 cancer epithelial in KS Panel A k=4 tetrahedron\n"
        "(archetypes from extended Panel A; cells projected, not refit)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out_fig = FIG / "Figure_4_ks_panelA_gse176078_sc_k4.png"
    fig.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Wrote", out_fig)
    print("Wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
