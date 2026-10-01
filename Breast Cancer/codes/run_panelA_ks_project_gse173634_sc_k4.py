#!/usr/bin/env python3
"""Project GSE173634 scRNA-seq into KS Panel A k=4 tetrahedron.

Uses archetypes from panel_a_ks_genelist_extendedk/ (k=4, p=0.018).
Reuses KS UMI cache from the k=3 sc projection when available.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/run_panelA_ks_project_gse173634_sc_k4.py"
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
CACHE = BREAST / "results" / "panel_a_ks_genelist_sc_gse173634" / "ks_counts_cache.npz"
OUT = BREAST / "results" / "panel_a_ks_genelist_sc_gse173634_k4"
FIG = BREAST / "figures"

SAHOO = ROOT / "Sahoo"
SERIES_PATH = SAHOO / "GSE173634_series_matrix"
BIOMART = ROOT / (
    "Hausser_Original /Universal cancer tasks, evolutionary tradeoffs, "
    " and the functions of driver mutations/hs_gene_coordinates_biomart.tsv"
)
UMI_PATH = SAHOO / "GSE173634_RAW_UMI_counts.txt"

N_GENES_FILE = 47096
N_CELLS_FILE = 35276
NNZ = 114572483
MIN_UMI = 500
K = 4
ARC_COLORS = ["#4C78A8", "#F58518", "#E45756", "#72B7B2"]
SUBTYPE_COLORS = {
    "LA": "#4C78A8",
    "LB": "#F58518",
    "H": "#E45756",
    "TNA": "#72B7B2",
    "TNB": "#54A24B",
    "Basal-like": "#B279A2",
    "NA": "#CCCCCC",
}


def parse_series_subtypes(path: Path) -> dict[str, str]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    titles, subtypes = None, None
    for line in lines:
        if line.startswith("!Sample_title"):
            titles = [t.strip().strip('"') for t in line.split("\t")[1:]]
        if line.startswith("!Sample_characteristics_ch1") and "breast cancer subtype:" in line:
            subtypes = [
                t.strip().strip('"').replace("breast cancer subtype: ", "")
                for t in line.split("\t")[1:]
            ]
    out = {}
    if not titles or not subtypes:
        return out
    for t, s in zip(titles, subtypes):
        key = t.rsplit("_", 1)[0] if t[-2:] in ("_1", "_2") else t
        out[key.upper().replace("-", "")] = s
        out[t.upper().replace("-", "")] = s
    return out


def ensembl_to_symbol_map(biomart: Path, needed_symbols: set[str]) -> dict[str, str]:
    df = pd.read_csv(biomart, sep="\t")
    df = df.dropna(subset=["ID", "name"])
    want = {s.upper(): s for s in needed_symbols}
    mapping = {}
    for ens, sym in zip(df["ID"].astype(str), df["name"].astype(str)):
        if sym.upper() in want:
            mapping.setdefault(ens, want[sym.upper()])
    return mapping


def read_umi_tail_labels(path: Path):
    with open(path, "rb") as f:
        f.seek(0, 2)
        size = f.tell()
        f.seek(-min(size, 50_000_000), 2)
        data = f.read().decode("utf-8", errors="replace")
    lines = data.splitlines()
    barcodes = lines[-(N_CELLS_FILE + N_GENES_FILE) : -N_GENES_FILE]
    genes = lines[-N_GENES_FILE:]
    return barcodes, genes


def load_ks_umi_counts(path: Path, keep_gene_idx: list[int], n_genes_keep: int):
    row_of = {g + 1: i for i, g in enumerate(keep_gene_idx)}
    data = np.zeros((n_genes_keep, N_CELLS_FILE), dtype=np.float32)
    total_umi = np.zeros(N_CELLS_FILE, dtype=np.float64)
    print(f"Streaming UMI matrix for {n_genes_keep} KS genes …")
    with open(path) as f:
        assert f.readline().startswith("%%MatrixMarket")
        dims = f.readline().split()
        assert int(dims[0]) == N_GENES_FILE and int(dims[1]) == N_CELLS_FILE
        for k, line in enumerate(f):
            if k >= NNZ:
                break
            a, b, c = line.split()
            gi, ci, v = int(a), int(b), float(c)
            total_umi[ci - 1] += v
            r = row_of.get(gi)
            if r is not None:
                data[r, ci - 1] = v
            if (k + 1) % 20_000_000 == 0:
                print(f"  … {k+1:,}/{NNZ:,}")
    return data, total_umi


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


def lookup_subtype(cl, subtypes):
    for key in (cl, cl.replace("MDAMB", "MDA-MB")):
        if key in subtypes:
            return subtypes[key]
    for k, v in subtypes.items():
        if k.replace("-", "") == cl:
            return v
    return "NA"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    arcs_path = PANEL_A / f"archetypes_k{K}_parti.npy"
    if not arcs_path.is_file():
        print(f"Missing {arcs_path}")
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

    # Load / build KS sc counts
    if CACHE.is_file():
        print(f"Loading KS count cache: {CACHE}")
        z = np.load(CACHE, allow_pickle=True)
        counts = z["counts"]
        total_umi = z["total_umi"]
        uniq_sym = list(z["uniq_sym"])
        barcodes = list(z["barcodes"])
    else:
        print("No cache; streaming UMI file …")
        barcodes, ens_genes = read_umi_tail_labels(UMI_PATH)
        ens_to_sym = ensembl_to_symbol_map(BIOMART, set(expr_cl.index.astype(str)))
        uniq_idx, uniq_sym, seen = [], [], {}
        for i, ens in enumerate(ens_genes):
            sym = ens_to_sym.get(ens)
            if sym is not None and sym in expr_cl.index and sym not in seen:
                seen[sym] = i
                uniq_idx.append(i)
                uniq_sym.append(sym)
        counts, total_umi = load_ks_umi_counts(UMI_PATH, uniq_idx, len(uniq_sym))
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            CACHE,
            counts=counts,
            total_umi=total_umi,
            uniq_sym=np.array(uniq_sym),
            barcodes=np.array(barcodes),
        )

    n_ks = (counts > 0).sum(axis=0)
    keep = (total_umi >= MIN_UMI) & (n_ks >= 10)
    print(f"Cells after QC: {int(keep.sum())}/{len(keep)}")

    log_sc = log_normalize_with_lib(counts[:, keep], total_umi[keep])
    barcodes_k = np.asarray(barcodes)[keep]
    cell_line = (
        pd.Series(barcodes_k).str.split("_").str[0].str.upper().str.replace("-", "", regex=False)
    )

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

    subtypes = parse_series_subtypes(SERIES_PATH)
    subtype = cell_line.map(lambda x: lookup_subtype(x, subtypes))

    meta = pd.DataFrame(
        {
            "barcode": barcodes_k,
            "cell_line": cell_line.values,
            "subtype": subtype.values,
            "PC1": sc_pc[:, 0],
            "PC2": sc_pc[:, 1],
            "PC3": sc_pc[:, 2],
            "inside_simplex": inside,
            "nearest_archetype": w.argmax(axis=1) + 1,
        }
    )
    for i in range(K):
        meta[f"w_arc{i+1}"] = w[:, i]
    meta.to_csv(OUT / "sc_cells_in_panelA_k4_space.csv", index=False)
    pd.DataFrame(
        scores_lines[:, :3],
        index=expr_cl.columns.astype(str),
        columns=["PC1", "PC2", "PC3"],
    ).to_csv(OUT / "cellline_pc_scores.csv")
    pd.DataFrame(
        arcs_pc[:, :3],
        index=[f"arc{i+1}" for i in range(K)],
        columns=["PC1", "PC2", "PC3"],
    ).to_csv(OUT / "archetypes_pc123.csv")
    pd.Series(
        {
            "k": K,
            "n_shared_genes": len(shared),
            "n_cells": int(sc_pc.shape[0]),
            "frac_inside": float(inside.mean()),
            "panel_a_p": 0.018,
            "n_pcs_panel_a": n_pcs,
        }
    ).to_json(OUT / "projection_report.json")

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
        label="GSE173634 sc",
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
        idx = np.where(subtype.values == st)[0]
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
    ax2.set_title("sc cells by subtype in k=4 space")
    ax2.legend(loc="upper left", fontsize=6)

    fig.suptitle(
        "Hausser Fig. 4–style: GSE173634 sc in KS Panel A k=4 tetrahedron\n"
        "(archetypes from extended Panel A; cells projected, not refit)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out_fig = FIG / "Figure_4_ks_panelA_gse173634_sc_k4.png"
    fig.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Wrote", out_fig)
    print("Wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
