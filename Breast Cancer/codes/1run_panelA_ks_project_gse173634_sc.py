#!/usr/bin/env python3
"""Project GSE173634 breast cell-line scRNA-seq into KS Panel A k=3 space.

Question (Hausser Fig. 4 style):
  We fitted a k=3 triangle on DepMap breast cell lines (KS genes, Panel A,
  t-ratio p=0.022). Do single cells from Gambardella et al. (GSE173634)
  fall inside / near that same triangle?

MAGIC (Sahoo et al. 2024 protocol, from their GitHub Rmagic code):
  lib-size normalize + sqrt on the FULL gene x QC-cell count matrix, then
  magic(..., genes=<KS ENSG IDs>, solver="approximate"). Downstream
  batch-match / combined PCA / archetype projection is unchanged.

Set RUN_MAGIC = False to reproduce the pre-MAGIC (log1p CP10k) pipeline.
MAGIC and non-MAGIC runs write to separate results/figures paths.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/1run_panelA_ks_project_gse173634_sc.py"
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
sys.path.insert(0, str(HERE))  # magic_utils.py alongside this script

from src.io import load_expression_csv
from src.pca import align_pca_signs, fit_pca, inverse_transform_scores
from magic_utils import stream_full_sparse_mtx, run_magic_on_ks_genes

CL_MATRIX = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
PANEL_A = BREAST / "results" / "panel_a_ks_genelist"
FIG = BREAST / "figures"

SAHOO = ROOT / "Sahoo"
UMI_PATH = SAHOO / "GSE173634_RAW_UMI_counts.txt"
SERIES_PATH = SAHOO / "GSE173634_series_matrix"
BIOMART = ROOT / (
    "Hausser_Original /Universal cancer tasks, evolutionary tradeoffs, "
    " and the functions of driver mutations/hs_gene_coordinates_biomart.tsv"
)

N_GENES_FILE = 47096
N_CELLS_FILE = 35276
NNZ = 114572483
MIN_UMI = 500

# --- MAGIC (Sahoo / Rmagic defaults) ---
RUN_MAGIC = True
MAGIC_KNN = 5
MAGIC_T = "auto"
MAGIC_N_PCA = 100
MAGIC_SOLVER = "approximate"  # as in Sahoo's Rmagic calls

SUFFIX = "_magic" if RUN_MAGIC else ""
OUT = BREAST / "results" / f"panel_a_ks_genelist_sc_gse173634{SUFFIX}"


def parse_series_subtypes(path: Path) -> dict[str, str]:
    """cell_line -> subtype code (LA/LB/H/TNA/TNB/Basal-like)."""
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
    if not titles or not subtypes:
        return {}
    out = {}
    for t, s in zip(titles, subtypes):
        key = t.rsplit("_", 1)[0] if t[-2:] in ("_1", "_2") else t
        out[key.upper().replace("-", "")] = s
        out[t.upper().replace("-", "")] = s
    return out


def ensembl_to_symbol_map(biomart: Path, needed_symbols: set[str]) -> dict[str, str]:
    """ENSG → KS symbol for genes we care about."""
    df = pd.read_csv(biomart, sep="\t")
    df = df.dropna(subset=["ID", "name"])
    df["ID"] = df["ID"].astype(str)
    df["name"] = df["name"].astype(str)
    want = {s.upper(): s for s in needed_symbols}
    mapping = {}
    for ens, sym in zip(df["ID"], df["name"]):
        if sym.upper() in want:
            mapping.setdefault(ens, want[sym.upper()])
    return mapping


def read_umi_tail_labels(path: Path):
    """Barcodes then gene ENSG IDs appended after MatrixMarket body."""
    with open(path, "rb") as f:
        f.seek(0, 2)
        size = f.tell()
        f.seek(-min(size, 50_000_000), 2)
        data = f.read().decode("utf-8", errors="replace")
    lines = data.splitlines()
    barcodes = lines[-(N_CELLS_FILE + N_GENES_FILE) : -N_GENES_FILE]
    genes = lines[-N_GENES_FILE:]
    if len(barcodes) != N_CELLS_FILE or len(genes) != N_GENES_FILE:
        raise RuntimeError(
            f"Label block size mismatch: {len(barcodes)} barcodes, {len(genes)} genes"
        )
    return barcodes, genes


def load_ks_umi_counts(path: Path, keep_gene_idx: list[int], n_genes_keep: int):
    """Stream MatrixMarket; return dense KS genes × cells and total UMI per cell."""
    row_of = {g + 1: i for i, g in enumerate(keep_gene_idx)}  # MTX is 1-based
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
                print(f"  … {k+1:,}/{NNZ:,} entries")
    print(f"  nonzero KS hits stored: {n_hit:,}")
    return data, total_umi


def log_normalize_with_lib(counts: np.ndarray, lib: np.ndarray, scale: float = 1e4) -> np.ndarray:
    """counts: genes × cells; lib: total UMI per cell → log1p(CP10k)."""
    lib = np.asarray(lib, dtype=float).copy()
    lib[lib <= 0] = np.nan
    return np.log1p(counts / lib * scale)


def draw_triangle(ax, arcs, color="#333333", lw=1.5):
    pts = list(arcs) + [arcs[0]]
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ax.plot(xs, ys, color=color, lw=lw, zorder=5)
    ax.scatter(arcs[:, 0], arcs[:, 1], s=120, c=["#4C78A8", "#F58518", "#E45756"],
               edgecolors="white", linewidths=0.8, zorder=6)


def barycentric(points, vertices):
    a = np.vstack([vertices.T, np.ones((1, vertices.shape[0]))])
    b = np.vstack([points.T, np.ones((1, points.shape[0]))])
    w, *_ = np.linalg.lstsq(a, b, rcond=None)
    return w.T


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    if not UMI_PATH.is_file():
        print(f"Missing {UMI_PATH}")
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
    print(f"MAGIC imputation: {RUN_MAGIC}")

    pca_cl, scores_cl = fit_pca(expr_cl.T.values, n_components=n_pcs)
    pca_cl, scores_cl, _ = align_pca_signs(scores_cl, saved, pca_cl)
    print(f"Rebuilt PCA; max |score diff| vs saved = {np.max(np.abs(scores_cl - saved)):.4g}")

    print("\n=== Map KS symbols ↔ GSE173634 ENSG ===")
    barcodes, ens_genes = read_umi_tail_labels(UMI_PATH)
    ens_to_sym = ensembl_to_symbol_map(BIOMART, set(expr_cl.index.astype(str)))
    keep_idx = []
    keep_symbols = []
    for i, ens in enumerate(ens_genes):
        sym = ens_to_sym.get(ens)
        if sym is not None and sym in expr_cl.index:
            keep_idx.append(i)
            keep_symbols.append(sym)
    seen = {}
    uniq_idx, uniq_sym = [], []
    for i, s in zip(keep_idx, keep_symbols):
        if s not in seen:
            seen[s] = i
            uniq_idx.append(i)
            uniq_sym.append(s)
    uniq_ensg = [ens_genes[i] for i in uniq_idx]
    print(f"KS genes mappable in sc atlas: {len(uniq_sym)}/{expr_cl.shape[0]}")
    missing = [g for g in expr_cl.index.astype(str) if g not in seen]
    if missing:
        print(f"  missing from sc ({len(missing)}): {missing[:12]}{'…' if len(missing)>12 else ''}")
    if len(uniq_sym) < 50:
        print("Too few overlapping genes.")
        return 1

    cache = OUT / "ks_counts_cache.npz"
    if cache.is_file():
        print(f"Loading cached KS counts: {cache}")
        z = np.load(cache, allow_pickle=True)
        counts = z["counts"]
        total_umi = z["total_umi"]
        uniq_sym = list(z["uniq_sym"])
        barcodes = list(z["barcodes"])
        if "uniq_ensg" in z.files:
            uniq_ensg = list(z["uniq_ensg"])
        else:
            # older cache without ENSG list — rebuild from current mapping
            uniq_ensg = [ens_genes[seen[s]] for s in uniq_sym]
    else:
        counts, total_umi = load_ks_umi_counts(UMI_PATH, uniq_idx, len(uniq_sym))
        np.savez_compressed(
            cache,
            counts=counts,
            total_umi=total_umi,
            uniq_sym=np.array(uniq_sym),
            uniq_ensg=np.array(uniq_ensg),
            barcodes=np.array(barcodes),
        )
        print("Cached KS counts to", cache)

    n_ks_genes = (counts > 0).sum(axis=0)
    keep_cells = (total_umi >= MIN_UMI) & (n_ks_genes >= 10)
    print(
        f"Cells after QC (total UMI>={MIN_UMI}, KS genes>=10): "
        f"{int(keep_cells.sum())}/{len(keep_cells)}"
    )
    print(
        f"  total UMI median={np.median(total_umi):.0f}; "
        f"KS genes detected median={np.median(n_ks_genes):.0f}"
    )
    barcodes_k = np.asarray(barcodes)[keep_cells]
    cell_line = (
        pd.Series(barcodes_k)
        .str.split("_")
        .str[0]
        .str.upper()
        .str.replace("-", "", regex=False)
    )

    shared = [g for g in expr_cl.index.astype(str) if g in set(uniq_sym)]
    cl_genes = expr_cl.loc[shared]
    # ENSG ids for shared KS genes (Panel A order)
    sym_to_ensg = dict(zip(uniq_sym, uniq_ensg))
    shared_ensg = [sym_to_ensg[g] for g in shared]

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
            full_mat = stream_full_sparse_mtx(
                UMI_PATH, N_GENES_FILE, N_CELLS_FILE, NNZ, keep_cell_idx
            )  # cells x genes (ENSG columns by position)
            imputed = run_magic_on_ks_genes(
                full_mat,
                all_gene_names=list(ens_genes),
                ks_gene_names=shared_ensg,
                n_pca=MAGIC_N_PCA,
                knn=MAGIC_KNN,
                t=MAGIC_T,
                solver=MAGIC_SOLVER,
            )  # cells x KS-ENSG
            imputed.columns = shared  # rename ENSG → symbol
            imputed.index = barcodes_k
            np.savez_compressed(
                magic_cache,
                imputed=imputed.T.values,  # genes x cells
                genes=np.array(shared),
                barcodes=np.array(barcodes_k),
            )
            print("Cached MAGIC output to", magic_cache)
            sc_genes = imputed.T.loc[shared]
    else:
        log_sc = log_normalize_with_lib(counts[:, keep_cells], total_umi[keep_cells])
        sc_genes = pd.DataFrame(log_sc, index=uniq_sym, columns=barcodes_k).loc[shared]

    print("\n=== Match sc gene mean/SD to Panel A cell lines (batch correction) ===")
    # .copy(): pandas ≥2.0 CoW can hand back a read-only .values view
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
    print(f"Projected {sc_pc.shape[0]} cells into {sc_pc.shape[1]}-D batch-matched PC space")

    w = barycentric(sc_pc[:, :2], arcs_pc[:, :2])
    inside = (w >= -1e-6).all(axis=1)
    print(f"Cells inside triangle: {int(inside.sum())}/{len(inside)} ({100*inside.mean():.1f}%)")

    subtypes = parse_series_subtypes(SERIES_PATH)

    def lookup_subtype(cl):
        for key in (cl, cl.replace("MDAMB", "MDA-MB"), cl):
            if key in subtypes:
                return subtypes[key]
        for kk, v in subtypes.items():
            if kk.replace("-", "") == cl:
                return v
        return "NA"

    subtype = cell_line.map(lookup_subtype)

    meta = pd.DataFrame(
        {
            "barcode": barcodes_k,
            "cell_line": cell_line.values,
            "subtype": subtype.values,
            "PC1": sc_pc[:, 0],
            "PC2": sc_pc[:, 1] if sc_pc.shape[1] > 1 else 0.0,
            "inside_simplex": inside,
            "w_arc1": w[:, 0],
            "w_arc2": w[:, 1],
            "w_arc3": w[:, 2],
            "nearest_archetype": w.argmax(axis=1) + 1,
        }
    )
    meta.to_csv(OUT / "sc_cells_in_panelA_space.csv", index=False)
    pd.DataFrame(scores_cl_share, index=expr_cl.columns.astype(str), columns=["PC1", "PC2"]).to_csv(
        OUT / "cellline_pc_scores_shared_genes.csv"
    )
    pd.DataFrame(arcs_pc, columns=[f"PC{i+1}" for i in range(arcs_pc.shape[1])],
                 index=[f"arc{i+1}" for i in range(k)]).to_csv(OUT / "archetypes_in_shared_pc.csv")
    pd.Series(shared).to_csv(OUT / "shared_ks_genes.csv", index=False, header=["gene"])

    report = {
        "dataset": "GSE173634",
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

    subtype_colors = {
        "LA": "#4C78A8",
        "LB": "#F58518",
        "H": "#E45756",
        "TNA": "#72B7B2",
        "TNB": "#54A24B",
        "Basal-like": "#B279A2",
        "NA": "#CCCCCC",
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    ax = axes[0]
    ax.scatter(scores_cl_share[:, 0], scores_cl_share[:, 1], s=28, c="#888888",
               alpha=0.9, label="DepMap cell lines (Panel A)", zorder=3, edgecolors="white", linewidths=0.3)
    rng = np.random.default_rng(0)
    take = rng.choice(sc_pc.shape[0], size=min(8000, sc_pc.shape[0]), replace=False)
    ax.scatter(sc_pc[take, 0], sc_pc[take, 1], s=2, c="#4C78A8", alpha=0.15, label="GSE173634 sc cells", zorder=2)
    draw_triangle(ax, arcs_pc[:, :2])
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"KS Panel A k=3 triangle\n{int(inside.sum())}/{len(inside)} cells inside ({100*inside.mean():.1f}%)")
    ax.legend(loc="best", fontsize=7, markerscale=1.2)
    ax.set_aspect("equal", adjustable="datalim")

    ax = axes[1]
    for st, col in subtype_colors.items():
        idx = np.where(subtype.values == st)[0]
        if len(idx) == 0:
            continue
        idx_plot = idx if len(idx) <= 3000 else rng.choice(idx, 3000, replace=False)
        ax.scatter(
            sc_pc[idx_plot, 0],
            sc_pc[idx_plot, 1],
            s=3,
            c=col,
            alpha=0.35,
            label=f"{st} (n={len(idx)})",
        )
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
    ax.set_title("Single cells colored by cell-line subtype (GSE173634)")
    ax.legend(loc="best", fontsize=6, markerscale=2)
    ax.set_aspect("equal", adjustable="datalim")

    magic_note = " (MAGIC: libsize+sqrt)" if RUN_MAGIC else " (no MAGIC)"
    fig.suptitle(
        f"Hausser Fig. 4–style check: GSE173634 scRNA-seq in KS Panel A k=3 space{magic_note}\n"
        "(archetypes fitted on DepMap breast lines; cells projected, not refit)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out_fig = FIG / f"Figure_4_ks_panelA_gse173634_sc{SUFFIX}.png"
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
