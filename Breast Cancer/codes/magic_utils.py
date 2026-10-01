"""MAGIC imputation for Panel A scRNA-seq projection scripts.

Matches Sahoo et al. 2024 (iScience) / their GitHub Rmagic usage:

    expr <- t(counts)                         # cells x genes
    expr <- library.size.normalize(expr)
    expr <- sqrt(expr)
    magic(expr, genes=<genes of interest>, solver="approximate")

Paper Methods only say count matrices were MAGIC-imputed before scoring;
their deposited code (Single_cell_RNA_Seq_cellline/code_MAGIC_AUCell.R,
Spatial_data_analysis/spatial_analysis_pilot.R, receptor_ligand_analysis/
code_version_R.R) is the concrete protocol: lib-size normalize + sqrt on the
FULL gene matrix, then Rmagic with solver='approximate'. Graph PCA uses the
default n_pca=100; genes= only controls which genes are returned.

We mirror that in Python via scprep + magic-impute. Downstream Panel A
batch-matching absorbs the sqrt-vs-log1p scale difference.

Install (Python 3.12 venvs often need --no-deps because scprep still pins
pandas<2.1):

    pip install --no-deps magic-impute graphtools scprep pygsp tasklogger future decorator
    pip install deprecated   # graphtools runtime dep skipped by --no-deps
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def stream_full_sparse_mtx(
    path: Path,
    n_genes_file: int,
    n_cells_file: int,
    nnz: int,
    keep_cell_idx: np.ndarray,
) -> sp.csr_matrix:
    """Stream MatrixMarket (genes x cells); keep all genes, only keep_cell_idx.

    Returns cells x genes CSR in the same cell order as sorted(keep_cell_idx),
    which matches ``np.where(keep_cells)[0]`` / ``barcodes[keep_cells]``.
    """
    keep_cell_idx = np.asarray(keep_cell_idx, dtype=np.int64)
    # 1-based MTX column ids → dense position among kept cells (ascending)
    col_remap = {int(c) + 1: i for i, c in enumerate(np.sort(keep_cell_idx))}
    n_keep = len(col_remap)

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    print(
        f"Streaming FULL matrix, keeping {n_keep}/{n_cells_file} cells "
        f"x {n_genes_file} genes …"
    )
    with open(path) as f:
        assert f.readline().startswith("%%MatrixMarket")
        dims = f.readline().split()
        assert int(dims[0]) == n_genes_file and int(dims[1]) == n_cells_file
        for k, line in enumerate(f):
            if k >= nnz:
                break
            a, b, c = line.split()
            gi, ci, v = int(a), int(b), float(c)
            j = col_remap.get(ci)
            if j is not None:
                rows.append(gi - 1)
                cols.append(j)
                vals.append(v)
            if (k + 1) % 20_000_000 == 0:
                print(f"  … {k + 1:,}/{nnz:,} entries")
    mat_genes_x_cells = sp.csr_matrix(
        (vals, (rows, cols)), shape=(n_genes_file, n_keep), dtype=np.float32
    )
    print(f"  kept {mat_genes_x_cells.nnz:,} nonzeros")
    return mat_genes_x_cells.T.tocsr()  # cells x genes


def run_magic_on_ks_genes(
    counts_cells_x_genes: sp.csr_matrix,
    all_gene_names: list[str],
    ks_gene_names: list[str],
    *,
    n_pca: int = 100,
    knn: int = 5,
    t: str | int = "auto",
    solver: str = "approximate",
    random_state: int = 0,
    n_jobs: int = -1,
    rescale: str | float = "median",
) -> pd.DataFrame:
    """Lib-size normalize + sqrt the full matrix, MAGIC-impute, return KS genes.

    Parameters mirror Sahoo's Rmagic calls. ``rescale="median"`` matches R
    ``library.size.normalize`` (scprep's default 10000 is also fine; median is
    closer to their R pipeline).

    Returns
    -------
    DataFrame, shape (n_cells, n_ks_genes_found), columns = ks gene names
    in the order they appear in ``ks_gene_names``.
    """
    import magic
    import scprep

    if len(all_gene_names) != counts_cells_x_genes.shape[1]:
        raise ValueError(
            f"all_gene_names length {len(all_gene_names)} != "
            f"n_genes {counts_cells_x_genes.shape[1]}"
        )

    name_to_idx = {g: i for i, g in enumerate(all_gene_names)}
    missing = [g for g in ks_gene_names if g not in name_to_idx]
    if missing:
        raise ValueError(
            f"{len(missing)} KS genes absent from MAGIC matrix "
            f"(e.g. {missing[:5]})"
        )
    ks_idx = [name_to_idx[g] for g in ks_gene_names]

    print(
        f"Library-size normalize (rescale={rescale!r}) + sqrt on "
        f"{counts_cells_x_genes.shape[0]} cells x {counts_cells_x_genes.shape[1]} genes …"
    )
    data = scprep.normalize.library_size_normalize(
        counts_cells_x_genes, rescale=rescale
    )
    data = scprep.transform.sqrt(data)

    print(
        f"Running MAGIC (knn={knn}, t={t}, n_pca={n_pca}, solver={solver!r}) "
        f"returning {len(ks_idx)} KS genes …"
    )
    magic_op = magic.MAGIC(
        knn=knn,
        t=t,
        n_pca=n_pca,
        solver=solver,
        random_state=random_state,
        n_jobs=n_jobs,
        verbose=True,
    )
    imputed = magic_op.fit_transform(data, genes=ks_idx)
    arr = np.asarray(imputed, dtype=np.float64)
    if arr.ndim != 2:
        raise RuntimeError(f"Unexpected MAGIC output shape: {getattr(arr, 'shape', None)}")
    return pd.DataFrame(arr, columns=list(ks_gene_names))
