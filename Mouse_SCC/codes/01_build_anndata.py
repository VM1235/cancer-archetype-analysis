import re
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

P = "Mouse_SCC/data/processed/"

print("Reading counts (this takes a few minutes)...", flush=True)
m = scipy.io.mmread(P + "counts.mtx").tocsr()          # genes x cells
genes = [g.strip() for g in open(P + "genes.txt")]
barcodes = [b.strip() for b in open(P + "barcodes.txt")]
meta = pd.read_csv(P + "cell_metadata.csv").set_index("barcode")

assert m.shape == (len(genes), len(barcodes)), "matrix shape does not match gene/barcode lists"
assert list(meta.index) == barcodes, "metadata order does not match barcodes"
assert len(set(barcodes)) == len(barcodes), "duplicate barcodes"

adata = ad.AnnData(X=m.T.tocsr().astype(np.float32), obs=meta,
                   var=pd.DataFrame(index=genes))
adata.var_names_make_unique()
adata.layers["counts"] = adata.X.copy()
print("cells x genes:", adata.shape)

# gene filter: detected in at least 1% of cells, drop mito and ribosomal genes
det = np.asarray((adata.X > 0).sum(0)).ravel() >= 0.01 * adata.n_obs
bad = adata.var_names.str.match(r"^(mt-|Mt-)") | adata.var_names.str.match(r"^(Rpl|Rps|Mrpl|Mrps)\d")
keep = det & ~np.asarray(bad)
print("genes detected in >=1% of cells:", int(det.sum()))
print("mito/ribo genes removed among those:", int((det & bad).sum()))
adata = adata[:, keep].copy()
print("genes kept:", adata.n_vars)

# normalise to 10,000 per cell, then log(x+1)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.write_h5ad(P + "expr_log_filtered.h5ad")

print("\nCells per tumour and EMT state:")
print(pd.crosstab(adata.obs["sample"], adata.obs["cluster_annotation"], margins=True))
print("\nCells per assay:")
print(adata.obs["assay"].value_counts())
print("\nFraction of cells with a label:", adata.obs["cluster_annotation"].notna().mean())
print("Saved", P + "expr_log_filtered.h5ad")
