#!/usr/bin/env python3
"""Step 1e: depth-matched copy of the cells (a control for sequencing depth).

Takes the cells that pass the >= 10,000 UMI cleaning, re-reads their raw counts, and thins every cell to the same
total (TARGET UMIs) by binomial downsampling of each gene count (p = TARGET / cell total). Then the same
normalisation as step 1 (CP10k, log1p) on the same genes. All labels (arm, position, KS_EMT, states, phase) are copied
from expr_tumour.h5ad. Every cell then has the same depth by construction (nCount_RNA ~ TARGET).
Output: data/processed/expr_tumour_depthmatched.h5ad  (use with EXPR_NAME=expr_tumour_depthmatched.h5ad)
"""
import importlib.util, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc
HERE = Path(__file__).resolve().parent.parent
P = HERE / "data" / "processed"
TARGET = 10000
spec = importlib.util.spec_from_file_location("b", Path(__file__).parent / "01_build_anndata_funcs.py")
b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
ref = ad.read_h5ad(P / "expr_tumour.h5ad")
keep = ref.obs_names[ref.obs["nCount_RNA"].values >= 10000]
ref = ref[keep].copy()
rng = np.random.default_rng(0)
parts = []
for line, prefix, meta, header in [("MCF10A", "pseudospace", "GSE114687_pseudospace_metadata.tsv.gz", False),
                                   ("HuMEC", "HuMEC", "GSE114687_HuMEC_pseudospace_metadata.tsv.gz", True)]:
    raw = b.load_line(line, prefix, meta, header)
    raw = raw[[c for c in raw.obs_names if c in set(ref.obs_names)]].copy()
    X = raw.X.tocsr().astype(np.int64)
    tot = np.asarray(X.sum(1)).ravel()
    p = np.repeat(TARGET / tot, np.diff(X.indptr))
    X.data = rng.binomial(X.data, np.clip(p, 0, 1))
    X.eliminate_zeros()
    gi = [raw.var_names.get_loc(g) for g in ref.var_names]
    parts.append(ad.AnnData(X=X[:, gi].astype(np.float32), obs=pd.DataFrame(index=raw.obs_names), var=pd.DataFrame(index=ref.var_names)))
    print(line, parts[-1].n_obs, "cells thinned to", TARGET, flush=True)
a = ad.concat(parts, join="inner")
a = a[ref.obs_names].copy()
a.obs = ref.obs.copy()
a.obs["nCount_orig"] = ref.obs["nCount_RNA"].values
a.obs["nCount_RNA"] = np.asarray(a.X.sum(1)).ravel()
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
a.write_h5ad(P / "expr_tumour_depthmatched.h5ad")
print(a.shape, "median thinned UMI", int(np.median(a.obs["nCount_RNA"])))
