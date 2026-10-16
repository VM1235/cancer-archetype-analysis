#!/usr/bin/env python3
"""Step 8c (adapted from Cook 08): is a triangle needed beyond a 1-D path?

The shuffle test asks whether data are more triangular than a structureless blob. A smooth path
(epithelial -> mesenchymal) is structure too. Here the null is a PATH:
  1. take the first two PCs (x) of a fit and a one-dimensional coordinate p for each cell
  2. fit a cubic f(p) to each PC; residuals r = x - f(p); path-null x* = f(p) + r[permuted cells]
  3. refit PCHA k=3 on x*, record the t-ratio; compare with the observed t-ratio
If the observed t-ratio is not above the path-null, a path explains the shape.
Cook used its own pseudotime; 4T1 has none, so p = KS_EMT (mouse Tan-2014 mesenchymal minus epithelial z-score),
rank-scaled. This is an EMT-only path: a shape that is not along EMT would still beat it.

Units: a clone (cleaned, decontaminated within-clone fit, results/within_unit/<gs>_<clone>_clean_dc) or
       'pooled' (cleaned pooled fit, results/<gs>_clean, fit cells only; needs EXPR_NAME=expr_tumour_clean.h5ad).
Usage (from Project_1/):
  .venv/bin/python -u "Breast Cancer/GSE114687_pseudospace/codes/08c_path_vs_triangle.py" --unit EM2 --geneset ks
"""
import argparse, json, time
import numpy as np, pandas as pd
from scipy.spatial import ConvexHull, QhullError
from scipy.stats import rankdata
import anndata as ad
import scc_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--unit", required=True)
ap.add_argument("--geneset", default="all")
ap.add_argument("--n-perm", type=int, default=200)
ap.add_argument("--suffix", default="_clean_dc")
args = ap.parse_args()
t0 = time.time()
k = 3
if args.unit == "pooled":
    d = C.RESULTS / f"{args.geneset}_clean"
    z = np.load(d / "pca.npz")
    rows = z["fit_idx"]; X = z["scores_all"][rows][:, :2]
    obs = ad.read_h5ad(C.EXPR, backed="r").obs
    arcs0 = None
else:
    d = C.RESULTS / "within_unit" / f"{args.geneset}_{args.unit}{args.suffix}"
    z = np.load(d / "pca.npz")
    rows = z["rows"]; X = z["scores"][:, :2]
    obs = ad.read_h5ad(C.PROC / "expr_tumour.h5ad", backed="r").obs
p = obs["KS_EMT"].values[rows].astype(float)
u = (rankdata(p) - 0.5) / len(p)
V = np.vander(u, 4)
coef, *_ = np.linalg.lstsq(V, X, rcond=None)
F = V @ coef; R = X - F
r2 = 1 - R.var(0) / X.var(0)

def t_ratio(s, arcs):
    return float(C.simplex_volume(arcs) / ConvexHull(s[:, : k - 1]).volume)

fit = C.fit_simplex(X, k)
obs_t = fit["t_ratio"]
rng = np.random.default_rng(1)
null = []
for i in range(args.n_perm):
    Xs = F + R[rng.permutation(len(R))]
    try:
        np.random.seed(20_000 + i)
        arcs, *_ = C.fit_pcha_best(Xs, k, n_init=C.N_INIT_NULL, delta=0.0)
        null.append(t_ratio(Xs, arcs))
    except (QhullError, ValueError, RuntimeError):
        null.append(np.nan)
nl = np.asarray(null); fin = nl[np.isfinite(nl)]
pv = (1 + np.sum(fin >= obs_t)) / (1 + len(fin))
res = dict(unit=args.unit, geneset=args.geneset, n_cells=int(len(p)), r2_pc1_by_KS_EMT=float(r2[0]),
           r2_pc2_by_KS_EMT=float(r2[1]), t_ratio=obs_t, path_null_median=float(np.median(fin)),
           path_null_95=float(np.quantile(fin, .95)), n_null=int(len(fin)), p_vs_path=float(pv))
out = d / "path_vs_triangle.json"
out.write_text(json.dumps(res, indent=2))
print(json.dumps(res), f"({time.time() - t0:.0f}s)")
