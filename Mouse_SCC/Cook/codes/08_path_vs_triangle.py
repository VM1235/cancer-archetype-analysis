#!/usr/bin/env python3
"""Step 8 (Cook): is a triangle needed beyond a 1-D path?

The shuffle test asks whether data are more triangular than a structureless blob. But a smooth path
(epithelial -> mesenchymal) is structure too, and a curved path can look triangular. Here the null is
a PATH:
  1. take the cleaned cells of one line, their first two PCs (x) and the dataset's own pseudotime (p)
  2. fit a smooth curve f(p) to each PC (cubic polynomial in rank-scaled pseudotime)
  3. residuals r = x - f(p); path-null data: x* = f(p) + r[permuted cells]  (whole residual vectors permuted)
  4. refit PCHA k=3 on x*, record the t-ratio; compare with the observed t-ratio.
If the observed t-ratio is not above the path-null, a path explains the shape and no extra corner is needed.
Also reports how much of PC1 and PC2 the path explains (R^2).

Usage (from Project_1/):
  .venv/bin/python -u Mouse_SCC/Cook/codes/08_path_vs_triangle.py --line A549 --geneset all
"""
import argparse, json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.spatial import ConvexHull, QhullError
from scipy.stats import rankdata
import cook_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--line", required=True)
ap.add_argument("--geneset", default="all")
ap.add_argument("--variant", default="clean")
ap.add_argument("--n-perm", type=int, default=200)
args = ap.parse_args()
t0 = time.time()
d = C.RESULTS / "within_line" / f"{args.geneset}_{args.line}_{args.variant}"
z = np.load(d / "pca.npz")
a = C.load_adata()
p = a.obs["pseudotime"].values[z["rows"]].astype(float)
X = z["scores"][:, :2]
u = (rankdata(p) - 0.5) / len(p)                     # rank-scaled pseudotime in (0,1)
V = np.vander(u, 4)                                  # cubic
coef, *_ = np.linalg.lstsq(V, X, rcond=None)
F = V @ coef
R = X - F
r2 = 1 - R.var(0) / X.var(0)
k = 3
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
res = dict(line=args.line, geneset=args.geneset, n_cells=int(len(p)), r2_pc1=float(r2[0]), r2_pc2=float(r2[1]),
           t_ratio=obs_t, path_null_median=float(np.median(fin)), path_null_95=float(np.quantile(fin, .95)),
           n_null=int(len(fin)), p_vs_path=float(pv))
(d / "path_vs_triangle.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res), f"({time.time() - t0:.0f}s)")
