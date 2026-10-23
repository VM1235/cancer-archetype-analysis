#!/usr/bin/env python3
"""Step 8d: fit and test a simplex inside ONE tumour, for any gene set (all or KS), raw expression.

Fills the gaps left by 08_within_tumour.py, which ran only SC15, SC32 and SC46 on all genes.
Same method as 08_within_tumour.py --raw: PCA on that tumour's cells only (exact, genes centred,
10 PCs kept), PCHA on the first k-1 PCs (15 starts, delta 0), permutation test (each PC shuffled
independently, 5 starts per shuffle, p = (1 + hits) / (1 + shuffles)). Output goes to
results/within_tumour_gs/<geneset>_<tumour>/ (fit_k*.json, archetypes_k*.npy, null_k*.npy, pca.npz).

Usage (from Project_1/):
  .venv/bin/python -u Mouse_SCC/codes/08d_within_tumour_geneset.py --tumour SC32 --geneset ks --ks 3 4
"""
import argparse
import json
import time

import numpy as np
from scipy.spatial import ConvexHull, QhullError

import scc_common as C


def t_ratio(sc, arcs):
    k = arcs.shape[0]
    return float(C.simplex_volume(arcs) / ConvexHull(sc[:, : k - 1]).volume)


ap = argparse.ArgumentParser()
ap.add_argument("--tumour", required=True)
ap.add_argument("--geneset", choices=["all", "ks", "hallmark"], required=True)
ap.add_argument("--ks", type=int, nargs="+", default=[3, 4])
ap.add_argument("--n-perm", type=int, default=200)
args = ap.parse_args()

a = C.load_adata()
gi = C.gene_index(a, args.geneset)
rows = np.where(a.obs["sample"].astype(str).values == args.tumour)[0]
out = C.RESULTS / "within_tumour_gs" / f"{args.geneset}_{args.tumour}"
out.mkdir(parents=True, exist_ok=True)
Xd = a.X[rows][:, gi].toarray().astype(np.float64)
pca = C.fit_pca(Xd, n_pcs=10)
sc = pca["scores_fit"]
np.savez_compressed(out / "pca.npz", rows=rows, scores=sc, mean=pca["mean"],
                    components=pca["components"], evr=pca["evr"][:50], gene_idx=gi)
for k in args.ks:
    t0 = time.time()
    fit = C.fit_simplex(sc, k)
    obs_t = fit["t_ratio"]
    rng = np.random.default_rng(C.SEED + 1000 + k)
    null = []
    for i in range(args.n_perm):
        sh = C._shuffle_columns(sc[:, : k - 1], rng)
        try:
            np.random.seed(C.SEED + 10_000 * k + i)
            arcs, *_ = C.fit_pcha_best(sh, k, n_init=C.N_INIT_NULL, delta=0.0)
            null.append(t_ratio(sh, arcs))
        except (QhullError, ValueError, RuntimeError):
            null.append(np.nan)
    nl = np.asarray(null)
    fin = nl[np.isfinite(nl)]
    p = (1 + np.sum(fin >= obs_t)) / (1 + len(fin))
    np.save(out / f"null_k{k}.npy", nl)
    np.save(out / f"archetypes_k{k}.npy", fit["archetypes"])
    res = dict(tumour=args.tumour, geneset=args.geneset, k=k, n_cells=int(len(rows)),
               n_genes=int(len(gi)), t_ratio=obs_t, null_median=float(np.median(fin)),
               null_95=float(np.quantile(fin, 0.95)), null_max=float(fin.max()),
               n_null=int(len(fin)), p_value=float(p),
               evr_pc1_pc3=[float(x) for x in pca["evr"][:3]])
    (out / f"fit_k{k}.json").write_text(json.dumps(res, indent=2))
    print(args.geneset, args.tumour, "k", k, "t", round(obs_t, 4), "p", round(float(p), 4),
          f"({time.time() - t0:.0f}s)", flush=True)
