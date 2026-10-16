#!/usr/bin/env python3
"""Step 8: fit a simplex INSIDE one library (or one mouse) and test it with shuffles.

Question (same as the within-tumour step in SCC): is there a triangle inside a single
sample, even though the pooled fit shows none? PCA is done on that unit's cells only
(centred, genes not scaled), PCHA on the first k-1 PCs (15 starts, delta 0), then the
usual shuffle test (each PC shuffled on its own, 5 starts). No expression correction.

Units: a library name (M2-Met, ...), a mouse (M1, M2 = its in-vivo libraries, at most
1000 cells per library), or 'invitro' (both cultured libraries, a negative-ish control).
A unit larger than --max-cells is randomly subsampled (seed 0) to keep runs fast.

Usage (from Project_1/):
  .venv/bin/python -u "Breast Cancer/Mouse_4T1_GSE303505"/codes/08_within_unit.py --unit M2-Met --geneset all --k 3
Output: results/within_unit/<geneset>_<unit>/ (fit_k<k>.json, archetypes, pca.npz, null)
"""
import argparse, json, time
import numpy as np, pandas as pd
from scipy.spatial import ConvexHull, QhullError
import scc_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--unit", required=True)
ap.add_argument("--geneset", choices=["all", "ks"], default="all")
ap.add_argument("--k", type=int, default=3)
ap.add_argument("--n-perm", type=int, default=200)
ap.add_argument("--max-cells", type=int, default=4000)
ap.add_argument("--clean", action="store_true",
                help="drop cells with < 10,000 UMI and technical genes (hepatocyte ambient RNA, "
                     "haemoglobin, Malat1/Neat1/Gm26870/Gm42418/AY036118); outputs get suffix _clean")
ap.add_argument("--decontam", action="store_true",
                help="also drop cells whose immune/stromal panel score is > 0.3 (use with --clean); suffix _clean_dc")
args_min_umi = 10000
PANEL = "Ptprc Lyz2 C1qa C1qb Cd3e Cd79a Itgam Adgre1 Col1a1 Col3a1 Dcn Lum Pecam1".split()   # not expressed by 4T1 tumour cells
DC_CUT = 0.3
TECH = ("Alb Apoa1 Apoa2 Apoh Fga Fgb Fgg Hp Saa1 Saa2 Serpina1a Serpina1b Serpina1c Serpina1d "
        "Serpina1e Hbb-bs Hbb-bt Hba-a1 Hba-a2 Malat1 Neat1 Gm26870 Gm42418 AY036118").split()
args = ap.parse_args()
t0 = time.time()
out = C.RESULTS / "within_unit" / f"{args.geneset}_{args.unit}{'_clean' if args.clean else ''}{'_dc' if args.decontam else ''}"
out.mkdir(parents=True, exist_ok=True)

a = C.load_adata()
obs = a.obs
lib = obs["sample"].astype(str).values
if args.unit in ("M1", "M2"):
    pool = np.where(np.char.startswith(lib.astype(str), args.unit + "-"))[0]
    sub = obs.iloc[pool]
    rows = pool[C.balanced_fit_cells(sub, max_per_tumour=1000)]
elif args.unit == "invivo":          # both mice, all in-vivo libraries, balanced per library
    pool = np.where(obs["assay"].values == "in_vivo")[0]
    sub = obs.iloc[pool]
    rows = pool[C.balanced_fit_cells(sub, max_per_tumour=1000)]
elif args.unit == "invitro":
    rows = np.where(obs["assay"].values == "in_vitro")[0]
else:
    rows = np.where(lib == args.unit)[0]
if args.clean:
    rows = rows[obs["nCount_RNA"].values[rows] >= args_min_umi]
if args.decontam:
    pi = [a.var_names.get_loc(g) for g in PANEL if g in a.var_names]
    score = np.asarray(a.X[rows][:, pi].mean(1)).ravel()
    print(f"decontam: dropping {(score > DC_CUT).sum()} of {len(rows)} cells with panel score > {DC_CUT}", flush=True)
    rows = rows[score <= DC_CUT]
if len(rows) > args.max_cells:
    rows = np.sort(np.random.default_rng(0).choice(rows, args.max_cells, replace=False))
gi = C.gene_index(a, args.geneset)
if args.clean:
    bad = set(TECH)
    gi = np.array([i for i in gi if a.var_names[i] not in bad])
Xd = a.X[rows][:, gi].toarray().astype(np.float64)
genes = a.var_names[gi].tolist()
pca = C.fit_pca(Xd, n_pcs=10)
sc = pca["scores_fit"]
k = args.k
fit = C.fit_simplex(sc, k)
obs_t = fit["t_ratio"]
print(f"{args.unit}/{args.geneset} k={k}: {len(rows)} cells, t={obs_t:.4f}", flush=True)

def t_ratio(s, arcs):
    return float(C.simplex_volume(arcs) / ConvexHull(s[:, : k - 1]).volume)

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
nl = np.asarray(null); fin = nl[np.isfinite(nl)]
p = (1 + np.sum(fin >= obs_t)) / (1 + len(fin))
np.save(out / f"null_k{k}.npy", nl)
np.save(out / f"archetypes_k{k}.npy", fit["archetypes"])
np.savez_compressed(out / "pca.npz", rows=rows, scores=sc, mean=pca["mean"],
                    components=pca["components"], evr=pca["evr"][:50], genes=np.array(genes))
res = dict(unit=args.unit, geneset=args.geneset, k=k, n_cells=int(len(rows)), t_ratio=obs_t,
           null_median=float(np.median(fin)), null_95=float(np.quantile(fin, 0.95)),
           null_max=float(fin.max()), n_null=int(len(fin)), p_value=float(p),
           evr_pc1_pc3=[float(x) for x in pca["evr"][:3]])
(out / f"fit_k{k}.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res), f"({time.time() - t0:.0f}s)")
