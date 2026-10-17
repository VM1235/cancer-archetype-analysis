#!/usr/bin/env python3
"""Step 7 (Cook revisit): per-line simplex fits after cleaning, optionally balanced by time.

Why: the first Cook fits used every forward-arm cell. The PDAC work showed that low-depth cells
(and technical genes) can hide a within-sample simplex, and Cook's mito fraction explains 0.28-0.31
of PC1 in A549 and DU145. Late time points are also scarce (MCF7 7d: 22 cells), so a late corner
carries little weight.

Variants
  clean        forward arm; drop cells with < 10,000 UMI or mito fraction > 0.07
  clean_tbal   as 'clean', then at most --tcap (default 300) cells per time point
Everything else as in the first run: exact PCA on those cells (centred), PCHA on the first k-1 PCs
(15 starts, delta 0), shuffle each PC on its own (5 starts), 200 shuffles.

Usage (from Project_1/):
  .venv/bin/python -u Mouse_SCC/Cook/codes/07_within_line_clean.py --line A549 --geneset all --k 3 --variant clean
Output: results/within_line/<geneset>_<line>_<variant>/
"""
import argparse, json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
from scipy.spatial import ConvexHull, QhullError
import cook_common as C

ap = argparse.ArgumentParser()
ap.add_argument("--line", required=True)
ap.add_argument("--geneset", choices=["all", "ks"], default="all")
ap.add_argument("--k", type=int, default=3)
ap.add_argument("--variant", choices=["clean", "clean_tbal"], default="clean")
ap.add_argument("--tcap", type=int, default=300)
ap.add_argument("--n-perm", type=int, default=200)
ap.add_argument("--drop-human-tech", action="store_true",
                help="also drop human technical genes (MALAT1, NEAT1, haemoglobin, liver ambient); output suffix _ht")
TECH_HUMAN = "MALAT1 NEAT1 HBB HBA1 HBA2".split()   # human symbols (cultured cells: no liver ambient genes); the mouse list above matches nothing in human data
args = ap.parse_args()
t0 = time.time()
out = C.RESULTS / "within_line" / f"{args.geneset}_{args.line}_{args.variant}{'_ht' if args.drop_human_tech else ''}"
out.mkdir(parents=True, exist_ok=True)

a = C.load_adata()
o = a.obs
ok = ((o["cell_line"].values == args.line) & (~o["removal_arm"].values)
      & (o["nCount_RNA"].values >= 10000) & (o["pct_mito"].values <= 0.07))
rows = np.where(ok)[0]
if args.variant == "clean_tbal":
    df = pd.DataFrame({"i": rows, "t": o["time"].astype(str).values[rows]})
    rows = np.sort(pd.concat([g.sample(min(len(g), args.tcap), random_state=0)
                              for _, g in df.groupby("t")])["i"].values)
gi = C.gene_index(a, args.geneset)
if args.drop_human_tech:
    bad = set(TECH_HUMAN)
    print('technical genes dropped:', sorted(g for g in a.var_names[gi] if g in bad), flush=True)
    gi = np.array([i for i in gi if a.var_names[i] not in bad])
genes = a.var_names[gi].tolist()
Xd = a.X[rows][:, gi].toarray().astype(np.float64)
pca = C.fit_pca(Xd, n_pcs=10)
sc = pca["scores_fit"]
k = args.k
fit = C.fit_simplex(sc, k)
obs_t = fit["t_ratio"]
print(f"{args.line}/{args.geneset}/{args.variant} k={k}: {len(rows)} cells, t={obs_t:.4f}", flush=True)

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
np.save(out / f"null_k{k}.npy", nl); np.save(out / f"archetypes_k{k}.npy", fit["archetypes"])
np.savez_compressed(out / "pca.npz", rows=rows, scores=sc, mean=pca["mean"], components=pca["components"],
                    evr=pca["evr"][:50], genes=np.array(genes))
res = dict(line=args.line, geneset=args.geneset, variant=args.variant, k=k, n_cells=int(len(rows)),
           t_ratio=obs_t, null_median=float(np.median(fin)), null_95=float(np.quantile(fin, 0.95)),
           n_null=int(len(fin)), p_value=float(p))
(out / f"fit_k{k}.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res), f"({time.time() - t0:.0f}s)")
