#!/usr/bin/env python3
"""Step 8 (follow-up 2): fit a simplex inside single tumours, compare corners by genes.

Question: does the epithelial / hybrid / late-EMT triangle from step 1 (V1 score
removed) exist INSIDE individual tumours, or does it only appear when different
tumours are pooled?

Data: V1-corrected expression (step 1, 07_regress_v1.py): for every cell and gene
  x_res = x - beta_g * (score - score_mean),  beta and score from all_scrna_noV1/pca.npz.
Means of x_res over any set of cells are computed without densifying:
  mean(x_res[set]) = mean(x[set]) - (mean(score[set]) - score_mean) * beta

Modes
  --fit --tumour SC32 --k 3 --n-perm 200
      PCA on that tumour's cells only (exact, centred), PCHA (15 starts, delta 0)
      on the first k-1 PCs, permutation test (shuffle each PC, 5 starts), saved to
      results/within_tumour/<tumour>/.
  --compare
      For each tumour fit (k=3) and the global step-1 fit (k=3, scRNA cells):
      corner signature = mean x_res in the 10% of cells closest to the corner minus
      mean x_res in the other cells (all 14,502 genes). Pearson correlation between
      all corner signatures; state composition of each corner's specialists.

--raw: uncorrected expression (a per-tumour V1 shift is removed by centring anyway);
       outputs in results/within_tumour_raw/, pooled reference = the six-tumour fit (all_scrna6).

Usage (from Project_1/):
  .venv/bin/python -u Mouse_SCC/codes/08_within_tumour.py --fit --tumour SC32 --k 3
  .venv/bin/python -u Mouse_SCC/codes/08_within_tumour.py --compare
"""
import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, QhullError

import scc_common as C

TUMOURS = ["SC15", "SC32", "SC46"]
OUT = C.RESULTS / "within_tumour"
NOV1 = C.RESULTS / "all_scrna_noV1"
GLOBAL = NOV1          # pooled reference fit whose corners the tumour corners are compared with
RAW = False            # --raw: no V1 correction; reference = six-tumour fit (all_scrna6)


def load():
    a = C.load_adata()
    if RAW:
        return a, np.zeros(a.n_vars), np.zeros(a.n_obs), 0.0
    z = np.load(NOV1 / "pca.npz")
    assert list(a.obs_names) == list(z["barcodes"])
    s_mean = float(json.loads((NOV1 / "pca_info.json").read_text())["score_mean_fit"])
    return a, z["beta"], z["v1_score"], s_mean


def resid_dense(X, rows, beta, score, s_mean):
    return X[rows].toarray().astype(np.float64) - np.outer(score[rows] - s_mean, beta)


def resid_mean(X, rows, beta, score, s_mean):
    return np.asarray(X[rows].mean(0)).ravel() - (score[rows].mean() - s_mean) * beta


def signatures(X, rows, pts, arcs, beta, score, s_mean, frac=0.10):
    """Corner signature per corner: mean residual expression, closest 10% minus the rest."""
    sig = []
    d = np.linalg.norm(pts[:, None, :] - arcs[None, :, :], axis=2)
    for j in range(arcs.shape[0]):
        close = d[:, j] <= np.quantile(d[:, j], frac)
        sig.append(resid_mean(X, rows[close], beta, score, s_mean)
                   - resid_mean(X, rows[~close], beta, score, s_mean))
    return np.vstack(sig)


def t_ratio(sc, arcs):
    k = arcs.shape[0]
    return float(C.simplex_volume(arcs) / ConvexHull(sc[:, : k - 1]).volume)


def fit_mode(args):
    t0 = time.time()
    a, beta, score, s_mean = load()
    rows = np.where(a.obs["sample"].astype(str).values == args.tumour)[0]
    out = OUT / args.tumour
    out.mkdir(parents=True, exist_ok=True)
    Xd = resid_dense(a.X, rows, beta, score, s_mean)
    pca = C.fit_pca(Xd, n_pcs=10)
    sc = pca["scores_fit"]
    fit = C.fit_simplex(sc, args.k)
    obs_t = fit["t_ratio"]
    print(f"{args.tumour} k={args.k}: {len(rows)} cells, t={obs_t:.4f}", flush=True)
    rng = np.random.default_rng(C.SEED + 1000 + args.k)
    null = []
    for i in range(args.n_perm):
        sh = C._shuffle_columns(sc[:, : args.k - 1], rng)
        try:
            np.random.seed(C.SEED + 10_000 * args.k + i)
            arcs, *_ = C.fit_pcha_best(sh, args.k, n_init=C.N_INIT_NULL, delta=0.0)
            null.append(t_ratio(sh, arcs))
        except (QhullError, ValueError, RuntimeError):
            null.append(np.nan)
    nl = np.asarray(null)
    fin = nl[np.isfinite(nl)]
    p = (1 + np.sum(fin >= obs_t)) / (1 + len(fin))
    np.save(out / f"null_k{args.k}.npy", nl)
    np.save(out / f"archetypes_k{args.k}.npy", fit["archetypes"])
    np.savez_compressed(out / "pca.npz", rows=rows, scores=sc, mean=pca["mean"],
                        components=pca["components"], evr=pca["evr"][:50])
    res = dict(tumour=args.tumour, k=args.k, n_cells=int(len(rows)), t_ratio=obs_t,
               null_median=float(np.median(fin)), null_95=float(np.quantile(fin, 0.95)),
               null_max=float(fin.max()), n_null=int(len(fin)), p_value=float(p),
               evr_pc1_pc3=[float(x) for x in pca["evr"][:3]])
    (out / f"fit_k{args.k}.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res), f"({time.time() - t0:.0f}s)")


def compare_mode():
    a, beta, score, s_mean = load()
    obs = a.obs
    X = a.X
    state = obs["cluster_annotation"].astype(str).values
    tum = obs["sample"].astype(str).values
    names, sigs, comp_rows = [], [], []

    # global step-1 corners (k=3), scRNA cells only
    zg = np.load(GLOBAL / "pca.npz")
    ag = np.load(GLOBAL / "archetypes_k3.npy")
    ref = np.where(obs["assay"].values == "scRNA")[0]
    gpts = zg["scores_all"][ref, :2]
    gsig = signatures(X, ref, gpts, ag, beta, score, s_mean)
    # name each pooled corner by the top EMT state among its scRNA specialists
    gw = C.barycentric(gpts, ag)
    gpur, gnear = C.purity_and_nearest(gw)
    short = {"Epithelial": "epi", "Early-Hybrid": "early", "Hybrid": "hybrid", "Late-EMT": "late"}
    glab = {}
    for j in range(3):
        st = pd.Series(state[ref][(gpur >= C.PURITY_CUT) & (gnear == j)]).value_counts()
        nm = "G_" + short[st.index[0]]
        glab[j] = nm if nm not in glab.values() else f"{nm}{j + 1}"
    print("pooled reference corners:", glab)
    for j in range(3):
        names.append(glab[j]); sigs.append(gsig[j])

    for t in TUMOURS:
        d = OUT / t
        z = np.load(d / "pca.npz")
        rows, sc = z["rows"], z["scores"]
        arcs = np.load(d / "archetypes_k3.npy")
        w = C.barycentric(sc[:, :2], arcs)
        pur, near = C.purity_and_nearest(w)
        sig = signatures(X, rows, sc[:, :2], arcs, beta, score, s_mean)
        for j in range(3):
            m = (pur >= C.PURITY_CUT) & (near == j)
            st = pd.Series(state[rows][m]).value_counts(normalize=True)
            ph = obs["phase"].astype(str).values[rows][m]
            names.append(f"{t}_V{j + 1}"); sigs.append(sig[j])
            comp_rows.append(dict(corner=f"{t}_V{j + 1}", n_specialists=int(m.sum()),
                                  **{s: round(100 * st.get(s, 0), 1) for s in C.STATE_ORDER},
                                  pct_cycling=round(100 * np.isin(ph, ["S", "G2M"]).mean(), 1),
                                  top_genes=", ".join(a.var_names[np.argsort(sig[j])[::-1][:8]])))
        res = json.loads((d / "fit_k3.json").read_text())
        comp_rows.append(dict(corner=f"{t}: t={res['t_ratio']:.3f}, p={res['p_value']:.3f}",
                              n_specialists=None))
    S = np.vstack(sigs)
    R = pd.DataFrame(np.corrcoef(S), index=names, columns=names).round(2)
    R.to_csv(OUT / "corner_signature_correlation.csv")
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(OUT / "corner_composition.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 80)
    print("Corner composition (specialists = purity gap >= 0.35):\n" + comp.to_string(index=False))
    print("\nCorrelation of corner gene signatures (global G_* vs within-tumour corners):")
    print(R.loc[[n for n in names if not n.startswith("G_")], list(glab.values())].to_string())
    print("\nWithin-tumour corner vs corner (all):\n" + R.to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", action="store_true")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--tumour")
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--raw", action="store_true", help="no V1 correction; outputs in within_tumour_raw/")
    a = ap.parse_args()
    if a.raw:
        RAW = True
        OUT = C.RESULTS / "within_tumour_raw"
        GLOBAL = C.RESULTS / "all_scrna6"
    if a.fit:
        fit_mode(a)
    if a.compare:
        compare_mode()
