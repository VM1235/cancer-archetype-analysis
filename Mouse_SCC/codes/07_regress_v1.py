#!/usr/bin/env python3
"""Step 7 (follow-up 1): remove the V1 "long-gene" score, then redo version B's PCA.

Why: in version B (scRNA fit cells, all genes) corner V1 tracks a sample-level
signal (all multiome samples, SC15, partly SC46) that looks technical. If we remove
that signal and refit, does V1 disappear, and what k does the data then support?

What it does
  1. Score per cell = mean log-expression of the top-20 V1 marker genes of the
     version B k=3 fit (results/all_scrna/describe_k3/vertex_marker_genes.csv, V1_up).
  2. On the version-B fit cells (same 6,400 scRNA cells, same seed), for every gene g:
        beta_g = cov(x_g, score) / var(score)
     and every cell (fit and projected) gets  x_g - beta_g * (score - mean_score).
     This removes the part of each gene that moves linearly with the score.
  3. Exact PCA on the residual fit cells; project all 23,959 cells.
  4. Saves pca.npz in the same format as 02_pca.py, so 03_fit_and_test.py and
     04_describe_simplex.py run unchanged with --tag scrna_noV1.
  5. Diagnostics: share of PCs 1-10 explained by state / tumour / assay / depth /
     the V1 score, before (version B) and after.

Usage (from Project_1/):
  .venv/bin/python -u Mouse_SCC/codes/07_regress_v1.py
"""
import time

import numpy as np
import pandas as pd

import scc_common as C

TAG_IN, TAG_OUT, N_MARKERS = "scrna", "scrna_noV1", 20


def drivers(scores, o, score, evr):
    rows = []
    state = o["cluster_annotation"].astype(str).values
    for j in range(10):
        s = scores[:, j]
        r = s - pd.Series(s).groupby(state).transform("mean").values
        rows.append(dict(PC=j + 1, pct_variance=100 * evr[j],
                         state=C.eta_squared(s, state), clusterU=C.eta_squared(s, o["clusterU"]),
                         tumour=C.eta_squared(s, o["sample"]),
                         tumour_in_state=C.eta_squared(r, o["sample"]),
                         phase=C.eta_squared(s, o["phase"]),
                         log_depth=float(np.corrcoef(s, np.log10(o["nCount_RNA"].values))[0, 1] ** 2),
                         v1_score=float(np.corrcoef(s, score)[0, 1] ** 2)))
    return pd.DataFrame(rows)


def main():
    t0 = time.time()
    out = C.run_dir("all", TAG_OUT)
    a = C.load_adata()
    obs = a.obs.copy()
    X = a.X.tocsr()
    genes = a.var_names.tolist()

    mk = pd.read_csv(C.RESULTS / "all_scrna" / "describe_k3" / "vertex_marker_genes.csv")["V1_up"]
    markers = mk.tolist()[:N_MARKERS]
    idx = [a.var_names.get_loc(g) for g in markers]
    score = np.asarray(X[:, idx].mean(1)).ravel()
    del a
    fit_idx = C.balanced_fit_cells(obs, only_assay="scRNA")
    ref_fit = np.load(C.RESULTS / "all_scrna" / "pca.npz")["fit_idx"]
    assert np.array_equal(fit_idx, ref_fit), "fit cells differ from version B"
    print(f"markers: {markers}\nfit cells {len(fit_idx)} (identical to version B)", flush=True)

    Xfit = X[fit_idx].toarray().astype(np.float64)
    s_fit = score[fit_idx]
    s_mean = s_fit.mean()
    sc = s_fit - s_mean
    beta = (sc @ (Xfit - Xfit.mean(0))) / (sc @ sc)              # one slope per gene
    Xfit -= np.outer(sc, beta)
    pca = C.fit_pca(Xfit)
    del Xfit
    # all cells: (x - beta*(s - s_mean)) projected  ==  x@V' - mean@V' - (s - s_mean) * (beta@V')
    V = pca["components"]
    scores_all = C.project(X, pca) - np.outer(score - s_mean, beta @ V.T)
    err = np.abs(scores_all[fit_idx] - pca["scores_fit"]).max()
    assert err < 1e-6 * max(1.0, np.abs(pca["scores_fit"]).max()), f"projection mismatch {err}"
    print(f"PCA done ({time.time() - t0:.0f}s)", flush=True)

    # diagnostics: before (version B) vs after
    o = obs.iloc[fit_idx]
    zb = np.load(C.RESULTS / "all_scrna" / "pca.npz")
    before = drivers(zb["scores_all"][fit_idx], o, s_fit, zb["evr"])
    after = drivers(scores_all[fit_idx], o, s_fit, pca["evr"])
    before.to_csv(out / "pc_drivers_before.csv", index=False)
    after.to_csv(out / "pc_drivers.csv", index=False)
    pd.set_option("display.width", 200)
    print("\nBEFORE (version B):\n" + before.round(3).to_string(index=False))
    print("\nAFTER (V1 score regressed out):\n" + after.round(3).to_string(index=False))

    top = {}
    for j in range(4):
        v = V[j]
        top[f"PC{j + 1}+"] = [genes[i] for i in np.argsort(v)[::-1][:15]]
        top[f"PC{j + 1}-"] = [genes[i] for i in np.argsort(v)[:15]]
    top = pd.DataFrame(top)
    top.to_csv(out / "pc_top_genes.csv", index=False)
    print("\nTop genes after regression:\n" + top.head(8).to_string(index=False))

    np.savez_compressed(out / "pca.npz", scores_all=scores_all, fit_idx=fit_idx, mean=pca["mean"],
                        components=V, evr=pca["evr"][:200], genes=np.array(genes),
                        barcodes=np.array(obs.index, dtype=str), beta=beta, v1_score=score)
    C.save_json(dict(geneset="all", tag=TAG_OUT, source_tag=TAG_IN, markers=markers,
                     n_fit_cells=int(len(fit_idx)), score_mean_fit=float(s_mean),
                     note="V1 marker score regressed out of every gene (linear, fit-cell slopes)"),
                out / "pca_info.json")
    print(f"\nwrote {out} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
