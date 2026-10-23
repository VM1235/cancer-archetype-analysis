#!/usr/bin/env python3
"""Step 2: PCA on the library-balanced fit cells, then diagnostics.

What it does
  1. Picks the fit cells: at most 800 random cells from each of the 12 tumours,
     so no single big tumour decides the axes.
  2. Exact PCA (full SVD) on those cells. Genes centred, not scaled.
  3. Projects all all cells onto the same PCs.
  4. Diagnostics: for PCs 1-10, how much of each PC is explained by tumour,
     assay, EMT state, finer cluster, cell-cycle phase and sequencing depth.
     If tumour or assay explain more than state, the axes are batch-driven.
  5. PC-count rule (Groves 2022, single cell): keep PCs until the next one adds
     less than 0.1% of the variance. Reported, used for the ESV curve.
  6. Checks a fast (randomized) PCA against the exact one on PCs 1-10.

Usage (from Project_1/):
  .venv/bin/python -u "Breast Cancer/GSE114687_pseudospace"/codes/02_pca.py --geneset all
  .venv/bin/python -u "Breast Cancer/GSE114687_pseudospace"/codes/02_pca.py --geneset ks
Optional: --exclude-tumours SC15 SC46 (fit cells from the other tumours only; use with --tag)
          --remove-phase  (subtracts per-phase gene means before PCA; saved as tag 'nophase')
          --assay in_vivo     (fit cells from in_vivo tumours only; saved as tag 'in_vivo').
                            Multiome cells are still projected, as guests.
"""
import argparse
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scc_common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geneset", choices=["all", "ks", "hallmark"], default="all")
    ap.add_argument("--remove-phase", action="store_true")
    ap.add_argument("--pca-method", choices=["gram", "svd"], default="gram")
    ap.add_argument("--assay", choices=["in_vivo", "in_vitro"], default=None)
    ap.add_argument("--exclude-tumours", nargs="*", default=[],
                    help="tumours left out of the fit cells (still projected afterwards)")
    ap.add_argument("--tag", default=None, help="output folder suffix (default main/nophase/in_vivo)")
    args = ap.parse_args()
    tag = args.tag or ("nophase" if args.remove_phase else
                       args.assay.lower() if args.assay else "main")
    out = C.run_dir(args.geneset, tag)
    t0 = time.time()

    a = C.load_adata()
    gi = C.gene_index(a, args.geneset)
    genes = a.var_names[gi].tolist()
    obs = a.obs.copy()
    X = a.X[:, gi].tocsr()
    del a
    fit_idx = C.balanced_fit_cells(obs, only_assay=args.assay, exclude_tumours=args.exclude_tumours)
    print(f"gene set '{args.geneset}': {len(genes)} genes; fit cells {len(fit_idx)} "
          f"of {X.shape[0]}", flush=True)

    Xfit = X[fit_idx].toarray().astype(np.float64)
    shift_all = None
    if args.remove_phase:
        pm = C.phase_means(Xfit, obs["phase"].values[fit_idx])
        shift_fit = np.vstack([pm[p] for p in obs["phase"].values[fit_idx]])
        Xfit = Xfit - shift_fit
        shift_all = np.vstack([pm[p] for p in obs["phase"].values])
        print("removed per-phase gene means (cell-cycle regression)", flush=True)

    pca = C.fit_pca(Xfit, method=args.pca_method)
    print(f"exact PCA ({args.pca_method}) done ({time.time() - t0:.0f}s)", flush=True)
    scores_all = C.project(X, pca, shift=shift_all)
    # sanity: projecting the fit cells must reproduce the SVD scores
    err = np.abs(scores_all[fit_idx] - pca["scores_fit"]).max()
    assert err < 1e-6 * max(1.0, np.abs(pca["scores_fit"]).max()), f"projection mismatch {err}"

    # ---- PC-count rule
    evr = pca["evr"]
    below = np.where(evr < 0.001)[0]
    n_rule = int(below[0]) if len(below) else len(evr)
    cum = np.cumsum(evr)
    print(f"PC1..PC10 variance %: {np.round(100 * evr[:10], 2).tolist()}")
    print(f"0.1% rule keeps {n_rule} PCs, which together explain {100 * cum[n_rule - 1]:.1f}%")

    # ---- diagnostics on the fit cells
    o = obs.iloc[fit_idx]
    logdepth = np.log10(o["nCount_RNA"].values)
    rows = []
    for j in range(10):
        s = pca["scores_fit"][:, j]
        rows.append(dict(
            PC=j + 1, pct_variance=100 * evr[j],
            tumour=C.eta_squared(s, o["sample"]), assay=C.eta_squared(s, o["assay"]),
            state=C.eta_squared(s, o["cluster_annotation"]),
            clusterU=C.eta_squared(s, o["clusterU"]), phase=C.eta_squared(s, o["phase"]),
            position=C.eta_squared(s, o["position"]), arm=C.eta_squared(s, o["arm"]),
            log_depth=float(np.corrcoef(s, logdepth)[0, 1] ** 2),
            pct_mito=float(np.corrcoef(s, o["pct_counts_mt"].values)[0, 1] ** 2),
        ))
    diag = pd.DataFrame(rows)
    diag.to_csv(out / "pc_drivers.csv", index=False)
    print("\nShare of each PC explained by each factor (0-1; fit cells):")
    print(diag.round(3).to_string(index=False))

    # ---- fast vs exact PCA on the first 10 PCs
    from sklearn.decomposition import PCA
    rp = PCA(n_components=10, svd_solver="randomized", iterated_power=7, random_state=0)
    rs = rp.fit_transform(Xfit)
    agree = [abs(np.corrcoef(rs[:, j], pca["scores_fit"][:, j])[0, 1]) for j in range(10)]
    print("\nrandomized vs exact |corr| PC1..10:", np.round(agree, 4).tolist())

    # ---- top genes per PC
    top = {}
    for j in range(5):
        v = pca["components"][j]
        top[f"PC{j + 1}+"] = [genes[i] for i in np.argsort(v)[::-1][:15]]
        top[f"PC{j + 1}-"] = [genes[i] for i in np.argsort(v)[:15]]
    pd.DataFrame(top).to_csv(out / "pc_top_genes.csv", index=False)

    # ---- save
    np.savez_compressed(out / "pca.npz", scores_all=scores_all, fit_idx=fit_idx,
                        mean=pca["mean"], components=pca["components"], evr=evr[:200],
                        genes=np.array(genes), barcodes=np.array(obs.index, dtype=str))
    C.save_json(dict(geneset=args.geneset, tag=tag, n_genes=len(genes),
                     n_fit_cells=int(len(fit_idx)), n_cells=int(X.shape[0]), fit_assay=args.assay, excluded_tumours=args.exclude_tumours,
                     max_per_tumour=C.MAX_PER_TUMOUR, seed=C.SEED, n_pcs_rule=n_rule,
                     cumvar_at_rule=float(cum[n_rule - 1]),
                     randomized_vs_exact_abs_corr=agree, remove_phase=args.remove_phase),
                out / "pca_info.json")

    # ---- figures
    fig, ax = plt.subplots(1, 1, figsize=(5, 3.4))
    ax.plot(np.arange(1, 51), 100 * evr[:50], "o-", ms=3)
    ax.axhline(0.1, ls="--", c="0.5", lw=0.8)
    ax.axvline(n_rule + 0.5, ls=":", c="0.5", lw=0.8)
    ax.set_xlabel("PC"); ax.set_ylabel("% variance"); ax.set_yscale("log")
    ax.set_title(f"Scree ({args.geneset}, {tag}); 0.1% rule -> {n_rule} PCs", fontsize=9)
    fig.savefig(out / "scree.png", dpi=200, bbox_inches="tight"); plt.close(fig)

    s = pca["scores_fit"]
    fig, axes = plt.subplots(1, 5, figsize=(22, 4.2))
    for ax, col in zip(axes, ["cluster_annotation", "sample", "position", "assay", "phase"]):
        lab = o[col].astype(str).values
        for x in sorted(np.unique(lab)):
            m = lab == x
            ax.scatter(s[m, 0], s[m, 1], s=2, alpha=0.5, label=x, rasterized=True)
        ax.set_title(col); ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        ax.legend(markerscale=5, fontsize=6, loc="best")
    fig.savefig(out / "pc1_pc2_by_label.png", dpi=180, bbox_inches="tight"); plt.close(fig)
    print(f"\nwrote {out}  ({time.time() - t0:.0f}s total)")


if __name__ == "__main__":
    main()
