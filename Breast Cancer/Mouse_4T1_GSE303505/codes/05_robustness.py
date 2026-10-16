#!/usr/bin/env python3
"""Step 5: do the vertices survive when the fit cells change?

For each variant, the whole fit is redone from scratch (new PCA + new PCHA),
all all cells are placed in the new simplex, and the result is compared
with the main fit cell by cell.

Variants
  loto     leave one library out, for each library
  lomo     leave one mouse out (all libraries of M1, then of M2)
  in_vivo    fit on in_vivo cells only (8 tumours); removes the in_vitro assay
  nophase  per-phase gene means removed before PCA (cell cycle regressed out)
  seed     a different random library-balanced subsample (seeds 1 and 2)

Comparison (vertices from two fits live in different PC spaces, so they are
compared through the cells, not coordinates)
  agreement_pct   % of cells whose nearest vertex matches, after the best
                  one-to-one matching of new vertices to main vertices
  per vertex      the matched vertex's top EMT state and clusterU among its
                  specialists, and the share of its largest tumour
  heldout_*       (loto) where the left-out tumour's cells land

Results are appended after every variant run, so a stopped run resumes.

Usage (from Project_1/):
  .venv/bin/python -u "Breast Cancer/Mouse_4T1_GSE303505"/codes/05_robustness.py --geneset all --k 3 --variant loto
  .venv/bin/python -u "Breast Cancer/Mouse_4T1_GSE303505"/codes/05_robustness.py --geneset all --tag in_vivo --k 3 --variant loto
"""
import argparse
import time

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

import scc_common as C


def place(scores, arcs, k):
    w = C.barycentric(scores[:, : k - 1], arcs)
    pur, near = C.purity_and_nearest(w)
    return near, pur >= C.PURITY_CUT


def match(near_main, near_new, k):
    tab = np.zeros((k, k))
    np.add.at(tab, (near_main, near_new), 1)
    r, c = linear_sum_assignment(-tab)
    new_to_main = dict(zip(c, r))
    return np.array([new_to_main[x] for x in near_new]), tab[r, c].sum() / tab.sum()


def vertex_profile(obs, near, spec, k, prefix):
    out = {}
    for v in range(k):
        m = spec & (near == v)
        o = obs[m]
        out[f"{prefix}V{v + 1}_n_spec"] = int(m.sum())
        if m.sum() == 0:
            continue
        st = o["cluster_annotation"].value_counts(normalize=True)
        cu = o["clusterU"].value_counts(normalize=True)
        tu = o["sample"].value_counts(normalize=True)
        out[f"{prefix}V{v + 1}_state"] = f"{st.index[0]} {100 * st.iloc[0]:.0f}%"
        out[f"{prefix}V{v + 1}_clusterU"] = f"{cu.index[0]} {100 * cu.iloc[0]:.0f}%"
        out[f"{prefix}V{v + 1}_top_tumour"] = f"{tu.index[0]} {100 * tu.iloc[0]:.0f}%"
    return out


def refit(X, obs, fit_idx, k, genes_idx, remove_phase=False):
    Xfit = X[fit_idx].toarray().astype(np.float64)
    shift = None
    if remove_phase:
        pm = C.phase_means(Xfit, obs["phase"].values[fit_idx])
        Xfit -= np.vstack([pm[p] for p in obs["phase"].values[fit_idx]])
        shift = np.vstack([pm[p] for p in obs["phase"].values])
    pca = C.fit_pca(Xfit, n_pcs=max(k - 1, 10))
    del Xfit
    scores = C.project(X, pca, shift=shift)
    fit = C.fit_simplex(scores[fit_idx], k)
    return scores, fit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geneset", choices=["all", "ks"], default="all")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--variant", choices=["loto", "nophase", "seed"], required=True)
    ap.add_argument("--tag", default="main", help="which main fit to test (main or in_vivo)")
    args = ap.parse_args()
    k = args.k
    run = C.run_dir(args.geneset, args.tag)
    # the in_vivo-only fit keeps its own rule: fit cells come from in_vivo tumours only
    base_sel = dict(only_assay="in_vivo") if args.tag == "in_vivo" else {}
    if args.tag == "in_vivo" and args.variant == "in_vivo":
        raise SystemExit("variant 'in_vivo' is the in_vivo fit itself; nothing to test")
    out = run / "robustness"
    out.mkdir(exist_ok=True)
    res_path = out / f"k{k}_{args.variant}.csv"
    done = set(pd.read_csv(res_path)["variant"]) if res_path.exists() else set()

    z = np.load(run / "pca.npz")
    arcs_main = np.load(run / f"archetypes_k{k}.npy")
    near_main, spec_main = place(z["scores_all"], arcs_main, k)

    a = C.load_adata()
    assert list(a.obs_names) == list(z["barcodes"])
    gi = C.gene_index(a, args.geneset)
    X = a.X[:, gi].tocsr()
    obs = a.obs.copy()
    for c in ["sample", "assay", "cluster_annotation", "clusterU", "phase"]:
        obs[c] = obs[c].astype(str)
    del a

    if args.variant == "loto":
        tumours = sorted(obs.loc[obs["assay"] == "in_vivo", "sample"].unique()
                         if args.tag == "in_vivo" else obs["sample"].unique())
        jobs = [(f"without_{t}", dict(exclude_tumours=[t], **base_sel), False, t)
                for t in tumours]
    elif args.variant == "lomo":      # leave one mouse out: M1 libraries out, then M2 libraries out
        libs = obs.loc[obs["assay"] == "in_vivo", "sample"].unique() if args.tag == "in_vivo" else obs["sample"].unique()
        jobs = [(f"without_{m}", dict(exclude_tumours=[l for l in libs if l.startswith(m + "-")], **base_sel), False, m)
                for m in ("M1", "M2")]
    elif args.variant == "in_vivo":
        jobs = [("in_vivo_only", dict(only_assay="in_vivo"), False, None)]
    elif args.variant == "nophase":
        jobs = [("phase_regressed", dict(base_sel), True, None)]
    else:
        jobs = [(f"seed{s}", dict(seed=s, **base_sel), False, None) for s in (1, 2)]

    base = dict(geneset=args.geneset, k=k)
    base.update(vertex_profile(obs, near_main, spec_main, k, "main_"))
    for name, sel, rmphase, heldout in jobs:
        if name in done:
            print(f"skip {name} (done)")
            continue
        t0 = time.time()
        fit_idx = C.balanced_fit_cells(obs, **sel)
        scores, fit = refit(X, obs, fit_idx, k, gi, remove_phase=rmphase)
        near_new, spec_new = place(scores, fit["archetypes"], k)
        near_new, agree = match(near_main, near_new, k)
        both_spec = spec_main & spec_new
        agree_spec = float((near_main[both_spec] == near_new[both_spec]).mean())
        row = dict(variant=name, n_fit_cells=len(fit_idx), t_ratio=fit["t_ratio"],
                   agreement_pct=100 * agree, agreement_specialists_pct=100 * agree_spec,
                   **vertex_profile(obs, near_new, spec_new, k, ""))
        if heldout is not None:
            m = obs["sample"].str.startswith(heldout + "-").values if args.variant == "lomo" else (obs["sample"] == heldout).values
            row["heldout_n"] = int(m.sum())
            row["heldout_agreement_pct"] = 100 * float((near_main[m] == near_new[m]).mean())
            g = np.where(spec_new[m], near_new[m], -1)
            for v in range(k):
                row[f"heldout_V{v + 1}_pct"] = 100 * float((g == v).mean())
            row["heldout_interior_pct"] = 100 * float((g == -1).mean())
        row.update(base)
        df = pd.DataFrame([row])
        df.to_csv(res_path, mode="a", header=not res_path.exists(), index=False)
        print(f"{name}: agreement {100 * agree:.1f}% (specialists {100 * agree_spec:.1f}%), "
              f"t={fit['t_ratio']:.3f}  [{time.time() - t0:.0f}s]", flush=True)
        for v in range(k):
            print(f"   V{v + 1}: main {base.get(f'main_V{v + 1}_state')} | "
                  f"new {row.get(f'V{v + 1}_state')}, {row.get(f'V{v + 1}_top_tumour')}")
    print("wrote", res_path)


if __name__ == "__main__":
    main()
