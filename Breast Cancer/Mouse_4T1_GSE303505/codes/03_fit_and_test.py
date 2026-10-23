#!/usr/bin/env python3
"""Step 3: fit simplices for k = 3..7 and test them with a permutation null.

What it does
  1. ESV curve: how much variance (in the first 10 PCs) a k-vertex simplex
     explains, k = 2..10. The "elbow" is where adding a vertex stops helping.
     10 PCs is a fixed working space; it holds every simplex up to k = 11.
  2. Real fit for each k: PCHA on the first k-1 PCs of the fit cells,
     15 starts, keep the largest-volume simplex (repo engine, same as Pal).
     t-ratio = simplex volume / volume of the data's convex hull.
  3. Permutation test: shuffle each of the first k-1 PCs independently
     (this keeps each PC's spread but destroys the shape), refit with 5 starts,
     and record the t-ratio. p = (1 + #null >= observed) / (1 + #shuffles).
     A real simplex has a t-ratio closer to 1 (tighter fit) than shuffled data.
     Saved after every shuffle, so a stopped run resumes where it left off.

Each k can run as its own process (they are independent), then --summarize
collects the per-k results into one table and figure. run_fits.sh does this.

Usage (from Project_1/):
  .venv/bin/python -u "Breast Cancer/Mouse_4T1_GSE303505"/codes/03_fit_and_test.py --geneset all --n-perm 3 --pilot
  .venv/bin/python -u "Breast Cancer/Mouse_4T1_GSE303505"/codes/03_fit_and_test.py --geneset all --esv
  .venv/bin/python -u "Breast Cancer/Mouse_4T1_GSE303505"/codes/03_fit_and_test.py --geneset all --k 3 --n-perm 200
  .venv/bin/python -u "Breast Cancer/Mouse_4T1_GSE303505"/codes/03_fit_and_test.py --geneset all --summarize
"""
import json
import argparse
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, QhullError

import scc_common as C
from src.archetypes import fit_pcha

ESV_DIMS = 10
ESV_K = range(2, 11)


def esv_curve(scores):
    rows = []
    X = scores[:, :ESV_DIMS]
    for k in ESV_K:
        np.random.seed(C.SEED)
        best = max(fit_pcha(X, k, delta=0.0)[2] for _ in range(5))
        rows.append(dict(k=k, esv=best))
    df = pd.DataFrame(rows)
    df["gain"] = df["esv"].diff()
    return df


def elbow(esv):
    """ParTI DimensionFinder: k farthest from the straight line joining the curve ends."""
    k = esv["k"].values.astype(float)
    y = esv["esv"].values
    p0, p1 = np.array([k[0], y[0]]), np.array([k[-1], y[-1]])
    d = np.abs(np.cross(p1 - p0, np.column_stack([k, y]) - p0)) / np.linalg.norm(p1 - p0)
    return int(k[np.argmax(d)])


def t_ratio(scores, arcs):
    k = arcs.shape[0]
    return float(C.simplex_volume(arcs) / ConvexHull(scores[:, : k - 1]).volume)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geneset", choices=["all", "ks", "hallmark"], default="all")
    ap.add_argument("--tag", default="main")
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--pilot", action="store_true", help="time a few shuffles; write nothing final")
    ap.add_argument("--k", type=int, nargs="*", default=list(C.K_VALUES))
    ap.add_argument("--esv", action="store_true", help="only compute the ESV curve")
    ap.add_argument("--summarize", action="store_true", help="only collect per-k results")
    args = ap.parse_args()
    out = C.run_dir(args.geneset, args.tag)
    z = np.load(out / "pca.npz")
    scores = z["scores_all"][z["fit_idx"]]
    print(f"{args.geneset}/{args.tag}: {scores.shape[0]} fit cells", flush=True)

    if args.esv:
        t0 = time.time()
        esv = esv_curve(scores)
        esv.to_csv(out / "esv_curve.csv", index=False)
        print(esv.round(4).to_string(index=False))
        print(f"ESV elbow at k={elbow(esv)} ({time.time() - t0:.0f}s)", flush=True)
        return
    if args.summarize:
        summarize(out, args)
        return

    rows = []
    for k in args.k:
        t0 = time.time()
        fit = C.fit_simplex(scores, k)
        np.save(out / f"archetypes_k{k}.npy", fit["archetypes"])
        obs_t = fit["t_ratio"]
        print(f"\nk={k}: t-ratio {obs_t:.4f}, {fit['n_ok']}/{C.N_INIT} starts ok "
              f"({time.time() - t0:.1f}s)", flush=True)

        null_path = out / f"null_k{k}.npy"
        null = [] if args.pilot or not null_path.exists() else list(np.load(null_path))
        rng = np.random.default_rng(C.SEED + 1000 + k)
        for _ in range(len(null)):                     # replay RNG to the resume point
            C._shuffle_columns(scores[:, : k - 1], rng)
        t1 = time.time()
        start = len(null)
        for i in range(start, args.n_perm):
            sh = C._shuffle_columns(scores[:, : k - 1], rng)
            try:
                np.random.seed(C.SEED + 10_000 * k + i)
                arcs, *_ = C.fit_pcha_best(sh, k, n_init=C.N_INIT_NULL, delta=0.0)
                null.append(t_ratio(sh, arcs))
            except (QhullError, ValueError, RuntimeError) as e:
                print(f"   shuffle {i + 1} failed: {e}", flush=True)
                null.append(np.nan)
            if not args.pilot:
                np.save(null_path, np.asarray(null))
            if (i + 1) % 25 == 0 or args.pilot:
                print(f"   {i + 1}/{args.n_perm} shuffles, null t max so far "
                      f"{np.nanmax(null):.4f}", flush=True)
        per = (time.time() - t1) / max(1, args.n_perm - start)
        nl = np.asarray(null, float)
        nl = nl[np.isfinite(nl)]
        p = (1 + np.sum(nl >= obs_t)) / (1 + len(nl))
        row = dict(k=k, t_ratio=obs_t, null_median=float(np.median(nl)),
                   null_95pct=float(np.quantile(nl, 0.95)), null_max=float(nl.max()),
                   n_null=len(nl), p_value=float(p), esv_in_k_minus_1_pcs=fit["esv"],
                   n_starts_ok=fit["n_ok"], sec_per_shuffle=per)
        rows.append(row)
        if not args.pilot:
            (out / f"fit_k{k}.json").write_text(json.dumps(row, indent=2))
        print(f"k={k}: p = {p:.4f} from {len(nl)} shuffles; {per:.1f}s per shuffle", flush=True)
    print("\n" + pd.DataFrame(rows).round(4).to_string(index=False))


def summarize(out, args):
    rows = [json.loads((out / f"fit_k{k}.json").read_text())
            for k in C.K_VALUES if (out / f"fit_k{k}.json").exists()]
    tab = pd.DataFrame(rows)
    tab.to_csv(out / "fit_summary.csv", index=False)
    esv = pd.read_csv(out / "esv_curve.csv")
    k_elbow = elbow(esv)
    print(tab.round(4).to_string(index=False))
    print(f"ESV elbow k={k_elbow}")
    sig = tab[tab["p_value"] < 0.01]["k"].tolist()
    print("k with p < 0.01:", sig if sig else "none")

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    axes[0].plot(esv["k"], esv["esv"], "o-")
    axes[0].axvline(k_elbow, ls="--", c="0.5")
    axes[0].set_xlabel("k (vertices)"); axes[0].set_ylabel(f"ESV in first {ESV_DIMS} PCs")
    axes[0].set_title(f"ESV curve, elbow k={k_elbow}")
    for k in tab["k"]:
        nl = np.load(out / f"null_k{k}.npy")
        axes[1].scatter(np.full(len(nl), k) + np.random.uniform(-.15, .15, len(nl)), nl,
                        s=4, c="0.6")
    axes[1].plot(tab["k"], tab["t_ratio"], "rD", label="observed")
    axes[1].set_xlabel("k"); axes[1].set_ylabel("t-ratio"); axes[1].legend()
    axes[1].set_title("Observed t-ratio vs shuffled null")
    fig.suptitle(f"{args.geneset}/{args.tag}", fontsize=9)
    fig.savefig(out / "fit_and_null.png", dpi=200, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
