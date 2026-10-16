#!/usr/bin/env python3
"""Step 9 (follow-up 3): project every tumour into one tumour's frozen triangle.

Question: is a corner found inside one tumour (the host) shared by the other
tumours? If late-EMT cells of other tumours land at the host's late corner, the
corner is shared biology; if they land elsewhere or outside, it is not.

Host: the within-tumour k=3 fit of --host (default SC32, the largest mixed tumour)
from 08_within_tumour.py. Guests: all other cells, V1-corrected the same way
(step 1), centred on the HOST's gene means and projected on the host's PCs.
No per-tumour rescaling, so each tumour's fingerprint is kept (on purpose).

Outputs (results/cross_projection_<host>/):
  cells.csv.gz         host-space PC1/PC2, weights, corner/interior, inside flag
  summary.md           per tumour x state: % at each host corner, interior, inside
  guests_in_host.png   one panel per tumour, cells coloured by EMT state

Usage (from Project_1/):
  .venv/bin/python -u Mouse_SCC/codes/09_cross_projection.py --host SC32
"""
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scc_common as C

NOV1 = C.RESULTS / "all_scrna_noV1"
WT = C.RESULTS / "within_tumour"
SCOL = {"Epithelial": "#1f77b4", "Early-Hybrid": "#ff7f0e", "Hybrid": "#2ca02c", "Late-EMT": "#d62728"}
ORDER = ["SC32", "SC15", "SC46", "SC33", "SC2b1", "SC2b3", "SC35", "SC43", "MU3", "MU4", "MU7", "MU8"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="SC32")
    ap.add_argument("--raw", action="store_true", help="no V1 correction; uses within_tumour_raw/")
    args = ap.parse_args()
    wt = C.RESULTS / "within_tumour_raw" if args.raw else WT
    out = C.RESULTS / (f"cross_projection_{args.host}" + ("_raw" if args.raw else ""))
    out.mkdir(exist_ok=True)

    a = C.load_adata()
    if args.raw:
        beta, score, s_mean = np.zeros(a.n_vars), np.zeros(a.n_obs), 0.0
    else:
        z = np.load(NOV1 / "pca.npz")
        assert list(a.obs_names) == list(z["barcodes"])
        beta, score = z["beta"], z["v1_score"]
        s_mean = float(json.loads((NOV1 / "pca_info.json").read_text())["score_mean_fit"])
    h = np.load(wt / args.host / "pca.npz")
    arcs = np.load(wt / args.host / "archetypes_k3.npy")
    V = h["components"][:2]

    # (x_res - host_mean) @ V'  with  x_res = x - beta*(score - s_mean)
    pts = np.asarray(a.X @ V.T) - np.outer(score - s_mean, beta @ V.T) - h["mean"] @ V.T
    hr = h["rows"]
    assert np.abs(pts[hr] - h["scores"][:, :2]).max() < 1e-6 * np.abs(h["scores"]).max(), "host mismatch"

    w = C.barycentric(pts, arcs)
    pur, near = C.purity_and_nearest(w)
    Vn = ["V1", "V2", "V3"]
    comp = pd.read_csv(wt / "corner_composition.csv")
    names = {}
    for j, v in enumerate(Vn):
        r = comp[comp["corner"] == f"{args.host}_{v}"].iloc[0]
        top = max(C.STATE_ORDER, key=lambda s: r[s])
        names[v] = f"{v} ({top} {r[top]:.0f}%)"
    cells = pd.DataFrame({"PC1": pts[:, 0], "PC2": pts[:, 1],
                          "w1": w[:, 0], "w2": w[:, 1], "w3": w[:, 2],
                          "inside": (w >= -1e-6).all(1),
                          "group": np.where(pur >= C.PURITY_CUT, np.array(Vn)[near], "interior"),
                          "sample": a.obs["sample"].astype(str).values,
                          "assay": a.obs["assay"].astype(str).values,
                          "state": a.obs["cluster_annotation"].astype(str).values},
                         index=a.obs_names)
    cells["is_host"] = cells["sample"] == args.host
    cells.to_csv(out / "cells.csv.gz")

    lines = [f"# Cross-projection into the {args.host} triangle (k=3, {'uncorrected' if args.raw else 'V1-corrected'})", "",
             "Host corners (from 08_within_tumour.py --compare): " + "; ".join(names.values()), "",
             "Rows: tumour x state. Columns: % of those cells at each host corner (purity gap >= 0.35),",
             "% interior, and % inside the triangle (all weights >= 0).", ""]
    rows = []
    for t in ORDER:
        for s in C.STATE_ORDER:
            m = (cells["sample"] == t) & (cells["state"] == s)
            if m.sum() < 20:
                continue
            g = cells.loc[m, "group"].value_counts(normalize=True)
            rows.append(dict(tumour=t + (" (host)" if t == args.host else ""), state=s, n=int(m.sum()),
                             **{names[v]: round(100 * g.get(v, 0), 1) for v in Vn},
                             interior=round(100 * g.get("interior", 0), 1),
                             inside=round(100 * cells.loc[m, "inside"].mean(), 1)))
    tab = pd.DataFrame(rows)
    tab.to_csv(out / "tumour_state_in_host.csv", index=False)
    pd.set_option("display.width", 250)
    lines.append(tab.to_string(index=False))
    # pooled over guests by state
    G = cells[~cells["is_host"] & (cells["assay"] == "scRNA")]
    pool = (100 * pd.crosstab(G["state"], G["group"], normalize="index")).round(1).reindex(C.STATE_ORDER)
    pool["inside"] = (100 * G.groupby("state")["inside"].mean()).round(1)
    lines += ["", "## All scRNA guest cells pooled, by state (%)", "", pool.to_string()]
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    fig, axes = plt.subplots(3, 4, figsize=(18, 13), sharex=True, sharey=True)
    tri = np.vstack([arcs[:, :2], arcs[:1, :2]])
    for ax, t in zip(axes.ravel(), ORDER):
        m = (cells["sample"] == t).values
        sub = cells[m].sample(frac=1, random_state=0)
        ax.scatter(cells["PC1"], cells["PC2"], s=0.5, c="#e6e6e6", rasterized=True)
        ax.scatter(sub["PC1"], sub["PC2"], s=2, c=[SCOL[x] for x in sub["state"]], rasterized=True)
        ax.plot(tri[:, 0], tri[:, 1], "k-", lw=1)
        for j, v in enumerate(Vn):
            ax.annotate(names[v], arcs[j, :2], fontsize=7, fontweight="bold")
        ax.set_title(f"{t}{' (host)' if t == args.host else ''}  n={m.sum()}  "
                     f"inside {100 * sub['inside'].mean():.0f}%", fontsize=10)
    for s, c in SCOL.items():
        axes[0, 0].scatter([], [], c=c, s=20, label=s)
    axes[0, 0].legend(fontsize=7, loc="lower left")
    fig.suptitle(f"Every tumour projected into the {args.host} within-tumour triangle "
                 f"({'uncorrected' if args.raw else 'V1-corrected'} expression; host PC1 vs PC2; grey = all cells)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out / "guests_in_host.png", dpi=150)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
