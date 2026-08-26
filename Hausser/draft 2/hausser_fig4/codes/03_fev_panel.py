#!/usr/bin/env python3
"""Fig. 4d-style FEV: do METABRIC PCs explain single-cell variance?

Port of scBC.R FEV block (~786–874) / paper Methods:
  First 5 METABRIC PCs explain ~25.4% of the variance explained by the first
  5 single-cell PCs; shuffled METABRIC PCs explain ~1.5%.

Usage:
  .venv/bin/python -u "Hausser/draft 2/hausser_fig4/codes/03_fev_panel.py"
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy import stats

from _paths import ANALYSIS, METABRIC_ORIG, SC_COUNTS

RESULTS = ANALYSIS / "results" / "projection"
FIGURES = ANALYSIS / "figures"
LN2 = float(np.log(2))


def log_normalize(counts, scale=1e5):
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log1p(counts.div(lib, axis=1) * scale)


def fev(X, loadings_k):
    """Fraction of variance in X (cells × genes) explained by loadings (genes × k)."""
    tot = float(np.var(X))
    if tot <= 0:
        return np.nan
    proj = X @ loadings_k
    recon = proj @ loadings_k.T
    res = float(np.var(X - recon))
    return 1.0 - res / tot


def main():
    meta_path = RESULTS / "sc_cell_metadata.csv"
    if not meta_path.is_file():
        print("Run 01_project_sc_onto_metabric_shape.py first.")
        return 1

    meta = pd.read_csv(meta_path)
    cell_ids = meta["cell_id"].tolist()
    tumor_ids = meta["tumor_id"].tolist()

    # Reload METABRIC genes + expression
    gene_names = pd.read_csv(METABRIC_ORIG / "geneListExp.list", header=None)[0].astype(str)
    exp = pd.read_csv(METABRIC_ORIG / "expMatrix.csv", header=None)
    exp.columns = gene_names
    filt = pd.read_csv(METABRIC_ORIG / "geneNamesAfterExprFiltering.list", header=None)[0].astype(str)
    keep = [g for g in filt if g in exp.columns]
    exp = exp.loc[:, keep]

    # Reload sc cancer cells (same IDs as projection)
    sc = pd.read_csv(SC_COUNTS, sep="\t")
    sc = sc.rename(columns={sc.columns[0]: "gene"}).set_index("gene")
    sc.index = sc.index.astype(str)
    sc = sc.loc[~sc.index.duplicated(keep="first")]
    missing = [c for c in cell_ids if c not in sc.columns]
    if missing:
        print(f"Missing {len(missing)} cells in counts; abort.")
        return 1
    sc = sc.loc[:, cell_ids]
    log_sc = log_normalize(sc)

    frac = (log_sc > 0).mean(axis=1)
    sc_genes = frac[frac >= 0.5].index
    common = [g for g in keep if g in set(sc_genes)]
    print(f"FEV common genes: {len(common)}; cells: {len(cell_ids)}")

    # Centered log2 sc matrix: cells × genes
    sc_mat = log_sc.loc[common].values.astype(float)
    sc0 = (sc_mat - sc_mat.mean(axis=1, keepdims=True)) / LN2
    Xsc = sc0.T  # cells × genes

    # METABRIC PCA on common genes (center, no scale)
    Xmb = exp.loc[:, common].values.astype(float)
    Xmb0 = Xmb - Xmb.mean(axis=0)
    n_pcs = min(50, Xmb0.shape[0] - 1, Xmb0.shape[1])
    pca_mb = PCA(n_components=n_pcs, svd_solver="full", random_state=0).fit(Xmb0)
    Lmb = pca_mb.components_.T  # genes × pcs

    # Shuffled METABRIC: permute each gene across tumors
    rng = np.random.default_rng(0)
    Xrnd = Xmb0.copy()
    for j in range(Xrnd.shape[1]):
        Xrnd[:, j] = rng.permutation(Xrnd[:, j])
    pca_rnd = PCA(n_components=n_pcs, svd_solver="full", random_state=1).fit(Xrnd)
    Lrnd = pca_rnd.components_.T

    # Per-tumor FEV spectra (tumors with ≥50 cells, as in scBC.R)
    from collections import Counter

    counts = Counter(tumor_ids)
    big_tumors = [t for t, n in counts.items() if n >= 50]
    print("Tumors with ≥50 cells:", big_tumors)

    rows = []
    for t in big_tumors:
        idx = [i for i, tid in enumerate(tumor_ids) if tid == t]
        Xt = Xsc[idx, :]
        # self PCA
        n_self = min(100, Xt.shape[0] - 1, Xt.shape[1])
        pca_self = PCA(n_components=n_self, svd_solver="full", random_state=0).fit(Xt)
        Lself = pca_self.components_.T
        for k in range(1, min(21, n_self, n_pcs) + 1):
            rows.append(
                {
                    "tumor": t,
                    "nPCs": k,
                    "sc": fev(Xt, Lself[:, :k]),
                    "mb": fev(Xt, Lmb[:, :k]),
                    "rnd": fev(Xt, Lrnd[:, :k]),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "fev_by_tumor.csv", index=False)

    # Paper summary at 5 PCs: mb / sc * 100
    at5 = df[df["nPCs"] == 5].copy()
    at5["pct_of_sc"] = 100 * at5["mb"] / at5["sc"]
    at5["pct_rnd_of_sc"] = 100 * at5["rnd"] / at5["sc"]
    mb_mean, mb_sd = float(at5["pct_of_sc"].mean()), float(at5["pct_of_sc"].std())
    rnd_mean, rnd_sd = float(at5["pct_rnd_of_sc"].mean()), float(at5["pct_rnd_of_sc"].std())
    ttest = stats.ttest_ind(at5["pct_of_sc"], at5["pct_rnd_of_sc"], equal_var=False)
    print(f"At 5 PCs: METABRIC explains {mb_mean:.1f}±{mb_sd:.1f}% of sc-PC variance")
    print(f"          shuffled explains {rnd_mean:.1f}±{rnd_sd:.1f}%")
    print(f"          paper: 25.4±3% vs 1.5±0.4%; t-test p={ttest.pvalue:.4g}")

    summary = {
        "n_tumors": len(big_tumors),
        "mb_pct_of_sc_mean": mb_mean,
        "mb_pct_of_sc_sd": mb_sd,
        "rnd_pct_of_sc_mean": rnd_mean,
        "rnd_pct_of_sc_sd": rnd_sd,
        "ttest_pvalue": float(ttest.pvalue),
        "paper_mb": 25.4,
        "paper_rnd": 1.5,
    }
    (RESULTS / "fev_summary.json").write_text(json.dumps(summary, indent=2))

    FIGURES.mkdir(parents=True, exist_ok=True)
    # Panel like Fig 4d: box/strip of % variance at 5 PCs
    fig, ax = plt.subplots(figsize=(3.2, 3.5))
    plot_df = at5.melt(
        id_vars=["tumor"],
        value_vars=["pct_of_sc", "pct_rnd_of_sc"],
        var_name="source",
        value_name="pct",
    )
    plot_df["source"] = plot_df["source"].map(
        {"pct_of_sc": "METABRIC PCs", "pct_rnd_of_sc": "shuffled PCs"}
    )
    positions = {"METABRIC PCs": 0, "shuffled PCs": 1}
    for src, color in [("METABRIC PCs", "#4C78A8"), ("shuffled PCs", "#B0B0B0")]:
        vals = plot_df.loc[plot_df["source"] == src, "pct"].values
        x = np.full(len(vals), positions[src], dtype=float)
        x = x + (np.random.default_rng(0).random(len(vals)) - 0.5) * 0.08
        ax.scatter(x, vals, s=35, c=color, zorder=3, edgecolors="white", linewidths=0.4)
        ax.hlines(np.mean(vals), positions[src] - 0.15, positions[src] + 0.15, colors="black", lw=2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["METABRIC PCs", "shuffled PCs"])
    ax.set_ylabel("% of single-cell PC variance explained")
    ax.set_title(f"Fig. 4d style\nMETABRIC {mb_mean:.1f}% vs shuffle {rnd_mean:.1f}%\n(p={ttest.pvalue:.2g})")
    ax.set_ylim(0, max(40, mb_mean + 10))
    fig.tight_layout()
    out = FIGURES / "Figure_4d_fev.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print("Wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
