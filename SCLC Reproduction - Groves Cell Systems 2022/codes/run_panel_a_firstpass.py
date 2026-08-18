#!/usr/bin/env python3
"""First-pass Panel A: PCA, PCHA ESV curve, and 100-shuffle t-ratio tests."""

from pathlib import Path
import sys

SCLC = Path(__file__).resolve().parents[1]
ROOT = SCLC.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.io import load_expression_csv
from src.pca import fit_pca, cumulative_variance
from src.archetypes import esv_curve, permutation_t_ratio

OUT = SCLC / "results" / "panel_a"
OUT.mkdir(parents=True, exist_ok=True)

PAPER = {
    4: {"p": 0.059, "t_ratio": 0.247},
    5: {"p": 0.034, "t_ratio": 0.107},
    6: {"p": 0.016, "t_ratio": 0.043},
}

N_PCS = 12
K_ESV = range(2, 16)
K_PERM = (4, 5, 6)
N_PERM = 100
DELTA = 0.1
SEED = 0


def main():
    expr = load_expression_csv()
    # samples × genes
    X = expr.T.values
    print(f"Loaded expression: {expr.shape[0]} genes × {expr.shape[1]} samples")

    pca12, scores12 = fit_pca(X, n_components=N_PCS)
    cumvar = cumulative_variance(pca12)
    print(
        f"PCA: top {N_PCS} PCs explain {100 * cumvar[-1]:.1f}% variance "
        f"(paper ~50%)"
    )
    np.save(OUT / "pc_scores_12.npy", scores12)
    pd.DataFrame(
        {
            "pc": np.arange(1, N_PCS + 1),
            "explained_variance_ratio": pca12.explained_variance_ratio_,
            "cumulative": cumvar,
        }
    ).to_csv(OUT / "pca_variance.csv", index=False)

    print("Fitting PCHA for ESV curve, k=2..15 ...")
    fits = esv_curve(scores12, K_ESV, delta=DELTA, seed=SEED)
    esv = pd.DataFrame([{"k": f["k"], "esv": f["esv"]} for f in fits])
    esv["delta_esv"] = esv["esv"].diff()
    esv.to_csv(OUT / "esv_curve.csv", index=False)
    print(esv.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(esv["k"], esv["esv"], marker="o")
    axes[0].set_xlabel("k archetypes")
    axes[0].set_ylabel("Explained sample variance")
    axes[0].set_title("ESV vs k")
    axes[1].plot(esv["k"].iloc[1:], esv["delta_esv"].iloc[1:], marker="o")
    axes[1].set_xlabel("k archetypes")
    axes[1].set_ylabel("ESV(k) − ESV(k−1)")
    axes[1].set_title("ESV elbow")
    fig.tight_layout()
    fig.savefig(OUT / "esv_curve.png", dpi=150)
    plt.close(fig)

    perm_rows = []
    for k in K_PERM:
        print(f"t-ratio permutation test: k={k}, n_perm={N_PERM} ...")
        fit_k = next(f for f in fits if f["k"] == k)
        result = permutation_t_ratio(
            scores12,
            k=k,
            n_perm=N_PERM,
            delta=DELTA,
            seed=SEED + k,
            observed_archetypes=fit_k["archetypes"],
        )
        np.save(OUT / f"archetypes_k{k}.npy", result["archetypes"])
        np.save(OUT / f"S_k{k}.npy", fit_k["S"])
        np.save(OUT / f"null_t_ratios_k{k}.npy", result["null_t_ratios"])
        paper = PAPER[k]
        perm_rows.append(
            {
                "k": k,
                "t_ratio": result["t_ratio"],
                "p_value": result["p_value"],
                "n_perm": result["n_perm"],
                "n_success": result["n_success"],
                "paper_p": paper["p"],
                "paper_t_ratio": paper["t_ratio"],
            }
        )
        print(
            f"  k={k}: t={result['t_ratio']:.4f} (paper {paper['t_ratio']:.3f}), "
            f"p={result['p_value']:.3f} from {result['n_success']}/{N_PERM} "
            f"(paper p={paper['p']:.3f})"
        )

    table = pd.DataFrame(perm_rows)
    table.to_csv(OUT / "t_ratio_firstpass.csv", index=False)

    # 2-D visualization PCA (separate from the 12-PC fitting space)
    pca2, scores2 = fit_pca(X, n_components=2)
    arcs5 = next(f for f in fits if f["k"] == 5)["archetypes"]
    # map 12-PC archetypes into the 2-PC visualization via the gene-space reconstruction
    # X ≈ scores12 @ pca12.components_  (centered); archetype gene space:
    gene_space_arcs = arcs5 @ pca12.components_ + pca12.mean_
    arcs2 = pca2.transform(gene_space_arcs)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(scores2[:, 0], scores2[:, 1], s=18, alpha=0.7, c="#4C78A8", label="cell lines")
    ax.scatter(arcs2[:, 0], arcs2[:, 1], s=80, c="#E45756", zorder=3, label="archetypes")
    for i in range(5):
        for j in range(i + 1, 5):
            ax.plot(
                [arcs2[i, 0], arcs2[j, 0]],
                [arcs2[i, 1], arcs2[j, 1]],
                color="#E45756",
                lw=1,
                zorder=2,
            )
    ax.set_xlabel("PC1 (visualization)")
    ax.set_ylabel("PC2 (visualization)")
    ax.set_title("Panel A (first pass): k=5 polytope")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "polytope_k5_2d.png", dpi=150)
    plt.close(fig)

    print(f"\nWrote outputs to {OUT}")
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
