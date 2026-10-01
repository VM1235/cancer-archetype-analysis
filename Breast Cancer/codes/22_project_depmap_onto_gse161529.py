#!/usr/bin/env python3
"""Project 63 DepMap KS lines into Pal GSE161529 k=3 space (reverse Hausser Fig. 4).

Pal PCA + k=3 vertices frozen. DepMap genes mean/SD-matched to Pal MAGIC KS,
then pca.transform. No PCHA refit.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/22_project_depmap_onto_gse161529.py"
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.io import load_expression_csv
from src.pca import align_pca_signs, fit_pca

PAL_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse161529_magic" / "ks_genes_magic_imputed.npz"
)
PAL_FIT = BREAST / "results" / "panel_a_ks_genelist_gse161529"
CL_MATRIX = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
PAM50 = BREAST / "results" / "panel_b" / "pam50_labels_panelA.csv"
OUT = BREAST / "results" / "depmap_on_gse161529"
FIG = BREAST / "figures"

K = 3
ARC_COLORS = ["#4C78A8", "#F58518", "#E45756"]
PAM50_COLORS = {
    "Basal": "#54A24B",
    "Her2": "#E45756",
    "LumB": "#F58518",
    "LumA": "#4C78A8",
    "Normal": "#B0B0B0",
}


def barycentric(points, vertices):
    a = np.vstack([vertices.T, np.ones((1, vertices.shape[0]))])
    b = np.vstack([points.T, np.ones((1, points.shape[0]))])
    w, *_ = np.linalg.lstsq(a, b, rcond=None)
    return w.T


def draw_triangle(ax, arcs):
    pts = list(arcs) + [arcs[0]]
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#333333", lw=1.5, zorder=5)
    ax.scatter(
        arcs[:, 0], arcs[:, 1], s=120, c=ARC_COLORS[: arcs.shape[0]],
        edgecolors="white", linewidths=0.8, zorder=6,
    )
    for i in range(arcs.shape[0]):
        ax.annotate(
            str(i + 1), arcs[i], textcoords="offset points", xytext=(6, 6),
            fontsize=8, fontweight="bold",
        )


def match_to_reference(guest: pd.DataFrame, ref: pd.DataFrame) -> np.ndarray:
    ref_mean = ref.mean(axis=1).values.copy()
    ref_sd = ref.std(axis=1, ddof=1).values.copy()
    ref_sd[ref_sd < 1e-8] = 1.0
    g_mean = guest.mean(axis=1).values.copy()
    g_sd = guest.std(axis=1, ddof=1).values.copy()
    g_sd[g_sd < 1e-8] = 1.0
    return ((guest.values - g_mean[:, None]) / g_sd[:, None]) * ref_sd[:, None] + ref_mean[:, None]


def style():
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.bbox": "tight",
            "savefig.dpi": 240,
        }
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    for p in (PAL_MAGIC, PAL_FIT / "pc_scores.npy", CL_MATRIX, PAM50):
        if not Path(p).is_file():
            print(f"Missing {p}")
            return 1

    n_pcs = int((PAL_FIT / "n_pcs.txt").read_text().strip().splitlines()[0])
    saved = np.load(PAL_FIT / "pc_scores.npy")
    arcs = np.load(PAL_FIT / f"archetypes_k{K}_parti.npy")

    print("=== Pal MAGIC KS + frozen PCA ===")
    zpal = np.load(PAL_MAGIC, allow_pickle=True)
    pal_genes = [str(g) for g in zpal["genes"]]
    pal_expr = pd.DataFrame(zpal["imputed"], index=pal_genes, columns=list(zpal["barcodes"]))
    pca, pal_scores = fit_pca(pal_expr.T.values, n_components=n_pcs)
    pca, pal_scores, _ = align_pca_signs(pal_scores, saved, pca)
    print(f"Pal {pal_expr.shape[1]} cells; PCA rebuild max|diff|={np.max(np.abs(pal_scores - saved)):.4g}")

    print("=== DepMap KS Panel A matrix ===")
    cl = load_expression_csv(CL_MATRIX)
    cl.index = cl.index.astype(str)
    print(f"DepMap: {cl.shape[0]} genes x {cl.shape[1]} lines")
    missing = [g for g in pal_genes if g not in set(cl.index)]
    extra = [g for g in cl.index if g not in set(pal_genes)]
    print(f"DepMap missing vs Pal: {missing or 'none'}; extra: {len(extra)}")
    cl_aligned = cl.reindex(pal_genes)
    pal_mean = pal_expr.mean(axis=1)
    for g in missing:
        cl_aligned.loc[g] = pal_mean.loc[g]
        print(f"  filled {g} with Pal mean {float(pal_mean.loc[g]):.4g}")

    print("=== Mean/SD match DepMap → Pal, transform ===")
    cl_matched = match_to_reference(cl_aligned, pal_expr)
    cl_pc = pca.transform(cl_matched.T)
    w = barycentric(cl_pc[:, : arcs.shape[1]], arcs)
    inside = (w >= -1e-6).all(axis=1)
    nearest = w.argmax(axis=1) + 1
    print(f"DepMap inside Pal k=3: {int(inside.sum())}/{len(inside)} ({100 * inside.mean():.1f}%)")

    pal_w = barycentric(pal_scores[:, : arcs.shape[1]], arcs)
    pal_inside = (pal_w >= -1e-6).all(axis=1)
    print(f"Pal self (barycentric): {int(pal_inside.sum())}/{len(pal_inside)} ({100 * pal_inside.mean():.1f}%)")

    pam = pd.read_csv(PAM50)
    pam["cell_line"] = pam["cell_line"].astype(str)
    pam = pam.set_index("cell_line").reindex(cl.columns.astype(str))
    subtype = pam["pam50_subtype"].fillna("NA").astype(str).values
    print("PAM50 counts:\n", pd.Series(subtype).value_counts().to_string())

    out = pd.DataFrame(
        {
            "cell_line": cl.columns.astype(str),
            "pam50": subtype,
            "PC1": cl_pc[:, 0],
            "PC2": cl_pc[:, 1],
            "inside_simplex": inside,
            "nearest_archetype": nearest,
            "w_arc1": w[:, 0],
            "w_arc2": w[:, 1],
            "w_arc3": w[:, 2],
        }
    )
    out.to_csv(OUT / "depmap_in_gse161529_k3.csv", index=False)
    ct = pd.crosstab(out["pam50"], out["nearest_archetype"])
    ct.to_csv(OUT / "pam50_by_nearest_arc_counts.csv")
    (ct.div(ct.sum(axis=1), axis=0)).to_csv(OUT / "pam50_by_nearest_arc_rowfrac.csv")
    print("\nPAM50 vs nearest Pal arc (counts):")
    print(ct.to_string())
    print("\nInside by PAM50:")
    for st, g in out.groupby("pam50"):
        print(f"  {st}: {g.inside_simplex.mean():.1%} (n={len(g)})")

    report = {
        "host": "GSE161529",
        "guest": "DepMap_KS_63",
        "k": K,
        "n_lines": int(len(inside)),
        "n_inside": int(inside.sum()),
        "frac_inside": float(inside.mean()),
        "frac_inside_pal_self": float(pal_inside.mean()),
        "filled_genes": missing,
        "batch_match": "gene_mean_sd_to_pal",
        "inside_by_pam50": {
            st: {"frac": float(g.inside_simplex.mean()), "n": int(len(g))}
            for st, g in out.groupby("pam50")
        },
    }
    (OUT / "projection_report.json").write_text(json.dumps(report, indent=2) + "\n")

    style()
    rng = np.random.default_rng(0)
    take_pal = rng.choice(pal_scores.shape[0], size=min(8000, pal_scores.shape[0]), replace=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    ax = axes[0]
    ax.scatter(
        pal_scores[take_pal, 0], pal_scores[take_pal, 1],
        s=3, c="#B0B0B0", alpha=0.2, linewidths=0, label="GSE161529 Pal", zorder=1,
    )
    ax.scatter(
        cl_pc[:, 0], cl_pc[:, 1],
        s=42, c="#4C78A8", alpha=0.95, edgecolors="white", linewidths=0.4,
        label="DepMap lines", zorder=3,
    )
    draw_triangle(ax, arcs[:, :2])
    ax.set_xlabel("PC1 (Pal)")
    ax.set_ylabel("PC2 (Pal)")
    ax.set_title(
        f"DepMap in Pal k=3 triangle\n"
        f"{int(inside.sum())}/{len(inside)} inside ({100 * inside.mean():.1f}%)"
    )
    ax.legend(loc="best", fontsize=7)
    ax.set_aspect("equal", adjustable="datalim")

    ax = axes[1]
    for st, col in PAM50_COLORS.items():
        idx = np.where(subtype == st)[0]
        if len(idx) == 0:
            continue
        ax.scatter(
            cl_pc[idx, 0], cl_pc[idx, 1],
            s=48, c=col, alpha=0.95, edgecolors="white", linewidths=0.4,
            label=st, zorder=3,
        )
    draw_triangle(ax, arcs[:, :2])
    ax.set_xlabel("PC1 (Pal)")
    ax.set_ylabel("PC2 (Pal)")
    ax.set_title("DepMap lines by PAM50")
    ax.legend(loc="best", fontsize=7)
    ax.set_aspect("equal", adjustable="datalim")
    fig.suptitle(
        "DepMap KS lines projected into GSE161529-fit k=3 (no PCHA refit)",
        y=1.03,
    )
    fig.savefig(FIG / "Figure_4_depmap_on_gse161529_k3.png")
    plt.close(fig)
    print("Wrote", FIG / "Figure_4_depmap_on_gse161529_k3.png")
    print("Wrote", OUT / "projection_report.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
