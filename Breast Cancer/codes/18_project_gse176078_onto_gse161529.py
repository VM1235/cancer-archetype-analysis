#!/usr/bin/env python3
"""Project GSE176078 (Wu) MAGIC KS cells into GSE161529 (Pal) self-fit space.

Pal PCA + PCHA vertices are frozen (17.9_run_panelA_ks_genelist_gse161529.py).
Wu cells are gene-mean/SD matched to Pal, then transformed with Pal's PCA.
No PCHA refit. Containment is barycentric in (k-1)-D at k=3 and k=4.

Wu MAGIC cache is missing LHFPL6 (205/206 KS genes). That gene is filled
with Pal's mean so the frozen 206-gene PCA stays valid.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/18_project_gse176078_onto_gse161529.py"
"""

from __future__ import annotations

from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

from src.pca import align_pca_signs, fit_pca

PAL_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse161529_magic" / "ks_genes_magic_imputed.npz"
)
WU_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse176078_magic" / "ks_genes_magic_imputed.npz"
)
PAL_FIT = BREAST / "results" / "panel_a_ks_genelist_gse161529"
META_PATH = ROOT / "Sahoo" / "GSE176078_sc" / "Wu_etal_2021_BRCA_scRNASeq" / "metadata.csv"
OUT = BREAST / "results" / "gse176078_on_gse161529"
FIG = BREAST / "figures"

K_VALUES = (3, 4)
ARC_COLORS = ["#4C78A8", "#F58518", "#E45756", "#72B7B2"]
SUBTYPE_COLORS = {"ER+": "#4C78A8", "HER2+": "#E45756", "TNBC": "#54A24B"}


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
        ax.annotate(str(i + 1), arcs[i], textcoords="offset points", xytext=(6, 6),
                    fontsize=8, fontweight="bold")


def draw_tetrahedron(ax, arcs):
    k = arcs.shape[0]
    for i in range(k):
        for j in range(i + 1, k):
            ax.plot(
                [arcs[i, 0], arcs[j, 0]],
                [arcs[i, 1], arcs[j, 1]],
                [arcs[i, 2], arcs[j, 2]],
                color="#666666", lw=1.2, alpha=0.9,
            )
    ax.scatter(
        arcs[:, 0], arcs[:, 1], arcs[:, 2], s=160, c=ARC_COLORS[:k],
        edgecolors="white", linewidths=0.8, depthshade=False, zorder=10,
    )


def fit_3d_view(ax, pts, pad=0.08):
    spans = np.ptp(pts, axis=0)
    spans = np.maximum(spans, 1e-6)
    for i, setter in enumerate([ax.set_xlim, ax.set_ylim, ax.set_zlim]):
        lo, hi = pts[:, i].min(), pts[:, i].max()
        p = spans[i] * pad
        setter(lo - p, hi + p)
    try:
        ax.set_box_aspect(spans / spans.max())
    except AttributeError:
        pass


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


def match_to_reference(guest: pd.DataFrame, ref: pd.DataFrame) -> np.ndarray:
    """Gene-wise mean/SD match guest onto reference (genes x cells)."""
    ref_mean = ref.mean(axis=1).values
    ref_sd = ref.std(axis=1, ddof=1).values
    ref_sd[ref_sd < 1e-8] = 1.0
    g_mean = guest.mean(axis=1).values
    g_sd = guest.std(axis=1, ddof=1).values
    g_sd[g_sd < 1e-8] = 1.0
    return ((guest.values - g_mean[:, None]) / g_sd[:, None]) * ref_sd[:, None] + ref_mean[:, None]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    for p in (PAL_MAGIC, WU_MAGIC, PAL_FIT / "pc_scores.npy", META_PATH):
        if not p.is_file():
            print(f"Missing {p}")
            return 1

    n_pcs = int((PAL_FIT / "n_pcs.txt").read_text().strip().splitlines()[0])
    saved = np.load(PAL_FIT / "pc_scores.npy")
    k_elbow = int((PAL_FIT / "k_from_finder.txt").read_text().strip().splitlines()[0])

    print("=== Pal GSE161529 MAGIC KS + frozen PCA ===")
    zpal = np.load(PAL_MAGIC, allow_pickle=True)
    pal_genes = [str(g) for g in zpal["genes"]]
    pal_expr = pd.DataFrame(zpal["imputed"], index=pal_genes, columns=list(zpal["barcodes"]))
    print(f"Pal: {pal_expr.shape[0]} genes x {pal_expr.shape[1]} cells")

    pca, pal_scores = fit_pca(pal_expr.T.values, n_components=n_pcs)
    pca, pal_scores, _ = align_pca_signs(pal_scores, saved, pca)
    print(f"Rebuilt Pal PCA dim={n_pcs}; max |score diff| vs saved = {np.max(np.abs(pal_scores - saved)):.4g}")

    print("=== Wu GSE176078 MAGIC KS ===")
    zwu = np.load(WU_MAGIC, allow_pickle=True)
    wu_genes = [str(g) for g in zwu["genes"]]
    wu_expr = pd.DataFrame(zwu["imputed"], index=wu_genes, columns=list(zwu["barcodes"]))
    print(f"Wu: {wu_expr.shape[0]} genes x {wu_expr.shape[1]} cells")

    missing = [g for g in pal_genes if g not in set(wu_genes)]
    extra = [g for g in wu_genes if g not in set(pal_genes)]
    print(f"Wu missing vs Pal: {missing or 'none'}; Wu extra: {extra or 'none'}")
    wu_aligned = wu_expr.reindex(pal_genes)
    pal_mean = pal_expr.mean(axis=1)
    for g in missing:
        wu_aligned.loc[g] = pal_mean.loc[g]
        print(f"  filled {g} with Pal mean {float(pal_mean.loc[g]):.4g}")

    print("=== Gene-wise mean/SD match Wu → Pal, transform with Pal PCA ===")
    wu_matched = match_to_reference(wu_aligned, pal_expr)
    wu_pc = pca.transform(wu_matched.T)
    print(f"Projected {wu_pc.shape[0]} Wu cells into Pal {wu_pc.shape[1]}-D PC space")

    meta = pd.read_csv(META_PATH)
    if "Unnamed: 0" in meta.columns:
        meta = meta.rename(columns={"Unnamed: 0": "barcode"})
    meta = meta.set_index("barcode").reindex(wu_expr.columns)
    subtype = meta["subtype"].astype(str).values if "subtype" in meta.columns else np.array(["NA"] * wu_pc.shape[0])
    patient = (
        meta["orig.ident"].astype(str).values
        if "orig.ident" in meta.columns
        else np.array(["NA"] * wu_pc.shape[0])
    )

    report = {
        "host": "GSE161529",
        "guest": "GSE176078",
        "n_pal_cells": int(pal_expr.shape[1]),
        "n_wu_cells": int(wu_pc.shape[0]),
        "n_genes_pca": len(pal_genes),
        "n_genes_wu_native": len(wu_genes),
        "filled_genes": missing,
        "n_pcs": n_pcs,
        "k_elbow": k_elbow,
        "batch_match": "gene_mean_sd_to_pal",
        "p_values": None,
        "by_k": {},
    }

    rng = np.random.default_rng(0)
    take_pal = rng.choice(pal_scores.shape[0], size=min(8000, pal_scores.shape[0]), replace=False)
    take_wu = rng.choice(wu_pc.shape[0], size=min(8000, wu_pc.shape[0]), replace=False)

    for k in K_VALUES:
        arcs = np.load(PAL_FIT / f"archetypes_k{k}_parti.npy")
        n_vol = arcs.shape[1]
        pal_w = barycentric(pal_scores[:, :n_vol], arcs)
        wu_w = barycentric(wu_pc[:, :n_vol], arcs)
        pal_inside = (pal_w >= -1e-6).all(axis=1)
        wu_inside = (wu_w >= -1e-6).all(axis=1)

        print(f"\n--- k={k} ({n_vol}-D) ---")
        print(
            f"Pal (self, barycentric): {int(pal_inside.sum())}/{len(pal_inside)} "
            f"({100 * pal_inside.mean():.1f}%)"
        )
        print(
            f"Wu inside Pal simplex: {int(wu_inside.sum())}/{len(wu_inside)} "
            f"({100 * wu_inside.mean():.1f}%)"
        )

        by_st = {}
        print("Wu inside by subtype:")
        for st in sorted(set(subtype)):
            idx = subtype == st
            frac = float(wu_inside[idx].mean()) if idx.any() else float("nan")
            by_st[st] = {"frac_inside": frac, "n": int(idx.sum())}
            print(f"  {st}: {frac:.1%} (n={int(idx.sum())})")

        cols = {
            "barcode": list(wu_expr.columns),
            "patient": patient,
            "subtype": subtype,
            "inside_simplex": wu_inside,
            "nearest_archetype": wu_w.argmax(axis=1) + 1,
        }
        for j in range(min(3, wu_pc.shape[1])):
            cols[f"PC{j + 1}"] = wu_pc[:, j]
        for j in range(k):
            cols[f"w_arc{j + 1}"] = wu_w[:, j]
        out_df = pd.DataFrame(cols)
        out_df.to_csv(OUT / f"gse176078_in_gse161529_k{k}.csv", index=False)

        report["by_k"][str(k)] = {
            "n_vol": n_vol,
            "frac_inside_pal_self": float(pal_inside.mean()),
            "n_inside_pal_self": int(pal_inside.sum()),
            "frac_inside_wu": float(wu_inside.mean()),
            "n_inside_wu": int(wu_inside.sum()),
            "n_wu": int(len(wu_inside)),
            "by_subtype": by_st,
        }

        style()
        if k == 3:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
            ax = axes[0]
            ax.scatter(
                pal_scores[take_pal, 0], pal_scores[take_pal, 1],
                s=3, c="#B0B0B0", alpha=0.25, linewidths=0, label="GSE161529 Pal (host)", zorder=1,
            )
            ax.scatter(
                wu_pc[take_wu, 0], wu_pc[take_wu, 1],
                s=3, c="#4C78A8", alpha=0.2, linewidths=0, label="GSE176078 Wu (guest)", zorder=2,
            )
            draw_triangle(ax, arcs[:, :2])
            ax.set_xlabel("PC1 (Pal)")
            ax.set_ylabel("PC2 (Pal)")
            ax.set_title(
                f"Wu in Pal k=3 triangle\n"
                f"{int(wu_inside.sum())}/{len(wu_inside)} inside ({100 * wu_inside.mean():.1f}%)"
            )
            ax.legend(loc="best", fontsize=7)
            ax.set_aspect("equal", adjustable="datalim")

            ax = axes[1]
            for st, col in SUBTYPE_COLORS.items():
                idx = np.where(subtype == st)[0]
                if len(idx) == 0:
                    continue
                idx_plot = idx if len(idx) <= 4000 else rng.choice(idx, 4000, replace=False)
                ax.scatter(
                    wu_pc[idx_plot, 0], wu_pc[idx_plot, 1],
                    s=3, c=col, alpha=0.25, linewidths=0, label=st, zorder=2,
                )
            draw_triangle(ax, arcs[:, :2])
            ax.set_xlabel("PC1 (Pal)")
            ax.set_ylabel("PC2 (Pal)")
            ax.set_title("Wu cells by clinical subtype")
            ax.legend(loc="best", fontsize=7)
            ax.set_aspect("equal", adjustable="datalim")
            fig.suptitle(
                "GSE176078 projected into GSE161529-fit k=3 (point estimate, no p-value)",
                y=1.03,
            )
            fig.savefig(FIG / "Figure_4_gse176078_on_gse161529_k3.png")
            plt.close(fig)
            print("Wrote", FIG / "Figure_4_gse176078_on_gse161529_k3.png")
        else:
            fig = plt.figure(figsize=(12, 5.6))
            ax = fig.add_subplot(1, 2, 1, projection="3d")
            ax.scatter(
                pal_scores[take_pal, 0], pal_scores[take_pal, 1], pal_scores[take_pal, 2],
                s=2, c="#B0B0B0", alpha=0.08, linewidths=0,
            )
            ax.scatter(
                wu_pc[take_wu, 0], wu_pc[take_wu, 1], wu_pc[take_wu, 2],
                s=2, c="#4C78A8", alpha=0.12, linewidths=0,
            )
            draw_tetrahedron(ax, arcs[:, :3])
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.set_zlabel("PC3")
            ax.set_title(
                f"Wu in Pal k=4 tetrahedron\n"
                f"{int(wu_inside.sum())}/{len(wu_inside)} inside ({100 * wu_inside.mean():.1f}%)"
            )
            fit_3d_view(ax, np.vstack([arcs[:, :3], wu_pc[take_wu, :3], pal_scores[take_pal, :3]]))

            ax = fig.add_subplot(1, 2, 2, projection="3d")
            for st, col in SUBTYPE_COLORS.items():
                idx = np.where(subtype == st)[0]
                if len(idx) == 0:
                    continue
                idx_plot = idx if len(idx) <= 4000 else rng.choice(idx, 4000, replace=False)
                ax.scatter(
                    wu_pc[idx_plot, 0], wu_pc[idx_plot, 1], wu_pc[idx_plot, 2],
                    s=2, c=col, alpha=0.15, linewidths=0, label=st,
                )
            draw_tetrahedron(ax, arcs[:, :3])
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.set_zlabel("PC3")
            ax.set_title("Wu by clinical subtype")
            ax.legend(loc="best", fontsize=7)
            fit_3d_view(ax, np.vstack([arcs[:, :3], wu_pc[take_wu, :3]]))
            fig.suptitle(
                "GSE176078 projected into GSE161529-fit k=4 (elbow; no p-value)",
                y=1.02,
            )
            fig.savefig(FIG / "Figure_4_gse176078_on_gse161529_k4.png")
            plt.close(fig)
            print("Wrote", FIG / "Figure_4_gse176078_on_gse161529_k4.png")

    (OUT / "projection_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print("\nWrote", OUT / "projection_report.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
