#!/usr/bin/env python3
"""Project GSE173634 (Gambardella cell-line sc) into GSE161529 (Pal) self-fit space.

Same protocol as 18_project_gse176078_onto_gse161529.py: Pal PCA/PCHA frozen,
guest gene-mean/SD matched to Pal, pca.transform, barycentric k=3 and k=4.
No PCHA refit. GSE173634 MAGIC cache also lacks LHFPL6; filled with Pal mean.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/19_project_gse173634_onto_gse161529.py"
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
GUEST_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse173634_magic" / "ks_genes_magic_imputed.npz"
)
PAL_FIT = BREAST / "results" / "panel_a_ks_genelist_gse161529"
SERIES_PATH = ROOT / "Sahoo" / "GSE173634_series_matrix"
OUT = BREAST / "results" / "gse173634_on_gse161529"
FIG = BREAST / "figures"

K_VALUES = (3, 4)
ARC_COLORS = ["#4C78A8", "#F58518", "#E45756", "#72B7B2"]
SUBTYPE_COLORS = {
    "LA": "#4C78A8",
    "LB": "#F58518",
    "H": "#E45756",
    "TNA": "#72B7B2",
    "TNB": "#54A24B",
    "Basal-like": "#B279A2",
    "NA": "#CCCCCC",
}


def parse_series_subtypes(path: Path) -> dict[str, str]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    titles, subtypes = None, None
    for line in lines:
        if line.startswith("!Sample_title"):
            titles = [t.strip().strip('"') for t in line.split("\t")[1:]]
        if line.startswith("!Sample_characteristics_ch1") and "breast cancer subtype:" in line:
            subtypes = [
                t.strip().strip('"').replace("breast cancer subtype: ", "")
                for t in line.split("\t")[1:]
            ]
    if not titles or not subtypes:
        return {}
    out = {}
    for t, s in zip(titles, subtypes):
        key = t.rsplit("_", 1)[0] if t[-2:] in ("_1", "_2") else t
        out[key.upper().replace("-", "")] = s
        out[t.upper().replace("-", "")] = s
    return out


def lookup_subtype(cl: str, subtypes: dict[str, str]) -> str:
    for key in (cl, cl.replace("MDAMB", "MDA-MB"), cl):
        if key in subtypes:
            return subtypes[key]
    for kk, v in subtypes.items():
        if kk.replace("-", "") == cl:
            return v
    return "NA"


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
    ref_mean = ref.mean(axis=1).values.copy()
    ref_sd = ref.std(axis=1, ddof=1).values.copy()
    ref_sd[ref_sd < 1e-8] = 1.0
    g_mean = guest.mean(axis=1).values.copy()
    g_sd = guest.std(axis=1, ddof=1).values.copy()
    g_sd[g_sd < 1e-8] = 1.0
    return ((guest.values - g_mean[:, None]) / g_sd[:, None]) * ref_sd[:, None] + ref_mean[:, None]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    for p in (PAL_MAGIC, GUEST_MAGIC, PAL_FIT / "pc_scores.npy", SERIES_PATH):
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

    print("=== GSE173634 MAGIC KS ===")
    zg = np.load(GUEST_MAGIC, allow_pickle=True)
    guest_genes = [str(g) for g in zg["genes"]]
    guest_expr = pd.DataFrame(zg["imputed"], index=guest_genes, columns=list(zg["barcodes"]))
    print(f"GSE173634: {guest_expr.shape[0]} genes x {guest_expr.shape[1]} cells")

    missing = [g for g in pal_genes if g not in set(guest_genes)]
    extra = [g for g in guest_genes if g not in set(pal_genes)]
    print(f"Guest missing vs Pal: {missing or 'none'}; extra: {extra or 'none'}")
    guest_aligned = guest_expr.reindex(pal_genes)
    pal_mean = pal_expr.mean(axis=1)
    for g in missing:
        guest_aligned.loc[g] = pal_mean.loc[g]
        print(f"  filled {g} with Pal mean {float(pal_mean.loc[g]):.4g}")

    print("=== Gene-wise mean/SD match GSE173634 → Pal, transform with Pal PCA ===")
    guest_matched = match_to_reference(guest_aligned, pal_expr)
    guest_pc = pca.transform(guest_matched.T)
    print(f"Projected {guest_pc.shape[0]} GSE173634 cells into Pal {guest_pc.shape[1]}-D PC space")

    barcodes = pd.Index(guest_expr.columns.astype(str))
    cell_line = barcodes.str.split("_").str[0].str.upper().str.replace("-", "", regex=False)
    subtypes = parse_series_subtypes(SERIES_PATH)
    subtype = np.array([lookup_subtype(cl, subtypes) for cl in cell_line])

    report = {
        "host": "GSE161529",
        "guest": "GSE173634",
        "n_pal_cells": int(pal_expr.shape[1]),
        "n_guest_cells": int(guest_pc.shape[0]),
        "n_genes_pca": len(pal_genes),
        "n_genes_guest_native": len(guest_genes),
        "filled_genes": missing,
        "n_pcs": n_pcs,
        "k_elbow": k_elbow,
        "batch_match": "gene_mean_sd_to_pal",
        "p_values": None,
        "by_k": {},
    }

    rng = np.random.default_rng(0)
    take_pal = rng.choice(pal_scores.shape[0], size=min(8000, pal_scores.shape[0]), replace=False)
    take_g = rng.choice(guest_pc.shape[0], size=min(8000, guest_pc.shape[0]), replace=False)

    for k in K_VALUES:
        arcs = np.load(PAL_FIT / f"archetypes_k{k}_parti.npy")
        n_vol = arcs.shape[1]
        pal_w = barycentric(pal_scores[:, :n_vol], arcs)
        guest_w = barycentric(guest_pc[:, :n_vol], arcs)
        pal_inside = (pal_w >= -1e-6).all(axis=1)
        guest_inside = (guest_w >= -1e-6).all(axis=1)

        print(f"\n--- k={k} ({n_vol}-D) ---")
        print(
            f"Pal (self, barycentric): {int(pal_inside.sum())}/{len(pal_inside)} "
            f"({100 * pal_inside.mean():.1f}%)"
        )
        print(
            f"GSE173634 inside Pal simplex: {int(guest_inside.sum())}/{len(guest_inside)} "
            f"({100 * guest_inside.mean():.1f}%)"
        )

        by_st = {}
        print("GSE173634 inside by cell-line subtype:")
        for st in sorted(set(subtype)):
            idx = subtype == st
            frac = float(guest_inside[idx].mean()) if idx.any() else float("nan")
            by_st[st] = {"frac_inside": frac, "n": int(idx.sum())}
            print(f"  {st}: {frac:.1%} (n={int(idx.sum())})")

        cols = {
            "barcode": list(barcodes),
            "cell_line": cell_line.to_numpy(),
            "subtype": subtype,
            "inside_simplex": guest_inside,
            "nearest_archetype": guest_w.argmax(axis=1) + 1,
        }
        for j in range(min(3, guest_pc.shape[1])):
            cols[f"PC{j + 1}"] = guest_pc[:, j]
        for j in range(k):
            cols[f"w_arc{j + 1}"] = guest_w[:, j]
        pd.DataFrame(cols).to_csv(OUT / f"gse173634_in_gse161529_k{k}.csv", index=False)

        report["by_k"][str(k)] = {
            "n_vol": n_vol,
            "frac_inside_pal_self": float(pal_inside.mean()),
            "n_inside_pal_self": int(pal_inside.sum()),
            "frac_inside_guest": float(guest_inside.mean()),
            "n_inside_guest": int(guest_inside.sum()),
            "n_guest": int(len(guest_inside)),
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
                guest_pc[take_g, 0], guest_pc[take_g, 1],
                s=3, c="#4C78A8", alpha=0.2, linewidths=0, label="GSE173634 sc (guest)", zorder=2,
            )
            draw_triangle(ax, arcs[:, :2])
            ax.set_xlabel("PC1 (Pal)")
            ax.set_ylabel("PC2 (Pal)")
            ax.set_title(
                f"GSE173634 in Pal k=3 triangle\n"
                f"{int(guest_inside.sum())}/{len(guest_inside)} inside "
                f"({100 * guest_inside.mean():.1f}%)"
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
                    guest_pc[idx_plot, 0], guest_pc[idx_plot, 1],
                    s=3, c=col, alpha=0.25, linewidths=0, label=st, zorder=2,
                )
            draw_triangle(ax, arcs[:, :2])
            ax.set_xlabel("PC1 (Pal)")
            ax.set_ylabel("PC2 (Pal)")
            ax.set_title("GSE173634 cells by cell-line subtype")
            ax.legend(loc="best", fontsize=7)
            ax.set_aspect("equal", adjustable="datalim")
            fig.suptitle(
                "GSE173634 projected into GSE161529-fit k=3 (point estimate, no p-value)",
                y=1.03,
            )
            fig.savefig(FIG / "Figure_4_gse173634_on_gse161529_k3.png")
            plt.close(fig)
            print("Wrote", FIG / "Figure_4_gse173634_on_gse161529_k3.png")
        else:
            fig = plt.figure(figsize=(12, 5.6))
            ax = fig.add_subplot(1, 2, 1, projection="3d")
            ax.scatter(
                pal_scores[take_pal, 0], pal_scores[take_pal, 1], pal_scores[take_pal, 2],
                s=2, c="#B0B0B0", alpha=0.08, linewidths=0,
            )
            ax.scatter(
                guest_pc[take_g, 0], guest_pc[take_g, 1], guest_pc[take_g, 2],
                s=2, c="#4C78A8", alpha=0.12, linewidths=0,
            )
            draw_tetrahedron(ax, arcs[:, :3])
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.set_zlabel("PC3")
            ax.set_title(
                f"GSE173634 in Pal k=4 tetrahedron\n"
                f"{int(guest_inside.sum())}/{len(guest_inside)} inside "
                f"({100 * guest_inside.mean():.1f}%)"
            )
            fit_3d_view(ax, np.vstack([arcs[:, :3], guest_pc[take_g, :3], pal_scores[take_pal, :3]]))

            ax = fig.add_subplot(1, 2, 2, projection="3d")
            for st, col in SUBTYPE_COLORS.items():
                idx = np.where(subtype == st)[0]
                if len(idx) == 0:
                    continue
                idx_plot = idx if len(idx) <= 4000 else rng.choice(idx, 4000, replace=False)
                ax.scatter(
                    guest_pc[idx_plot, 0], guest_pc[idx_plot, 1], guest_pc[idx_plot, 2],
                    s=2, c=col, alpha=0.15, linewidths=0, label=st,
                )
            draw_tetrahedron(ax, arcs[:, :3])
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
            ax.set_zlabel("PC3")
            ax.set_title("GSE173634 by cell-line subtype")
            ax.legend(loc="best", fontsize=7)
            fit_3d_view(ax, np.vstack([arcs[:, :3], guest_pc[take_g, :3]]))
            fig.suptitle(
                "GSE173634 projected into GSE161529-fit k=4 (elbow; no p-value)",
                y=1.02,
            )
            fig.savefig(FIG / "Figure_4_gse173634_on_gse161529_k4.png")
            plt.close(fig)
            print("Wrote", FIG / "Figure_4_gse173634_on_gse161529_k4.png")

    (OUT / "projection_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print("\nWrote", OUT / "projection_report.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
