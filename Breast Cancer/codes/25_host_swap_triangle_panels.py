#!/usr/bin/env python3
"""Side-by-side k=3 triangles: DepMap host (top) vs Pal host (bottom).

Each column is one dataset. Points use the saved projection coordinates.
The percent in each title is the saved inside-simplex fraction (not a
recount from the PC1–PC2 drawing). The wireframe is the PC1–PC2 view
of that fit, matching the earlier single-guest figures.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/25_host_swap_triangle_panels.py"
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
RES = BREAST / "results"
FIG = BREAST / "figures"

INSIDE_C = "#4C78A8"
OUTSIDE_C = "#E45756"
HOST_C = "#C8C8C8"
ARC_C = ["#4C78A8", "#F58518", "#E45756"]
N_SHOW = 6000
SEED = 0


def style():
    plt.rcParams.update(
        {
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.bbox": "tight",
            "savefig.dpi": 220,
        }
    )


def take(df: pd.DataFrame, n: int, rng: np.random.Generator) -> pd.DataFrame:
    if len(df) <= n:
        return df
    idx = rng.choice(len(df), size=n, replace=False)
    return df.iloc[np.sort(idx)]


def draw_triangle(ax, arcs: np.ndarray):
    pts = np.vstack([arcs[:, :2], arcs[:1, :2]])
    ax.plot(pts[:, 0], pts[:, 1], color="#222222", lw=1.3, zorder=5)
    ax.scatter(
        arcs[:, 0], arcs[:, 1], s=36, c=ARC_C[: len(arcs)],
        edgecolors="white", linewidths=0.6, zorder=6,
    )
    for i, (x, y) in enumerate(arcs[:, :2]):
        ax.annotate(
            str(i + 1), (x, y), textcoords="offset points", xytext=(4, 4),
            fontsize=7, fontweight="bold", zorder=7,
        )


def scatter_in_out(ax, df: pd.DataFrame, s: float):
    out = df.loc[~df["inside"].to_numpy()]
    inn = df.loc[df["inside"].to_numpy()]
    if len(out):
        ax.scatter(
            out["PC1"], out["PC2"], s=s, c=OUTSIDE_C, alpha=0.28,
            linewidths=0, zorder=2, label="outside",
        )
    if len(inn):
        ax.scatter(
            inn["PC1"], inn["PC2"], s=s, c=INSIDE_C, alpha=0.35,
            linewidths=0, zorder=3, label="inside",
        )


def limits(ax, xs, ys, arcs):
    pts = np.column_stack([np.asarray(xs, float), np.asarray(ys, float)])
    if len(pts) > 20:
        lo = np.quantile(pts, 0.005, axis=0)
        hi = np.quantile(pts, 0.995, axis=0)
    else:
        lo = pts.min(axis=0)
        hi = pts.max(axis=0)
    lo = np.minimum(lo, arcs[:, :2].min(axis=0))
    hi = np.maximum(hi, arcs[:, :2].max(axis=0))
    pad = 0.08 * np.maximum(hi - lo, 1e-3)
    ax.set_xlim(lo[0] - pad[0], hi[0] + pad[0])
    ax.set_ylim(lo[1] - pad[1], hi[1] + pad[1])
    ax.set_aspect("equal", adjustable="box")


def panel(ax, df, arcs, host_xy, title, point_s, host_s):
    if host_xy is not None and len(host_xy):
        ax.scatter(
            host_xy[:, 0], host_xy[:, 1], s=host_s, c=HOST_C, alpha=0.45,
            linewidths=0, zorder=1,
        )
    scatter_in_out(ax, df, point_s)
    draw_triangle(ax, arcs)
    xs = np.concatenate([df["PC1"].to_numpy(), arcs[:, 0]])
    ys = np.concatenate([df["PC2"].to_numpy(), arcs[:, 1]])
    if host_xy is not None and len(host_xy):
        xs = np.concatenate([xs, host_xy[:, 0]])
        ys = np.concatenate([ys, host_xy[:, 1]])
    limits(ax, xs, ys, arcs)
    ax.set_title(title)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")


def from_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=lambda c: c in {"PC1", "PC2", "inside_simplex"})
    df = df.rename(columns={"inside_simplex": "inside"})
    df["inside"] = df["inside"].astype(bool)
    return df


def arcs_csv(path: Path) -> np.ndarray:
    a = pd.read_csv(path)
    return a[["PC1", "PC2"]].to_numpy(float)


def main():
    style()
    rng = np.random.default_rng(SEED)
    FIG.mkdir(parents=True, exist_ok=True)

    depmap_pc = np.load(RES / "panel_a_ks_genelist" / "pc_scores.npy")
    depmap_arcs = np.load(RES / "panel_a_ks_genelist" / "archetypes_k3_parti.npy")[:, :2]
    # Self-containment in the Panel A plane (k=3 uses the first 2 PCs).
    a = np.vstack([depmap_arcs.T, np.ones((1, 3))])
    b = np.vstack([depmap_pc[:, :2].T, np.ones((1, depmap_pc.shape[0]))])
    w, *_ = np.linalg.lstsq(a, b, rcond=None)
    depmap_inside = (w.T >= -1e-6).all(axis=1)
    depmap_self = pd.DataFrame(
        {"PC1": depmap_pc[:, 0], "PC2": depmap_pc[:, 1], "inside": depmap_inside}
    )

    pal_pc = np.load(RES / "panel_a_ks_genelist_gse161529" / "pc_scores.npy")
    pal_arcs = np.load(RES / "panel_a_ks_genelist_gse161529" / "archetypes_k3_parti.npy")[:, :2]
    pal_w_csv = RES / "panel_b_gse161529" / "distances_subtype.csv"
    # inside flag for the Pal self cloud: barycentric in the 2 PCs the
    # wireframe uses is not the 6-D test. Use the published self fraction
    # from the host-swap table by recomputing with the saved archetype
    # matrix in its native dimension.
    pal_arcs_full = np.load(RES / "panel_a_ks_genelist_gse161529" / "archetypes_k3_parti.npy")
    n_vol = pal_arcs_full.shape[1]
    aa = np.vstack([pal_arcs_full.T, np.ones((1, pal_arcs_full.shape[0]))])
    bb = np.vstack([pal_pc[:, :n_vol].T, np.ones((1, pal_pc.shape[0]))])
    ww, *_ = np.linalg.lstsq(aa, bb, rcond=None)
    pal_inside = (ww.T >= -1e-6).all(axis=1)
    pal_self = pd.DataFrame(
        {"PC1": pal_pc[:, 0], "PC2": pal_pc[:, 1], "inside": pal_inside}
    )

    guests_depmap = {
        "GSE161529\n(Pal)": (
            from_csv(RES / "panel_a_ks_genelist_sc_gse161529_magic" / "sc_cells_in_panelA_space.csv"),
            arcs_csv(RES / "panel_a_ks_genelist_sc_gse161529_magic" / "archetypes_in_shared_pc.csv"),
            pd.read_csv(
                RES / "panel_a_ks_genelist_sc_gse161529_magic" / "cellline_pc_scores_shared_genes.csv"
            )[["PC1", "PC2"]].to_numpy(float),
        ),
        "GSE176078\n(Wu)": (
            from_csv(RES / "panel_a_ks_genelist_sc_gse176078_magic" / "sc_cells_in_panelA_space.csv"),
            arcs_csv(RES / "panel_a_ks_genelist_sc_gse176078_magic" / "archetypes_in_shared_pc.csv"),
            pd.read_csv(
                RES / "panel_a_ks_genelist_sc_gse176078_magic" / "cellline_pc_scores_shared_genes.csv"
            )[["PC1", "PC2"]].to_numpy(float),
        ),
        "GSE173634": (
            from_csv(RES / "panel_a_ks_genelist_sc_gse173634_magic" / "sc_cells_in_panelA_space.csv"),
            arcs_csv(RES / "panel_a_ks_genelist_sc_gse173634_magic" / "archetypes_in_shared_pc.csv"),
            pd.read_csv(
                RES / "panel_a_ks_genelist_sc_gse173634_magic" / "cellline_pc_scores_shared_genes.csv"
            )[["PC1", "PC2"]].to_numpy(float),
        ),
    }
    guests_pal = {
        "GSE176078\n(Wu)": from_csv(RES / "gse176078_on_gse161529" / "gse176078_in_gse161529_k3.csv"),
        "GSE173634": from_csv(RES / "gse173634_on_gse161529" / "gse173634_in_gse161529_k3.csv"),
        "DepMap\n(63 lines)": from_csv(RES / "depmap_on_gse161529" / "depmap_in_gse161529_k3.csv"),
    }

    fig, axes = plt.subplots(2, 4, figsize=(13.6, 7.4))
    pal_bg = take(pal_self, N_SHOW, rng)[["PC1", "PC2"]].to_numpy(float)

    # Top row: DepMap host. Last column is the lines themselves.
    order_top = ["GSE161529\n(Pal)", "GSE176078\n(Wu)", "GSE173634"]
    for j, name in enumerate(order_top):
        df, arcs, host = guests_depmap[name]
        shown = take(df, N_SHOW, rng)
        frac = float(df["inside"].mean())
        panel(
            axes[0, j], shown, arcs, host,
            f"{name}\n{100 * frac:.1f}% inside",
            point_s=3, host_s=8,
        )
    shown_lines = depmap_self
    panel(
        axes[0, 3], shown_lines, depmap_arcs, None,
        f"DepMap\n(63 lines)\n{100 * depmap_self['inside'].mean():.1f}% inside",
        point_s=16, host_s=8,
    )
    axes[0, 0].set_ylabel("Host: DepMap KS k=3\nPC2")

    # Bottom row: Pal host.
    shown_pal = take(pal_self, N_SHOW, rng)
    panel(
        axes[1, 0], shown_pal, pal_arcs, None,
        f"GSE161529\n(Pal)\n{100 * pal_self['inside'].mean():.1f}% inside",
        point_s=3, host_s=8,
    )
    for j, name in enumerate(["GSE176078\n(Wu)", "GSE173634", "DepMap\n(63 lines)"], start=1):
        df = guests_pal[name]
        shown = take(df, N_SHOW, rng)
        frac = float(df["inside"].mean())
        s = 18 if "DepMap" in name else 3
        panel(
            axes[1, j], shown, pal_arcs, pal_bg,
            f"{name}\n{100 * frac:.1f}% inside",
            point_s=s, host_s=2,
        )
    axes[1, 0].set_ylabel("Host: Pal GSE161529 k=3\nPC2")

    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=INSIDE_C, markersize=6, label="inside simplex"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=OUTSIDE_C, markersize=6, label="outside simplex"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=HOST_C, markersize=6, label="host points"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(
        "Same four datasets in two k=3 triangles\n"
        "Top: archetypes fit on DepMap lines. Bottom: archetypes fit on Pal cells. "
        "PC axes are not shared across panels.",
        fontsize=11, y=1.08,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG / "Figure_host_swap_triangles_k3.png"
    fig.savefig(out)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)
    print("DepMap self inside", f"{depmap_self['inside'].mean():.3%}", "n", len(depmap_self))
    print("Pal self inside", f"{pal_self['inside'].mean():.3%}", "n", len(pal_self), "n_vol", n_vol)
    print("unused", pal_w_csv.name)
    print("Wrote", out)


if __name__ == "__main__":
    main()
