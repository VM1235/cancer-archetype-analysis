#!/usr/bin/env python3
"""Pal Panel B reframed on KS EMT programs (epi / mes / hybrid), IHC as overlay.

Defines per-cell program states from Tan KS Mes−Epi and joint high epi+mes (hybrid).
Runs the same distance-bin hypergeometric engine as clinical subtype Panel B.
Summarizes barycentric "edge" cells (two dominant weights) vs vertex cells.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/24_panelB_emt_hybrid_gse161529.py"
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.enrichment import hypergeometric_enrichment

PAL_FIT = BREAST / "results" / "panel_a_ks_genelist_gse161529"
PANEL_B = BREAST / "results" / "panel_b_gse161529"
OUT = PANEL_B / "emt_hybrid"
FIG = BREAST / "figures"

K = 3
N_BINS = 5
FDR = 0.1
PROGRAM_LEVELS = ("epithelial", "mesenchymal", "hybrid", "intermediate")
ARC_COLORS = {0: "#4C78A8", 1: "#F58518", 2: "#E45756"}
PROGRAM_COLORS = {
    "epithelial": "#4C78A8",
    "mesenchymal": "#E45756",
    "hybrid": "#9C755F",
    "intermediate": "#B0B0B0",
}
SUBTYPE_COLORS = {"ER+": "#4C78A8", "HER2+": "#E45756", "TNBC": "#54A24B"}


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


def classify_program(epi: np.ndarray, mes: np.ndarray, emt: np.ndarray) -> np.ndarray:
    """Mutually exclusive program labels (hybrid takes priority)."""
    epi_q75, mes_q75 = np.quantile(epi, 0.75), np.quantile(mes, 0.75)
    emt_lo, emt_hi = np.quantile(emt, [1 / 3, 2 / 3])
    hybrid = (epi >= epi_q75) & (mes >= mes_q75)
    out = np.full(len(emt), "intermediate", dtype=object)
    out[hybrid] = "hybrid"
    rest = ~hybrid
    out[rest & (emt <= emt_lo)] = "epithelial"
    out[rest & (emt >= emt_hi)] = "mesenchymal"
    return out


def barycentric_edge_label(w: np.ndarray) -> np.ndarray:
    """Vertex vs edge between top-2 archetypes (k=3)."""
    order = np.argsort(-w, axis=1)
    top = order[:, 0]
    second = order[:, 1]
    w_sorted = np.sort(w, axis=1)
    purity = w_sorted[:, -1] - w_sorted[:, -2]
    labels = np.empty(w.shape[0], dtype=object)
    for i in range(w.shape[0]):
        if purity[i] >= 0.35:
            labels[i] = f"vertex_arc{top[i] + 1}"
        else:
            a, b = sorted((int(top[i]), int(second[i])))
            labels[i] = f"edge_arc{a + 1}_arc{b + 1}"
    return labels


def plot_program_enrichment(table: pd.DataFrame, path: Path, title: str):
    fig, axes = plt.subplots(1, len(PROGRAM_LEVELS), figsize=(11.5, 3.5), sharey=True)
    for ax, program in zip(axes, PROGRAM_LEVELS):
        sub = table[table["subtype"] == program]
        for arc in range(K):
            arc_tab = sub[sub["archetype"] == arc].sort_values("bin")
            ax.plot(
                arc_tab["bin"],
                arc_tab["fold_enrichment"],
                marker="o",
                color=ARC_COLORS[arc],
                label=f"Arc {arc + 1}",
            )
            sig = arc_tab[arc_tab["sig_peak_at_bin0"]]
            if len(sig):
                ax.scatter(
                    sig["bin"],
                    sig["fold_enrichment"],
                    s=70,
                    facecolors="none",
                    edgecolors=ARC_COLORS[arc],
                    linewidths=1.6,
                    zorder=3,
                )
        ax.axhline(1.0, color="0.7", lw=1, ls="--")
        ax.set_title(program)
        ax.set_xlabel("distance bin (0 = closest)")
        ax.set_xticks(range(N_BINS))
    axes[0].set_ylabel("fold enrichment")
    axes[-1].legend(loc="upper right", fontsize=7, frameon=False)
    fig.suptitle(title, y=1.06, fontsize=11)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def draw_triangle(ax, arcs):
    pts = list(arcs) + [arcs[0]]
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color="#333333", lw=1.4, zorder=5)
    ax.scatter(
        arcs[:, 0],
        arcs[:, 1],
        s=110,
        c=[ARC_COLORS[i] for i in range(arcs.shape[0])],
        edgecolors="white",
        linewidths=0.8,
        zorder=6,
    )
    pole_names = ["hybrid", "epi", "mes"]
    for i in range(arcs.shape[0]):
        ax.annotate(
            f"{i + 1}\n({pole_names[i]})",
            arcs[i],
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=7,
            fontweight="bold",
        )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    style()

    cells = pd.read_csv(PANEL_B / "pal_cells_named.csv")
    bins = pd.read_csv(PANEL_B / "bins_subtype.csv")
    vertex = pd.read_csv(PANEL_B / "vertex_ks_scores.csv")
    cells = cells.merge(bins, on="barcode", how="inner", validate="one_to_one")

    w = np.load(PAL_FIT / "S_k3_parti.npy")
    arcs = np.load(PAL_FIT / f"archetypes_k{K}_parti.npy")
    scores = np.load(PAL_FIT / "pc_scores.npy")
    if w.shape[0] != len(cells):
        raise ValueError("weight / cell count mismatch")

    epi = cells["epi_score"].values
    mes = cells["mes_score"].values
    emt = cells["emt_score"].values
    program = classify_program(epi, mes, emt)
    cells["program_state"] = program
    cells["bary_edge"] = barycentric_edge_label(w)

    bin_ids = cells[[f"arc{i+1}" for i in range(K)]].values.astype(int)
    print("Program state counts:\n", pd.Series(program).value_counts().to_string())
    print("\nVertex KS names:\n", vertex.to_string(index=False))

    table = hypergeometric_enrichment(
        cells["program_state"], bin_ids, fdr=FDR, subtype_levels=PROGRAM_LEVELS
    )
    table = table.rename(columns={"subtype": "program"})
    table.to_csv(OUT / "enrichment_program_state.csv", index=False)

    hits = table[table["sig_peak_at_bin0"]]
    print("\nSignificant bin-0 peaks (program states):")
    if hits.empty:
        print("  none")
    else:
        print(
            hits[["program", "archetype", "fold_enrichment", "q_value"]].to_string(
                index=False
            )
        )

    # IHC overlay on same bins (re-run for reference in this folder)
    ihc = hypergeometric_enrichment(
        cells["subtype"], bin_ids, fdr=FDR, subtype_levels=("ER+", "HER2+", "TNBC")
    )
    ihc = ihc.rename(columns={"subtype": "clinical_subtype"})
    ihc.to_csv(OUT / "enrichment_clinical_subtype.csv", index=False)

    # Mean scores by distance bin per archetype
    rows = []
    for j in range(K):
        for b in range(N_BINS):
            m = bin_ids[:, j] == b
            rows.append(
                {
                    "archetype": j + 1,
                    "bin": b,
                    "n": int(m.sum()),
                    "epi_mean": float(epi[m].mean()),
                    "mes_mean": float(mes[m].mean()),
                    "emt_mean": float(emt[m].mean()),
                    "hybrid_frac": float((program[m] == "hybrid").mean()),
                }
            )
    by_bin = pd.DataFrame(rows)
    by_bin.to_csv(OUT / "mean_ks_by_arc_distance_bin.csv", index=False)

    edge_summary = (
        cells.groupby(["bary_edge", "program_state"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    edge_summary.to_csv(OUT / "program_by_barycentric_region.csv")
    edge_frac = edge_summary.div(edge_summary.sum(axis=1), axis=0)
    edge_frac.to_csv(OUT / "program_by_barycentric_region_rowfrac.csv")

    hybrid_on_13 = cells["bary_edge"] == "edge_arc1_arc3"
    print(
        f"\nHybrid program on Arc1–Arc3 edge: "
        f"{(cells.loc[hybrid_on_13, 'program_state'] == 'hybrid').mean():.1%} "
        f"(n={hybrid_on_13.sum()})"
    )
    print(
        f"Hybrid at Arc1 vertex: "
        f"{(cells.loc[cells.bary_edge == 'vertex_arc1', 'program_state'] == 'hybrid').mean():.1%}"
    )

    plot_program_enrichment(
        table.rename(columns={"program": "subtype"}),
        FIG / "Figure_emt_panelB_program_enrichment.png",
        "Pal k=3 — KS program enrichment vs distance to vertex",
    )

    rng = np.random.default_rng(0)
    take = rng.choice(len(cells), size=min(10000, len(cells)), replace=False)

    fig, axes = plt.subplots(2, 2, figsize=(11.2, 9.2))
    ax = axes[0, 0]
    for prog in PROGRAM_LEVELS:
        idx = np.where(program == prog)[0]
        if len(idx) == 0:
            continue
        idx_plot = idx if len(idx) <= 4000 else rng.choice(idx, 4000, replace=False)
        ax.scatter(
            scores[idx_plot, 0],
            scores[idx_plot, 1],
            s=3,
            c=PROGRAM_COLORS[prog],
            alpha=0.35,
            linewidths=0,
            label=prog,
        )
    draw_triangle(ax, arcs[:, :2])
    ax.set_title("KS program state (not IHC)")
    ax.legend(fontsize=7, markerscale=2)
    ax.set_aspect("equal", adjustable="datalim")

    ax = axes[0, 1]
    for st in ("ER+", "HER2+", "TNBC"):
        idx = np.where(cells["subtype"].values == st)[0]
        if len(idx) == 0:
            continue
        idx_plot = idx if len(idx) <= 4000 else rng.choice(idx, 4000, replace=False)
        ax.scatter(
            scores[idx_plot, 0],
            scores[idx_plot, 1],
            s=3,
            c=SUBTYPE_COLORS[st],
            alpha=0.35,
            linewidths=0,
            label=st,
        )
    draw_triangle(ax, arcs[:, :2])
    ax.set_title("Clinical subtype (overlay)")
    ax.legend(fontsize=7, markerscale=2)
    ax.set_aspect("equal", adjustable="datalim")

    ax = axes[1, 0]
    sc = ax.scatter(
        epi[take],
        mes[take],
        s=4,
        c=emt[take],
        cmap="coolwarm",
        vmin=-1.5,
        vmax=1.5,
        alpha=0.45,
        linewidths=0,
    )
    ax.axhline(np.quantile(mes, 0.75), color="0.5", ls=":", lw=1)
    ax.axvline(np.quantile(epi, 0.75), color="0.5", ls=":", lw=1)
    ax.set_xlabel("cell epi score (mean KS Epi z)")
    ax.set_ylabel("cell mes score (mean KS Mes z)")
    ax.set_title("Hybrid quadrant: both scores ≥ 75th pctile")
    fig.colorbar(sc, ax=ax, shrink=0.8, label="EMT (mes − epi)")

    ax = axes[1, 1]
    means = (
        cells.groupby("nearest_archetype", observed=True)[["epi_score", "mes_score", "emt_score"]]
        .mean()
        .reindex([1, 2, 3])
    )
    x = np.arange(3)
    width = 0.25
    ax.bar(x - width, means["epi_score"], width, color="#4C78A8", label="Epi")
    ax.bar(x, means["mes_score"], width, color="#E45756", label="Mes")
    ax.bar(x + width, means["emt_score"], width, color="#9C755F", label="EMT (mes−epi)")
    ax.axhline(0, color="0.7", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(["Arc 1\nhybrid pole", "Arc 2\nepi pole", "Arc 3\nmes pole"])
    ax.set_ylabel("mean cell score")
    ax.set_title("Nearest-vertex cells")
    ax.legend(fontsize=7)

    fig.suptitle(
        "Pal GSE161529 k=3 as an EMT simplex (Tan KS); TNBC maps partly to hybrid pole",
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(FIG / "Figure_emt_hybrid_pal_triangle.png")
    plt.close(fig)

    # Barycentric edge vs hybrid
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    order = edge_frac.index.tolist()
    bottom = np.zeros(len(order))
    for prog in PROGRAM_LEVELS:
        if prog not in edge_frac.columns:
            continue
        vals = edge_frac.loc[order, prog].values
        ax.bar(order, vals, bottom=bottom, label=prog, color=PROGRAM_COLORS[prog], width=0.7)
        bottom += vals
    ax.set_ylabel("fraction of cells in region")
    ax.set_title("Program composition: triangle vertices vs edges (barycentric)")
    ax.legend(loc="upper right", fontsize=7)
    plt.xticks(rotation=25, ha="right")
    fig.tight_layout()
    fig.savefig(FIG / "Figure_emt_hybrid_barycentric_regions.png")
    plt.close(fig)

    pole_map = {1: "hybrid", 2: "epithelial", 3: "mesenchymal"}
    summary = {
        "program_definition": {
            "hybrid": "epi and mes cell scores both >= 75th percentile",
            "epithelial": "not hybrid and EMT score <= 33rd percentile",
            "mesenchymal": "not hybrid and EMT score >= 67th percentile",
            "intermediate": "remaining cells",
        },
        "vertex_ks_names": dict(zip(vertex["archetype"], vertex["ks_name"])),
        "program_bin0_peaks": [],
        "clinical_bin0_peaks": [],
        "interpretation": (
            "Arc 2 = epithelial pole, Arc 3 = mesenchymal pole, Arc 1 = hybrid/mixed pole "
            "(high both at vertex); clinical TNBC enriches at Arc 1 but program hybrid is the "
            "primary KS alignment."
        ),
    }
    for _, r in hits.iterrows():
        summary["program_bin0_peaks"].append(
            {
                "program": r["program"],
                "archetype": int(r["archetype"]) + 1,
                "pole": pole_map.get(int(r["archetype"]) + 1),
                "fold_enrichment": float(r["fold_enrichment"]),
                "q_value": float(r["q_value"]),
            }
        )
    ihc_hits = ihc[ihc["sig_peak_at_bin0"]]
    for _, r in ihc_hits.iterrows():
        summary["clinical_bin0_peaks"].append(
            {
                "subtype": r["clinical_subtype"],
                "archetype": int(r["archetype"]) + 1,
                "fold_enrichment": float(r["fold_enrichment"]),
            }
        )
    (OUT / "emt_hybrid_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    cells.to_csv(OUT / "pal_cells_with_program_state.csv", index=False)

    print("\nWrote", OUT)
    print("Figures:", FIG / "Figure_emt_hybrid_pal_triangle.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
