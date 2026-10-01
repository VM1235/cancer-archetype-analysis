#!/usr/bin/env python3
"""Host-swap containment table/figure + Pal Arc-1 (bin-0) patient audit.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/23_host_swap_comparison_and_arc1_audit.py"
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

PANEL_B = BREAST / "results" / "panel_b_gse161529"
OUT_TABLE = BREAST / "results" / "host_swap_comparison"
FIG = BREAST / "figures"

# Precomputed from projection_report.json + DepMap self barycentric on Panel A fit.
ROWS = ("GSE161529 (Pal)", "GSE176078 (Wu)", "GSE173634", "DepMap (63 lines)")
DEPMAP_HOST = (0.1325855181, 0.4553431007, 0.5358793343, 0.5238095238)
PAL_HOST = (0.7575589230, 0.7838341197, 0.9608460208, 0.9047619048)


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


def arc1_patient_audit() -> dict:
    meta = pd.read_csv(PANEL_B / "pal_cells_named.csv")
    bins = pd.read_csv(PANEL_B / "bins_subtype.csv")
    df = meta.merge(bins, on="barcode", how="inner", validate="one_to_one")
    # arc1 column: distance bin to Arc 1 (1-based label); bin 0 = closest to vertex.
    bin0 = df.loc[df["arc1"] == 0].copy()
    n_bin0 = len(bin0)
    by_patient = (
        bin0.groupby(["patient", "subtype"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    by_patient["frac_of_arc1_bin0"] = by_patient["n_cells"] / n_bin0
    by_patient.to_csv(PANEL_B / "pal_arc1_bin0_by_patient.csv", index=False)

    tnbc_bin0 = bin0.loc[bin0["subtype"] == "TNBC"]
    tnbc_by_patient = (
        tnbc_bin0.groupby("patient", observed=True)
        .size()
        .sort_values(ascending=False)
    )
    top_patient = tnbc_by_patient.index[0] if len(tnbc_by_patient) else None
    top_n = int(tnbc_by_patient.iloc[0]) if len(tnbc_by_patient) else 0
    top_frac = float(top_n / len(tnbc_bin0)) if len(tnbc_bin0) else float("nan")

    subtype_counts = bin0["subtype"].value_counts()
    nearest1 = df.loc[df["nearest_archetype"] == 1].copy()
    near_by_patient = (
        nearest1.groupby(["patient", "subtype"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )
    near_by_patient["frac_of_nearest_arc1"] = near_by_patient["n_cells"] / len(nearest1)
    near_by_patient.to_csv(PANEL_B / "pal_nearest_arc1_by_patient.csv", index=False)
    tnbc_near = nearest1.loc[nearest1["subtype"] == "TNBC"]
    top_near = (
        tnbc_near.groupby("patient", observed=True).size().sort_values(ascending=False)
    )

    summary = {
        "definition": "Pal cells with distance bin 0 to Arc 1 (IHC: TNBC peak at bin 0)",
        "n_cells_arc1_bin0": n_bin0,
        "n_patients_arc1_bin0": int(bin0["patient"].nunique()),
        "subtype_counts": subtype_counts.astype(int).to_dict(),
        "tnbc_in_bin0": int(len(tnbc_bin0)),
        "tnbc_frac_of_bin0": float(len(tnbc_bin0) / n_bin0) if n_bin0 else None,
        "n_tnbc_patients_in_bin0": int(tnbc_bin0["patient"].nunique()),
        "top_tnbc_patient_in_bin0": top_patient,
        "top_tnbc_patient_n_cells": top_n,
        "top_tnbc_patient_frac_of_tnbc_bin0": top_frac,
        "nearest_arc1_cells": int(len(nearest1)),
        "nearest_arc1_n_patients": int(nearest1["patient"].nunique()),
        "top_tnbc_patient_nearest_arc1": (
            str(top_near.index[0]) if len(top_near) else None
        ),
        "top_tnbc_patient_frac_nearest_arc1": (
            float(top_near.iloc[0] / len(tnbc_near)) if len(tnbc_near) else None
        ),
    }
    (PANEL_B / "pal_arc1_bin0_audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    return summary


def plot_host_swap():
    OUT_TABLE.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    table = pd.DataFrame(
        {
            "guest": ROWS,
            "frac_inside_depmap_host": DEPMAP_HOST,
            "frac_inside_pal_host": PAL_HOST,
        }
    )
    table["pct_depmap_host"] = (100 * table["frac_inside_depmap_host"]).round(1)
    table["pct_pal_host"] = (100 * table["frac_inside_pal_host"]).round(1)
    table.to_csv(OUT_TABLE / "inside_fractions_k3.csv", index=False)

    x = np.arange(len(ROWS))
    w = 0.36
    style()
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    b1 = ax.bar(x - w / 2, table["pct_depmap_host"], w, label="Host: DepMap KS (k=3)", color="#4C78A8")
    b2 = ax.bar(x + w / 2, table["pct_pal_host"], w, label="Host: GSE161529 Pal (k=3)", color="#F58518")
    ax.set_xticks(x, ROWS, rotation=15, ha="right")
    ax.set_ylabel("% inside simplex (barycentric weights ≥ 0)")
    ax.set_ylim(0, 105)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title(
        "Asymmetric containment: tumor Pal simplex holds guests;\n"
        "DepMap simplex does not hold Pal (MAGIC sc → DepMap PCA; Pal host: mean/SD match)"
    )

    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                h + 1.2,
                f"{h:.1f}%",
                ha="center",
                va="bottom",
                fontsize=7,
            )

    fig.savefig(FIG / "Figure_host_swap_inside_fractions_k3.png")
    fig.savefig(FIG / "Figure_host_swap_inside_fractions_k3.pdf")
    plt.close(fig)

    fig2, ax2 = plt.subplots(figsize=(7.2, 2.4))
    ax2.axis("off")
    cell_text = [
        [f"{a:.1f}%", f"{b:.1f}%"]
        for a, b in zip(table["pct_depmap_host"], table["pct_pal_host"])
    ]
    tbl = ax2.table(
        cellText=cell_text,
        rowLabels=ROWS,
        colLabels=["DepMap host", "Pal host"],
        loc="center",
        cellLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.1, 1.4)
    ax2.set_title("Fraction inside k=3 simplex by guest dataset", pad=12)
    fig2.savefig(FIG / "Figure_host_swap_inside_table_k3.png")
    plt.close(fig2)


def main() -> int:
    summary = arc1_patient_audit()
    print("=== Pal Arc 1 bin-0 patient audit ===")
    print(json.dumps(summary, indent=2))
    plot_host_swap()
    print("Wrote", PANEL_B / "pal_arc1_bin0_by_patient.csv")
    print("Wrote", PANEL_B / "pal_arc1_bin0_audit_summary.json")
    print("Wrote", OUT_TABLE / "inside_fractions_k3.csv")
    print("Wrote", FIG / "Figure_host_swap_inside_fractions_k3.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
