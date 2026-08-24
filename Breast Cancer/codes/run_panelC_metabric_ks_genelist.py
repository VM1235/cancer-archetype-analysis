#!/usr/bin/env python3
"""METABRIC Panel C on KS (Tan 2014) genelist-restricted archetypes.

Reads prepared METABRIC tables from panel_c_metabric_ks_genelist/ (inputs only).
Reuses k=3 archetypes from panel_a_ks_genelist/ without refitting Panel A.
"""

from pathlib import Path
import argparse
import subprocess
import sys

HERE = Path(__file__).resolve().parent
BREAST = HERE.parent
ROOT = BREAST.parent
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.utils import shuffle

from src.combat import combat as combat_py
from src.gene_lists import rename_tumor_genes_to_cell_line
from src.io import load_expression_csv
from src.pca import align_pca_signs, fit_pca, inverse_transform_scores

CL_MATRIX = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
SCORES_PATH = BREAST / "results" / "panel_a_ks_genelist" / "pc_scores.npy"
N_PCS_PATH = BREAST / "results" / "panel_a_ks_genelist" / "n_pcs.txt"
SUGGESTED = BREAST / "results" / "panel_a_ks_genelist" / "suggested_k.txt"
TUMOR_EXPR = BREAST / "results" / "panel_c_metabric_ks_genelist" / "metabric_log_shared_genes.csv"
TUMOR_META = BREAST / "results" / "panel_c_metabric_ks_genelist" / "metabric_clinical_metadata.csv"
OUT = BREAST / "results" / "panel_c_metabric_ks_genelist"
FIG = BREAST / "figures"
COMBAT_R = HERE / "combat_cellline_tumor.R"

N_PCS_VAR = 20
SEED = 0

ER_COL = "ER_STATUS"
HER2_COL = "HER2_STATUS"
CLAUDIN_COL = "CLAUDIN_SUBTYPE"

IHC_COLOR = {
    "ER+/HER2-": "#3B6FA0",
    "ER+/HER2+": "#D989B5",
    "ER-/HER2+": "#E07A3D",
    "ER-/HER2-": "#54A24B",
}
CLAUDIN_COLOR = {
    "LumA": "#4C78A8",
    "LumB": "#F58518",
    "Her2": "#E45756",
    "Basal": "#72B7B2",
    "claudin-low": "#B279A2",
    "Normal": "#9D9D9D",
    "NC": "#CCCCCC",
}
ARC_COLORS = ["#4C78A8", "#F58518", "#E45756", "#72B7B2"]


def cumulative(values):
    return np.cumsum(np.asarray(values, dtype=float))


def tumor_variance_on_pcs(tumor_X, pca, total_var):
    scores = pca.transform(tumor_X)
    return 100.0 * scores.var(axis=0, ddof=1) / total_var


def ihc_group(meta):
    ok = meta[ER_COL].isin(["Positive", "Negative"]) & meta[HER2_COL].isin(
        ["Positive", "Negative"]
    )
    group = pd.Series(pd.NA, index=meta.index, dtype="object")
    group.loc[ok] = (
        "ER"
        + meta.loc[ok, ER_COL].map({"Positive": "+", "Negative": "-"})
        + "/HER2"
        + meta.loc[ok, HER2_COL].map({"Positive": "+", "Negative": "-"})
    )
    return ok, group


def barycentric_weights(points, vertices):
    k = vertices.shape[0]
    a = np.vstack([vertices.T, np.ones((1, k))])
    b = np.vstack([points.T, np.ones((1, points.shape[0]))])
    w, *_ = np.linalg.lstsq(a, b, rcond=None)
    return w.T


def breakdown_tables(table, group_col, label, k):
    """Nearest-archetype counts and mean distance to nearest archetype by group."""
    sub = table.loc[table[group_col].notna() & (table[group_col].astype(str) != "")]
    ct = pd.crosstab(sub[group_col], sub["nearest_archetype"], margins=True)
    ct.to_csv(OUT / f"{label}_by_nearest_archetype_counts.csv")
    nearest_dist = sub[[f"dist_arc{i+1}" for i in range(k)]].to_numpy()
    nearest_idx = sub["nearest_archetype"].to_numpy().astype(int) - 1
    dist_to_nearest = nearest_dist[np.arange(len(sub)), nearest_idx]
    mean_df = (
        pd.DataFrame({"group": sub[group_col], "dist_to_nearest": dist_to_nearest})
        .groupby("group", dropna=False)["dist_to_nearest"]
        .agg(["count", "mean", "std"])
        .round(4)
    )
    mean_df.to_csv(OUT / f"{label}_mean_distance_to_nearest_archetype.csv")
    return ct, mean_df


def run_combat(merged, batch, force, prefer_r=True):
    merged_path = OUT / "merged_uncorrected_shared_genes.csv"
    batch_path = OUT / "combat_batch.csv"
    out_path = OUT / "CCLE_METABRIC_COMBAT.csv"
    if out_path.exists() and not force:
        print("Using existing", out_path)
        return pd.read_csv(out_path, index_col=0)

    merged.to_csv(merged_path)
    pd.DataFrame({"sample": merged.columns.astype(str), "batch": batch}).to_csv(
        batch_path, index=False
    )
    rlib = BREAST / "rlib"
    sva_ok = (rlib / "sva").exists()
    if prefer_r and sva_ok:
        print("ComBat via sva::ComBat (Groves bc2: mod=~1, ref.batch=cell_line)")
        cmd = [
            "Rscript",
            str(COMBAT_R),
            str(merged_path),
            str(batch_path),
            str(out_path),
            "cell_line",
            str(rlib),
        ]
        subprocess.check_call(cmd)
        return pd.read_csv(out_path, index_col=0)

    print("sva not available; parametric ComBat in Python (same model: ~1, ref=cell_line)")
    adj = combat_py(merged.values, np.asarray(batch), ref_batch="cell_line")
    out = pd.DataFrame(adj, index=merged.index, columns=merged.columns)
    out.to_csv(out_path)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-shuffle", type=int, default=20)
    parser.add_argument("--force-combat", action="store_true")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    n_pcs = int(N_PCS_PATH.read_text().strip().splitlines()[0])
    k_star = int(SUGGESTED.read_text().strip().splitlines()[0])
    ARCS_PATH = BREAST / "results" / "panel_a_ks_genelist" / f"archetypes_k{k_star}_parti.npy"
    expr_cl = load_expression_csv(CL_MATRIX)
    saved = np.load(SCORES_PATH)
    arcs = np.load(ARCS_PATH)
    n_vol = arcs.shape[1]
    k = arcs.shape[0]
    print(f"Cell-line matrix: {expr_cl.shape[0]} genes × {expr_cl.shape[1]} lines")
    print(f"Panel A: {n_pcs} PCs; archetypes k={k} in {n_vol}-D")
    print("PCHA is NOT refit on METABRIC.")

    pca_cl, scores_cl = fit_pca(expr_cl.T.values, n_components=n_pcs)
    pca_cl, scores_cl, _ = align_pca_signs(scores_cl, saved, pca_cl)
    print(f"Rebuilt cell-line PCA; max |score diff| vs saved = {np.max(np.abs(scores_cl - saved)):.4g}")

    arcs_full = np.zeros((k, n_pcs))
    arcs_full[:, :n_vol] = arcs
    gene_arcs = pd.DataFrame(
        inverse_transform_scores(pca_cl, arcs_full).T,
        index=expr_cl.index.astype(str),
        columns=[f"arc{i+1}" for i in range(k)],
    )

    tumors = pd.read_csv(TUMOR_EXPR, index_col=0)
    tumors.index = tumors.index.astype(str)
    tumors.columns = tumors.columns.astype(str)
    rename, unmatched_cl = rename_tumor_genes_to_cell_line(tumors.index, expr_cl.index)
    if rename:
        tumors = tumors.rename(index=rename)
        print("Renamed tumor genes to cell-line symbols:", rename)
    if unmatched_cl:
        print("Cell-line restricted genes still absent from tumors:", unmatched_cl)

    meta = pd.read_csv(TUMOR_META, index_col=0)
    meta.index = meta.index.astype(str)
    tumors = tumors.loc[:, [c for c in tumors.columns if c in meta.index]]
    meta = meta.reindex(tumors.columns)

    shared = expr_cl.index.astype(str).intersection(tumors.index)
    print(f"Shared genes for ComBat merge: {len(shared)}")
    cl_s = expr_cl.loc[shared]
    cl_s.columns = cl_s.columns.astype(str)
    tu_s = tumors.loc[shared]
    merged = pd.concat([cl_s, tu_s], axis=1)
    batch = (["cell_line"] * cl_s.shape[1]) + (["tumor"] * tu_s.shape[1])
    print(f"Uncorrected merge: {merged.shape[0]} genes × {merged.shape[1]} samples")

    combined = run_combat(merged, batch, force=args.force_combat)
    combined.columns = combined.columns.astype(str)
    combined = combined.loc[shared]
    gene_arcs = gene_arcs.loc[shared]
    is_cell = np.array([True] * cl_s.shape[1] + [False] * tu_s.shape[1])
    is_tumor = ~is_cell
    print(f"ComBat matrix: {combined.shape[0]} genes × {combined.shape[1]} samples")

    n_comp = min(N_PCS_VAR, combined.shape[1] - 1, combined.shape[0])
    pca_all, scores_all = fit_pca(combined.T.values, n_components=n_comp)
    arc_scores = pca_all.transform(gene_arcs.T.values)
    print("Combined PCA fit; archetypes projected (not refit).")

    tumor_X = combined.loc[:, is_tumor].T.values
    tumor_total_var = float(tumor_X.var(axis=0, ddof=1).sum())
    pca_tumor, _ = fit_pca(tumor_X, n_components=n_comp)
    ev_tumor_only = 100.0 * pca_tumor.explained_variance_ratio_
    ev_combined = tumor_variance_on_pcs(tumor_X, pca_all, tumor_total_var)

    score_cols = [f"PC{i+1}" for i in range(n_comp)]
    pd.DataFrame(scores_all, index=combined.columns, columns=score_cols).to_csv(
        OUT / "combined_pca_scores.csv"
    )
    pd.DataFrame(arc_scores, index=gene_arcs.columns, columns=score_cols).to_csv(
        OUT / "archetypes_in_combined_pca.csv"
    )
    gene_arcs.to_csv(OUT / "archetypes_gene_space_shared.csv")
    pd.Series(
        np.where(is_cell, "cell_line", "tumor"), index=combined.columns, name="source"
    ).to_csv(OUT / "sample_sources.csv")

    t3 = scores_all[is_tumor][:, :n_vol]
    a3 = arc_scores[:, :n_vol]
    dist = np.sqrt(((t3[:, None, :] - a3[None, :, :]) ** 2).sum(axis=2))
    nearest = dist.argmin(axis=1) + 1
    weights = barycentric_weights(t3, a3)
    inside = (weights >= -1e-6).all(axis=1)
    frac_inside = float(inside.sum()) / inside.size
    print(f"\nTumors inside simplex (all barycentric weights >= 0): {int(inside.sum())}/{inside.size} ({100*frac_inside:.1f}%)")
    print(f"Tumors outside simplex: {int((~inside).sum())}/{inside.size} ({100*(1-frac_inside):.1f}%)")

    ok, group = ihc_group(meta)
    tumor_ids = combined.columns[is_tumor]
    table = meta.copy()
    table["IHC_group"] = group
    for i, c in enumerate(score_cols[:8]):
        table[c] = scores_all[is_tumor, i]
    dist_df = pd.DataFrame(
        dist, index=tumor_ids, columns=[f"dist_arc{i+1}" for i in range(k)]
    )
    dist_df["nearest_archetype"] = nearest
    dist_df["inside_simplex"] = inside
    for i in range(k):
        dist_df[f"w_arc{i+1}"] = weights[:, i]
    table = table.join(dist_df)
    table.to_csv(OUT / "metabric_sample_subtype_archetype.csv")

    # CLAUDIN_SUBTYPE breakdown
    claudin_ct, claudin_mean = breakdown_tables(table, CLAUDIN_COL, "claudin_subtype", k)
    print("\n=== Nearest archetype by CLAUDIN_SUBTYPE ===")
    print(claudin_ct.to_string())
    print("\nMean distance to nearest archetype by CLAUDIN_SUBTYPE:")
    print(claudin_mean.to_string())

    # ER/HER2 breakdown
    ihc = table.loc[ok].copy()
    ihc_ct, ihc_mean = breakdown_tables(ihc, "IHC_group", "ihc", k)
    print(f"\n=== Nearest archetype by ER_STATUS/HER2_STATUS (n={len(ihc)}) ===")
    print(ihc_ct.to_string())
    print("\nMean distance to nearest archetype by IHC group:")
    print(ihc_mean.to_string())

    inside_ct = pd.crosstab(
        table[CLAUDIN_COL].replace("", pd.NA),
        table["inside_simplex"],
        margins=True,
    )
    inside_ct.to_csv(OUT / "claudin_inside_simplex_counts.csv")
    ihc_inside = pd.crosstab(ihc["IHC_group"], ihc["inside_simplex"], margins=True)
    ihc_inside.to_csv(OUT / "ihc_inside_simplex_counts.csv")

    n_shuffle = int(args.n_shuffle)
    ev_null = None
    ratio = ev_combined / np.maximum(ev_tumor_only, 1e-12)
    cum_t = cumulative(ev_tumor_only)
    cum_c = cumulative(ev_combined)
    ratio_cum = cum_c / np.maximum(cum_t, 1e-12)
    if n_shuffle > 0:
        print(f"\nVariance null: {n_shuffle} gene-row shuffles of the ComBat matrix")
        rng = np.random.default_rng(SEED)
        rows = []
        for i in range(n_shuffle):
            ran = shuffle(combined, random_state=int(rng.integers(0, 2**31 - 1)))
            ran.index = combined.index
            pca_rand, _ = fit_pca(ran.T.values, n_components=n_comp)
            rows.append(tumor_variance_on_pcs(tumor_X, pca_rand, tumor_total_var))
            print(f"  shuffle {i+1}/{n_shuffle}", flush=True)
        ev_null = np.vstack(rows)
        np.save(OUT / "ev_null_shuffles.npy", ev_null)

    pd.DataFrame(
        {
            "n_components": np.arange(1, n_comp + 1),
            "ev_tumor_only": ev_tumor_only,
            "ev_combined_on_tumors": ev_combined,
            "cum_tumor_only": cum_t,
            "cum_combined_on_tumors": cum_c,
            "ratio_combined_over_tumor": ratio_cum,
        }
    ).to_csv(OUT / "variance_explained_combined_pca.csv", index=False)

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.3))
    ax = axes[0]
    ax.scatter(
        scores_all[is_cell, 0],
        scores_all[is_cell, 1],
        s=22,
        c="#B8B8B8",
        zorder=1,
        label=f"cell lines (n={int(is_cell.sum())})",
    )
    unlabeled = table[table[CLAUDIN_COL].isna() | (table[CLAUDIN_COL].astype(str) == "")]
    if len(unlabeled):
        ax.scatter(
            unlabeled["PC1"],
            unlabeled["PC2"],
            s=10,
            c="#DDDDDD",
            alpha=0.55,
            zorder=2,
            label=f"tumor, CLAUDIN unknown (n={len(unlabeled)})",
        )
    for gname, color in CLAUDIN_COLOR.items():
        sub = table[table[CLAUDIN_COL] == gname]
        if len(sub):
            ax.scatter(
                sub["PC1"],
                sub["PC2"],
                s=12,
                c=color,
                alpha=0.8,
                zorder=3,
                label=f"{gname} (n={len(sub)})",
            )
    for i in range(k):
        for j in range(i + 1, k):
            ax.plot(
                [arc_scores[i, 0], arc_scores[j, 0]],
                [arc_scores[i, 1], arc_scores[j, 1]],
                color="#555555",
                lw=0.9,
                zorder=4,
            )
        ax.scatter(
            arc_scores[i, 0],
            arc_scores[i, 1],
            s=90,
            c=ARC_COLORS[i],
            edgecolors="k",
            linewidths=0.4,
            zorder=5,
        )
        ax.annotate(
            f"Arc {i+1}",
            (arc_scores[i, 0], arc_scores[i, 1]),
            textcoords="offset points",
            xytext=(5, 5),
            fontsize=8,
            fontweight="bold",
        )
    ax.set_xlabel("PC1 (combined, ComBat)")
    ax.set_ylabel("PC2 (combined, ComBat)")
    ax.set_title("METABRIC tumors in cell-line archetype space")
    ax.legend(frameon=False, fontsize=7, loc="best")

    ax = axes[1]
    xs = np.arange(1, n_comp + 1)
    if ev_null is not None:
        for row in ev_null:
            ax.plot(xs, cumulative(row), color="0.75", lw=0.8, alpha=0.8)
    ax.plot(xs, cum_t, marker="o", color="#4C78A8", label="tumor-only PCA (ceiling)")
    ax.plot(xs, cum_c, marker="o", color="#F58518", label="combined PCA on tumors")
    ax.axvline(k, color="0.5", ls="--", lw=1)
    ax.set_xlabel("number of components")
    ax.set_ylabel("cumulative % tumor variance")
    pct = 100.0 * float(ratio_cum[k - 1])
    ax.set_title(f"{k} combined PCs = {pct:.0f}% of tumor-only ceiling")
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle("Figure 1C  —  METABRIC-BRCA, KS-restricted archetypes", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "Figure_1C_metabric_ks_genelist.png", dpi=200)
    plt.close(fig)
    print("Wrote", FIG / "Figure_1C_metabric_ks_genelist.png")

    print("\n=== Summary ===")
    print(f"cell lines={int(is_cell.sum())}  tumors={int(is_tumor.sum())}  genes={combined.shape[0]}")
    print(f"inside simplex={int(inside.sum())}/{inside.size} ({100*frac_inside:.1f}%)")
    print(f"outside simplex={int((~inside).sum())}/{inside.size} ({100*(1-frac_inside):.1f}%)")
    print("Wrote outputs under", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
