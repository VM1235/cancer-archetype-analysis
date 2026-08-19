#!/usr/bin/env python3
"""Panel C for hypothesis-driven k=2: map the two cell-line poles into tumors.

Do not refit PCHA on tumors. Map k=2 archetypes to gene space via the saved
12-PC PCA, reuse existing ComBat (cell line vs tumor), fit a NEW combined PCA,
project the two poles. Inside the simplex = between the two ends.
"""

from pathlib import Path
import argparse
import subprocess
import sys

HERE = Path(__file__).resolve().parent
GBM = HERE.parent
ROOT = GBM.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.combat import combat as combat_py
from src.io import load_expression_csv
from src.pca import align_pca_signs, fit_pca, inverse_transform_scores
from src.paths import BREAST
from k2_utils import (
    SUBTYPE_COLOR,
    corr_table,
    line_barycentric,
    signature_scores,
    style_mpl,
)

CL_MATRIX = GBM / "data" / "processed" / "input_panelA_glioblastoma_ccle_logtpm_filtered.csv"
SCORES_PATH = GBM / "results" / "panel_a" / "pc_scores.npy"
N_PCS_PATH = GBM / "results" / "panel_a" / "n_pcs.txt"
ARCS_PATH = GBM / "results" / "panel_a_k2" / "archetypes_k2_parti.npy"
POLES = GBM / "results" / "panel_a_k2" / "pole_assignment.txt"
TCGA_EXPR = GBM / "results" / "panel_c_tcga" / "tcga_primary_log_shared_genes.csv"
TCGA_META = GBM / "results" / "panel_c_tcga" / "tcga_primary_metadata.csv"
COMBAT_EXISTING = GBM / "results" / "panel_c_tcga" / "CCLE_TCGA_COMBAT.csv"
OUT = GBM / "results" / "panel_c_k2"
FIG = GBM / "figures"
COMBAT_R = HERE / "combat_cellline_tumor.R"

K = 2
N_PCS_VAR = 20
SEED = 0
SUB_COL = "GeneExp_Subtype"


def read_poles():
    mes_idx, pn_idx = 1, 0
    for line in POLES.read_text().splitlines():
        if line.startswith("mes_idx="):
            mes_idx = int(line.split("=", 1)[1])
        if line.startswith("pn_idx="):
            pn_idx = int(line.split("=", 1)[1])
    return mes_idx, pn_idx


def run_combat(merged, batch, force, prefer_r=True):
    if COMBAT_EXISTING.exists() and not force:
        print("Using existing ComBat matrix", COMBAT_EXISTING)
        return pd.read_csv(COMBAT_EXISTING, index_col=0)

    merged_path = OUT / "merged_uncorrected_shared_genes.csv"
    batch_path = OUT / "combat_batch.csv"
    out_path = OUT / "CCLE_TCGA_COMBAT.csv"
    print("--force-combat: writing ComBat under panel_c_k2 (not overwriting panel_c_tcga)")
    merged.to_csv(merged_path)
    pd.DataFrame({"sample": merged.columns.astype(str), "batch": batch}).to_csv(
        batch_path, index=False
    )
    rlib_candidates = [GBM / "rlib", BREAST / "rlib"]
    rlib = next((p for p in rlib_candidates if (p / "sva").exists()), rlib_candidates[0])
    sva_ok = (rlib / "sva").exists()
    if prefer_r and sva_ok:
        print("ComBat via sva::ComBat (mod=~1, ref.batch=cell_line)")
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

    print("sva not available; parametric ComBat in Python")
    adj = combat_py(merged.values, np.asarray(batch), ref_batch="cell_line")
    out = pd.DataFrame(adj, index=merged.index, columns=merged.columns)
    out.to_csv(out_path)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--force-combat", action="store_true")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    n_pcs = int(N_PCS_PATH.read_text().strip().splitlines()[0])
    expr_cl = load_expression_csv(CL_MATRIX)
    saved = np.load(SCORES_PATH)
    arcs = np.load(ARCS_PATH)
    n_vol = arcs.shape[1]
    mes_idx, pn_idx = read_poles()
    print(f"Cell-line matrix: {expr_cl.shape[0]} genes × {expr_cl.shape[1]} lines")
    print(f"Saved PCA: {n_pcs} PCs; k=2 archetypes in {n_vol}-D (hypothesis-driven)")
    print("PCHA is NOT refit on TCGA.")
    print(f"Poles: Arc {mes_idx + 1}=MES, Arc {pn_idx + 1}=PN")

    pca_cl, scores_cl = fit_pca(expr_cl.T.values, n_components=n_pcs)
    pca_cl, scores_cl, _ = align_pca_signs(scores_cl, saved, pca_cl)
    print(f"Rebuilt cell-line PCA; max |score diff| vs saved = {np.max(np.abs(scores_cl - saved)):.4g}")

    arcs_full = np.zeros((K, n_pcs))
    arcs_full[:, :n_vol] = arcs
    gene_arcs = pd.DataFrame(
        inverse_transform_scores(pca_cl, arcs_full).T,
        index=expr_cl.index.astype(str),
        columns=[f"arc{i + 1}" for i in range(K)],
    )

    tumors = pd.read_csv(TCGA_EXPR, index_col=0)
    tumors.index = tumors.index.astype(str)
    tumors.columns = tumors.columns.astype(str)
    meta = pd.read_csv(TCGA_META, index_col=0)
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
    print("Combined PCA fit; two archetypes projected (not refit).")

    score_cols = [f"PC{i + 1}" for i in range(n_comp)]
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

    # Between the two poles: barycentric on the line in combined PCA.
    t_pts = scores_all[is_tumor]
    w_line = line_barycentric(t_pts, arc_scores[pn_idx], arc_scores[mes_idx])
    w_pn, w_mes = w_line[:, 0], w_line[:, 1]
    inside = (w_line >= -1e-6).all(axis=1)
    print(
        f"Tumors between the two poles (barycentric ≥ 0 on the PN–MES line, {n_comp}-D): "
        f"{int(inside.sum())}/{inside.size}"
    )

    # Also report 1-D Groves-style weights on combined PC1 only (k−1 = 1).
    w_pc1 = line_barycentric(
        scores_all[is_tumor][:, :1], arc_scores[pn_idx, :1], arc_scores[mes_idx, :1]
    )
    inside_pc1 = (w_pc1 >= -1e-6).all(axis=1)
    print(
        f"Tumors between poles on combined PC1 only: "
        f"{int(inside_pc1.sum())}/{inside_pc1.size}"
    )

    tumor_ids = combined.columns[is_tumor]
    table = meta.copy()
    if SUB_COL in table.columns:
        table[SUB_COL] = table[SUB_COL].astype(str)
    for i, c in enumerate(score_cols[:8]):
        table[c] = scores_all[is_tumor, i]
    table["w_PN"] = w_pn
    table["w_MES"] = w_mes
    table["inside_simplex"] = inside
    table["w_PN_pc1"] = w_pc1[:, 0]
    table["w_MES_pc1"] = w_pc1[:, 1]
    table["inside_simplex_pc1"] = inside_pc1
    dist = np.sqrt(((t_pts[:, None, :] - arc_scores[None, :, :]) ** 2).sum(axis=2))
    table["nearest_archetype"] = dist.argmin(axis=1) + 1
    table["dist_PN_pole"] = dist[:, pn_idx]
    table["dist_MES_pole"] = dist[:, mes_idx]
    table.to_csv(OUT / "tcga_sample_subtype_archetype.csv")

    labeled = (
        table[table[SUB_COL].isin(SUBTYPE_COLOR)] if SUB_COL in table.columns else table.iloc[0:0]
    )
    print(f"\nPrimary tumors with Verhaak GeneExp_Subtype: {len(labeled)}/{len(table)}")
    if len(labeled):
        print(labeled[SUB_COL].value_counts().to_string())
        ct = pd.crosstab(labeled[SUB_COL], labeled["nearest_archetype"], margins=True)
        ct.to_csv(OUT / "subtype_by_nearest_archetype_counts.csv")
        print("\nNearest pole (1=Arc1, 2=Arc2), rows=Verhaak:")
        print(ct.to_string())
        inside_ct = pd.crosstab(labeled[SUB_COL], labeled["inside_simplex"], margins=True)
        inside_ct.to_csv(OUT / "subtype_inside_simplex_counts.csv")
        print("\nBetween the two poles by Verhaak:")
        print(inside_ct.to_string())
        pn_mes = labeled[labeled[SUB_COL].isin(["Proneural", "Mesenchymal"])]
        print(f"\nPN vs MES tumors only: n={len(pn_mes)}")
        if len(pn_mes):
            print(
                pn_mes.groupby(SUB_COL)["w_MES"]
                .agg(["count", "mean", "median"])
                .to_string()
            )
            inside_pm = pd.crosstab(pn_mes[SUB_COL], pn_mes["inside_simplex"], margins=True)
            print("Between poles, PN vs MES:")
            print(inside_pm.to_string())

    # Signatures on tumors (z within TCGA matrix) and on ComBat tumors.
    tu_sig, _ = signature_scores(tu_s)
    tu_sig.to_csv(OUT / "signature_scores_tumors_uncorrected.csv")
    cb_tumors = combined.loc[:, is_tumor]
    cb_sig, _ = signature_scores(cb_tumors)
    cb_sig.to_csv(OUT / "signature_scores_tumors_combat.csv")
    w_df = pd.DataFrame({"w_MES": w_mes, "w_PN": w_pn}, index=tumor_ids.astype(str))
    tu_sig = tu_sig.reindex(w_df.index)
    cb_sig = cb_sig.reindex(w_df.index)
    cor_u = pd.concat(
        [corr_table(w_df["w_MES"], tu_sig, "w_MES"), corr_table(w_df["w_PN"], tu_sig, "w_PN")]
    )
    cor_u.insert(0, "matrix", "tcga_uncorrected")
    cor_c = pd.concat(
        [corr_table(w_df["w_MES"], cb_sig, "w_MES"), corr_table(w_df["w_PN"], cb_sig, "w_PN")]
    )
    cor_c.insert(0, "matrix", "combat_tumors")
    cor = pd.concat([cor_u, cor_c], ignore_index=True)
    cor.to_csv(OUT / "signature_vs_mixture_tumors.csv", index=False)
    print("\nTumor signature vs PCHA mixture (Pearson r):")
    print(cor[["matrix", "weight", "signature", "pearson_r", "pearson_p"]].to_string(index=False))

    style_mpl(plt)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.2))
    ax = axes[0]
    ax.scatter(
        scores_all[is_cell, 0],
        scores_all[is_cell, 1],
        s=22,
        c="#B8B8B8",
        zorder=1,
        label=f"cell lines (n={int(is_cell.sum())})",
    )
    unlabeled = table[~table.index.isin(labeled.index)] if len(labeled) else table
    if len(unlabeled):
        ax.scatter(
            unlabeled["PC1"],
            unlabeled["PC2"],
            s=10,
            c="#DDDDDD",
            alpha=0.55,
            zorder=2,
            label=f"tumor, subtype missing (n={len(unlabeled)})",
        )
    for gname, color in SUBTYPE_COLOR.items():
        sub = labeled[labeled[SUB_COL] == gname] if len(labeled) else labeled
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
    ax.plot(
        [arc_scores[pn_idx, 0], arc_scores[mes_idx, 0]],
        [arc_scores[pn_idx, 1], arc_scores[mes_idx, 1]],
        color="#555555",
        lw=1.2,
        zorder=4,
    )
    for idx, name in ((pn_idx, "PN pole"), (mes_idx, "MES pole")):
        color = SUBTYPE_COLOR["Proneural" if name.startswith("PN") else "Mesenchymal"]
        ax.scatter(
            arc_scores[idx, 0],
            arc_scores[idx, 1],
            s=110,
            c=color,
            edgecolors="k",
            linewidths=0.4,
            zorder=5,
        )
        ax.annotate(
            f"Arc {idx + 1}\n{name}",
            (arc_scores[idx, 0], arc_scores[idx, 1]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=8,
            fontweight="bold",
        )
    ax.set_xlabel("PC1 (combined, ComBat)")
    ax.set_ylabel("PC2 (combined, ComBat)")
    ax.set_title("TCGA tumors on the cell-line PN–MES axis")
    ax.legend(frameon=False, fontsize=7, loc="best")

    ax = axes[1]
    order = ["Proneural", "Neural", "Classical", "Mesenchymal"]
    data, colors, tick = [], [], []
    for gname in order:
        sub = labeled[labeled[SUB_COL] == gname] if len(labeled) else labeled
        if len(sub) == 0:
            continue
        data.append(sub["w_MES"].to_numpy(dtype=float))
        colors.append(SUBTYPE_COLOR[gname])
        tick.append(f"{gname}\n(n={len(sub)})")
    if data:
        bp = ax.boxplot(data, patch_artist=True, widths=0.55)
        ax.set_xticklabels(tick)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.45)
        rng = np.random.default_rng(SEED)
        for i, vals in enumerate(data, start=1):
            x = i + rng.uniform(-0.12, 0.12, size=len(vals))
            ax.scatter(x, vals, s=8, c="0.25", alpha=0.5, zorder=3)
    ax.axhline(0.0, color="0.7", lw=0.8, ls="--")
    ax.axhline(1.0, color="0.7", lw=0.8, ls="--")
    ax.set_ylabel("mixture weight toward MES pole")
    ax.set_title("Tumor position on the PN (0) – MES (1) line")
    ax.set_ylim(-0.35, 1.35)

    fig.suptitle(
        "Figure 1C (k=2)  —  TCGA-GBM on hypothesis-driven PN–MES axis",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(FIG / "Figure_1C_gbm_k2.png", dpi=200)
    plt.close(fig)
    print("Wrote", FIG / "Figure_1C_gbm_k2.png")
    print(f"cell lines={int(is_cell.sum())}  tumors={int(is_tumor.sum())}  genes={combined.shape[0]}")
    print(f"inside simplex (between poles)={int(inside.sum())}/{inside.size}")
    print("Wrote outputs under", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
