#!/usr/bin/env python3
"""Name Pal GSE161529 k=3 vertices: (1) ER/HER2/TNBC Panel B, (2) Tan KS Epi/Mes.

Uses frozen Pal PCA + archetypes from 17.9 (not refit). Distance-bin
hypergeometric enrichment is the same engine as DepMap Panel B. KS scores
are mean z-scored epithelial minus mesenchymal genes (Tan 2014 Table S1B).

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/21_panelB_gse161529.py"
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

from src.enrichment import distance_bins, hypergeometric_enrichment
from src.pca import align_pca_signs, fit_pca, inverse_transform_scores

PAL_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse161529_magic" / "ks_genes_magic_imputed.npz"
)
WU_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse176078_magic" / "ks_genes_magic_imputed.npz"
)
PAL_META = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse161529_magic" / "sc_cells_in_panelA_space.csv"
)
WU_PROJ = BREAST / "results" / "gse176078_on_gse161529" / "gse176078_in_gse161529_k3.csv"
PAL_FIT = BREAST / "results" / "panel_a_ks_genelist_gse161529"
GENELIST = BREAST / "data" / "processed" / "ks_cellline_signature_tan2014.csv"
OUT = BREAST / "results" / "panel_b_gse161529"
FIG = BREAST / "figures"

K = 3
N_BINS = 5
FDR = 0.1
KEEP = ("ER+", "HER2+", "TNBC")
ARC_COLORS = {0: "#4C78A8", 1: "#F58518", 2: "#E45756"}
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
    for i in range(arcs.shape[0]):
        ax.annotate(
            str(i + 1),
            arcs[i],
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=8,
            fontweight="bold",
        )


def plot_enrichment(table, path, title):
    fig, axes = plt.subplots(1, len(KEEP), figsize=(10.2, 3.4), sharey=True)
    for ax, subtype in zip(axes, KEEP):
        sub = table[table["subtype"] == subtype]
        for arc in sorted(sub["archetype"].unique()):
            arc_tab = sub[sub["archetype"] == arc].sort_values("bin")
            ax.plot(
                arc_tab["bin"],
                arc_tab["fold_enrichment"],
                marker="o",
                color=ARC_COLORS[int(arc)],
                label=f"Arc {int(arc) + 1}",
            )
            sig = arc_tab[arc_tab["sig_peak_at_bin0"]]
            if len(sig):
                ax.scatter(
                    sig["bin"],
                    sig["fold_enrichment"],
                    s=70,
                    facecolors="none",
                    edgecolors=ARC_COLORS[int(arc)],
                    linewidths=1.6,
                    zorder=3,
                )
        ax.axhline(1.0, color="0.7", lw=1, ls="--")
        ax.set_title(subtype)
        ax.set_xlabel("distance bin (0 = closest)")
        ax.set_xticks(range(N_BINS))
    axes[0].set_ylabel("fold enrichment")
    axes[-1].legend(loc="upper right", fontsize=7, frameon=False)
    fig.suptitle(title, y=1.06, fontsize=12)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def zscore_genes(expr: pd.DataFrame) -> pd.DataFrame:
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1, ddof=1).replace(0, np.nan).fillna(1.0)
    sd = sd.clip(lower=1e-8)
    return expr.sub(mu, axis=0).div(sd, axis=0)


def ks_cell_scores(expr: pd.DataFrame, epi: list[str], mes: list[str]) -> pd.DataFrame:
    z = zscore_genes(expr)
    epi_g = [g for g in epi if g in z.index]
    mes_g = [g for g in mes if g in z.index]
    epi_s = z.loc[epi_g].mean(axis=0) if epi_g else pd.Series(0.0, index=z.columns)
    mes_s = z.loc[mes_g].mean(axis=0) if mes_g else pd.Series(0.0, index=z.columns)
    out = pd.DataFrame({"epi_score": epi_s, "mes_score": mes_s})
    out["emt_score"] = out["mes_score"] - out["epi_score"]
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    style()

    scores = np.load(PAL_FIT / "pc_scores.npy")
    barcodes = pd.read_csv(PAL_FIT / "pc_scores_barcodes.csv")["barcode"].astype(str)
    arcs = np.load(PAL_FIT / f"archetypes_k{K}_parti.npy")
    n_pcs = int((PAL_FIT / "n_pcs.txt").read_text().strip().splitlines()[0])
    n_vol = arcs.shape[1]
    if len(barcodes) != scores.shape[0]:
        raise ValueError("barcode/score length mismatch")
    print(f"Pal scores {scores.shape[0]} cells, archetypes {arcs.shape}")

    meta = pd.read_csv(PAL_META, usecols=["barcode", "patient", "subtype"])
    meta["barcode"] = meta["barcode"].astype(str)
    meta = meta.drop_duplicates("barcode").set_index("barcode").reindex(barcodes)
    if meta["subtype"].isna().any():
        n_miss = int(meta["subtype"].isna().sum())
        print(f"Warning: {n_miss} barcodes missing subtype; filling NA")
        meta["subtype"] = meta["subtype"].fillna("NA")
        meta["patient"] = meta["patient"].fillna("NA")

    print("Subtype counts:\n", meta["subtype"].value_counts().to_string())

    print("\n=== (1) Clinical subtype Panel B ===")
    scores_u = scores[:, :n_vol]
    bin_ids, distances = distance_bins(scores_u, arcs, n_bins=N_BINS)
    for j in range(arcs.shape[0]):
        sizes = [(bin_ids[:, j] == b).sum() for b in range(N_BINS)]
        print(f"  Arc {j + 1} bin sizes: {sizes}")

    table = hypergeometric_enrichment(
        meta["subtype"], bin_ids, fdr=FDR, subtype_levels=KEEP
    )
    keep_cols = [
        "archetype",
        "subtype",
        "bin",
        "fold_enrichment",
        "p_value",
        "q_value",
        "sig_peak_at_bin0",
        "k_in_bin",
        "bin_size",
        "K_subtype",
        "significant",
        "peak_bin",
    ]
    table = table[keep_cols]
    table.to_csv(OUT / "enrichment_subtype.csv", index=False)
    pd.DataFrame(distances, index=barcodes, columns=[f"arc{i+1}" for i in range(K)]).to_csv(
        OUT / "distances_subtype.csv"
    )
    pd.DataFrame(bin_ids, index=barcodes, columns=[f"arc{i+1}" for i in range(K)]).to_csv(
        OUT / "bins_subtype.csv"
    )

    w = np.load(PAL_FIT / "S_k3_parti.npy")
    nearest = w.argmax(axis=1) + 1
    ct = pd.crosstab(meta["subtype"].values, nearest, rownames=["subtype"], colnames=["nearest_arc"])
    ct.to_csv(OUT / "nearest_archetype_by_subtype_counts.csv")
    (ct.div(ct.sum(axis=1), axis=0)).to_csv(OUT / "nearest_archetype_by_subtype_rowfrac.csv")
    print("\nNearest archetype vs subtype (counts):")
    print(ct.to_string())

    hits = table[table["sig_peak_at_bin0"]]
    print("\nSignificant bin-0 peaks (q < 0.1 and peak at bin 0):")
    if hits.empty:
        print("  none")
    else:
        print(
            hits[["subtype", "archetype", "fold_enrichment", "p_value", "q_value"]].to_string(
                index=False
            )
        )

    names_subtype = {}
    for arc in range(K):
        sub_hits = hits[hits["archetype"] == arc]
        if sub_hits.empty:
            names_subtype[arc + 1] = None
            print(f"  Arc {arc + 1}: no significant subtype bin-0 peak")
        else:
            labs = list(sub_hits.sort_values("fold_enrichment", ascending=False)["subtype"])
            names_subtype[arc + 1] = labs
            print(f"  Arc {arc + 1}: {labs}")

    plot_enrichment(
        table,
        FIG / "Figure_1B_gse161529_subtype.png",
        "Pal GSE161529 k=3 — subtype enrichment vs distance to vertex",
    )
    plot_enrichment(table, OUT / "enrichment_subtype.png", "Pal subtype enrichment")

    print("\n=== (2) Tan KS epithelial / mesenchymal programs ===")
    genelist = pd.read_csv(GENELIST)
    cat = dict(zip(genelist["gene_symbol"].astype(str), genelist["category"].astype(str)))
    zpal = np.load(PAL_MAGIC, allow_pickle=True)
    genes = [str(g) for g in zpal["genes"]]
    expr = pd.DataFrame(zpal["imputed"], index=genes, columns=list(zpal["barcodes"]))
    expr = expr.reindex(columns=list(barcodes))
    epi = [g for g in genes if cat.get(g) == "Epi"]
    mes = [g for g in genes if cat.get(g) == "Mes"]
    print(f"KS genes in Pal MAGIC: Epi={len(epi)}  Mes={len(mes)}")

    ks = ks_cell_scores(expr, epi, mes)
    ks = ks.reindex(barcodes)
    print(
        f"EMT score (Mes−Epi)  mean={ks['emt_score'].mean():.3f}  "
        f"sd={ks['emt_score'].std():.3f}"
    )

    pca, rebuilt = fit_pca(expr.T.values, n_components=n_pcs)
    pca, rebuilt, _ = align_pca_signs(rebuilt, scores, pca)
    print(f"PCA rebuild max |diff| vs saved scores: {np.max(np.abs(rebuilt - scores)):.4g}")
    arcs_full = np.zeros((K, n_pcs))
    arcs_full[:, :n_vol] = arcs
    gene_arcs = inverse_transform_scores(pca, arcs_full)
    gene_arc_df = pd.DataFrame(gene_arcs.T, index=genes, columns=[f"arc{i+1}" for i in range(K)])
    gene_arc_df.to_csv(OUT / "archetypes_gene_space.csv")

    mu = expr.mean(axis=1)
    sd = expr.std(axis=1, ddof=1).clip(lower=1e-8)
    z_arc = gene_arc_df.sub(mu, axis=0).div(sd, axis=0)
    vertex_ks = []
    top_genes = []
    for i in range(K):
        col = f"arc{i+1}"
        epi_v = float(z_arc.loc[epi, col].mean()) if epi else np.nan
        mes_v = float(z_arc.loc[mes, col].mean()) if mes else np.nan
        emt_v = mes_v - epi_v
        if emt_v > 0.2:
            ks_name = "mesenchymal"
        elif emt_v < -0.2:
            ks_name = "epithelial"
        else:
            ks_name = "mixed / high both"
        vertex_ks.append(
            {
                "archetype": i + 1,
                "epi_z": epi_v,
                "mes_z": mes_v,
                "emt_z": emt_v,
                "ks_name": ks_name,
            }
        )
        top = z_arc[col].sort_values(ascending=False).head(12)
        for gene, val in top.items():
            top_genes.append(
                {
                    "archetype": i + 1,
                    "gene": gene,
                    "category": cat.get(gene, ""),
                    "z_vs_pal_mean": float(val),
                }
            )
        print(
            f"  Arc {i + 1} reconstructed: Epi z={epi_v:.3f}  Mes z={mes_v:.3f}  "
            f"EMT={emt_v:.3f}  → {vertex_ks[-1]['ks_name']}"
        )
        print("    top genes:", ", ".join(f"{g}({cat.get(g, '?')})" for g in top.index[:8]))

    vertex_ks_df = pd.DataFrame(vertex_ks)
    vertex_ks_df.to_csv(OUT / "vertex_ks_scores.csv", index=False)
    pd.DataFrame(top_genes).to_csv(OUT / "vertex_top_genes.csv", index=False)

    cell_by_arc = (
        pd.DataFrame(
            {
                "nearest_arc": nearest,
                "epi_score": ks["epi_score"].values,
                "mes_score": ks["mes_score"].values,
                "emt_score": ks["emt_score"].values,
                "subtype": meta["subtype"].values,
            }
        )
        .groupby("nearest_arc")[["epi_score", "mes_score", "emt_score"]]
        .mean()
    )
    cell_by_arc.to_csv(OUT / "ks_scores_by_nearest_archetype.csv")
    print("\nMean cell KS scores by nearest archetype:")
    print(cell_by_arc.to_string())

    bin0_ks = []
    for j in range(K):
        m = bin_ids[:, j] == 0
        bin0_ks.append(
            {
                "archetype": j + 1,
                "n_bin0": int(m.sum()),
                "emt_mean_bin0": float(ks["emt_score"].values[m].mean()),
                "epi_mean_bin0": float(ks["epi_score"].values[m].mean()),
                "mes_mean_bin0": float(ks["mes_score"].values[m].mean()),
            }
        )
    pd.DataFrame(bin0_ks).to_csv(OUT / "ks_scores_bin0.csv", index=False)

    # Wu hold-out: same KS score, Pal nearest-arc from projection
    wu_note = None
    if WU_MAGIC.is_file() and WU_PROJ.is_file():
        zwu = np.load(WU_MAGIC, allow_pickle=True)
        wu_expr = pd.DataFrame(zwu["imputed"], index=list(map(str, zwu["genes"])), columns=list(zwu["barcodes"]))
        wu_ks = ks_cell_scores(wu_expr, epi, mes)
        wu_p = pd.read_csv(WU_PROJ)
        wu_p["barcode"] = wu_p["barcode"].astype(str)
        wu_ks = wu_ks.reindex(wu_p["barcode"])
        wu_p = wu_p.copy()
        wu_p["emt_score"] = wu_ks["emt_score"].values
        wu_sum = wu_p.groupby("nearest_archetype")["emt_score"].agg(["mean", "count"])
        wu_sum.to_csv(OUT / "wu_emt_by_pal_nearest_arc.csv")
        print("\nWu EMT (Mes−Epi) by Pal nearest archetype (hold-out):")
        print(wu_sum.to_string())
        wu_note = wu_sum["mean"].to_dict()

    cells = pd.DataFrame(
        {
            "barcode": barcodes,
            "patient": meta["patient"].values,
            "subtype": meta["subtype"].values,
            "PC1": scores[:, 0],
            "PC2": scores[:, 1],
            "nearest_archetype": nearest,
            "w_arc1": w[:, 0],
            "w_arc2": w[:, 1],
            "w_arc3": w[:, 2],
            "epi_score": ks["epi_score"].values,
            "mes_score": ks["mes_score"].values,
            "emt_score": ks["emt_score"].values,
        }
    )
    cells.to_csv(OUT / "pal_cells_named.csv", index=False)

    rng = np.random.default_rng(0)
    take = rng.choice(len(barcodes), size=min(12000, len(barcodes)), replace=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    ax = axes[0]
    for st, col in SUBTYPE_COLORS.items():
        idx = np.where(meta["subtype"].values == st)[0]
        if len(idx) == 0:
            continue
        idx_plot = idx if len(idx) <= 5000 else rng.choice(idx, 5000, replace=False)
        ax.scatter(
            scores[idx_plot, 0],
            scores[idx_plot, 1],
            s=3,
            c=col,
            alpha=0.25,
            linewidths=0,
            label=st,
        )
    draw_triangle(ax, arcs[:, :2])
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Pal cells by tumor subtype")
    ax.legend(fontsize=7)
    ax.set_aspect("equal", adjustable="datalim")

    ax = axes[1]
    sc = ax.scatter(
        scores[take, 0],
        scores[take, 1],
        s=3,
        c=ks["emt_score"].values[take],
        cmap="coolwarm",
        vmin=-1.5,
        vmax=1.5,
        alpha=0.5,
        linewidths=0,
    )
    draw_triangle(ax, arcs[:, :2])
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Pal cells by KS EMT score (Mes − Epi)")
    fig.colorbar(sc, ax=ax, shrink=0.75, label="EMT")
    ax.set_aspect("equal", adjustable="datalim")
    fig.suptitle("Pal GSE161529 k=3 vertices — subtype vs KS programs", y=1.02)
    fig.savefig(FIG / "Figure_1B_gse161529_named_scatter.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    x = np.arange(K)
    width = 0.35
    ax.bar(x - width / 2, vertex_ks_df["epi_z"], width, color="#4C78A8", label="Epi z")
    ax.bar(x + width / 2, vertex_ks_df["mes_z"], width, color="#E45756", label="Mes z")
    ax.axhline(0, color="0.7", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([f"Arc {i+1}\n({n})" for i, n in zip(range(K), vertex_ks_df["ks_name"])])
    ax.set_ylabel("mean gene z at vertex (vs Pal cells)")
    ax.legend(fontsize=8)
    ax.set_title("Reconstructed Pal vertices: KS Epi vs Mes")
    fig.savefig(FIG / "Figure_1B_gse161529_ks_vertices.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    data, labels = [], []
    for i in range(1, K + 1):
        data.append(ks["emt_score"].values[nearest == i])
        labels.append(f"Arc {i}")
    ax.boxplot(data, tick_labels=labels, showfliers=False)
    ax.axhline(0, color="0.7", lw=1, ls="--")
    ax.set_ylabel("cell EMT score (Mes − Epi)")
    ax.set_title("Pal cells nearest each vertex")
    fig.savefig(FIG / "Figure_1B_gse161529_emt_by_arc.png")
    plt.close(fig)

    summary = {
        "k": K,
        "n_cells": int(len(barcodes)),
        "n_bins": N_BINS,
        "fdr": FDR,
        "subtype_bin0_peaks": {str(k): v for k, v in names_subtype.items()},
        "vertex_ks": vertex_ks,
        "n_epi_genes": len(epi),
        "n_mes_genes": len(mes),
        "wu_emt_by_pal_arc": {str(k): float(v) for k, v in wu_note.items()} if wu_note else None,
    }
    (OUT / "naming_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("\nWrote", OUT)
    print("Figures:", FIG / "Figure_1B_gse161529_subtype.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
