#!/usr/bin/env python3
"""Refit the Pal KS k=3 triangle without patient TN-0135, then project guests.

TN-0135 is left out of the PCA and the PCHA fit. Those cells are then
projected into the new triangle in the same gene space, without a second
mean/SD match, because they come from the same MAGIC matrix. Wu, GSE173634,
and the 63 DepMap lines are mean/SD-matched to the training cells (Pal
without TN-0135) and transformed with that PCA. The training cells are
scored too, so the inside-fraction table has a self row.

k=3 only. Same PCHA settings as the full Pal fit: 15 starts, delta=0,
keep the largest volume. No permutation test.

Usage (repo root):
  .venv/bin/python -u "Breast Cancer/codes/28_pal_leaveout_tn0135.py"
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

from src.archetypes import fit_pcha_best, t_ratio
from src.io import load_expression_csv
from src.pca import fit_pca

PAL_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse161529_magic" / "ks_genes_magic_imputed.npz"
)
WU_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse176078_magic" / "ks_genes_magic_imputed.npz"
)
G173_MAGIC = (
    BREAST / "results" / "panel_a_ks_genelist_sc_gse173634_magic" / "ks_genes_magic_imputed.npz"
)
PAL_CELLS = (
    BREAST / "results" / "panel_b_gse161529" / "emt_hybrid" / "pal_cells_with_program_state.csv"
)
DEPMAP = BREAST / "data" / "processed" / "input_panelA_ks_genelist.csv"
OUT = BREAST / "results" / "pal_leave_tn0135"
FIG = BREAST / "figures"

HOLD_OUT = "TN-0135"
K = 3
N_INIT = 15
DELTA = 0.0
ARC_COLORS = ["#4C78A8", "#F58518", "#E45756"]


def barycentric(points, vertices):
    a = np.vstack([vertices.T, np.ones((1, vertices.shape[0]))])
    b = np.vstack([points.T, np.ones((1, points.shape[0]))])
    w, *_ = np.linalg.lstsq(a, b, rcond=None)
    return w.T


def match_to_reference(guest: pd.DataFrame, ref: pd.DataFrame) -> np.ndarray:
    ref_mean = ref.mean(axis=1).to_numpy()
    ref_sd = ref.std(axis=1, ddof=1).to_numpy()
    ref_sd[ref_sd < 1e-8] = 1.0
    g_mean = guest.mean(axis=1).to_numpy()
    g_sd = guest.std(axis=1, ddof=1).to_numpy()
    g_sd[g_sd < 1e-8] = 1.0
    return ((guest.to_numpy() - g_mean[:, None]) / g_sd[:, None]) * ref_sd[:, None] + ref_mean[:, None]


def prepare_guest(guest: pd.DataFrame, genes: list[str], ref: pd.DataFrame) -> np.ndarray:
    missing = [g for g in genes if g not in guest.index]
    aligned = guest.reindex(genes)
    ref_mean = ref.mean(axis=1)
    for g in missing:
        aligned.loc[g] = float(ref_mean.loc[g])
    if missing:
        print(f"  filled {missing} with training-set mean")
    return match_to_reference(aligned, ref)


def score_block(name: str, pc: np.ndarray, arcs: np.ndarray) -> dict:
    w = barycentric(pc[:, : arcs.shape[1]], arcs)
    inside = (w >= -1e-6).all(axis=1)
    nearest = w.argmax(axis=1) + 1
    print(
        f"{name}: {int(inside.sum())}/{len(inside)} inside "
        f"({100 * inside.mean():.1f}%)  nearest "
        + ", ".join(f"arc{i}={(nearest == i).mean():.1%}" for i in (1, 2, 3))
    )
    return {"w": w, "inside": inside, "nearest": nearest, "pc": pc}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    print("Loading Pal MAGIC")
    z = np.load(PAL_MAGIC, allow_pickle=True)
    genes = [str(g) for g in z["genes"]]
    barcodes = [str(b) for b in z["barcodes"]]
    expr = pd.DataFrame(z["imputed"], index=genes, columns=barcodes)

    meta = pd.read_csv(PAL_CELLS, usecols=["barcode", "patient", "subtype", "nearest_archetype"])
    meta["barcode"] = meta["barcode"].astype(str)
    meta = meta.drop_duplicates("barcode").set_index("barcode").reindex(barcodes)
    if meta["patient"].isna().any():
        raise SystemExit(f"{int(meta['patient'].isna().sum())} MAGIC barcodes missing a patient label")
    hold = meta["patient"].to_numpy() == HOLD_OUT
    print(f"Pal {expr.shape[1]} cells; holding out {HOLD_OUT}: {int(hold.sum())} cells")
    if hold.sum() < 1000:
        raise SystemExit("held-out patient is too small; check the patient id")

    train_expr = expr.loc[:, ~hold]
    hold_expr = expr.loc[:, hold]
    print(f"Training matrix: {train_expr.shape[0]} genes x {train_expr.shape[1]} cells")

    print("PCA on the training cells")
    pca, train_scores = fit_pca(train_expr.T.to_numpy(), n_components=6)
    print(
        "PC variance:",
        ", ".join(f"PC{i+1}={v:.1%}" for i, v in enumerate(pca.explained_variance_ratio_)),
    )

    print(f"PCHA k={K}, {N_INIT} starts, delta={DELTA}", flush=True)
    archetypes, weights, varexpl, vol, n_ok = fit_pcha_best(
        train_scores, K, n_init=N_INIT, delta=DELTA
    )
    t = t_ratio(train_scores, archetypes)
    print(f"inits_ok={n_ok}/{N_INIT}  ESV={varexpl:.3f}  volume={vol:.4g}  t-ratio={t:.3f}")
    np.save(OUT / "archetypes_k3.npy", archetypes)
    np.save(OUT / "train_pc_scores.npy", train_scores[:, :2])
    pd.Series(train_expr.columns, name="barcode").to_csv(OUT / "train_barcodes.csv", index=False)

    blocks = {}
    blocks["Pal without TN-0135"] = score_block("Pal without TN-0135", train_scores, archetypes)
    hold_pc = pca.transform(hold_expr.T.to_numpy())
    blocks[HOLD_OUT] = score_block(HOLD_OUT, hold_pc, archetypes)

    print("Projecting other datasets onto the training PCA")
    wu = np.load(WU_MAGIC, allow_pickle=True)
    wu_expr = pd.DataFrame(wu["imputed"], index=[str(g) for g in wu["genes"]], columns=[str(b) for b in wu["barcodes"]])
    wu_pc = pca.transform(prepare_guest(wu_expr, genes, train_expr).T)
    blocks["Wu GSE176078"] = score_block("Wu GSE176078", wu_pc, archetypes)

    g173 = np.load(G173_MAGIC, allow_pickle=True)
    g173_expr = pd.DataFrame(
        g173["imputed"], index=[str(g) for g in g173["genes"]], columns=[str(b) for b in g173["barcodes"]]
    )
    g173_pc = pca.transform(prepare_guest(g173_expr, genes, train_expr).T)
    blocks["GSE173634"] = score_block("GSE173634", g173_pc, archetypes)

    dep = load_expression_csv(DEPMAP)
    dep.index = dep.index.astype(str)
    dep_pc = pca.transform(prepare_guest(dep, genes, train_expr).T)
    blocks["DepMap"] = score_block("DepMap", dep_pc, archetypes)

    # Where the old corners went, using only training cells.
    old_nearest = meta.loc[train_expr.columns, "nearest_archetype"].to_numpy(int)
    new_nearest = blocks["Pal without TN-0135"]["nearest"]
    cross = pd.crosstab(
        pd.Series(old_nearest, name="old_nearest_arc"),
        pd.Series(new_nearest, name="new_nearest_arc"),
        normalize="index",
    )
    cross.to_csv(OUT / "old_arc_vs_new_arc_train.csv")
    print("\nTraining cells: fraction of each old arc landing on each new arc")
    print(cross.round(3).to_string())

    hold_nearest = blocks[HOLD_OUT]["nearest"]
    hold_w = blocks[HOLD_OUT]["w"]
    print(
        f"\n{HOLD_OUT} median weights "
        f"{np.median(hold_w, axis=0).round(3).tolist()}  "
        f"inside {blocks[HOLD_OUT]['inside'].mean():.1%}"
    )

    rows = []
    order = ["Pal without TN-0135", HOLD_OUT, "Wu GSE176078", "GSE173634", "DepMap"]
    # Host-swap inside rates on the original full-Pal triangle
    # (results/host_swap/inside_fractions_k3.csv). The saved Pal weights in
    # pal_cells_with_program_state.csv were clipped at 0, so they cannot
    # recover a per-patient inside rate. The published self rate is the
    # whole Pal cohort, 75.8%.
    old_frac = {
        "Pal without TN-0135": 0.7578,
        HOLD_OUT: 0.7578,
        "Wu GSE176078": 0.7838341197,
        "GSE173634": 0.9608460208,
        "DepMap": 0.9047619048,
    }
    old_note = {
        "Pal without TN-0135": (
            "Original 0.758 is the full Pal cohort on the old triangle, not this subset. "
            "Saved Pal weights were clipped, so a per-patient old inside rate is not available."
        ),
        HOLD_OUT: "Original 0.758 is the full Pal cohort on the old triangle, not this patient.",
        "Wu GSE176078": "Original rate is Wu inside the full-Pal triangle.",
        "GSE173634": "Original rate is GSE173634 inside the full-Pal triangle.",
        "DepMap": "Original rate is DepMap inside the full-Pal triangle.",
    }

    for name in order:
        b = blocks[name]
        rows.append(
            {
                "dataset": name,
                "n": int(len(b["inside"])),
                "frac_inside_new": float(b["inside"].mean()),
                "frac_inside_original_pal_triangle": old_frac[name],
                "frac_nearest_arc1": float((b["nearest"] == 1).mean()),
                "frac_nearest_arc2": float((b["nearest"] == 2).mean()),
                "frac_nearest_arc3": float((b["nearest"] == 3).mean()),
                "median_w1": float(np.median(b["w"][:, 0])),
                "median_w2": float(np.median(b["w"][:, 1])),
                "median_w3": float(np.median(b["w"][:, 2])),
                "original_rate_scope": old_note[name],
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "inside_fractions.csv", index=False)
    print("\n", summary.to_string(index=False))

    # KS program at the new vertices: mean score of the closest 5% of training cells.
    gl = pd.read_csv(BREAST / "data" / "processed" / "ks_cellline_signature_tan2014.csv")
    epi_genes = [g for g in gl.loc[gl["category"] == "Epi", "gene_symbol"] if g in train_expr.index]
    mes_genes = [g for g in gl.loc[gl["category"] == "Mes", "gene_symbol"] if g in train_expr.index]
    zscore = train_expr.sub(train_expr.mean(axis=1), axis=0)
    zscore = zscore.div(train_expr.std(axis=1, ddof=1).replace(0, np.nan).fillna(1.0), axis=0)
    epi_s = zscore.loc[epi_genes].mean(axis=0).to_numpy()
    mes_s = zscore.loc[mes_genes].mean(axis=0).to_numpy()
    dist = np.linalg.norm(train_scores[:, :2][:, None, :] - archetypes[None, :, :2], axis=2)
    n_take = int(np.ceil(0.05 * train_scores.shape[0]))
    vertex_rows = []
    for j in range(3):
        idx = np.argpartition(dist[:, j], n_take - 1)[:n_take]
        vertex_rows.append(
            {
                "new_arc": j + 1,
                "mean_epi": float(epi_s[idx].mean()),
                "mean_mes": float(mes_s[idx].mean()),
                "mean_emt": float((mes_s[idx] - epi_s[idx]).mean()),
            }
        )
    vertices = pd.DataFrame(vertex_rows)
    vertices.to_csv(OUT / "new_vertex_ks_means.csv", index=False)
    print("\nNew vertices, KS means on the closest 5% of training cells")
    print(vertices.round(3).to_string(index=False))

    report = {
        "held_out_patient": HOLD_OUT,
        "n_held_out": int(hold.sum()),
        "n_train": int((~hold).sum()),
        "k": K,
        "n_init": N_INIT,
        "n_init_ok": int(n_ok),
        "delta": DELTA,
        "esv": float(varexpl),
        "t_ratio": float(t),
        "note": "t-ratio is the observed volume ratio only. Shuffles were not rerun.",
        "projection": (
            "TN-0135 transformed with the training PCA in the shared MAGIC gene space. "
            "Wu, GSE173634, and DepMap mean/SD-matched to the training cells first."
        ),
    }
    (OUT / "fit_report.json").write_text(json.dumps(report, indent=2))

    # Small weights table for TN-0135 only, not the full cell tables.
    pd.DataFrame(
        {
            "barcode": hold_expr.columns,
            "patient": HOLD_OUT,
            "subtype": meta.loc[hold_expr.columns, "subtype"].to_numpy(),
            "inside_simplex": blocks[HOLD_OUT]["inside"],
            "nearest_archetype": hold_nearest,
            "PC1": hold_pc[:, 0],
            "PC2": hold_pc[:, 1],
            "w_arc1": hold_w[:, 0],
            "w_arc2": hold_w[:, 1],
            "w_arc3": hold_w[:, 2],
        }
    ).to_csv(OUT / "tn0135_in_new_triangle.csv", index=False)

    rng = np.random.default_rng(0)
    style_scatter(train_scores, hold_pc, archetypes, blocks, rng)
    print("Wrote", OUT)


def style_scatter(train_scores, hold_pc, arcs, blocks, rng):
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
            "savefig.dpi": 200,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))

    def draw(ax):
        pts = list(arcs[:, :2]) + [arcs[:1, :2]]
        xy = np.vstack(pts)
        ax.plot(xy[:, 0], xy[:, 1], color="#333333", lw=1.4, zorder=5)
        ax.scatter(
            arcs[:, 0], arcs[:, 1], s=90, c=ARC_COLORS, edgecolors="white", linewidths=0.6, zorder=6,
        )
        for i in range(3):
            ax.annotate(str(i + 1), arcs[i, :2], textcoords="offset points", xytext=(6, 6), fontsize=8)

    take = rng.choice(train_scores.shape[0], size=8000, replace=False)
    ax = axes[0]
    ax.scatter(train_scores[take, 0], train_scores[take, 1], s=2, c="#B0B0B0", linewidths=0, alpha=0.35)
    draw(ax)
    inn = blocks["Pal without TN-0135"]["inside"].mean()
    ax.set_title(f"Training cells, Pal without TN-0135\n{100 * inn:.1f}% inside")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

    ax = axes[1]
    ax.scatter(train_scores[take, 0], train_scores[take, 1], s=2, c="#E6E6E6", linewidths=0, alpha=0.25)
    htake = rng.choice(hold_pc.shape[0], size=min(8000, hold_pc.shape[0]), replace=False)
    ax.scatter(hold_pc[htake, 0], hold_pc[htake, 1], s=3, c="#4C78A8", linewidths=0, alpha=0.35)
    draw(ax)
    inn = blocks[HOLD_OUT]["inside"].mean()
    ax.set_title(f"TN-0135 projected in\n{100 * inn:.1f}% inside")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    fig.savefig(FIG / "Figure_pal_leaveout_tn0135.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
