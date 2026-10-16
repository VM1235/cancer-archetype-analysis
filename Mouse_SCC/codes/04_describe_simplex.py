#!/usr/bin/env python3
"""Step 4: place every cell in the fitted simplex and overlay the paper's labels.

What it does (for one chosen k)
  1. Barycentric weights for all 23,959 cells: how much of each vertex a cell is.
     Purity = largest weight - second largest. Specialist if purity >= 0.35,
     otherwise interior (a mix of vertices). Same rule as the Pal work.
  2. Tables, vertex by:
       EMT state (cluster_annotation), finer cluster (clusterU), tumour, assay,
       cell-cycle phase. Shown for specialists, since those define a vertex.
  3. Tumour check: what share of each vertex's specialists comes from its
     largest tumour, and how many tumours contribute at least 5%.
  4. Interior share per state: are early-hybrid cells more often interior?
  5. Distance-bin enrichment (ParTI style): split cells into 10 equal bins by
     distance to each vertex; test which labels are over-represented in the
     closest bin (hypergeometric, BH-corrected).
  6. Marker genes per vertex: mean log-expression in the closest 10% of cells
     minus the mean in all other cells (all 14,502 genes).
  7. KS EMT score (mean z of mouse KS mesenchymal genes - epithelial genes) and
     the authors' LKPR signature scores, averaged per vertex.
  8. Figures: cells in PC space with the simplex drawn, coloured by each label.

Vertex names are NOT assigned automatically. The script prints the evidence;
naming is a decision made after looking at the tumour breakdown.

Usage (from Project_1/):
  .venv/bin/python -u Mouse_SCC/codes/04_describe_simplex.py --geneset all --k 3
"""
import argparse
import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scc_common as C
from src.enrichment import distance_bins, hypergeometric_enrichment


def pct_table(rows, cols, row_order=None, col_order=None):
    t = pd.crosstab(rows, cols)
    if col_order:
        t = t.reindex(columns=[c for c in col_order if c in t.columns])
    if row_order:
        t = t.reindex([r for r in row_order if r in t.index])
    return t, (100 * t.div(t.sum(1), axis=0)).round(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geneset", choices=["all", "ks"], default="all")
    ap.add_argument("--tag", default="main")
    ap.add_argument("--k", type=int, required=True)
    args = ap.parse_args()
    k = args.k
    run = C.run_dir(args.geneset, args.tag)
    out = run / f"describe_k{k}"
    out.mkdir(exist_ok=True)

    z = np.load(run / "pca.npz")
    scores, fit_idx = z["scores_all"], z["fit_idx"]
    arcs = np.load(run / f"archetypes_k{k}.npy")
    a = C.load_adata()
    assert list(a.obs_names) == list(z["barcodes"]), "cell order differs from the PCA file"
    obs = a.obs.copy()
    V = [f"V{i + 1}" for i in range(k)]

    # ---- 1. weights
    w = C.barycentric(scores[:, : k - 1], arcs)
    comp = C.composition(w)
    pur, near = C.purity_and_nearest(w)
    cells = pd.DataFrame(comp, columns=[f"w_{v}" for v in V], index=obs.index)
    cells["inside"] = (w >= -1e-6).all(1)
    cells["nearest"] = [V[i] for i in near]
    cells["purity"] = pur
    cells["specialist"] = pur >= C.PURITY_CUT
    cells["group"] = np.where(cells["specialist"], cells["nearest"], "interior")
    cells["fit_cell"] = False
    cells.iloc[fit_idx, cells.columns.get_loc("fit_cell")] = True
    for c in ["sample", "assay", "cluster_annotation", "clusterU", "phase"]:
        cells[c] = obs[c].astype(str).values

    # ---- 7. scores per cell
    epi, mes = C.ks_genes_mouse(a.var_names)
    def zmean(genes):
        idx = [a.var_names.get_loc(g) for g in genes]
        M = a.X[:, idx].toarray()
        return ((M - M.mean(0)) / np.where(M.std(0) > 0, M.std(0), 1)).mean(1)
    cells["KS_EMT"] = zmean(mes) - zmean(epi)
    lkpr = [c for c in obs.columns if c.startswith("LKPR_")]
    for c in lkpr:
        cells[c] = obs[c].values
    cells.to_csv(out / "cells_in_simplex.csv.gz")

    sp = cells[cells["specialist"]]
    lines = [f"# Simplex description: gene set {args.geneset}, tag {args.tag}, k={k}", ""]
    lines.append(f"Cells: {len(cells)}; inside simplex {100 * cells['inside'].mean():.1f}%; "
                 f"specialists {100 * cells['specialist'].mean():.1f}%; "
                 f"interior {100 * (~cells['specialist']).mean():.1f}%")
    lines.append(f"Fit cells inside simplex: {100 * cells.loc[cells.fit_cell, 'inside'].mean():.1f}%")

    def add(title, df):
        lines.extend(["", f"## {title}", "", df.to_string()])
        print(f"\n== {title} ==\n{df.to_string()}")

    # ---- 2. label tables (specialists)
    t_state, p_state = pct_table(sp["nearest"], sp["cluster_annotation"], V, C.STATE_ORDER)
    add("Specialists: vertex x EMT state (counts)", t_state)
    add("Specialists: vertex x EMT state (% of row)", p_state)
    t_u, p_u = pct_table(sp["nearest"], sp["clusterU"], V, C.CLUSTERU_ORDER)
    add("Specialists: vertex x clusterU (% of row)", p_u)
    add("Specialists: vertex x phase (% of row)", pct_table(sp["nearest"], sp["phase"], V)[1])
    add("Specialists: vertex x assay (% of row)", pct_table(sp["nearest"], sp["assay"], V)[1])
    t_tum, p_tum = pct_table(sp["nearest"], sp["sample"], V)
    add("Specialists: vertex x tumour (% of row)", p_tum)
    # where does each state go (all cells, incl. interior)
    _, p_where = pct_table(cells["cluster_annotation"], cells["group"], C.STATE_ORDER,
                           V + ["interior"])
    add("Each EMT state: % of its cells at each vertex or interior", p_where)
    _, p_whereU = pct_table(cells["clusterU"], cells["group"], C.CLUSTERU_ORDER,
                            V + ["interior"])
    add("Each clusterU: % of its cells at each vertex or interior", p_whereU)
    _, p_whereT = pct_table(cells["sample"], cells["group"], None, V + ["interior"])
    add("Each tumour: % of its cells at each vertex or interior", p_whereT)

    # ---- 3. tumour check
    tc = []
    for v in V:
        s = sp.loc[sp["nearest"] == v, "sample"].value_counts(normalize=True)
        tc.append(dict(vertex=v, n_specialists=int((sp["nearest"] == v).sum()),
                       top_tumour=s.index[0] if len(s) else "", top_share_pct=100 * s.iloc[0]
                       if len(s) else np.nan, n_tumours_ge5pct=int((s >= 0.05).sum())))
    tc = pd.DataFrame(tc).set_index("vertex").round(1)
    add("Tumour check per vertex", tc)

    # ---- 7b. scores per vertex group
    sc = cells.groupby("group")[["KS_EMT"] + lkpr + ["purity"]].mean().round(3)
    add("Mean scores per vertex (specialists) and interior", sc)

    # ---- 5. distance-bin enrichment, all cells
    bins, _ = distance_bins(scores[:, : k - 1], arcs, n_bins=10)
    enr = []
    for col, levels in [("cluster_annotation", C.STATE_ORDER), ("clusterU", C.CLUSTERU_ORDER),
                        ("phase", None), ("sample", None), ("assay", None)]:
        lv = levels or sorted(cells[col].unique())
        e = hypergeometric_enrichment(cells[col].values, bins, subtype_levels=lv)
        e.insert(0, "label_type", col)
        enr.append(e)
    enr = pd.concat(enr)
    enr["archetype"] = enr["archetype"].map(lambda i: V[i])
    enr.to_csv(out / "distance_bin_enrichment.csv", index=False)
    b0 = enr[(enr["bin"] == 0) & enr["significant"] & (enr["fold_enrichment"] >= 1.5)]
    b0 = b0.sort_values(["archetype", "fold_enrichment"], ascending=[True, False])
    add("Closest 10% of cells to each vertex: labels enriched >= 1.5-fold (q < 0.1)",
        b0[["archetype", "label_type", "subtype", "k_in_bin", "fold_enrichment", "q_value"]]
        .set_index("archetype").round(4))

    # ---- 6. marker genes per vertex
    X = a.X
    mk = {}
    for j, v in enumerate(V):
        close = bins[:, j] == 0
        d = np.asarray(X[close].mean(0)).ravel() - np.asarray(X[~close].mean(0)).ravel()
        mk[f"{v}_up"] = a.var_names[np.argsort(d)[::-1][:30]].tolist()
        mk[f"{v}_up_diff"] = np.round(np.sort(d)[::-1][:30], 3)
    mk = pd.DataFrame(mk)
    mk.to_csv(out / "vertex_marker_genes.csv", index=False)
    add("Top 15 genes higher near each vertex (log-expression difference)",
        mk[[f"{v}_up" for v in V]].head(15))

    (out / "summary.md").write_text("\n".join(lines) + "\n")

    # ---- 8. figures
    rng = np.random.default_rng(0)
    show = rng.permutation(len(cells))
    pairs = [(0, 1)] if k == 3 else list(itertools.combinations(range(min(k - 1, 3)), 2))
    panels = [("cluster_annotation", C.STATE_ORDER), ("clusterU", C.CLUSTERU_ORDER),
              ("sample", None), ("assay", None), ("phase", None), ("KS_EMT", None),
              ("purity", None)]
    for (i, j) in pairs:
        fig, axes = plt.subplots(2, 4, figsize=(20, 9.5))
        for ax, (col, order) in zip(axes.ravel(), panels):
            x, y = scores[show, i], scores[show, j]
            if col in ("KS_EMT", "purity"):
                vals = cells[col].values[show]
                lo, hi = np.quantile(vals, [0.02, 0.98])
                im = ax.scatter(x, y, c=vals, s=1.5, cmap="coolwarm" if col == "KS_EMT"
                                else "viridis", vmin=lo, vmax=hi, rasterized=True)
                plt.colorbar(im, ax=ax, shrink=0.7)
            else:
                # draw all cells in one random order so no label hides another
                lab = cells[col].values[show]
                levels = [lv for lv in (order or sorted(np.unique(lab))) if lv in set(lab)]
                cmap = plt.get_cmap("tab10" if len(levels) <= 10 else "tab20")
                col_of = {lv: cmap(n) for n, lv in enumerate(levels)}
                ax.scatter(x, y, s=1.5, c=[col_of[v] for v in lab], rasterized=True)
                for lv in levels:
                    ax.scatter([], [], s=12, color=col_of[lv], label=lv)
                ax.legend(fontsize=6, loc="best", framealpha=0.7)
            for p, q in itertools.combinations(range(k), 2):
                ax.plot(arcs[[p, q], i], arcs[[p, q], j], c="k", lw=0.8)
            ax.scatter(arcs[:, i], arcs[:, j], s=60, c="k", zorder=5)
            for n in range(k):
                ax.annotate(V[n], (arcs[n, i], arcs[n, j]), fontsize=11, fontweight="bold",
                            xytext=(5, 5), textcoords="offset points")
            ax.set_title(col); ax.set_xlabel(f"PC{i + 1}"); ax.set_ylabel(f"PC{j + 1}")
        # last panel: where each state goes
        ax = axes.ravel()[-1]
        p_where.plot.barh(stacked=True, ax=ax, colormap="tab10", width=0.8)
        ax.set_xlabel("% of state"); ax.set_ylabel("")
        ax.set_title("Each state: vertex or interior")
        ax.legend(fontsize=7)
        fig.suptitle(f"{args.geneset}/{args.tag}  k={k}  (all {len(cells)} cells)")
        fig.tight_layout()
        fig.savefig(out / f"simplex_PC{i + 1}_PC{j + 1}.png", dpi=170)
        plt.close(fig)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
