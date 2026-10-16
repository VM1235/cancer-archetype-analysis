#!/usr/bin/env python3
"""Step 6: specialists vs generalists with Groves 2022's closest-5% rule.

Rule (Groves et al. 2022, Cell Systems, Figure 3E legend): for each archetype,
the 5% of cells closest to it are that archetype's specialists; every other
cell is a generalist. Distance = Euclidean distance in the simplex's own space
(the first k-1 PCs).

Which cells: the "reference" cells, i.e. the population the simplex was fitted
on (all cells for tag main; in_vivo cells for tag in_vivo). The 5% cut-off distance
is set on the reference cells. Other cells (e.g. in_vitro guests) are then
classified with the same cut-off distances and reported separately.

For comparison it also reports Groves' second rule (STAR Methods: weight > 0.95
on one archetype), and the purity-gap rule (>= 0.35) used in the Pal work.

Outputs (results/<geneset>_<tag>/groves5_k<k>/):
  cells_groves5.csv.gz   class of every cell
  summary.md             all tables

Usage (from Project_1/):
  .venv/bin/python -u "Breast Cancer/Mouse_4T1_GSE303505"/codes/06_groves_specialists.py --geneset all --tag in_vivo --k 3
"""
import argparse

import numpy as np
import pandas as pd

import scc_common as C

FRACTION = 0.05


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geneset", choices=["all", "ks"], default="all")
    ap.add_argument("--tag", default="in_vivo")
    ap.add_argument("--k", type=int, default=3)
    args = ap.parse_args()
    k = args.k
    run = C.run_dir(args.geneset, args.tag)
    out = run / f"groves5_k{k}"
    out.mkdir(exist_ok=True)

    z = np.load(run / "pca.npz")
    arcs = np.load(run / f"archetypes_k{k}.npy")
    meta = pd.read_csv(C.PROC / "cell_metadata.csv").set_index("barcode").loc[z["barcodes"]]   # cleaned fits use a subset of the cells
    assert list(meta.index) == list(z["barcodes"]), "cell order differs from the PCA file"
    pts = z["scores_all"][:, : k - 1]
    V = [f"V{i + 1}" for i in range(k)]

    ref = (meta["assay"] == "in_vivo").values if args.tag == "in_vivo" else np.ones(len(meta), bool)
    d = np.linalg.norm(pts[:, None, :] - arcs[None, :, : k - 1], axis=2)   # cells x vertices
    cut = np.quantile(d[ref], FRACTION, axis=0)                             # per-vertex cut-off
    close = d <= cut
    n_multi = close.sum(1)
    # a cell close to two vertices goes to the nearer one (relative to that vertex's cut-off)
    cls = np.where(n_multi == 0, "generalist",
                   np.array(V)[np.argmin(d / cut, axis=1)])

    w = C.composition(C.barycentric(pts, arcs))
    pur, near = C.purity_and_nearest(w)
    cells = pd.DataFrame({
        "groves5": cls,
        "weight095": np.where(w.max(1) > 0.95, np.array(V)[near], "generalist"),
        "gap035": np.where(pur >= C.PURITY_CUT, np.array(V)[near], "generalist"),
        "reference_cell": ref,
    }, index=meta.index)
    for c in ["sample", "assay", "cluster_annotation", "clusterU", "phase"]:
        cells[c] = meta[c].astype(str).values
    for j, v in enumerate(V):
        cells[f"dist_{v}"] = d[:, j]
    cells.to_csv(out / "cells_groves5.csv.gz")

    cols = V + ["generalist"]
    R = cells[ref]
    lines = [f"# Groves closest-5% specialists: {args.geneset}/{args.tag}, k={k}", "",
             f"Reference cells (5% cut-off set on these): {int(ref.sum())}",
             f"Cells within the 5% cut-off of two vertices: {int((n_multi[ref] > 1).sum())}",
             f"Specialists among reference cells: {100 * (R['groves5'] != 'generalist').mean():.1f}%"]

    def add(title, df, note=None):
        lines.extend(["", f"## {title}", ""] + ([note, ""] if note else []) + [df.to_string()])
        print(f"\n== {title} ==\n{df.to_string()}")

    def tabs(rows, order, label):
        t = pd.crosstab(R[rows], R["groves5"]).reindex(columns=cols, fill_value=0)
        t = t.reindex([o for o in order if o in t.index])
        t["total"] = t.sum(1)
        add(f"{label}: counts", t)
        row = (100 * t[cols].div(t["total"], axis=0)).round(1)
        add(f"{label}: % of each {label} in each class (rows sum to 100)", row,
            "Read across: where do the cells of this label go?")
        col = (100 * t[cols].div(t[cols].sum(0), axis=1)).round(1)
        base = (100 * t["total"] / t["total"].sum()).round(1).rename("all reference cells")
        add(f"Each class: % made up of each {label} (columns sum to 100)",
            pd.concat([col, base], axis=1),
            "Read down: what is this class made of? Compare with the last column (overall share).")

    tabs("cluster_annotation", C.STATE_ORDER, "state")
    tabs("clusterU", C.CLUSTERU_ORDER, "clusterU")

    tum = pd.crosstab(R["groves5"], R["sample"]).reindex(cols)
    add("Tumour of each class (% of class)", (100 * tum.div(tum.sum(1), axis=0)).round(1))
    ph = pd.crosstab(R["groves5"], R["phase"]).reindex(cols)
    add("Cell-cycle phase of each class (% of class)", (100 * ph.div(ph.sum(1), axis=0)).round(1))

    # rule comparison
    comp = pd.DataFrame({r: R[r].value_counts().reindex(cols, fill_value=0)
                         for r in ["groves5", "weight095", "gap035"]})
    add("Rule comparison, reference cells (counts)", comp,
        "groves5 = closest 5% (Fig 3E); weight095 = weight > 0.95 (STAR Methods); "
        "gap035 = purity gap >= 0.35 (Pal work)")
    sp = pd.DataFrame({r: (100 * (R[r] != "generalist").groupby(R["cluster_annotation"]).mean())
                       .reindex(C.STATE_ORDER) for r in ["groves5", "weight095", "gap035"]}).round(1)
    add("% of each state called specialist, by rule", sp)

    if (~ref).any():
        G = cells[~ref]
        g = pd.crosstab(G["cluster_annotation"], G["groves5"]).reindex(columns=cols, fill_value=0)
        g = g.reindex([o for o in C.STATE_ORDER if o in g.index])
        add("Non-reference cells (e.g. in_vitro), same cut-offs: % of each state",
            (100 * g.div(g.sum(1), axis=0)).round(1))

    (out / "summary.md").write_text("\n".join(lines) + "\n")
    if k == 3:
        plot_triangle(cells, pts, arcs, cut, ref, V, out, f"{args.geneset}/{args.tag}")
    print("\nwrote", out)


def plot_triangle(cells, pts, arcs, cut, ref, V, out, title):
    """Reference cells in PC1-PC2 with the triangle; Groves specialists coloured."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    R = cells[ref]
    P = pts[ref]
    order = np.random.default_rng(0).permutation(len(R))
    gen = (R["groves5"] == "generalist").values
    vcol = {"V1": "#7b3294", "V2": "#d7301f", "V3": "#1a9850"}
    scol = {"Epithelial": "#1f77b4", "Early-Hybrid": "#ff7f0e",
            "Late-Hybrid": "#2ca02c", "Mesenchymal": "#d62728"}
    tumours = sorted(R["sample"].unique())
    tcmap = plt.get_cmap("tab10")
    tcol = {t: tcmap(i) for i, t in enumerate(tumours)}

    def frame(ax):
        for i in range(3):
            j = (i + 1) % 3
            ax.plot(arcs[[i, j], 0], arcs[[i, j], 1], c="k", lw=1)
        for i, v in enumerate(V):
            ax.add_patch(plt.Circle(arcs[i, :2], cut[i], fill=False, ls="--", lw=0.8,
                                    color="0.3"))
            ax.annotate(v, arcs[i, :2], fontsize=12, fontweight="bold",
                        xytext=(6, 6), textcoords="offset points")
        ax.scatter(arcs[:, 0], arcs[:, 1], s=50, c="k", zorder=5)
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.set_aspect("equal", "datalim")

    def grey(ax):
        ax.scatter(P[gen, 0], P[gen, 1], s=1, c="#d0d0d0", rasterized=True)

    fig, axes = plt.subplots(1, 3, figsize=(20, 6.6))

    ax = axes[0]
    grey(ax)
    s = ~gen
    ax.scatter(P[s, 0], P[s, 1], s=3, c=[vcol[v] for v in R["groves5"].values[s]],
               rasterized=True)
    for v in V:
        ax.scatter([], [], s=25, c=vcol[v], label=f"{v} specialist (n={int((R['groves5'] == v).sum())})")
    ax.scatter([], [], s=25, c="#d0d0d0", label=f"generalist (n={int(gen.sum())})")
    ax.legend(fontsize=8, loc="best")
    ax.set_title("Groves closest-5% specialists per vertex")
    frame(ax)

    ax = axes[1]
    grey(ax)
    idx = np.where(~gen)[0]
    idx = idx[np.isin(idx, order)]
    st = R["cluster_annotation"].values
    ax.scatter(P[idx, 0], P[idx, 1], s=3, c=[scol[x] for x in st[idx]], rasterized=True)
    for x, c in scol.items():
        n = int(((~gen) & (st == x)).sum())
        ax.scatter([], [], s=25, c=c, label=f"{x} (n={n})")
    ax.scatter([], [], s=25, c="#d0d0d0", label="generalist")
    ax.legend(fontsize=8, loc="best", title="specialists coloured by state", title_fontsize=8)
    ax.set_title("Specialists by EMT state")
    frame(ax)

    ax = axes[2]
    grey(ax)
    tm = R["sample"].values
    ax.scatter(P[idx, 0], P[idx, 1], s=3, c=[tcol[x] for x in tm[idx]], rasterized=True)
    for t in tumours:
        n = int(((~gen) & (tm == t)).sum())
        ax.scatter([], [], s=25, color=tcol[t], label=f"{t} (n={n})")
    ax.legend(fontsize=8, loc="best", title="specialists coloured by tumour", title_fontsize=8)
    ax.set_title("Specialists by tumour")
    frame(ax)

    fig.suptitle(f"{title} k=3: {int(ref.sum())} reference cells. Dashed circles = 5% cut-off "
                 f"distance per vertex (PC1-PC2)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "triangle_groves5.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
