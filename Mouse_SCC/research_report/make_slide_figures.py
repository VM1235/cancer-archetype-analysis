#!/usr/bin/env python3
"""Slide-sized figures for the professor deck. Reads saved results only; refits nothing.

Writes Mouse_SCC/research_report/figures/slides/*.png
Usage (from Project_1/): .venv/bin/python Mouse_SCC/research_report/make_slide_figures.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"
OUT = HERE / "figures" / "slides"
OUT.mkdir(parents=True, exist_ok=True)
S = ["Epithelial", "Early-Hybrid", "Hybrid", "Late-EMT"]
SC = {"Epithelial": "#2b6cb0", "Early-Hybrid": "#e8913a", "Hybrid": "#3a9d5d", "Late-EMT": "#c0392b"}
plt.rcParams.update({"font.size": 15, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "#ffffff", "savefig.dpi": 160, "savefig.bbox": "tight"})


def triangle(run, labels, title, name):
    c = pd.read_csv(RES / run / "describe_k3" / "cells_in_simplex.csv.gz", index_col=0)
    z = np.load(RES / run / "pca.npz")
    arcs = np.load(RES / run / "archetypes_k3.npy")
    m = (c["assay"] == "scRNA").values
    pts = z["scores_all"][m, :2]
    st = c["cluster_annotation"].values[m]
    o = np.random.default_rng(0).permutation(len(pts))
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(pts[o, 0], pts[o, 1], s=3, c=[SC[x] for x in st[o]], rasterized=True)
    t = np.vstack([arcs, arcs[:1]])
    ax.plot(t[:, 0], t[:, 1], "k-", lw=2)
    ax.scatter(arcs[:, 0], arcs[:, 1], s=120, c="#1a1a1a", zorder=5)
    for j, lab in enumerate(labels):
        ax.annotate(lab, arcs[j], fontsize=15, fontweight="bold", xytext=(8, 8), textcoords="offset points",
                    bbox=dict(boxstyle="round,pad=0.3", fc="#ffffff", ec="#999999", alpha=0.9))
    for s in S:
        ax.scatter([], [], s=60, c=SC[s], label=s)
    ax.legend(loc="lower left", fontsize=12, frameon=True)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.set_title(title, fontsize=16)
    fig.savefig(OUT / name); plt.close(fig)


def nulls():
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, run, lab in [(axes[0], "all_scrna", "Before: with the V1 axis"),
                         (axes[1], "all_scrna_noV1", "After: V1 axis removed")]:
        nl = np.load(RES / run / "null_k3.npy")
        f = json.load(open(RES / run / "fit_k3.json"))
        ax.hist(nl, bins=25, color="#b8c2cc", edgecolor="white")
        ax.axvline(f["t_ratio"], color="#c0392b", lw=3)
        n_ge = int((nl >= f["t_ratio"]).sum())
        ax.set_title(f"{lab}\nreal t = {f['t_ratio']:.3f}; {n_ge}/200 shuffles as good; p = {f['p_value']:.3f}",
                     fontsize=14)
        ax.set_xlabel("triangle fit (t-ratio) of shuffled data")
    axes[0].set_ylabel("number of shuffles")
    fig.tight_layout(); fig.savefig(OUT / "nulls_before_after.png"); plt.close(fig)


def drivers():
    rows = []
    for run, lab in [("all_main", "All genes,\nall cells"), ("all_scrna", "All genes,\nscRNA only"),
                     ("ks_main", "KS EMT genes")]:
        a = pd.read_csv(RES / run / "pc_drivers.csv").set_index("PC").loc[1]
        b = pd.read_csv(RES / run / "pc_drivers_within_state.csv").set_index("PC").loc[1]
        rows.append((lab, a["state"], b["tumour_in_state"], a["assay"]))
    d = pd.DataFrame(rows, columns=["v", "EMT state", "Tumour identity\n(within a state)", "Assay"]).set_index("v")
    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(d)); w = 0.26
    for i, (col, c) in enumerate(zip(d.columns, ["#2b6cb0", "#c0392b", "#8a8a8a"])):
        ax.bar(x + (i - 1) * w, d[col], w, color=c, label=col)
        for xi, v in zip(x, d[col]):
            ax.text(xi + (i - 1) * w, v + 0.015, f"{v:.2f}", ha="center", fontsize=12)
    ax.set_xticks(x); ax.set_xticklabels(d.index)
    ax.set_ylabel("share of PC1 explained"); ax.set_ylim(0, 1)
    ax.legend(fontsize=12, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.12), frameon=False)
    fig.savefig(OUT / "pc1_drivers.png"); plt.close(fig)


def within():
    R = pd.read_csv(RES / "within_tumour" / "corner_signature_correlation.csv", index_col=0)
    comp = pd.read_csv(RES / "within_tumour" / "corner_composition.csv")
    rows = [r for r in R.index if not r.startswith("G_")]
    M = R.loc[rows, ["G_epi", "G_hybrid", "G_late"]]
    lab = []
    for r in rows:
        c = comp[comp["corner"] == r].iloc[0]
        top = max(S, key=lambda s: c[s])
        lab.append(f"{r.replace('_', ' corner ')}  ({top} {c[top]:.0f}%)")
    pv = {}
    for t in ["SC15", "SC32", "SC46"]:
        pv[t] = json.load(open(RES / "within_tumour" / t / "fit_k3.json"))["p_value"]
    fig, ax = plt.subplots(figsize=(11, 8))
    im = ax.imshow(M.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(3):
            v = M.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=14,
                    color="white" if abs(v) > 0.55 else "#1a1a1a", fontweight="bold" if v > 0.6 else "normal")
    ax.set_xticks(range(3)); ax.set_xticklabels(["pooled\nepithelial corner", "pooled\nhybrid corner",
                                                 "pooled\nlate-EMT corner"])
    ax.set_yticks(range(len(lab))); ax.set_yticklabels(lab, fontsize=13)
    for k, t in enumerate(["SC15", "SC32", "SC46"]):
        ax.axhline(3 * k - 0.5, color="#1a1a1a", lw=2)
        ax.text(2.62, 3 * k + 1, f"{t}\ntriangle p = {pv[t]:.3f}", va="center", fontsize=13)
    ax.xaxis.tick_top()
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.6, pad=0.04); cb.set_label("gene-signature correlation (1 = same genes up)")
    fig.savefig(OUT / "within_tumour_corners.png"); plt.close(fig)


if __name__ == "__main__":
    triangle("all_scrna", ["V1 (long-gene axis)", "V2 Late EMT", "V3 dividing epithelial"],
             "Version B, k = 3 (scRNA cells)", "triangle_before.png")
    triangle("all_scrna_noV1", ["Late EMT", "Proliferating hybrid", "Epithelial / early hybrid"],
             "After removing the V1 axis, k = 3 (scRNA cells)", "triangle_after.png")
    nulls(); drivers(); within()
    print("wrote", sorted(p.name for p in OUT.glob("*.png")))


# ---------------------------------------------------------------- deck v2: triangle-only main slides
TCOL = {"MU3": "#9e9ac8", "MU4": "#807dba", "MU7": "#6a51a3", "MU8": "#4a1486", "SC15": "#e6550d",
        "SC46": "#fd8d3c", "SC32": "#3182bd", "SC33": "#6baed6", "SC2b1": "#31a354", "SC2b3": "#74c476",
        "SC35": "#636363", "SC43": "#969696"}
ACOL = {"scRNA": "#3182bd", "multiome": "#6a51a3"}


def _tri(ax, pts, cols, arcs, labels, title):
    o = np.random.default_rng(0).permutation(len(pts))
    ax.scatter(pts[o, 0], pts[o, 1], s=2.5, c=[cols[i] for i in o], rasterized=True)
    t = np.vstack([arcs[:, :2], arcs[:1, :2]])
    ax.plot(t[:, 0], t[:, 1], "-", color="#1a1a1a", lw=2)
    ax.scatter(arcs[:, 0], arcs[:, 1], s=110, c="#1a1a1a", zorder=5)
    for j, lab in enumerate(labels):
        ax.annotate(lab, arcs[j, :2], fontsize=14, fontweight="bold", xytext=(8, 8), textcoords="offset points",
                    bbox=dict(boxstyle="round,pad=0.25", fc="#ffffff", ec="#aaaaaa", alpha=0.92))
    ax.set_title(title, fontsize=16)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2")


def _legend(ax, d, loc="lower left", ncol=1):
    for k, c in d.items():
        ax.scatter([], [], s=55, c=c, label=k)
    ax.legend(loc=loc, fontsize=11, ncol=ncol, frameon=True, framealpha=0.9)


def pair(run, labels, second, cells, name, titles):
    c = pd.read_csv(RES / run / "describe_k3" / "cells_in_simplex.csv.gz", index_col=0)
    z = np.load(RES / run / "pca.npz")
    arcs = np.load(RES / run / "archetypes_k3.npy")
    m = (c["assay"] == "scRNA").values if cells == "scRNA" else np.ones(len(c), bool)
    pts = z["scores_all"][m, :2]
    sub = c[m]
    fig, axes = plt.subplots(1, 2, figsize=(18, 7.2))
    _tri(axes[0], pts, [SC[x] for x in sub["cluster_annotation"]], arcs, labels, titles[0])
    _legend(axes[0], SC)
    cmap = TCOL if second == "sample" else ACOL
    _tri(axes[1], pts, [cmap[x] for x in sub[second]], arcs, labels, titles[1])
    _legend(axes[1], {k: v for k, v in cmap.items() if k in set(sub[second])}, ncol=2 if second == "sample" else 1)
    fig.tight_layout(); fig.savefig(OUT / name); plt.close(fig)


def within_triangles():
    a_obs = pd.read_csv(HERE.parent / "data" / "processed" / "cell_metadata.csv")
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))
    for ax, t in zip(axes, ["SC15", "SC32", "SC46"]):
        d = RES / "within_tumour" / t
        z = np.load(d / "pca.npz")
        arcs = np.load(d / "archetypes_k3.npy")
        comp = pd.read_csv(RES / "within_tumour" / "corner_composition.csv")
        # labels = what the corner's gene signature matches (08 --compare) + its state mix
        labs = {"SC15": ["Epithelial-like (early hyb. 73%)", "Late-like (late 48%, early hyb. 50%)",
                         "Pure epithelial (100%)"],
                "SC32": ["Late-like (late 57%)", "Epithelial (97%)", "Early hybrid (99%)"],
                "SC46": ["Late EMT, stress genes (86%)", "Late-like (late 96%)",
                         "Hybrid-like (hybrid 64%)"]}[t]
        st = a_obs["cluster_annotation"].values[z["rows"]]
        p = json.load(open(d / "fit_k3.json"))["p_value"]
        _tri(ax, z["scores"][:, :2], [SC[x] for x in st], arcs, labs,
             f"{t} alone ({len(st)} cells): triangle p = {p:.3f}")
    _legend(axes[1], SC, loc="upper left")
    fig.tight_layout(); fig.savefig(OUT / "within_triangles.png"); plt.close(fig)


def deck_v2():
    pair("all_main", ["V1 multiome (nuclei) corner", "V2 dividing epithelial", "V3 late EMT"], "assay", "all",
         "A_pair.png", ["Version A, all cells: coloured by EMT state", "Same triangle: coloured by technology"])
    pair("all_scrna", ["V1 long-gene axis", "V2 late EMT", "V3 dividing epithelial"], "sample", "scRNA",
         "B_pair.png", ["Version B, scRNA cells: coloured by EMT state", "Same triangle: coloured by tumour"])
    pair("ks_main", ["V1 epithelial", "V2 hybrid + late", "V3 early hybrid"], "sample", "all",
         "C_pair.png", ["Version C, KS EMT genes: coloured by EMT state", "Same triangle: coloured by tumour"])
    pair("all_scrna_noV1", ["Late EMT", "Proliferating hybrid", "Epithelial / early hybrid"], "sample", "scRNA",
         "noV1_pair.png", ["V1 axis removed: coloured by EMT state", "Same triangle: coloured by tumour"])
    within_triangles()


if __name__ == "__main__":
    deck_v2()
    print("deck v2 figures done")


def groves_slide():
    """Version B Groves closest-5% figure, legends placed below the panels (slide 4)."""
    run = RES / "all_scrna"
    c = pd.read_csv(run / "groves5_k3" / "cells_groves5.csv.gz", index_col=0)
    z = np.load(run / "pca.npz")
    arcs = np.load(run / "archetypes_k3.npy")
    ref = c["reference_cell"].values
    P = z["scores_all"][ref, :2]
    R = c[ref]
    cut = [R.loc[R["groves5"] == v, f"dist_{v}"].max() for v in ["V1", "V2", "V3"]]
    gen = (R["groves5"] == "generalist").values
    vcol = {"V1": "#7b3294", "V2": "#d7301f", "V3": "#1a9850"}
    names = {"V1": "V1 long-gene axis", "V2": "V2 late EMT", "V3": "V3 dividing epithelial"}
    fig, axes = plt.subplots(1, 3, figsize=(22, 6.6))
    panels = [("Specialists by corner", R["groves5"].values, vcol),
              ("Specialists by EMT state", R["cluster_annotation"].values, SC),
              ("Specialists by tumour", R["sample"].values, {k: v for k, v in TCOL.items() if k in set(R["sample"])})]
    for ax, (title, lab, cmap) in zip(axes, panels):
        ax.scatter(P[gen, 0], P[gen, 1], s=2, c="#d4d4d4", rasterized=True)
        s = ~gen
        ax.scatter(P[s, 0], P[s, 1], s=4, c=[cmap[x] for x in lab[s]], rasterized=True)
        t = np.vstack([arcs, arcs[:1]])
        ax.plot(t[:, 0], t[:, 1], "-", color="#1a1a1a", lw=1.8)
        for j, v in enumerate(["V1", "V2", "V3"]):
            ax.add_patch(plt.Circle(arcs[j, :2], cut[j], fill=False, ls="--", lw=1, color="#555555"))
            ax.annotate(names[v], arcs[j, :2], fontsize=13, fontweight="bold", xytext=(6, 6),
                        textcoords="offset points",
                        bbox=dict(boxstyle="round,pad=0.2", fc="#ffffff", ec="#aaaaaa", alpha=0.9))
        ax.set_xlim(-86, 26); ax.set_ylim(-36, 32)
        ax.set_aspect("equal", adjustable="box"); ax.set_title(title, fontsize=16)
        ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        keys = [k for k in cmap if (lab[s] == k).any()]
        for k in keys:
            n = int((lab[s] == k).sum())
            ax.scatter([], [], s=60, color=cmap[k], label=f"{k} ({n})")
        ax.scatter([], [], s=60, color="#d4d4d4", label=f"generalist ({int(gen.sum())})")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=3 if len(keys) <= 4 else 4,
                  fontsize=12, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "groves_slide.png")
    plt.close(fig)


if __name__ == "__main__":
    groves_slide()
    print("groves slide figure done")
