#!/usr/bin/env python3
"""Step 13: slide/e-mail figures for the 4T1 clone analysis -> results/mail_figures/.

fig1 clone_emt_spectrum        KS_EMT score per clone (cells, boxes)
fig2 pooled_permutation_p      permutation p per k for the four pooled fits
fig3 pooled_triangle_cleaned   cleaned all-genes k=3 triangle: clone, phase, EMT score
fig4 pooled_tetrahedron_lightQC  light-QC all-genes k=4, PC1-PC2: clone, phase
fig5 within_clone_p            within-clone p-values, light QC vs cleaned
Usage (from Project_1/): .venv/bin/python "Breast Cancer/Mouse_4T1_GSE303505/codes/13_mail_figures.py"
"""
import json
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scc_common as C

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12})
OUT = C.RESULTS / "mail_figures"; OUT.mkdir(exist_ok=True)
CL = ["E", "EM1", "EM2", "EM3", "M", "4T1"]
COL = {"E": "#2E5E8C", "EM1": "#5C9BC4", "EM2": "#8FB996", "EM3": "#E3A04F", "M": "#B5473A", "4T1": "#7A7A7A"}
PH = {"G1": "#9E9E9E", "S": "#E3A04F", "G2M": "#B5473A"}
full = ad.read_h5ad(C.PROC / "expr_tumour.h5ad", backed="r").obs
clean = ad.read_h5ad(C.PROC / "expr_tumour_clean.h5ad", backed="r").obs

def loadfit(tag, k, obs):
    d = C.RESULTS / f"all_{tag}"
    z = np.load(d / "pca.npz"); arcs = np.load(d / f"archetypes_k{k}.npy")
    o = obs.loc[z["barcodes"]]
    return z["scores_all"], arcs, o

def draw(ax, S, arcs, k, color):
    idx = np.random.default_rng(0).permutation(len(S))
    ax.scatter(S[idx, 0], S[idx, 1], s=3, c=np.asarray(color, dtype=object)[idx] if not isinstance(color, np.ndarray) or color.dtype == object else color[idx], rasterized=True)
    ord_ = [0, 1, 2] if k == 3 else None
    if k == 3:
        T = np.vstack([arcs[:, :2], arcs[:1, :2]]); ax.plot(T[:, 0], T[:, 1], "k-", lw=1.6, zorder=4)
    else:
        for i in range(k):
            for j in range(i + 1, k):
                ax.plot(arcs[[i, j], 0], arcs[[i, j], 1], "k-", lw=1.0, zorder=4)
    ax.scatter(arcs[:, 0], arcs[:, 1], s=55, c="k", zorder=6)
    for v in range(k):
        ax.annotate(f"V{v + 1}", arcs[v, :2], fontsize=13, fontweight="bold", xytext=(7, 7), textcoords="offset points",
                    bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="0.7"), zorder=7)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.spines[["top", "right"]].set_visible(False)

def legend(ax, d, title):
    for k_, c_ in d.items(): ax.scatter([], [], s=35, color=c_, label=k_)
    ax.legend(loc="best", fontsize=10, framealpha=.9, title=title, markerscale=1)

# ---- fig1
fig, ax = plt.subplots(figsize=(7.5, 4.8))
data = [full.loc[full["sample"] == c, "KS_EMT"].values for c in CL]
bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.6)
for b, c in zip(bp["boxes"], CL): b.set_facecolor(COL[c]); b.set_alpha(0.85)
for m in bp["medians"]: m.set_color("k")
ax.set_xticklabels([f"{c}\n(n={len(d)})" for c, d in zip(CL, data)])
ax.set_ylabel("EMT score (KS_EMT; higher = more mesenchymal)")
ax.set_title("4T1 clones along the EMT score: E is separate, the rest overlap", fontsize=12)
ax.spines[["top", "right"]].set_visible(False)
fig.savefig(OUT / "fig1_clone_emt_spectrum.png", dpi=170, bbox_inches="tight"); plt.close(fig)

# ---- fig2
rows = []
for gs, tag, lab in [("ks", "main", "KS genes, light QC"), ("all", "main", "All genes, light QC"),
                     ("ks", "clean", "KS genes, cleaned"), ("all", "clean", "All genes, cleaned")]:
    s = pd.read_csv(C.RESULTS / f"{gs}_{tag}" / "fit_summary.csv")
    for _, r in s.iterrows(): rows.append((lab, int(r["k"]), float(r["p_value"])))
t = pd.DataFrame(rows, columns=["fit", "k", "p"])
fig, ax = plt.subplots(figsize=(8.4, 5.2))
cols = ["#7A7A7A", "#2E5E8C", "#B8B8B8", "#B5473A"]
for (lab, g), c in zip(t.groupby("fit", sort=False), cols):
    ax.plot(g["k"] + (0.06 * cols.index(c) - 0.09), g["p"], "o-", color=c, label=lab, ms=7)
ax.axhline(0.01, ls="--", c="k", lw=1); ax.text(2.75, 0.0115, "p = 0.01", ha="left", fontsize=10)
ax.set_yscale("log"); ax.set_ylim(0.003, 1.3); ax.set_xticks(range(3, 8))
ax.set_xlabel("number of corners (k)"); ax.set_ylabel("permutation p (200 shuffles; floor 0.005)")
ax.set_title("Pooled 4T1 fits: is the shape better than shuffled data?", fontsize=12)
ax.legend(fontsize=10, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False); ax.spines[["top", "right"]].set_visible(False)
fig.savefig(OUT / "fig2_pooled_permutation_p.png", dpi=170, bbox_inches="tight"); plt.close(fig)

# ---- fig3
S, arcs, o = loadfit("clean", 3, clean)
fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), constrained_layout=True)
draw(axes[0], S, arcs, 3, np.array([COL[x] for x in o["sample"].astype(str)], dtype=object)); legend(axes[0], COL, "clone")
draw(axes[1], S, arcs, 3, np.array([PH[x] for x in o["phase"].astype(str)], dtype=object)); legend(axes[1], PH, "phase")
idx = np.random.default_rng(0).permutation(len(S))
sc = axes[2].scatter(S[idx, 0], S[idx, 1], s=3, c=o["KS_EMT"].values[idx], cmap="coolwarm", vmin=-1.2, vmax=1.2, rasterized=True)
T = np.vstack([arcs[:, :2], arcs[:1, :2]]); axes[2].plot(T[:, 0], T[:, 1], "k-", lw=1.6)
axes[2].scatter(arcs[:, 0], arcs[:, 1], s=55, c="k", zorder=6)
for v in range(3): axes[2].annotate(f"V{v + 1}", arcs[v, :2], fontsize=13, fontweight="bold", xytext=(7, 7), textcoords="offset points")
axes[2].set_xlabel("PC1"); axes[2].set_ylabel("PC2"); axes[2].spines[["top", "right"]].set_visible(False)
fig.colorbar(sc, ax=axes[2], label="EMT score")
axes[0].set_title("by clone"); axes[1].set_title("by cell-cycle phase"); axes[2].set_title("by EMT score")
fig.suptitle("Cleaned pooled fit (all genes, k=3, p=0.010): V1 cycling, V2 epithelial-leaning, V3 = clone M", fontsize=14, fontweight="bold")
fig.savefig(OUT / "fig3_pooled_triangle_cleaned.png", dpi=170); plt.close(fig)

# ---- fig4
S, arcs, o = loadfit("main", 4, full)
fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True)
draw(axes[0], S, arcs, 4, np.array([COL[x] for x in o["sample"].astype(str)], dtype=object)); legend(axes[0], COL, "clone")
draw(axes[1], S, arcs, 4, np.array([PH[x] for x in o["phase"].astype(str)], dtype=object)); legend(axes[1], PH, "phase")
axes[0].set_title("by clone"); axes[1].set_title("by cell-cycle phase")
fig.suptitle("Light-QC pooled fit (all genes, k=4, p=0.005): E, 4T1, M and a cycling corner. View: PC1-PC2 of a 3-D shape", fontsize=12, fontweight="bold")
fig.savefig(OUT / "fig4_pooled_tetrahedron_lightQC.png", dpi=170); plt.close(fig)

# ---- fig5
raw = pd.read_csv(C.RESULTS / "within_clone_summary.csv"); dc = pd.read_csv(C.RESULTS / "within_clone_clean_dc_summary.csv")
cols_ = [("all", 3), ("all", 4), ("ks", 3), ("ks", 4)]
lab = ["all\nk=3", "all\nk=4", "KS\nk=3", "KS\nk=4"]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
for ax, df, name in [(axes[0], raw, "light QC (all tumour cells)"), (axes[1], dc, "cleaned (>=10k UMI, no immune/stromal)")]:
    M = np.full((len(CL), 4), np.nan)
    for i, c in enumerate(CL):
        for j, (g, k) in enumerate(cols_):
            r = df[(df.unit == c) & (df.geneset == g) & (df.k == k)]
            if len(r): M[i, j] = r["p_value"].values[0]
    im = ax.imshow(M, cmap="RdYlBu", norm=LogNorm(vmin=0.005, vmax=1), aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(lab); ax.set_yticks(range(6)); ax.set_yticklabels(CL)
    for i in range(6):
        for j in range(4):
            if np.isnan(M[i, j]): ax.text(j, i, "n/a", ha="center", va="center", fontsize=9, color="0.4")
            else: ax.text(j, i, f"{M[i, j]:.3f}" + ("*" if M[i, j] < 0.05 else ""), ha="center", va="center", fontsize=10)
    ax.set_title(name, fontsize=11)
fig.colorbar(im, ax=axes, label="permutation p (log scale)", shrink=0.85)
fig.suptitle("Shape inside single clones: * = p < 0.05. Cleaning removes the light-QC passes (E not testable after cleaning)", fontsize=12, fontweight="bold")
fig.savefig(OUT / "fig5_within_clone_p.png", dpi=170); plt.close(fig)
print("ok")

# ---- fig6: the fitted k=3 triangle inside each clone
fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
for ax, c in zip(axes.ravel(), CL):
    if c == "E":
        d = C.RESULTS / "within_unit" / "all_E"; note = "light QC, not cleaned"
        pk = json.load(open(C.RESULTS / "within_unit" / "ks_E" / "fit_k3.json"))["p_value"]
    else:
        d = C.RESULTS / "within_unit" / f"all_{c}_clean_dc"; note = "cleaned"
        pk = json.load(open(C.RESULTS / "within_unit" / f"ks_{c}_clean_dc" / "fit_k3.json"))["p_value"]
    z = np.load(d / "pca.npz"); arcs = np.load(d / "archetypes_k3.npy"); Sx = z["scores"]
    pa = json.load(open(d / "fit_k3.json"))["p_value"]
    ph = full["phase"].astype(str).values[z["rows"]]
    idx = np.random.default_rng(0).permutation(len(ph))
    ax.scatter(Sx[idx, 0], Sx[idx, 1], s=5, c=[PH[x] for x in ph[idx]], rasterized=True)
    T = np.vstack([arcs[:, :2], arcs[:1, :2]]); ax.plot(T[:, 0], T[:, 1], "k-", lw=1.8, zorder=4)
    ax.scatter(arcs[:, 0], arcs[:, 1], s=55, c="k", zorder=6)
    for v in range(3):
        ax.annotate(f"V{v + 1}", arcs[v, :2], fontsize=13, fontweight="bold", xytext=(7, 7), textcoords="offset points",
                    bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="0.7"), zorder=7)
    ax.set_title(f"{c} (n={len(ph)}, {note})\nk=3: p={pa:.3f} all genes, p={pk:.2f} KS genes", fontsize=12)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.spines[["top", "right"]].set_visible(False)
legend(axes[0, 0], PH, "cell-cycle phase")
fig.suptitle("Each clone's own k = 3 triangle (all genes), drawn whether or not it is significant", fontsize=14, fontweight="bold")
fig.savefig(OUT / "fig6_within_clone_triangles.png", dpi=170); plt.close(fig)

# ---- fig7: specialists vs generalists (Groves closest-5% rule)
SC = {"V1": "#2E5E8C", "V2": "#8FB996", "V3": "#B5473A", "V4": "#E3A04F", "generalist": "#D9D9D9"}
ST = ["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"]
fig, axes = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True)
for r, (tag, k, obs, name) in enumerate([("main", 4, full, "Light QC, k=4"), ("clean", 3, clean, "Cleaned, k=3")]):
    S, arcs, o = loadfit(tag, k, obs)
    g = pd.read_csv(C.RESULTS / f"all_{tag}" / f"groves5_k{k}" / "cells_groves5.csv.gz", index_col=0).loc[o.index]
    cls = g["groves5"].values
    ax = axes[r, 0]
    gen = cls == "generalist"
    ax.scatter(S[gen, 0], S[gen, 1], s=2, c=SC["generalist"], rasterized=True)
    for v in [x for x in SC if x != "generalist"]:
        m = cls == v
        if m.any(): ax.scatter(S[m, 0], S[m, 1], s=5, c=SC[v], label=f"{v} ({m.sum()})", rasterized=True)
    ax.scatter(arcs[:, 0], arcs[:, 1], s=40, c="k", zorder=6)
    for v in range(k): ax.annotate(f"V{v + 1}", arcs[v, :2], fontsize=12, fontweight="bold", xytext=(6, 6), textcoords="offset points")
    ax.legend(markerscale=3, fontsize=9, loc="upper left", title="closest 5% = specialists"); ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    ax.set_title(f"{name}: {100 * (~gen).mean():.0f}% specialists, rest generalists (grey)")
    ax.spines[["top", "right"]].set_visible(False)
    classes = [f"V{i + 1}" for i in range(k)] + ["generalist"]
    st = pd.crosstab(o["cluster_annotation"].astype(str).values, cls, normalize="index").reindex(ST)[classes] * 100
    ax = axes[r, 1]; left = np.zeros(len(ST))
    for c_ in classes:
        ax.barh(range(len(ST)), st[c_].values, left=left, color=SC[c_], label=c_); left += st[c_].values
    ax.set_yticks(range(len(ST))); ax.set_yticklabels(ST); ax.invert_yaxis(); ax.set_xlabel("% of the state's cells")
    ax.set_title("Where each EMT state goes"); ax.spines[["top", "right"]].set_visible(False)
    cl_ = pd.crosstab(cls, o["sample"].astype(str).values, normalize="index").reindex(classes)[CL] * 100
    ax = axes[r, 2]; left = np.zeros(len(classes))
    for c_ in CL:
        ax.barh(range(len(classes)), cl_[c_].values, left=left, color=COL[c_], label=c_); left += cl_[c_].values
    ax.set_yticks(range(len(classes))); ax.set_yticklabels(classes); ax.invert_yaxis(); ax.set_xlabel("% of the class")
    ax.set_title("Which clone makes up each class"); ax.spines[["top", "right"]].set_visible(False)
fig.legend(handles=[plt.Line2D([], [], marker="s", ls="", color=COL[c], label=c) for c in CL], loc="outside lower center", ncol=6, frameon=False, title="clone (right-hand panels)")
fig.suptitle("Specialists vs generalists (Groves 2022 closest-5% rule, all genes)", fontsize=14, fontweight="bold")
fig.savefig(OUT / "fig7_specialists_generalists.png", dpi=170); plt.close(fig)
print("ok2")
