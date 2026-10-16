#!/usr/bin/env python3
"""Step 13: slide figures for GSE114687 -> results/mail_figures/.
fig1 EMT score by sample and colony position
fig2 pooled permutation p per k (KS vs all genes)
fig3 within-unit p-values, light QC vs cleaned
fig4 cleaned all-genes k=3 triangles per sample, coloured by colony position, corners named
fig5 generalist share per EMT state (gap and closest-5% rules)
fig6 Hallmark scores of every cleaned corner
Usage (from Project_1/): .venv/bin/python "Breast Cancer/GSE114687_pseudospace/codes/13_mail_figures.py"
"""
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scc_common as C

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12})
OUT = C.RESULTS / "mail_figures"; OUT.mkdir(exist_ok=True)
SAMPLES = ["MCF10A_Mock", "MCF10A_TGFB", "HuMEC_Mock", "HuMEC_TGFB"]
PCOL = {"inner": "#2E5E8C", "outer": "#B5473A"}
obs = ad.read_h5ad(C.PROC / "expr_tumour.h5ad", backed="r").obs

# fig1
fig, ax = plt.subplots(figsize=(8.6, 4.8)); pos = []; lab = []
for i, s in enumerate(SAMPLES):
    for j, p in enumerate(["inner", "outer"]):
        v = obs.loc[(obs["sample"] == s) & (obs["position"] == p), "KS_EMT"].values
        b = ax.boxplot(v, positions=[i * 3 + j], widths=0.8, patch_artist=True, showfliers=False)
        b["boxes"][0].set_facecolor(PCOL[p]); b["boxes"][0].set_alpha(0.85)
        for m in b["medians"]: m.set_color("k")
        pos.append(i * 3 + j + 0.0)
ax.set_xticks([i * 3 + 0.5 for i in range(4)]); ax.set_xticklabels(SAMPLES)
ax.set_ylabel("EMT score (KS_EMT)"); ax.spines[["top", "right"]].set_visible(False)
for p, c in PCOL.items(): ax.scatter([], [], marker="s", s=60, color=c, label=p + " colony")
ax.legend(frameon=False); ax.set_title("GSE114687: TGF-β shifts the EMT score; inner vs outer colony differs only a little")
fig.savefig(OUT / "fig1_emt_score_by_sample.png", dpi=170, bbox_inches="tight"); plt.close(fig)

# fig2
fig, ax = plt.subplots(figsize=(7.6, 4.8))
for gs, lab_, c in [("ks", "KS genes", "#7A7A7A"), ("all", "All genes", "#2E5E8C")]:
    t = pd.read_csv(C.RESULTS / f"{gs}_main" / "fit_summary.csv")
    ax.plot(t["k"], t["p_value"], "o-", color=c, ms=7, label=lab_)
ax.axhline(0.01, ls="--", c="k", lw=1); ax.set_yscale("log"); ax.set_ylim(0.003, 1.3); ax.set_xticks(range(3, 8))
ax.set_xlabel("number of corners (k)"); ax.set_ylabel("permutation p (200 shuffles; floor 0.005)")
ax.set_title("Pooled over the four samples: only k >= 5-6 pass, and the corners are samples", fontsize=11)
ax.legend(frameon=False); ax.spines[["top", "right"]].set_visible(False)
fig.savefig(OUT / "fig2_pooled_permutation_p.png", dpi=170, bbox_inches="tight"); plt.close(fig)

# fig3
sm = pd.read_csv(C.RESULTS / "within_unit_summary.csv")
units = ["MCF10A", "MCF10A_Mock", "MCF10A_TGFB", "HuMEC", "HuMEC_Mock", "HuMEC_TGFB"]
cols_ = [("all", 3), ("all", 4), ("ks", 3), ("ks", 4)]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.9), constrained_layout=True)
for ax, var, name in [(axes[0], "light", "light QC"), (axes[1], "clean", "cleaned (>= 10,000 UMI)")]:
    M = np.full((6, 4), np.nan)
    for i, u in enumerate(units):
        for j, (g, k) in enumerate(cols_):
            r = sm[(sm.variant == var) & (sm.unit == u) & (sm.geneset == g) & (sm.k == k)]
            if len(r): M[i, j] = r["p_value"].values[0]
    im = ax.imshow(M, cmap="RdYlBu", norm=LogNorm(vmin=0.005, vmax=1), aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(["all\nk=3", "all\nk=4", "KS\nk=3", "KS\nk=4"])
    ax.set_yticks(range(6)); ax.set_yticklabels(units)
    for i in range(6):
        for j in range(4):
            if not np.isnan(M[i, j]): ax.text(j, i, f"{M[i, j]:.3f}" + ("*" if M[i, j] < 0.05 else ""), ha="center", va="center", fontsize=10)
    ax.set_title(name, fontsize=11)
fig.colorbar(im, ax=axes, label="permutation p (log)", shrink=0.85)
fig.suptitle("Shape inside single samples: * = p < 0.05 (light QC 11 of 24, cleaned 18 of 24)", fontsize=12, fontweight="bold")
fig.savefig(OUT / "fig3_within_unit_p.png", dpi=170); plt.close(fig)

# fig4
cd = pd.read_csv(C.RESULTS / "within_unit" / "corner_description_clean.csv")
def names(unit):
    r = cd[cd.corner.str.startswith(unit + "|all|k3|")].reset_index(drop=True); r["V"] = r["corner"].str.split("|").str[-1]
    nm = {}
    for _, x in r.iterrows():
        if x["pct_cycling"] >= 50: nm[x["V"]] = "cycling"
    rest = r[~r["V"].isin(nm)].sort_values("KS_EMT")
    nm[rest.iloc[0]["V"]] = "differentiated epithelial"
    nm[rest.iloc[1]["V"]] = "EMT / matrix" if "TGFB" in unit else "basal / partial"
    return nm
fig, axes = plt.subplots(2, 2, figsize=(13, 10.5), constrained_layout=True)
for ax, s in zip(axes.ravel(), SAMPLES):
    d = C.RESULTS / "within_unit" / f"all_{s}_clean"
    z = np.load(d / "pca.npz"); arcs = np.load(d / "archetypes_k3.npy"); S = z["scores"]
    p = obs["position"].astype(str).values[z["rows"]]
    idx = np.random.default_rng(0).permutation(len(p))
    ax.scatter(S[idx, 0], S[idx, 1], s=4, c=[PCOL[x] for x in p[idx]], alpha=0.7, rasterized=True)
    T = np.vstack([arcs[:, :2], arcs[:1, :2]]); ax.plot(T[:, 0], T[:, 1], "k-", lw=1.8)
    nm = names(s)
    for v in range(3):
        ax.scatter(*arcs[v, :2], s=55, c="k", zorder=6)
        ax.annotate(f"V{v + 1}: {nm.get(f'V{v + 1}', '')}", arcs[v, :2], fontsize=11, fontweight="bold", xytext=(7, 7), textcoords="offset points",
                    bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="0.7"), zorder=7)
    pa = sm[(sm.variant == "clean") & (sm.unit == s) & (sm.geneset == "all") & (sm.k == 3)]["p_value"].values[0]
    ax.set_title(f"{s} (n={len(p)}), k=3, p={pa:.3f}"); ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.spines[["top", "right"]].set_visible(False)
for p, c in PCOL.items(): axes[0, 0].scatter([], [], s=30, color=c, label=p + " colony")
axes[0, 0].legend(loc="lower left", frameon=True, fontsize=10)
fig.suptitle("Cleaned, all genes: the same three corners recur in every sample", fontsize=14, fontweight="bold")
fig.savefig(OUT / "fig4_triangles_by_sample.png", dpi=170); plt.close(fig)

# fig5
g = pd.read_csv(C.RESULTS / "within_unit" / "generalist_share_by_state.csv")
ST = ["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"]
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), constrained_layout=True, sharey=True)
fits = list(g["fit"].unique()); w = 0.8 / len(fits); cm = plt.get_cmap("tab10")
for ax, rule, name in [(axes[0], "gap", "purity-gap rule (>= 0.35)"), (axes[1], "groves5", "closest-5% rule")]:
    for i, f in enumerate(fits):
        v = g[(g.fit == f) & (g.rule == rule)][ST].values[0]
        ax.bar(np.arange(4) + i * w - 0.4 + w / 2, v, w, color=cm(i), label=f)
    ax.set_xticks(range(4)); ax.set_xticklabels(ST); ax.set_title(name); ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("% of the state's cells that are generalists"); fig.legend(*axes[1].get_legend_handles_labels(), loc="outside lower center", ncol=4, fontsize=9, frameon=False)
fig.suptitle("Hybrid states are the most generalist in most fits (cleaned k=3 fits)", fontsize=13, fontweight="bold")
fig.savefig(OUT / "fig5_generalist_share_by_state.png", dpi=170); plt.close(fig)

# fig6
h = pd.read_csv(C.RESULTS / "within_unit" / "hallmark_corners_clean.csv")
sets = ["E2F", "G2-M", "EMT", "TGF-b", "Hypoxia", "ApicalJunction", "KRAS-up"]
M = h[sets].values
fig, ax = plt.subplots(figsize=(8.6, 9), constrained_layout=True)
im = ax.imshow(M, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
ax.set_xticks(range(len(sets))); ax.set_xticklabels(sets, rotation=30, ha="right")
ax.set_yticks(range(len(h))); ax.set_yticklabels([f"{a} {b} ({int(c)}% cycling)" for a, b, c in zip(h["fit"], h["corner"], h["cycling"])], fontsize=9)
for i in range(M.shape[0]):
    for j in range(M.shape[1]): ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=8)
fig.colorbar(im, ax=ax, label="Hallmark score: 10% nearest minus the rest (z)", shrink=0.6)
ax.set_title("Hallmark gene sets at each cleaned corner")
fig.savefig(OUT / "fig6_hallmark_corners.png", dpi=170); plt.close(fig)
print("ok")
