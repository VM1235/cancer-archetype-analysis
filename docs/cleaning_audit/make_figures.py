#!/usr/bin/env python3
"""Figures for the cleaning-audit section of the report. From Project_1/: .venv/bin/python docs/cleaning_audit/make_figures.py"""
import json, os, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
R = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
os.chdir(R); OUT = "docs/cleaning_audit/figures"
INK, MUTE, BLUE, ORANGE = "#1B2A3A", "#5B6770", "#2F6DB5", "#D9822B"
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": MUTE, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK, "font.size": 10})
COOK = "Mouse_SCC/Cook/results/within_line"; GSE = "Breast Cancer/GSE114687_pseudospace/results/within_unit"
T4 = "Breast Cancer/Mouse_4T1_GSE303505/results/within_unit"; PD = "Mouse_SCC/Mouse_PDAC_GSE173958/results/within_unit"
def P(path):
    try: return json.load(open(path))["p_value"]
    except Exception: return np.nan

# ---- Fig A1: share of cells kept by the 10,000-UMI cutoff (measured on the h5ad files)
d = [("SCC", 71, 30, "MU8"), ("Cook", 89, 83, "OVCA420"), ("PDAC", 80, 52, "iv-100cell"), ("4T1", 69, 21, "clone E"), ("GSE114687", 45, 35, "HuMEC Mock")]
d = sorted(d, key=lambda r: r[1])
fig, ax = plt.subplots(figsize=(7.6, 3.6)); y = np.arange(len(d))
ax.barh(y, [r[1] for r in d], height=0.5, color=BLUE, label="all cells in the dataset")
ax.scatter([r[2] for r in d], y, s=70, color=ORANGE, zorder=3, edgecolor="white", linewidth=1.5, label="lowest single sample")
for i, r in enumerate(d):
    ax.text(r[1] + 1.5, i, f"{r[1]}%", va="center", fontsize=9); ax.text(r[2] + 2.5, i + 0.27, f"{r[2]}% ({r[3]})", va="center", fontsize=8, color=INK)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in d]); ax.set_xlim(0, 112); ax.set_xlabel("% of cells with at least 10,000 UMI")
ax.set_title("The same 10,000-UMI cutoff keeps very different shares of cells", loc="left", fontsize=11, fontweight="bold")
ax.legend(frameon=False, loc="lower right", fontsize=9); [ax.spines[s].set_visible(False) for s in ("top", "right")]
fig.tight_layout(); fig.savefig(f"{OUT}/figA1_umi_keep_fraction.png", dpi=170); plt.close(fig)

# ---- Fig A2: before vs after the human technical-gene fix (Cook + GSE114687)
pts = []
for gs in ("all", "ks"):
    for l in ("A549", "DU145", "MCF7", "OVCA420"):
        for k in (3, 4): pts.append(("Cook", P(f"{COOK}/{gs}_{l}_clean/fit_k{k}.json"), P(f"{COOK}/{gs}_{l}_clean_ht/fit_k{k}.json")))
    for l in ("A549",): pts.append(("Cook", P(f"{COOK}/{gs}_{l}_clean_tbal/fit_k3.json"), P(f"{COOK}/{gs}_{l}_clean_tbal_ht/fit_k3.json")))
    for u in ("MCF10A_Mock", "HuMEC_Mock", "MCF10A_TGFB", "HuMEC_TGFB", "MCF10A", "HuMEC"):
        for k in (3, 4): pts.append(("GSE114687", P(f"{GSE}/{gs}_{u}_clean/fit_k{k}.json"), P(f"{GSE}/{gs}_{u}_clean_ht/fit_k{k}.json")))
df = pd.DataFrame(pts, columns=["ds", "old", "new"]).dropna()
fig, ax = plt.subplots(figsize=(5.2, 5.0)); lo = 0.0035
ax.plot([lo, 1.2], [lo, 1.2], color=MUTE, lw=1); ax.axvline(0.05, color=MUTE, lw=0.8, ls=":"); ax.axhline(0.05, color=MUTE, lw=0.8, ls=":")
rng = np.random.default_rng(1)
for ds, c in (("Cook", ORANGE), ("GSE114687", BLUE)):
    s = df[df.ds == ds]; ax.scatter(s.old * np.exp(rng.normal(0, 0.03, len(s))), s.new * np.exp(rng.normal(0, 0.03, len(s))), s=38, color=c, alpha=0.85, edgecolor="white", linewidth=0.8, label=f"{ds} ({len(s)} fits)")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, 1.25); ax.set_ylim(lo, 1.25)
ax.set_xlabel("p before (filter matched no gene)"); ax.set_ylabel("p after (MALAT1, NEAT1, haemoglobin dropped)")
ax.text(0.052, 0.0042, "p = 0.05", fontsize=8, color=MUTE); ax.legend(frameon=False, loc="upper left", fontsize=9)
ax.set_title("Dropping the human technical genes\nmoves no fit across p = 0.05", loc="left", fontsize=10.5, fontweight="bold")
[ax.spines[s].set_visible(False) for s in ("top", "right")]; fig.tight_layout(); fig.savefig(f"{OUT}/figA2_human_gene_fix.png", dpi=170); plt.close(fig)

def heat(ax, M, rows, cols, title, note):
    cmap = plt.get_cmap("Blues_r"); norm = LogNorm(vmin=0.005, vmax=1.0)
    ax.imshow(np.clip(M, 0.005, 1), cmap=cmap, norm=norm, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isnan(v): ax.text(j, i, "n/a", ha="center", va="center", fontsize=8, color=MUTE); continue
            ax.text(j, i, f"{v:.3f}" + ("*" if v < 0.05 else ""), ha="center", va="center", fontsize=8.5, color="white" if v < 0.1 else INK, fontweight="bold" if v < 0.05 else "normal")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols); ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=9)
    ax.xaxis.tick_top(); ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold", pad=26); ax.text(0, -0.09, note, transform=ax.transAxes, fontsize=8, color=MUTE)

# ---- Fig A3: GSE114687 UMI sweep
U = ("MCF10A_Mock", "HuMEC_Mock", "MCF10A_TGFB", "HuMEC_TGFB", "MCF10A", "HuMEC"); rows, M = [], []
for gs in ("all", "ks"):
    for u in U:
        rows.append(f"{gs} · {u}"); M.append([P(f"{GSE}/{gs}_{u}_clean_ht{s}/fit_k3.json") for s in ("_t5000", "", "_t15000")])
fig, ax = plt.subplots(figsize=(6.4, 6.2)); heat(ax, np.array(M), rows, ["5,000", "10,000", "15,000"], "GSE114687: within-sample triangle p (k = 3) by UMI cutoff", "* p < 0.05. Darker = smaller p. 15,000 keeps only 284 to 1,838 cells per fit.")
fig.tight_layout(); fig.savefig(f"{OUT}/figA3_gse_umi_sweep.png", dpi=170); plt.close(fig)

# ---- Fig A4: 4T1 contamination cutoff
rows, M = [], []
for gs in ("all", "ks"):
    for u in ("EM1", "EM2", "EM3", "M", "4T1"):
        rows.append(f"{gs} · {u}"); M.append([P(f"{T4}/{gs}_{u}_clean_dc{s}/fit_k3.json") for s in ("_c0.2", "", "_c0.5")])
fig, ax = plt.subplots(figsize=(6.0, 4.6)); heat(ax, np.array(M), rows, ["0.2", "0.3 (used)", "0.5"], "4T1: within-clone triangle p (k = 3) by contamination cutoff", "Cleaned (>= 10,000 UMI). * p < 0.05. Clone E is excluded (142 cells after cleaning).")
fig.tight_layout(); fig.savefig(f"{OUT}/figA4_4t1_cutoff.png", dpi=170); plt.close(fig)

# ---- Fig A5: Hallmark scores at corners, SCC and 4T1
h = pd.read_csv("docs/cleaning_audit/hallmark_scc_4t1.csv"); h = h[~h.fit.str.contains("light-QC")]
sets = ["Hypoxia", "EMT", "G2-M", "E2F", "TGF-b", "ApicalJunction", "KRAS-up"]; setlab = ["Hypoxia", "EMT", "G2-M", "E2F", "TGF-b", "Apical\njunction", "KRAS\nup"]
labels = [f"{r.fit.replace(' all genes', '').replace('cleaned ', '')}  {r.corner}" + (f" ({int(r.cycling_pct)}% cyc.)" if not pd.isna(r.cycling_pct) else "") for r in h.itertuples()]
V = h[sets].values
fig, ax = plt.subplots(figsize=(8.6, 0.27 * len(h) + 2.0)); im = ax.imshow(V, cmap="RdBu_r", norm=TwoSlopeNorm(0, -0.9, 0.9), aspect="auto")
for i in range(V.shape[0]):
    for j in range(V.shape[1]): ax.text(j, i, f"{V[i, j]:+.2f}", ha="center", va="center", fontsize=7.5, color="white" if abs(V[i, j]) > 0.5 else INK)
ax.set_xticks(range(len(sets))); ax.set_xticklabels(setlab, fontsize=9); ax.xaxis.tick_top(); ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=8)
ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
plt.colorbar(im, ax=ax, shrink=0.4, pad=0.02).set_label("Hallmark score (z, nearest 10% minus rest)")
ax.set_title("Hallmark check of corner names: mouse SCC and 4T1", loc="left", fontsize=10.5, fontweight="bold", pad=30)
fig.tight_layout(); fig.savefig(f"{OUT}/figA5_hallmark_scc_4t1.png", dpi=170, bbox_inches="tight"); plt.close(fig)
print("figures written:", sorted(os.listdir(OUT)))
