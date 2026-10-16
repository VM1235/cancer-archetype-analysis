#!/usr/bin/env python3
"""Do the corner names survive an independent gene list? Score every corner with MSigDB Hallmark sets.

Hallmark gene sets (human symbols, Enrichr copy of MSigDB_Hallmark_2020) were downloaded to
Mouse_SCC/hallmark_data/ and are only read as text. Sets used: Hypoxia, Epithelial Mesenchymal Transition,
G2-M Checkpoint, E2F Targets, TGF-beta Signaling. Mouse data: symbol capitalised (HIF1A -> Hif1a).
For each fitted simplex, the 10% of cells nearest each corner are scored (mean z of the set's genes,
z within the unit). My earlier name for each corner is compared with the Hallmark set that is highest there.
Usage (from Project_1/):  .venv/bin/python Mouse_SCC/hallmark_check/hallmark_corners.py
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = Path("/Users/apple/Desktop/thesis_iisc/Project_1/Mouse_SCC")
OUT = ROOT / "hallmark_check"
sets_all = {}
for line in open(ROOT / "hallmark_data" / "MSigDB_Hallmark_2020.txt"):
    f = line.rstrip("\n").split("\t"); sets_all[f[0].strip()] = [g for g in f[2:] if g]
USE = {"Hypoxia": "Hypoxia", "EMT": "Epithelial Mesenchymal Transition", "G2-M checkpoint": "G2-M Checkpoint",
       "E2F targets": "E2F Targets", "TGF-β signalling": "TGF-beta Signaling"}
def score(a, rows, mouse):
    out = {}
    X = a.X[rows]
    for short, full in USE.items():
        genes = [g.capitalize() if mouse else g for g in sets_all[full]]
        idx = [a.var_names.get_loc(g) for g in genes if g in a.var_names]
        M = X[:, idx].toarray(); M = (M - M.mean(0)) / np.where(M.std(0) > 0, M.std(0), 1)
        out[short] = M.mean(1)
    return out
def corners(a, rows, P, arcs, names, label, mouse, recs):
    S = score(a, rows, mouse); k = arcs.shape[0]
    for v in range(k):
        d = np.linalg.norm(P[:, :k - 1] - arcs[v, :k - 1], axis=1); c = d <= np.quantile(d, .10)
        r = {"fit": label, "corner": f"V{v + 1}", "my_name": names[v]}
        r.update({n: round(float(S[n][c].mean()), 2) for n in S}); recs.append(r)
recs = []
# ---- Simeonov (mouse symbols)
a = ad.read_h5ad(ROOT / "Mouse_PDAC_GSE173958/data/processed/expr_tumour.h5ad")
W = ROOT / "Mouse_PDAC_GSE173958/results/within_unit"
lib_roles = {"M1-PT": {"V3": "epithelial", "V1": "proliferation", "V2": "EMT signalling"}, "M1-SS": {"V2": "epithelial", "V1": "proliferation", "V3": "EMT signalling"},
 "M1-Liver": {"V2": "epithelial", "V1": "proliferation", "V3": "EMT signalling"}, "M2-Liver": {"V3": "epithelial", "V1": "proliferation", "V2": "hypoxia"},
 "M2-Met": {"V2": "epithelial", "V1": "proliferation", "V3": "hypoxia"}, "M2-PT": {"V2": "epithelial", "V1": "proliferation", "V3": "hypoxia"}}
for u, rm in lib_roles.items():
    z = np.load(W / f"all_{u}_clean/pca.npz"); arcs = np.load(W / f"all_{u}_clean/archetypes_k3.npy")
    corners(a, z["rows"], z["scores"][:, :2], arcs, [rm[f"V{i + 1}"] for i in range(3)], f"Simeonov {u}", True, recs)
extra = {"M2": (3, ["epithelial", "proliferation", "hypoxia"]), "M1": (3, ["proliferation", "epithelial", "EMT signalling"]),
         "invivo": (4, ["hypoxia-like", "proliferation", "EMT signalling", "epithelial"])}
for u, (k, nm) in extra.items():
    z = np.load(W / f"all_{u}_clean/pca.npz"); arcs = np.load(W / f"all_{u}_clean/archetypes_k{k}.npy")
    lab = {"M2": "Simeonov mouse 2", "M1": "Simeonov mouse 1", "invivo": "Simeonov both mice (k=4)"}[u]
    corners(a, z["rows"], z["scores"][:, :k - 1], arcs, nm, lab, True, recs)
del a
# ---- Cook (human symbols); my names from the earlier programme scores
sys.path.insert(0, str(ROOT / "Cook/codes")); import cook_common as C
a = C.load_adata(); prog = pd.read_csv(C.RESULTS / "within_line/corner_programmes.csv")
for L in C.LINES:
    z = np.load(C.RESULTS / f"within_line/all_{L}_clean/pca.npz"); arcs = np.load(C.RESULTS / f"within_line/all_{L}_clean/archetypes_k3.npy")
    p = prog[prog.line == L].reset_index(drop=True)
    pr = int(p["proliferation"].idxmax()); rest = [i for i in range(3) if i != pr]; tg = max(rest, key=lambda i: p.loc[i, "TGF-β response"])
    names = ["other"] * 3; names[pr] = "proliferation"; names[tg] = "TGF-β response"
    corners(a, z["rows"], z["scores"][:, :2], arcs, names, f"Cook {L}", False, recs)
df = pd.DataFrame(recs); df.to_csv(OUT / "hallmark_corner_scores.csv", index=False)
pd.set_option("display.width", 220); print(df.to_string(index=False))
# ---- agreement
expect = {"proliferation": ["G2-M checkpoint", "E2F targets"], "hypoxia": ["Hypoxia"], "hypoxia-like": ["Hypoxia"],
          "EMT signalling": ["EMT", "TGF-β signalling"], "TGF-β response": ["EMT", "TGF-β signalling"]}
sets = list(USE)
rows = []
for _, r in df.iterrows():
    if r.my_name in expect:
        top = max(sets, key=lambda s: r[s]); rows.append(dict(fit=r.fit, corner=r.corner, my_name=r.my_name, top_set=top, agree=top in expect[r.my_name],
                  expected_best=max(r[s] for s in expect[r.my_name])))
ag = pd.DataFrame(rows); ag.to_csv(OUT / "hallmark_agreement.csv", index=False)
print("\nAgreement (is the highest Hallmark set one of the expected ones?):"); print(ag.to_string(index=False))
print(ag.groupby("my_name").agree.agg(["sum", "count"]))
# ---- heatmap
labels = [f"{r.fit}  V{r.corner[1:]}  ({r.my_name})" for r in df.itertuples()]
H = df[sets].values
fig, ax = plt.subplots(figsize=(10, 0.27 * len(df) + 2.2)); im = ax.imshow(H, cmap="RdBu_r", vmin=-1.6, vmax=1.6, aspect="auto")
ax.set_xticks(range(len(sets))); ax.set_xticklabels([s.replace(" ", "\n") for s in sets], fontsize=10); ax.xaxis.tick_top()
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=8)
for i in range(H.shape[0]):
    for j in range(H.shape[1]): ax.text(j, i, f"{H[i, j]:+.1f}", ha="center", va="center", fontsize=7, color="white" if abs(H[i, j]) > 0.9 else "#1B2A3A")
ax.tick_params(length=0); [s.set_visible(False) for s in ax.spines.values()]
plt.colorbar(im, ax=ax, shrink=.4, pad=.02).set_label("Hallmark score (z)")
fig.suptitle("Corner names checked against MSigDB Hallmark gene sets", y=0.995, fontsize=11, fontweight="bold")
fig.savefig(OUT / "hallmark_corner_check.png", dpi=170, bbox_inches="tight")
