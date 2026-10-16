#!/usr/bin/env python3
"""Step 9: do the corner names survive an independent gene list? (as Mouse_SCC/hallmark_check, here human symbols)

For each cleaned k=3 within-unit fit, the 10% of the unit's cells nearest each corner are scored on MSigDB Hallmark
sets (mean z of the set's genes, z within the unit): Hypoxia, EMT, G2-M checkpoint, E2F targets, TGF-beta signalling,
plus Apical junction and KRAS signalling up (differentiation / inflammation controls), and the top set is reported.
Usage (from Project_1/): .venv/bin/python "Breast Cancer/GSE114687_pseudospace/codes/09_hallmark_corners.py"
"""
import numpy as np, pandas as pd, scc_common as C

H = {}
for line in open(C.ROOT / "Mouse_SCC" / "hallmark_data" / "MSigDB_Hallmark_2020.txt"):
    f = line.rstrip("\n").split("\t"); H[f[0].strip()] = [g for g in f[2:] if g]
USE = {"Hypoxia": "Hypoxia", "EMT": "Epithelial Mesenchymal Transition", "G2-M": "G2-M Checkpoint", "E2F": "E2F Targets",
       "TGF-b": "TGF-beta Signaling", "ApicalJunction": "Apical Junction", "KRAS-up": "KRAS Signaling Up"}
JOBS = [("all", "MCF10A"), ("all", "MCF10A_Mock"), ("all", "MCF10A_TGFB"), ("all", "HuMEC_Mock"), ("all", "HuMEC_TGFB"),
        ("ks", "MCF10A_Mock"), ("ks", "HuMEC_Mock"), ("ks", "HuMEC_TGFB")]
a = C.load_adata(); k = 3
rows = []
for gs, unit in JOBS:
    d = C.RESULTS / "within_unit" / f"{gs}_{unit}_clean"
    z = np.load(d / "pca.npz"); arcs = np.load(d / f"archetypes_k{k}.npy"); rws = z["rows"]
    pts = z["scores"][:, : k - 1]
    dist = np.linalg.norm(pts[:, None, :] - arcs[None, :, : k - 1], axis=2)
    X = a.X[rws]
    for v in range(k):
        close = dist[:, v] <= np.quantile(dist[:, v], 0.10)
        r = dict(fit=f"{gs}_{unit}", corner=f"V{v + 1}", cycling=round(100 * a.obs["phase"].astype(str).values[rws][close].__ne__("G1").mean(), 0))
        for short, full in USE.items():
            idx = [a.var_names.get_loc(g) for g in H[full] if g in a.var_names]
            M = X[:, idx].toarray(); M = (M - M.mean(0)) / np.where(M.std(0) > 0, M.std(0), 1)
            s = M.mean(1)
            r[short] = round(float(s[close].mean() - s[~close].mean()), 2)
        r["top_set"] = max(USE, key=lambda s_: r[s_])
        rows.append(r)
t = pd.DataFrame(rows)
t.to_csv(C.RESULTS / "within_unit" / "hallmark_corners_clean.csv", index=False)
pd.set_option("display.width", 220)
print(t.to_string(index=False))
