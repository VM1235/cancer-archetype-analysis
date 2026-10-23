#!/usr/bin/env python3
"""Annotate the corners of every best fit (best_fits.csv): what task / programme each corner represents, and how that was found.
For each fit and each corner j: S_j = the 5% of cells closest to the corner (closest-5% rule). Then
  * Hallmark score per cell = mean log-expression of each of the 50 MSigDB Hallmark gene sets; the corner's effect for a set =
    (mean over S_j - mean over all other cells) / SD over all cells of the fit; the top set (effect > 0.5) names the programme;
  * % of S_j in S or G2M phase (cycling); dominant EMT state, dominant sample (tumour / line / library / clone) and its share;
  * top marker genes = the genes with the largest mean(S_j) - mean(rest) among the genes used in the fit.
Output: corner_annotations.csv. Run from Project_1/:  .venv/bin/python docs/emt_single_cell/best_simplex/corners.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import analyze as A
ROOT = A.ROOT
EXPR = {"SCC": A.D["SCC"] / "data/processed/expr_log_filtered.h5ad", "COOK": A.D["COOK"] / "data/processed/expr_log_filtered.h5ad",
        "PDAC": A.D["PDAC"] / "data/processed/expr_tumour.h5ad", "T4": A.D["T4"] / "data/processed/expr_tumour.h5ad", "GSE": A.D["GSE"] / "data/processed/expr_tumour.h5ad"}
HUMAN = {"COOK", "GSE"}
SAMPLE = {"SCC": "sample", "COOK": "cell_line", "PDAC": "sample", "T4": "sample", "GSE": "sample"}
TASK = {"E2F Targets": "cycling", "G2-M Checkpoint": "cycling", "Mitotic Spindle": "cycling", "Epithelial Mesenchymal Transition": "EMT / matrix",
        "Hypoxia": "hypoxia / glycolysis", "Glycolysis": "hypoxia / glycolysis", "TNF-alpha Signaling via NF-kB": "inflammatory / stress response",
        "Inflammatory Response": "inflammatory / stress response", "Interferon Gamma Response": "interferon response", "Interferon Alpha Response": "interferon response",
        "p53 Pathway": "p53 / stress", "UV Response Up": "stress response", "TGF-beta Signaling": "TGF-β response", "KRAS Signaling Up": "KRAS-up (squamous / differentiation)",
        "Myogenesis": "muscle / matrix-like", "Myc Targets V1": "MYC / growth (translation)", "Myc Targets V2": "MYC / growth (translation)", "KRAS Signaling Dn": "KRAS-down (epithelial / differentiated)", "Coagulation": "coagulation / matrix-like", "Angiogenesis": "angiogenesis / matrix", "Bile Acid Metabolism": "metabolic", "Estrogen Response Early": "luminal / differentiation", "Estrogen Response Late": "luminal / differentiation",
        "Apical Junction": "epithelial junction", "Cholesterol Homeostasis": "metabolic", "Oxidative Phosphorylation": "metabolic", "MTORC1 Signaling": "metabolic / growth"}

def hallmark_sets(var_names, human):
    sets = {}
    for line in open(ROOT / "Mouse_SCC/hallmark_data/MSigDB_Hallmark_2020.txt"):
        p = line.rstrip("\n").split("\t"); g = [x for x in p[2:] if x]
        if not human: g = [x[0] + x[1:].lower() for x in g]
        pres = [x for x in g if x in set(var_names)]
        if len(pres) >= 5: sets[p[0]] = pres
    return sets

def main():
    B = pd.read_csv(HERE / "best_fits.csv"); obs = A.load_obs(); out = []
    for ds in ("SCC", "COOK", "PDAC", "T4", "GSE"):
        a = ad.read_h5ad(EXPR[ds]); X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X); var = np.array(a.var_names); bc = {b: i for i, b in enumerate(a.obs_names)}
        sets = hallmark_sets(var, ds in HUMAN); names = list(sets)
        gidx = {g: i for i, g in enumerate(var)}
        Mh = sp.lil_matrix((len(var), len(names)))
        for j, n in enumerate(names):
            for g in sets[n]: Mh[gidx[g], j] = 1.0 / len(sets[n])
        Mh = Mh.tocsr()
        for _, r in B[B.dataset == ds].iterrows():
            pts, V, o = A.cells_for(r, obs)
            rows = np.array([bc[b] for b in o.index])
            d = np.linalg.norm(pts[:, None, :] - V[None, :, :], axis=2); cut = np.quantile(d, 0.05, axis=0)
            Xr = X[rows]; S = np.asarray((Xr @ Mh).todense()); sd = S.std(0) + 1e-9
            gs_idx = None
            if r.geneset != "all":
                z = np.load(ROOT / r["dir"] / "pca.npz", allow_pickle=True)
                gs_idx = None
            state = o["cluster_annotation"].astype(str).values; phase = o["phase"].astype(str).values if "phase" in o else None
            samp = o[SAMPLE[ds]].astype(str).values
            for j in range(V.shape[0]):
                m = d[:, j] <= cut[j]
                if m.sum() < 5: continue
                eff = (S[m].mean(0) - S[~m].mean(0)) / sd
                if r.geneset == "hallmark":   # the fit is on the Hallmark EMT genes themselves, so that set is circular: name from the others
                    eff[names.index("Epithelial Mesenchymal Transition")] = -9
                top = np.argsort(-eff)[:2]
                dm = np.asarray(Xr[m].mean(0)).ravel() - np.asarray(Xr[~m].mean(0)).ravel()
                if r.geneset == "ks":
                    ks = pd.read_csv(ROOT / "Breast Cancer/data/processed/ks_cellline_signature_tan2014.csv")["gene_symbol"]
                    allow = set(ks) if ds in HUMAN else set(x[0] + x[1:].lower() for x in ks)
                    mask = np.array([g in allow for g in var]); dm = np.where(mask, dm, -9)
                elif r.geneset == "hallmark":
                    allow = set(sets["Epithelial Mesenchymal Transition"]); mask = np.array([g in allow for g in var]); dm = np.where(mask, dm, -9)
                tg = var[np.argsort(-dm)[:4]]
                sv = pd.Series(state[m]).value_counts(normalize=True); sm = pd.Series(samp[m]).value_counts(normalize=True)
                task = TASK.get(names[top[0]], names[top[0]]) if eff[top[0]] > 0.5 else "no Hallmark set above 0.5 SD"
                cyc = 100 * np.isin(phase[m], ["S", "G2M"]).mean() if phase is not None else np.nan
                if cyc >= 60 and task == "no Hallmark set above 0.5 SD": task = "cycling"
                out.append(dict(dataset=ds, level=r.level, unit=r.unit, geneset=r.geneset, arm=r.arm, k=int(r.k), p=r.p, corner=f"V{j + 1}", n_specialists=int(m.sum()),
                                top_set=names[top[0]], top_effect=round(float(eff[top[0]]), 2), second_set=names[top[1]], second_effect=round(float(eff[top[1]]), 2), task=task,
                                cycling_pct=round(cyc, 0), top_state=sv.index[0], top_state_pct=round(100 * sv.iloc[0], 0),
                                top_sample=sm.index[0], top_sample_pct=round(100 * sm.iloc[0], 0), top_genes=" ".join(tg)))
        print(ds, "annotated", flush=True)
    pd.DataFrame(out).to_csv(HERE / "corner_annotations.csv", index=False)

if __name__ == "__main__":
    main()
