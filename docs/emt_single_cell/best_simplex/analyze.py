#!/usr/bin/env python3
"""Best simplex per dataset / unit / gene set, and hybrid-vs-generalist numbers under three rules.

Reads every saved fit (fit_k*.json, archetypes_k*.npy, pca.npz) in the five dataset folders. Nothing is refitted.
Outputs (this folder): all_fits.csv (every candidate fit), best_fits.csv (one row per
dataset x level x unit x gene set x arm: lowest p, ties -> smallest k), rules_best.csv (specialist/generalist numbers).
Run from Project_1/:  .venv/bin/python docs/emt_single_cell/best_simplex/analyze.py
"""
import json, re, sys
from pathlib import Path
import numpy as np, pandas as pd
import anndata as ad

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
HYB = {"Early-Hybrid", "Hybrid", "Late-Hybrid"}
PUR = 0.35; FRAC = 0.05

D = dict(
    SCC=ROOT / "Mouse_SCC", COOK=ROOT / "Mouse_SCC/Cook", PDAC=ROOT / "Mouse_SCC/Mouse_PDAC_GSE173958",
    T4=ROOT / "Breast Cancer/Mouse_4T1_GSE303505", GSE=ROOT / "Breast Cancer/GSE114687_pseudospace")
NAME = dict(SCC="Mouse SCC", COOK="Cook", PDAC="Mouse PDAC", T4="4T1 clones", GSE="MCF10A/HuMEC")

# ------------------------------------------------------------------ labels
def load_obs():
    obs = {}
    m = pd.read_csv(D["SCC"] / "data/processed/cell_metadata.csv").set_index("barcode"); obs["SCC"] = m
    for key, p in (("PDAC", "data/processed/expr_tumour.h5ad"), ("T4", "data/processed/expr_tumour.h5ad"), ("GSE", "data/processed/expr_tumour.h5ad")):
        obs[key] = ad.read_h5ad(D[key] / p, backed="r").obs.copy()
    obs["COOK"] = cook_obs()
    return obs

def cook_obs():
    cache = OUT / "cook_ks_states.csv"
    if cache.exists():
        return pd.read_csv(cache, index_col=0)
    a = ad.read_h5ad(D["COOK"] / "data/processed/expr_log_filtered.h5ad")
    ks = pd.read_csv(ROOT / "Breast Cancer/data/processed/ks_cellline_signature_tan2014.csv")
    cat = ks["category"].str.lower().str[:3]; pres = set(a.var_names)
    epi = [g for g, c in zip(ks["gene_symbol"], cat) if c == "epi" and g in pres]
    mes = [g for g, c in zip(ks["gene_symbol"], cat) if c == "mes" and g in pres]
    def zm(gs):
        M = a[:, gs].X; M = M.toarray() if hasattr(M, "toarray") else np.asarray(M)
        s = M.std(0); return ((M - M.mean(0)) / np.where(s > 0, s, 1)).mean(1)
    o = a.obs.copy(); o["KS_EMT"] = zm(mes) - zm(epi)
    cuts = np.quantile(o["KS_EMT"], [.25, .5, .75])
    o["cluster_annotation"] = np.array(["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"])[np.searchsorted(cuts, o["KS_EMT"].values)]
    o.to_csv(cache); return o

# ------------------------------------------------------------------ enumerate fits
def fits():
    rows = []
    def add(ds, level, unit, gs, arm, d, filt=None):
        for f in sorted(Path(d).glob("fit_k*.json")):
            j = json.loads(f.read_text()); k = int(re.search(r"k(\d+)", f.name).group(1))
            el = None; lg = Path(d) / "log_03_esv.txt"
            if lg.exists():
                mm = re.search(r"ESV elbow at k=(\d+)", lg.read_text())
                el = int(mm.group(1)) if mm else None
            rows.append(dict(elbow=el, dataset=ds, level=level, unit=unit, geneset=gs, arm=arm, k=k, p=j["p_value"], t_ratio=j["t_ratio"],
                             null_median=j.get("null_median"), n_null=j.get("n_null", j.get("n_null")), n_cells_fit=j.get("n_cells"), dir=str(Path(d).relative_to(ROOT)), filt=filt))
    R = D["SCC"] / "results"
    for gs, tag, arm, filt in (("all", "all_main", "light QC", None), ("all", "all_scrna", "scRNA only", "scrna"),
                               ("all", "all_scrna6", "six clean tumours", "six"), ("ks", "ks_main", "light QC", None),
                               ("ks", "ks_scrna", "scRNA only", "scrna"), ("ks", "ks_scrna6", "six clean tumours", "six"),
                               ("hallmark", "hallmark_main", "light QC", None), ("hallmark", "hallmark_scrna", "scRNA only", "scrna"), ("hallmark", "hallmark_scrna6", "six clean tumours", "six")):
        add("SCC", "pooled", "pooled", gs, arm, R / tag, filt)
    for t in ("SC15", "SC32", "SC46"):
        add("SCC", "within", t, "all", "light QC", R / "within_tumour_raw" / t)
    for p in sorted((R / "within_tumour_gs").glob("*_*")) if (R / "within_tumour_gs").exists() else []:
        gs, t = p.name.split("_", 1); add("SCC", "within", t, gs, "light QC", p)
    C = D["COOK"] / "results"
    for gs in ("all", "ks"):
        add("COOK", "pooled", "pooled (4 lines)", gs, "light QC", C / f"{gs}_pooled", "cook_fwd")
        add("COOK", "pooled", "pooled (4 lines)", gs, "cleaned", C / f"{gs}_pooled_clean", "cook_fwd")
        for ln in ("A549", "DU145", "MCF7", "OVCA420"):
            add("COOK", "within", ln, gs, "light QC", C / f"{gs}_{ln}", f"cook_line:{ln}")
            add("COOK", "within", ln, gs, "cleaned", C / "within_line" / f"{gs}_{ln}_clean_ht", f"cook_line:{ln}")
    P = D["PDAC"] / "results"
    for gs in ("all", "ks"):
        add("PDAC", "pooled", "pooled (in vivo)", gs, "light QC", P / f"{gs}_in_vivo", "pdac_vivo")
        for u in ("M1-PT", "M1-SS", "M1-Liver", "M1-Met", "M2-PT", "M2-Liver", "M2-Met", "M1", "M2"):
            add("PDAC", "within", u, gs, "light QC", P / "within_unit" / f"{gs}_{u}")
            add("PDAC", "within", u, gs, "cleaned", P / "within_unit" / f"{gs}_{u}_clean")
    T = D["T4"] / "results"
    for gs in ("all", "ks"):
        add("T4", "pooled", "pooled (6 clones)", gs, "light QC", T / f"{gs}_main")
        add("T4", "pooled", "pooled (6 clones)", gs, "cleaned", T / f"{gs}_clean", "t4_clean")
        for u in ("E", "EM1", "EM2", "EM3", "M", "4T1"):
            add("T4", "within", u, gs, "light QC", T / "within_unit" / f"{gs}_{u}")
            add("T4", "within", u, gs, "cleaned", T / "within_unit" / f"{gs}_{u}_clean_dc")
    G = D["GSE"] / "results"
    for gs in ("all", "ks"):
        add("GSE", "pooled", "pooled (4 samples)", gs, "light QC", G / f"{gs}_main")
        for u in ("MCF10A", "MCF10A_Mock", "MCF10A_TGFB", "HuMEC", "HuMEC_Mock", "HuMEC_TGFB"):
            add("GSE", "within", u, gs, "light QC", G / "within_unit" / f"{gs}_{u}")
            add("GSE", "within", u, gs, "cleaned", G / "within_unit" / f"{gs}_{u}_clean_ht")
    # Hallmark EMT gene set (fitted for this analysis; pooled light QC, within-unit on the cleaned arm where one exists)
    for t in ("SC15", "SC32", "SC46", "SC33", "SC2b1", "SC2b3", "SC35", "SC43", "MU3", "MU4", "MU7", "MU8"):
        pass  # SCC hallmark within-tumour dirs are picked up by the within_tumour_gs glob above
    add("COOK", "pooled", "pooled (4 lines)", "hallmark", "light QC", C / "hallmark_pooled", "cook_fwd")
    add("COOK", "pooled", "pooled (4 lines)", "hallmark", "cleaned", C / "hallmark_pooled_clean", "cook_fwd")
    for ln in ("A549", "DU145", "MCF7", "OVCA420"):
        add("COOK", "within", ln, "hallmark", "cleaned", C / "within_line" / f"hallmark_{ln}_clean_ht", f"cook_line:{ln}")
    for ln in ("A549", "DU145", "MCF7", "OVCA420"):
        add("COOK", "within", ln, "hallmark", "light QC", C / f"hallmark_{ln}", f"cook_line:{ln}")
    add("PDAC", "pooled", "pooled (in vivo)", "hallmark", "light QC", P / "hallmark_in_vivo", "pdac_vivo")
    for g_ in ("all", "ks", "hallmark"):
        add("PDAC", "pooled", "pooled (in vivo)", g_, "cleaned", P / f"{g_}_in_vivo_clean", "pdac_vivo")
        add("GSE", "pooled", "pooled (4 samples)", g_, "cleaned", G / f"{g_}_clean")
        if g_ != "hallmark": pass
    for u in ("M1-PT", "M1-SS", "M1-Liver", "M1-Met", "M2-PT", "M2-Liver", "M2-Met", "M1", "M2"):
        add("PDAC", "within", u, "hallmark", "cleaned", P / "within_unit" / f"hallmark_{u}_clean")
        add("PDAC", "within", u, "hallmark", "light QC", P / "within_unit" / f"hallmark_{u}")
    add("T4", "pooled", "pooled (6 clones)", "hallmark", "light QC", T / "hallmark_main")
    add("T4", "pooled", "pooled (6 clones)", "hallmark", "cleaned", T / "hallmark_clean", "t4_clean")
    for u in ("EM1", "EM2", "EM3", "M", "4T1"):
        add("T4", "within", u, "hallmark", "cleaned", T / "within_unit" / f"hallmark_{u}_clean_dc")
        add("T4", "within", u, "hallmark", "light QC", T / "within_unit" / f"hallmark_{u}")
    add("T4", "within", "E", "hallmark", "light QC", T / "within_unit" / "hallmark_E")
    add("GSE", "pooled", "pooled (4 samples)", "hallmark", "light QC", G / "hallmark_main")
    for u in ("MCF10A", "MCF10A_Mock", "MCF10A_TGFB", "HuMEC", "HuMEC_Mock", "HuMEC_TGFB"):
        add("GSE", "within", u, "hallmark", "cleaned", G / "within_unit" / f"hallmark_{u}_clean_ht")
        add("GSE", "within", u, "hallmark", "light QC", G / "within_unit" / f"hallmark_{u}")
    return pd.DataFrame(rows)

# ------------------------------------------------------------------ weights and rules
def bary(pts, V):
    A = np.vstack([V.T, np.ones((1, V.shape[0]))]); B = np.vstack([pts.T, np.ones((1, len(pts)))])
    return np.linalg.solve(A, B).T

def rules(pts, V):
    """pts: n x (k-1); V: k x (k-1). Returns dict rule -> boolean 'generalist' per cell, plus inside flag."""
    k = V.shape[0]
    w = bary(pts, V); inside = (w >= -1e-6).all(1)
    c = np.clip(w, 0, None); c = c / np.where(c.sum(1, keepdims=True) == 0, 1, c.sum(1, keepdims=True))
    s = np.sort(c, 1); pur = s[:, -1] - s[:, -2]
    d = np.linalg.norm(pts[:, None, :] - V[None, :, :], axis=2)
    cut = np.quantile(d, FRAC, axis=0)
    gen = {"gap": pur < PUR, "groves5": (d <= cut).sum(1) == 0, "w95": c.max(1) <= 0.95}
    return gen, inside, c.argmax(1)

def cells_for(row, obs):
    d = ROOT / row["dir"]; z = np.load(d / "pca.npz", allow_pickle=True); k = int(row["k"])
    V = np.load(d / f"archetypes_k{k}.npy")[:, : k - 1]
    ds = row["dataset"]
    if "rows" in z.files:
        o = obs[ds].iloc[z["rows"]]; pts = z["scores"][:, : k - 1]
    else:
        o = obs[ds].reindex(z["barcodes"]); pts = z["scores_all"][:, : k - 1]
        f = row["filt"] if isinstance(row["filt"], str) else None; keep = np.ones(len(o), bool)
        if f == "scrna": keep = (o["assay"] == "scRNA").values
        elif f == "six": keep = o["sample"].isin(["SC2b1", "SC2b3", "SC32", "SC33", "SC35", "SC43"]).values
        elif f == "cook_fwd": keep = ~o["removal_arm"].astype(bool).values
        elif f and f.startswith("cook_line:"): keep = (~o["removal_arm"].astype(bool) & (o["cell_line"] == f.split(":")[1])).values
        elif f == "pdac_vivo": keep = (o["assay"] == "in_vivo").values
        elif f == "t4_clean": keep = o["cluster_annotation"].notna().values
        o = o[keep]; pts = pts[keep]
    return pts, V, o

def main():
    obs = load_obs()
    F = fits()
    F = F[F["dir"].map(lambda d: (ROOT / d / "pca.npz").exists())].copy()
    F.drop(columns="filt").to_csv(OUT / "all_fits.csv", index=False)
    keys = ["dataset", "level", "unit", "geneset", "arm"]
    # primary rule (the report's own): ignore k beyond the ESV elbow when an elbow was computed (pooled fits)
    ok = F[F["elbow"].isna() | (F["k"] <= F["elbow"])]
    best = ok.sort_values(["p", "k"]).groupby(keys, sort=False).head(1).reset_index(drop=True)
    anyk = F.sort_values(["p", "k"]).groupby(keys, sort=False).head(1).set_index(keys)[["k", "p"]].rename(columns={"k": "k_any", "p": "p_any"})
    best = best.merge(anyk.reset_index(), on=keys, how="left")
    recs = []
    for _, r in best.iterrows():
        pts, V, o = cells_for(r, obs)
        gen, inside, near = rules(pts, V)
        hyb = o["cluster_annotation"].astype(str).isin(HYB).values
        rec = r.to_dict(); rec.update(n_cells=len(o), n_hybrid=int(hyb.sum()), hybrid_pct=100 * hyb.mean(), inside_pct=100 * inside.mean())
        for name, g in gen.items():
            rec[f"gen_pct_{name}"] = 100 * g.mean()
            rec[f"hyb_gen_{name}"] = 100 * (hyb & g).sum() / max(hyb.sum(), 1)       # % of hybrids that are generalists
            rec[f"gen_hyb_{name}"] = 100 * (hyb & g).sum() / max(g.sum(), 1)         # % of generalists that are hybrids
        recs.append(rec)
    B = pd.DataFrame(recs)
    B.to_csv(OUT / "best_fits.csv", index=False)
    print(B[keys + ["k", "p", "n_cells", "hybrid_pct", "hyb_gen_gap", "hyb_gen_groves5", "hyb_gen_w95"]].round(3).to_string())

if __name__ == "__main__":
    main()
