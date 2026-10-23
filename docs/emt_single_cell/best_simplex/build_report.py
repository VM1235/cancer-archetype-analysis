#!/usr/bin/env python3
"""Write REPORT_best_simplex.md (methods + results per figure) from best_fits.csv. Run from Project_1/."""
from pathlib import Path
import pandas as pd, numpy as np
HERE = Path(__file__).resolve().parent
B = pd.read_csv(HERE / "best_fits.csv"); F = pd.read_csv(HERE / "all_fits.csv")
NAME = dict(SCC="Mouse skin SCC (GSE288137)", COOK="Cook TGF-β cell lines (GSE147405)", PDAC="Mouse PDAC (GSE173958)", T4="4T1 clones (GSE303505)", GSE="MCF10A / HuMEC pseudospace (GSE114687)")
sig = B[B.p <= 0.05]
GSN = {"all": "all genes", "ks": "KS genes", "hallmark": "Hallmark EMT genes"}
def md(df, cols, heads, fmt=None):
    out = ["| " + " | ".join(heads) + " |", "|" + "---|" * len(heads)]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            cells.append(f"{v:.3g}" if c == "p" else (f"{v:.0f}" if isinstance(v, (float, np.floating)) else str(v)))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)
COLS = ["unit", "geneset", "arm", "k", "p", "n_cells", "hybrid_pct", "inside_pct", "hyb_gen_gap", "hyb_gen_groves5", "hyb_gen_w95", "gen_hyb_gap", "gen_hyb_groves5", "gen_hyb_w95"]
HEADS = ["unit", "genes", "arm", "k", "p", "cells", "hybrid %", "inside %", "hyb→gen gap", "hyb→gen 5%", "hyb→gen .95", "gen→hyb gap", "gen→hyb 5%", "gen→hyb .95"]
L = []
w = L.append
w("# Best simplex per dataset and hybrid / generalist analysis: methods and results\n")
w("Thesis project, IISc · 21 October 2026 · companion to `Best_simplex_hybrid_generalists.pptx`. All numbers come from `best_fits.csv`, `all_fits.csv` and the saved fit folders; no cell-level number is typed by hand.\n")
w("## 1. The question and the definitions\n")
w("For each of the five datasets, for the **all-genes** fit and the **KS EMT-gene** fit, **pooled** and **within each unit** of the dataset, the questions were: (1) which simplex is the most significant, and at what k; (2) in that simplex, what share of hybrid cells are generalists and what share of generalists are hybrid, under all three specialist rules.\n")
w("| Term | Definition used here | Where in the code |\n|---|---|---|\n"
  "| Simplex fit | PCA on the unit's cells (genes centred, not scaled), then PCHA (archetypal analysis) on the first k−1 PCs: 15 starts, delta 0, keep the largest-volume simplex | `scc_common.fit_pca`, `fit_simplex`, `03_fit_and_test.py` / `08_within_unit.py` in each dataset's `codes/` |\n"
  "| t-ratio | volume of the fitted simplex ÷ volume of the convex hull of the cells in the same k−1 PCs | `src/archetypes.py`, `simplex_volume`; `08_within_tumour.py::t_ratio` |\n"
  "| p-value | 200 shuffles; each of the k−1 PCs is permuted independently, refitted (5 starts), t-ratio recorded; p = (1 + shuffles ≥ real) / (1 + shuffles). The floor is 1/201 = 0.00498 | same scripts |\n"
  "| Hybrid | Early-Hybrid, Hybrid (SCC) or Late-Hybrid state. SCC: authors' labels; PDAC, 4T1, GSE114687: KS-score bins (PDAC 40/75/95%, others quartiles); Cook: KS-score quartiles computed here (Section 3.5) | `analyze.py::HYB`, `cook_obs` |\n"
  "| Barycentric weights | the k weights (sum 1) that express a cell as a mix of the k corners; cells outside the simplex have negative weights, clipped to 0 and renormalised | `analyze.py::bary`, `rules` |\n"
  "| Rule 1: gap | specialist if (largest − second-largest weight) ≥ 0.35 | `rules()`, `PUR = 0.35` |\n"
  "| Rule 2: closest 5% | for each corner, the 5% of cells nearest to it (Euclidean, in the k−1 PCs) are its specialists; ties go to the corner for which the cell is nearest relative to that corner's cut-off | `rules()`, `FRAC = 0.05` |\n"
  "| Rule 3: weight > 0.95 | specialist if one weight exceeds 0.95 | `rules()` |\n"
  "| Hybrids → generalists | 100 × (hybrid and generalist) ÷ (all hybrid) | `hyb_gen_<rule>` |\n"
  "| Generalists → hybrids | 100 × (hybrid and generalist) ÷ (all generalist); compare with the hybrid share of all cells (`hybrid_pct`) | `gen_hyb_<rule>` |\n")
w("\n## 2. What was run\n")
w("Nothing was refitted for the main question: every fit already saved `pca.npz` (PCA scores), `archetypes_k*.npy` (corners) and `fit_k*.json` (t-ratio, null, p). `analyze.py` walks these folders (`all_fits.csv` lists all %d candidate fits). The one thing that was missing was the SCC within-tumour fits, so these were run now with a new script (Section 3.6)." % len(F))
w("\n**Candidate fits counted** (light QC and cleaned arms; UMI-cutoff sweeps, depth-matched, regression and contamination-cutoff variants excluded; the withdrawn SCC V1-regression fit excluded):\n")
w("| Dataset | Pooled fits | Within-unit fits |\n|---|---|---|\n"
  "| SCC | all genes, KS genes and Hallmark: all cells, scRNA only, six clean tumours (KS and Hallmark scRNA-only / six-tumour arms and the Hallmark all-cells fit are new) | 12 tumours, all genes and KS genes, raw expression (SC15, SC32, SC46 all-genes from the original run; the rest from `08d`) |\n"
  "| Cook | all genes and KS genes, forward arm, 4 lines | 4 lines × light QC / cleaned (human technical-gene fix `_ht`) |\n"
  "| PDAC | all genes and KS genes, in vivo libraries | 9 units (7 libraries, mouse 1, mouse 2) × light QC / cleaned |\n"
  "| 4T1 | all genes and KS genes, light QC and cleaned | 6 clones × light QC / cleaned (clone E has no cleaned arm: 142 cells) |\n"
  "| GSE114687 | all genes and KS genes, light QC | 6 units × light QC / cleaned (`_ht`) |\n")
w("\n## 3. Methods in detail\n")
w("### 3.1 Choosing the best simplex of a fit\n")
w("For every dataset × level × unit × gene set × arm the candidates are the saved k values (pooled: 3 to 7; within unit: 3 and 4). The **most significant** is the lowest p; ties (very common at the 0.00498 floor) go to the **smallest k**, the parsimonious choice and the project's rule (\"smallest k with p < 0.01\"). For pooled fits k **beyond the ESV elbow is ignored** (the report's rule, because k = 6 and 7 pass on almost any data); the unrestricted best is kept in `k_any` / `p_any` for transparency. A fit is called significant at p ≤ 0.05.\n")
w("```python\nok = F[F['elbow'].isna() | (F['k'] <= F['elbow'])]\nbest = ok.sort_values(['p','k']).groupby(keys).head(1)\n```\n")
w("**Caveats.** The best of several k and arms is an optimistic p; with 200 shuffles p cannot go below 0.00498, so \"most significant\" among many floor p-values is decided by k alone.\n")
w("### 3.2 Cells analysed in each fit\n")
w("The rules are applied to the same cells the fit describes. Pooled: SCC all cells (23,959), scRNA only (14,059) or the six clean tumours; Cook forward-arm cells only (the removal arm was projected as guests and is excluded); PDAC in vivo tumour cells; 4T1 and GSE114687 all pooled cells (cleaned 4T1: 13,256). Within unit: the cells that entered the fit; **fits are capped at 4,000 random cells**, as in the original within-unit scripts, so large units are subsamples.\n")
w("### 3.3 Barycentric weights and the three rules\n")
w("```python\nw = bary(pts, V)                        # k weights per cell\nc = clip(w, 0) / row-sum                # outside cells: negative weights set to 0, renormalised\npur = top1(c) - top2(c)                 # purity\nd = ||pts - V_j||                       # distance of every cell to every corner\ncut = quantile(d, 0.05, axis=0)         # per-corner 5% cut-off\ngeneralist = {'gap': pur < 0.35, 'groves5': (d <= cut).sum(1) == 0, 'w95': c.max(1) <= 0.95}\n```\n")
w("The three rules are those in each dataset's `06_groves_specialists.py` (gap, Groves closest-5%, Groves weight > 0.95). Because PCHA is run with delta = 0 the corners lie inside the data cloud, so almost no cell has a weight above 0.95; the third rule therefore calls ≈100% of cells generalist in every fit and carries little information. The closest-5% rule caps specialists at about 5% per corner (≈15% for k = 3), so it is the most generous to generalists. The gap rule is the strictest of the useful rules.\n")
w("### 3.4 Hybrid labels per dataset\n")
w("SCC: the authors' `cluster_annotation` (Early-Hybrid and Hybrid are hybrid; Late-EMT is not). PDAC, 4T1, GSE114687: `cluster_annotation` from the KS score (Early-Hybrid and Late-Hybrid are hybrid). These states are **cut from the same KS genes** whose fit is tested, so KS-gene results are partly circular (Gap 4 in the report).\n")
w("### 3.5 Cook hybrids (new, assumption chosen with you)\n")
w("Cook has no state labels, only time after TGF-β. Following your choice, hybrids are defined from the KS score as in GSE114687: `KS_EMT = mean z(KS mesenchymal genes) − mean z(KS epithelial genes)` (163 human symbols, z over all 9,084 cells), cut at the quartiles into Epithelial / Early-Hybrid / Late-Hybrid / Mesenchymal; the middle two are hybrid (`analyze.py::cook_obs`, cached in `cook_ks_states.csv`). **Caveat:** the cut-offs are pooled over the four cell lines, so a line's hybrid share reflects its line offset (A549 ≈ 71% hybrid, MCF7 ≈ 11%, OVCA420 ≈ 14%).\n")
w("### 3.6 New SCC within-tumour fits\n")
w("`Mouse_SCC/codes/08d_within_tumour_geneset.py` repeats `08_within_tumour.py --raw` for any gene set: PCA on one tumour's cells only (exact, 10 PCs), PCHA 15 starts, 200 shuffles, k = 3 and 4, outputs in `results/within_tumour_gs/<geneset>_<tumour>/`. Run for KS genes on all 12 tumours and for all genes on the nine tumours not covered before. Those nine are mostly one-state tumours (e.g. SC35 is 99% hybrid or late EMT), so a significant triangle there says little about EMT states.\n")
w("### 3.7 The Hallmark EMT gene set (new fits for this analysis)\n")
w("Gene set: MSigDB Hallmark 2020 \"Epithelial Mesenchymal Transition\" (200 human genes; `Mouse_SCC/hallmark_data/MSigDB_Hallmark_2020.txt`). Mouse datasets map symbols by capitalisation (COL3A1 → Col3a1), as for the KS genes; human datasets use the symbols unchanged; genes absent from a dataset are dropped (`hallmark_emt_genes()` in each `scc_common.py` / `cook_common.py`, selected with `--geneset hallmark` in `02_pca.py`, `03_fit_and_test.py`, `08_within_unit.py`, `07_within_line_clean.py`, `08d_within_tumour_geneset.py`). Everything else is identical to the KS and all-genes fits (same PCA, PCHA, 200 shuffles, seeds). **Arms:** pooled fits on the light-QC data (4T1 also cleaned); within-unit fits on the cleaned arm where a cleaning arm exists (Cook with the human technical-gene fix, PDAC ≥10,000 UMI, 4T1 ≥10,000 UMI plus contamination filter, GSE114687 ≥10,000 UMI with the human fix) and on raw expression for SCC; 4T1 clone E has no cleaned arm (142 cells) and uses light QC. The Hallmark fits were run on one arm only to limit compute, so their \"arm\" is not always the one with the best p among arms as for all genes and KS genes.\n")
w("### 3.8 Reproducing\n")
w("```\n.venv/bin/python docs/emt_single_cell/best_simplex/analyze.py     # all_fits.csv, best_fits.csv\n.venv/bin/python docs/emt_single_cell/best_simplex/figures.py     # figures/\n.venv/bin/python docs/emt_single_cell/best_simplex/build_deck.py  # the pptx\n.venv/bin/python docs/emt_single_cell/best_simplex/build_report.py\n```\n")
w("\n## 4. Results, figure by figure\n")
w("### Figure 1. How often a simplex passes the permutation test (`figures/summary_pass.png`)\n![](figures/summary_pass.png)\n")
tot = B.groupby(["level", "geneset"]).apply(lambda x: f"{int((x.p<=0.05).sum())}/{len(x)}", include_groups=False) if "include_groups" in pd.DataFrame.groupby.__code__.co_varnames else B.groupby(["level", "geneset"]).apply(lambda x: f"{int((x.p<=0.05).sum())}/{len(x)}")
w("**What it shows.** For each dataset the fraction of fits (best k, light QC and cleaned counted separately) with p ≤ 0.05, for pooled fits and within-unit fits and for the two gene sets. **How.** `figures.py::fig_summary` counts `best_fits.csv` rows. **Totals:** " + "; ".join(f"{lv} {('all genes' if g=='all' else 'KS genes')} {v}" for (lv, g), v in tot.items()) + ". **Reading.** All-genes fits pass much more often than KS-gene fits, within units more often than pooled except where the pooled shape is sample identity (Cook, GSE114687 pooled). Pooled passes in Cook and the all-genes pooled fits are identity corners (cell lines, clones, samples), see the report.\n")
fig_notes = {
"SCC": "Units are the 12 tumours (the nine tumours beyond SC15, SC32, SC46 are new). Pooled arms: all cells (version A), scRNA only (B), six clean tumours, KS genes all cells (C). Pooled k for KS genes is 5 (k = 6 and 7 pass but lie beyond the elbow 5).",
"COOK": "Per-line fits use the light-QC fit (`results/all_<line>`, PCA on the line's forward-arm cells) and the cleaned fit with the human technical-gene fix (`within_line/*_clean_ht`).",
"PDAC": "Units are the seven libraries plus mouse 1 and mouse 2, each with light-QC and cleaned arms. Pooled in vivo uses the 44,266 in vivo tumour cells.",
"T4": "Clone E has no cleaned arm (142 cells after the UMI rule). Pooled arms: light QC (19,930 cells) and cleaned (13,256 cells).",
"GSE": "Six units (MCF10A and HuMEC, each both arms and per arm) with light-QC and cleaned (`_ht`) arms; pooled over the four samples, light QC only.",
}
for n, ds in enumerate(("SCC", "COOK", "PDAC", "T4", "GSE"), 2):
    d = B[B.dataset == ds]; sg = d[d.p <= 0.05]
    w(f"### Figure {n}. {NAME[ds]}\n")
    w(f"**{n}a. Best p and k per unit (`figures/p_heat_{ds}.png`).**\n![](figures/p_heat_{ds}.png)\n")
    w(f"Rows are units, columns are gene set × arm; each cell is the best simplex of that fit (k and p; bold = p ≤ 0.05; \"not run\" = no fit exists). {fig_notes[ds]} **{len(sg)} of {len(d)} fits reach p ≤ 0.05** (all genes {int((sg.geneset=='all').sum())}, KS genes {int((sg.geneset=='ks').sum())}). Generated by `figures.py::fig_pheat`.\n")
    w(f"**{n}b. Hybrids and generalists under three rules (`figures/rules_{ds}.png`).**\n![](figures/rules_{ds}.png)\n")
    w("Left: % of hybrid cells that are generalists (three coloured dots = three rules). Right: % of generalist cells that are hybrid, with the white bar = hybrid share of all cells in that fit (if a dot sits on the bar, generalists are no more hybrid than average). Significant fits first, then (ns). `figures.py::fig_rules`.\n")
    t = sg.sort_values(["p", "k"]).head(18)
    w(f"**{n}c. Significant fits in detail** (top 18 by p; all rows in `best_fits.csv`):\n")
    w(md(t, COLS, HEADS) if len(t) else "No fit reaches p ≤ 0.05.")
    w("")
    if len(sg):
        w(f"**Summary for {NAME[ds]} (significant fits).** Hybrid→generalist medians: gap {sg.hyb_gen_gap.median():.0f}%, closest-5% {sg.hyb_gen_groves5.median():.0f}%, weight>0.95 {sg.hyb_gen_w95.median():.0f}%. Generalist→hybrid medians: gap {sg.gen_hyb_gap.median():.0f}%, closest-5% {sg.gen_hyb_groves5.median():.0f}%, weight>0.95 {sg.gen_hyb_w95.median():.0f}%, against a median hybrid share of {sg.hybrid_pct.median():.0f}% of cells.\n")
w("### Figure 7. All datasets together (`figures/cross_rules.png`)\n![](figures/cross_rules.png)\n")
w(f"One dot per significant fit (p ≤ 0.05), medians in white, for the three rules. Across all {len(sig)} significant fits, the median share of hybrids that are generalists is **{sig.hyb_gen_gap.median():.0f}% (gap), {sig.hyb_gen_groves5.median():.0f}% (closest 5%), {sig.hyb_gen_w95.median():.0f}% (weight > 0.95)**; the median share of generalists that are hybrid is {sig.gen_hyb_gap.median():.0f}%, {sig.gen_hyb_groves5.median():.0f}%, {sig.gen_hyb_w95.median():.0f}% against a hybrid share of {sig.hybrid_pct.median():.0f}%. `figures.py::fig_cross_rules`.\n")
w("### Figure 8. Where hybrid cells sit (`figures/triangles.png`)\n![](figures/triangles.png)\n")
w("For each dataset the most significant all-genes k = 3 fit with at least 200 hybrid cells and a hybrid share of 10 to 90% (within-unit preferred). Cells are placed in PC1 vs PC2 of the unit, the simplex is drawn in blue, hybrids are orange (generalist) or red (specialist) under the gap rule (top row) or the closest-5% rule (bottom row), non-hybrids grey. This is a picture of Section 3.3 on real cells. `figures.py::fig_triangles` (uses `analyze.cells_for` and `analyze.rules`).\n")
w("### Figure book: one figure per analysis (`figures/book/`)\n")
w("**What it shows.** For every dataset and every gene set (**A. entire gene set, B. KS EMT genes, C. Hallmark EMT genes**) there is a set of pages, three units per page: the pooled fit first, then every tumour / cell line / library / clone. Each row is the unit's best simplex (lowest p among the arms; ties go to the cleaned arm, then to the smaller k; the row title is red and says NOT significant when p > 0.05). The four panels per row are the same cells on PC1 vs PC2 of that unit: (1) coloured by EMT state, hybrids being the two middle states; (2) to (4) specialists coloured by the corner they belong to and generalists grey, under the gap, closest-5% and weight > 0.95 rules. Under each panel: for state, the hybrid share of cells and the share inside the simplex; for the rules, the % of hybrids that are generalists and the % of generalists that are hybrid. **How it is made.** `figurebook.py` selects the fit (`pick`), loads the saved PCA scores and corners (`analyze.cells_for`), computes the rules on all cells of the fit (`analyze.rules`) and draws at most 3,500 randomly chosen cells per panel (statistics use all cells). For k ≥ 4 the simplex cannot be drawn in 2D: the corners are projected on PC1/PC2 and joined by their convex hull, so some corners may fall inside the outline. \n")
for ds_ in ("SCC", "COOK", "PDAC", "T4", "GSE"):
    for gs_, gt in (("all", "A. Entire gene set"), ("ks", "B. KS EMT genes"), ("hallmark", "C. Hallmark EMT genes")):
        pgs = sorted((HERE / "figures" / "book").glob(f"{ds_}_{gs_}_p*.png"))
        if not pgs: continue
        w(f"#### {NAME[ds_]}: {gt}\n")
        for f in pgs: w(f"![]({f.relative_to(HERE).as_posix()})\n")
w("\n## 5. What the numbers say\n")
gm = lambda c: sig[c].median()
w(f"1. **Hybrid cells are mostly generalists under every rule, and the ordering of the rules is stable.** Over the {len(sig)} significant fits (of {len(B)}), the % of hybrids that are generalists has median {gm('hyb_gen_gap'):.0f}% (range {sig.hyb_gen_gap.min():.0f} to {sig.hyb_gen_gap.max():.0f}) under the gap rule, {gm('hyb_gen_groves5'):.0f}% ({sig.hyb_gen_groves5.min():.0f} to {sig.hyb_gen_groves5.max():.0f}) under closest-5% and {gm('hyb_gen_w95'):.0f}% ({sig.hyb_gen_w95.min():.0f} to {sig.hyb_gen_w95.max():.0f}) under weight > 0.95. The gap rule is the most variable: its low values (a few %) come from fits whose corners are cell lines or clones (for example the Cook pooled KS fit, where hybrid-state cells sit at the line corners).\n")
d_ = (sig.gen_hyb_gap - sig.hybrid_pct)
w(f"2. **Generalists are barely enriched for hybrids.** The % of generalists that are hybrid is, under the gap rule, a median of {d_.median():+.1f} percentage points from the hybrid share of all cells (middle half of fits: {d_.quantile(.25):+.1f} to {d_.quantile(.75):+.1f}); under closest-5% and weight > 0.95 it is almost equal to the baseline. Hybrids are generalists mainly because most cells are generalists.\n")
gs_med = sig.groupby("geneset")[["hyb_gen_gap", "hyb_gen_groves5"]].median().round(0)
w("3. **Gene sets.** Fits passing p ≤ 0.05: " + ", ".join(f"{GSN[g]} {int((B[B.geneset==g].p<=0.05).sum())}/{len(B[B.geneset==g])}" for g in ("all","ks","hallmark")) + ". The hybrid/generalist medians are similar across gene sets (gap rule: " + ", ".join(f"{GSN[g]} {gs_med.loc[g,'hyb_gen_gap']:.0f}%" for g in gs_med.index) + ").\n")
w("4. **The three rules are not three independent confirmations**: they are one construction at three stringencies. The informative contrast is gap vs closest-5%.\n")
w("\n### 5.5 Hallmark EMT genes vs KS EMT genes, pooled fits (the gating comparison)\n")
pp = B[B.level == "pooled"].copy()
rows_ = []
for (ds_, arm_), g_ in pp.groupby(["dataset", "arm"]):
    r_ = {gg: g_[g_.geneset == gg] for gg in ("all", "ks", "hallmark")}
    if all(len(v) for v in r_.values()):
        f_ = lambda v: f"k={int(v.k.iloc[0])}, p={v.p.iloc[0]:.3g}{chr(8224) if v.n_null.iloc[0] < 150 else ''}, inside {v.inside_pct.iloc[0]:.0f}%"
        rows_.append(f"| {NAME[ds_]} | {arm_} | {f_(r_['all'])} | {f_(r_['ks'])} | {f_(r_['hallmark'])} |")
w("| Dataset | Arm | All genes | KS genes | Hallmark EMT genes |\n|---|---|---|---|---|\n" + "\n".join(rows_) + "\n")
w("† = 100 shuffles (p floor 0.0099); fits without † used 200. Hallmark and the cleaned-pooled / SCC-extra fits used k = 3 to 5 only. **Reading.** Hallmark does not behave very differently from KS: both mostly fail pooled (PDAC, SCC all cells), both pass where the pooled shape is cell/sample identity (Cook; GSE114687 after cleaning), and they disagree only in single arms (GSE114687 light QC: KS p = 0.119, Hallmark p = 0.010; 4T1 light QC: KS 0.174, Hallmark 0.020; 4T1 cleaned: KS 0.030, Hallmark 0.119). There is no consistent large gain from Hallmark, so on the agreed rule (\"only if the difference is large\") the Hallmark simplex does not need a separate follow-up beyond the figures already made.\n")
w("\n## 6. Limitations\n")
w("- Best-of-many selection and the p floor (Section 3.1). Many fits were run per dataset; no multiplicity correction.\n- Within-unit fits are 4,000-cell subsamples.\n- Hybrid labels are circular with KS genes (PDAC, 4T1, GSE114687, Cook) and pooled-bin based in Cook.\n- Pooled all-genes simplices are often sample identity (cell lines, clones, mice, samples), so their hybrid/generalist numbers describe identity corners, not EMT corners.\n- Light-QC and cleaned arms are shown side by side and not harmonised across datasets (see the cleaning audit).\n- Hallmark, the cleaned pooled arms (Cook, PDAC, GSE114687) and the SCC KS scRNA-only / six-tumour arms were fitted for this analysis with k = 3 to 5 and 100 shuffles for the pooled fits (marked † in tables); the Hallmark within-unit fits and everything else use 200 shuffles. Hallmark pooled fits could not be examined beyond k = 5.\n- Hallmark within-unit fits were run on both arms (cleaned and light QC) for Cook, PDAC, 4T1 and GSE114687 and on raw expression for SCC.\n")
(HERE / "REPORT_best_simplex.md").write_text("\n".join(L)); print("report written", len("\n".join(L)))
