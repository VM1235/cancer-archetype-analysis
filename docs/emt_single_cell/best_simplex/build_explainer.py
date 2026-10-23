#!/usr/bin/env python3
"""Explainer document: for each dataset, gene set (A entire / B KS / C Hallmark EMT) and level (pooled / best within-unit fit), the simplex
figure beside the explanation: significance, what the corners are (tasks) and how that was found, and the three-rule hybrid/generalist numbers.
Run from Project_1/ after analyze.py and corners.py:  .venv/bin/python docs/emt_single_cell/best_simplex/build_explainer.py
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import analyze as A, figures as FG, figurebook as FB
OUT = HERE / "figures" / "explainer"; OUT.mkdir(parents=True, exist_ok=True)   # the beige deck points OUT elsewhere
B = pd.read_csv(HERE / "best_fits.csv"); CA = pd.read_csv(HERE / "corner_annotations.csv")
BG, FGC, MUT = FG.BG, FG.FG, FG.MUT
GEN_COL, OUTL_COL, SPINE = "#3a414c", "#e8e8e8", "#2a2f37"   # theme: overridden by the beige deck builder
GS = {"all": "A. Entire gene set", "ks": "B. KS EMT genes", "hallmark": "C. Hallmark EMT genes"}
NAME = FG.NAME
DS_MEANING = {
 "SCC": "Mouse skin squamous cell carcinoma: 12 real tumours (4 multiome, 8 scRNA), 23,959 cells, with the authors' four EMT states. Tumour, assay and EMT state are tangled.",
 "COOK": "Four human cancer cell lines given TGF-β and sampled over time. Time is an EMT coordinate independent of our score, but each line has its own identity.",
 "PDAC": "Mouse pancreatic tumour cells from two mice at several body sites (primary, liver, metastasis). EMT states are ours (KS-score bins), so partly circular.",
 "T4": "Five single-cell clones of mouse 4T1 breast cancer chosen along an epithelial–mesenchymal spectrum, plus the parental line, one tumour each.",
 "GSE": "Human mammary cells (MCF10A, HuMEC) grown as colonies so EMT rises with distance from the colony edge; Mock vs TGF-β. Induced EMT in culture, no replicates."}

def pick_fit(ds, gs, level):
    d = B[(B.dataset == ds) & (B.geneset == gs) & (B.level == level)].copy()
    d["pri"] = d.arm.map(lambda a: FB.ARM_PRI.get(a, 9))
    d = d.sort_values(["p", "pri", "n_cells"], ascending=[True, True, False])
    return d.iloc[0], d

def figure(ds, gs, level, row, obs):
    pts, V, o = A.cells_for(row, obs); gen, inside, near = A.rules(pts, V)
    order = FB.STATES.get(ds, ["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"])
    st = o["cluster_annotation"].astype(str).values; sidx = np.array([order.index(x) if x in order else 0 for x in st])
    ann = CA[(CA.dataset == ds) & (CA.unit == row.unit) & (CA.geneset == gs) & (CA.arm == row.arm) & (CA.k == row.k)]
    fig, axs = plt.subplots(1, 2, figsize=(9.2, 4.0), gridspec_kw=dict(left=.01, right=.99, top=.9, bottom=.2, wspace=.05))
    rng = np.random.default_rng(0); idx = rng.permutation(len(o))[:4000]
    for j in range(4):
        m = sidx[idx] == j; axs[0].scatter(pts[idx][m, 0], pts[idx][m, 1], s=3, color=FB.SCOL[j], alpha=.75, lw=0, label=order[j])
    g = gen["groves5"][idx]
    axs[1].scatter(pts[idx][g, 0], pts[idx][g, 1], s=2, color=GEN_COL, alpha=.7, lw=0, label="generalist")
    for j in range(V.shape[0]):
        m = (~g) & (near[idx] == j)
        lab = f"V{j + 1}"
        a_ = ann[ann.corner == lab]
        if len(a_):
            c_ = a_.iloc[0]
            lab += (f": {c_.top_sample} (identity)" if (level.startswith("pooled") and c_.top_sample_pct >= 70) else f": {c_.task}")
        axs[1].scatter(pts[idx][m, 0], pts[idx][m, 1], s=5, color=FB.VCOL[j % 7], alpha=.95, lw=0, label=lab)
    for ax in axs:
        FB.outline(ax, V[:, :2], OUTL_COL); ax.set_xticks([]); ax.set_yticks([])
        for s_ in ax.spines.values(): s_.set_color(SPINE)
    axs[0].set_title("cells by EMT state", fontsize=9, color=MUT); axs[1].set_title("closest-5% specialists by corner (grey = generalist)", fontsize=9, color=MUT)
    axs[0].legend(loc="upper center", bbox_to_anchor=(.5, -.02), ncol=2, frameon=False, fontsize=7.5, labelcolor=FGC, markerscale=3)
    axs[1].legend(loc="upper center", bbox_to_anchor=(.5, -.02), ncol=2, frameon=False, fontsize=7, labelcolor=FGC, markerscale=3)
    fig.suptitle(f"{NAME[ds]} | {GS[gs]} | {row.unit} ({row.arm}), k={int(row.k)}, p={row.p:.3g}", fontsize=9.5, color=FGC)
    fn = OUT / (f"{ds}_{gs}_{level}".replace("/", "-") + ".png"); fig.savefig(fn, dpi=150, facecolor=fig.get_facecolor()); plt.close(fig)
    ann = ann.copy(); ann["multiome_pct"] = np.nan
    if ds == "SCC" and "assay" in o and level == "pooled":   # within one multiome tumour every cell is multiome: no flag
        d = np.linalg.norm(pts[:, None, :] - V[None, :, :], axis=2); cut = np.quantile(d, 0.05, axis=0)
        for j in range(V.shape[0]):
            m = d[:, j] <= cut[j]; ann.loc[ann.corner == f"V{j + 1}", "multiome_pct"] = 100 * (o["assay"].values[m] == "multiome").mean()
    return fn, ann

def corner_text(ann, level):
    items = []
    for _, c in ann.iterrows():
        ev = f"{c.top_set} (+{c.top_effect} SD)" if c.top_effect > 0.5 else "no Hallmark set above 0.5 SD"
        extra = []
        if not np.isnan(c.cycling_pct): extra.append(f"{c.cycling_pct:.0f}% S/G2M")
        extra.append(f"{c.top_state_pct:.0f}% {c.top_state}")
        if c.top_sample_pct >= 60: extra.append(f"{c.top_sample_pct:.0f}% from {c.top_sample}")
        if not np.isnan(c.get("multiome_pct", np.nan)) and (c.multiome_pct >= 70): extra.append(f"<b>{c.multiome_pct:.0f}% multiome cells: technical / assay axis</b>")
        items.append(f"<li><b>{c.corner}: {c.task}</b> <span class='m'>({ev}; {'; '.join(extra)}; genes {c.top_genes})</span></li>")
    return "<ul>" + "".join(items) + "</ul>"

def reading(ds, gs, level, row, d_all, ann):
    sig = row.p <= 0.05
    n_sig = int((d_all.p <= 0.05).sum()); n = len(d_all)
    t = []
    if level == "pooled":
        t.append(f"<b>Is there a simplex?</b> {'Yes' if sig else 'No'}: the best pooled fit is k={int(row.k)}, p={row.p:.3g} ({'passes' if sig else 'fails'} the permutation test).")
        ident = ann[ann.top_sample_pct >= 70]
        if sig and len(ident) >= max(1, len(ann) - 1):
            t.append(f"<b>But</b> {len(ident)} of {len(ann)} corners are dominated by one sample/tumour/line ({', '.join(ident.top_sample.astype(str).unique())}): the simplex mostly reflects <b>sample identity</b>, not a shared EMT trade-off.")
    else:
        t.append(f"<b>Is there a simplex?</b> {n_sig} of {n} within-unit fits reach p ≤ 0.05 for this gene set. Shown: the most significant, {row.unit} ({row.arm}), k={int(row.k)}, p={row.p:.3g}{' (not significant)' if not sig else ''}.")
    tasks = ann.task.value_counts()
    note = {"hallmark": " The Hallmark EMT set itself is left out when naming, because these corners were fitted on those genes.", "ks": " The KS genes overlap the EMT set, so an 'EMT / matrix' label is partly expected here."}.get(gs, "")
    if len(tasks): t.append("<b>Corners as tasks:</b> " + "; ".join(f"{k}" for k in ann.task) + ". <span class='m'>(Named from the Hallmark gene set that is highest at each corner, with phase, state and marker genes as support.)" + note + "</span>")
    return "".join(f"<p>{x}</p>" for x in t)

def rules_html(row):
    rows = ""
    for nm, k in (("Rule 1: gap (top weight − second ≥ 0.35)", "gap"), ("Rule 2: closest 5% of cells to each corner", "groves5"), ("Rule 3: weight > 0.95 on one corner", "w95")):
        rows += f"<tr><td>{nm}</td><td>{row[f'hyb_gen_{k}']:.0f}%</td><td>{row[f'gen_hyb_{k}']:.0f}%</td></tr>"
    return (f"<table class='r'><tr><th>rule</th><th>hybrids that are generalists</th><th>generalists that are hybrid</th></tr>{rows}</table>"
            f"<p class='m'>Hybrids are {row.hybrid_pct:.0f}% of the cells in this fit ({int(row.n_cells):,} cells), so a generalist group with ~{row.hybrid_pct:.0f}% hybrids means no enrichment. {row.inside_pct:.0f}% of cells lie inside the simplex.</p>")

CSS = """<style>:root{color-scheme:dark}body{background:#14171c;color:#e8e8e8;font:13.5px/1.5 -apple-system,Helvetica,Arial,sans-serif;margin:0;-webkit-print-color-adjust:exact;print-color-adjust:exact}
main{max-width:1250px;margin:0 auto;padding:24px 18px}h1{font-size:26px}h2{font-size:20px;border-bottom:1px solid #444;padding-bottom:4px;margin-top:34px}h3{font-size:15px;margin:6px 0}
.blk{display:grid;grid-template-columns:46% 54%;gap:16px;background:#1b1f26;border:1px solid #2a2f37;border-radius:10px;padding:12px;margin:14px 0;page-break-inside:avoid;align-items:start}
.blk img{width:100%;border-radius:6px}.m{color:#9aa0a6}p{margin:4px 0}ul{margin:4px 0 4px 18px;padding:0}li{margin:3px 0}
table.r{border-collapse:collapse;width:100%;font-size:12px;margin:6px 0}table.r td,table.r th{border:1px solid #333;padding:3px 6px;text-align:left}table.r th{background:#2b4a6f}
.tag{display:inline-block;padding:1px 8px;border-radius:9px;font-size:11px;font-weight:700;margin-right:6px}.ok{background:#2f855a;color:#0b0d10}.no{background:#c53030;color:#fff}
.def{background:#1b1f26;border:1px solid #2a2f37;border-radius:10px;padding:10px 14px;margin:10px 0}@page{size:A3 landscape;margin:10mm;background:#14171c}</style>"""

def main():
    obs = A.load_obs(); H = [f"<!doctype html><meta charset=utf-8><title>Explainer</title>{CSS}<main>"]
    H.append("<h1>Simplex results explained: five datasets, three gene sets, pooled and within unit</h1><p class='m'>Every block: the simplex figure on the left, the explanation on the right. Numbers come from <code>best_fits.csv</code> and <code>corner_annotations.csv</code>.</p>")
    H.append("""<div class='def'><b>How to read a block.</b> (1) <b>Simplex</b>: PCA, then archetype fit with k corners; significant when the real simplex is more triangle-like than 100 to 200 shuffled copies (p ≤ 0.05). (2) <b>Corners as tasks</b>: the 5% of cells nearest each corner are scored on the 50 Hallmark gene sets; the top set names the programme. A corner that is >70% one sample is an <i>identity</i> corner, not a task. (3) <b>Hybrid vs generalist</b>: a cell is a specialist if a corner claims it and a generalist otherwise, under three rules. Left figure panel: cells coloured by EMT state; right panel: specialists coloured by corner.
<br><b>The three rules.</b> Rule 1 gap: top barycentric weight minus second ≥ 0.35 (used in the Pal work). Rule 2 closest 5%: the 5% of cells nearest each corner (Groves 2022, main text). Rule 3 weight > 0.95 on one corner (Groves 2022, methods section). In the code all three come from earlier work (see <code>06_groves_specialists.py</code>); Rule 3 is the strictest and calls about 100% of cells generalist when the corners sit inside the cloud (delta = 0).
<br><b>Two numbers per rule.</b> % of hybrid cells that are generalists, and % of generalist cells that are hybrid (compare with the hybrid share of all cells).</div>""")
    rows = ""
    for ds in ("SCC", "COOK", "PDAC", "T4", "GSE"):
        cells = ""
        for gs in ("all", "ks", "hallmark"):
            pr, _ = pick_fit(ds, gs, "pooled"); _, dw = pick_fit(ds, gs, "within")
            cells += f"<td>{'<b>' if pr.p <= .05 else ''}k={int(pr.k)}, p={pr.p:.3g}{'</b>' if pr.p <= .05 else ''}</td><td>{int((dw.p <= .05).sum())}/{len(dw)}</td>"
        rows += f"<tr><td>{NAME[ds]}</td>{cells}</tr>"
    H.append("<div class='def'><b>At a glance.</b> Best pooled simplex (bold = significant) and the number of within-unit fits with p \u2264 0.05 (of all within-unit fits, both arms), per gene set.<table class='r'><tr><th rowspan=2>dataset</th><th colspan=2>A. Entire gene set</th><th colspan=2>B. KS genes</th><th colspan=2>C. Hallmark EMT genes</th></tr><tr><th>pooled</th><th>within</th><th>pooled</th><th>within</th><th>pooled</th><th>within</th></tr>" + rows + "</table></div>")
    for ds in ("SCC", "COOK", "PDAC", "T4", "GSE"):
        H.append(f"<h2>{NAME[ds]}</h2><p>{DS_MEANING[ds]}</p>")
        for gs in ("all", "ks", "hallmark"):
            for level, lab in (("pooled", "Pooled fit"), ("within", "Within unit")):
                row, dall = pick_fit(ds, gs, level)
                fn, ann = figure(ds, gs, level, row, obs)
                tag = f"<span class='tag {'ok' if row.p <= 0.05 else 'no'}'>{'significant' if row.p <= 0.05 else 'not significant'}</span>"
                H.append(f"<div class='blk'><div><img src='figures/explainer/{fn.name}'></div><div><h3>{GS[gs]} — {lab} {tag}</h3>{reading(ds, gs, level, row, dall, ann)}"
                         f"{corner_text(ann, level)}{rules_html(row)}</div></div>")
        print(ds, "done", flush=True)
    sig = B[B.p <= 0.05]; Cs = CA.merge(sig[["dataset", "level", "unit", "geneset", "arm", "k"]], on=["dataset", "level", "unit", "geneset", "arm", "k"])
    tk = []
    for gs in ("all", "ks", "hallmark"):
        b = B[B.geneset == gs]; tk.append(f"<li><b>{GS[gs]}:</b> {int((b.p <= .05).sum())} of {len(b)} fits significant (pooled {int(((b.level == 'pooled') & (b.p <= .05)).sum())}/{int((b.level == 'pooled').sum())}, within unit {int(((b.level == 'within') & (b.p <= .05)).sum())}/{int((b.level == 'within').sum())}).</li>")
    pc = Cs[Cs.level == "pooled"]; ident = (pc.top_sample_pct >= 70).mean() * 100
    wc = Cs[Cs.level == "within"].task.value_counts().head(6)
    tk.append(f"<li><b>Pooled corners are mostly identity.</b> {ident:.0f}% of the corners of significant pooled simplices are \u2265 70% one sample, tumour, line or clone: pooling mixes samples, and the simplex then separates them.</li>")
    tk.append("<li><b>Within a unit the corners are programmes, not epithelial / hybrid / mesenchymal states.</b> Most frequent corner names across significant within-unit fits: " + "; ".join(f"{k} ({v})" for k, v in wc.items()) + f". Of the {len(Cs)} corners in significant fits, {int(((Cs.top_state.isin(['Early-Hybrid', 'Hybrid', 'Late-Hybrid'])) & (Cs.top_state_pct >= 60)).sum())} ({100 * ((Cs.top_state.isin(['Early-Hybrid', 'Hybrid', 'Late-Hybrid'])) & (Cs.top_state_pct >= 60)).mean():.0f}%) have a hybrid state as >= 60% of their specialists, and those are mainly cycling / MYC-growth corners (proliferating hybrid-state cells; hybrids are about half the cells, so this is near the base rate). Corners named 'EMT / matrix' are mostly mesenchymal or late-EMT ({int(Cs[Cs.task == 'EMT / matrix'].top_state.isin(['Mesenchymal', 'Late-EMT']).sum())} of {int((Cs.task == 'EMT / matrix').sum())}), not hybrid.</li>")
    sg = sig
    tk.append(f"<li><b>Hybrids are mostly generalists.</b> Median % of hybrids that are generalists across {len(sg)} significant fits: Rule 1 (gap) {sg.hyb_gen_gap.median():.0f}%, Rule 2 (closest 5%) {sg.hyb_gen_groves5.median():.0f}%, Rule 3 (weight > 0.95) {sg.hyb_gen_w95.median():.0f}%. Median % of generalists that are hybrid: {sg.gen_hyb_gap.median():.0f}%, {sg.gen_hyb_groves5.median():.0f}%, {sg.gen_hyb_w95.median():.0f}%, against a hybrid share of {sg.hybrid_pct.median():.0f}% of cells: generalists are barely enriched for hybrids.</li>")
    tk.append("<li><b>The three rules differ in strictness, not in answer.</b> Rule 3 is near 100% by construction; Rule 2 is capped by taking only 5% of cells per corner; Rule 1 is the most informative and the most variable (its lowest values come from fits whose corners are cell lines or clones).</li>")
    H.append("<h2>What it all says</h2><div class='def'><ul>" + "".join(tk) + "</ul><p class='m'>Caveats: best-of-many selection across arms and k; Hallmark, cleaned-pooled and SCC scRNA/six-tumour pooled fits used 100 shuffles and k = 3 to 5; hybrid labels for PDAC, 4T1, GSE114687 and Cook are cut from the same KS score; within-unit fits are capped at 4,000 cells.</p></div>")
    H.append("</main>")
    (HERE / "EXPLAINER_results.html").write_text("".join(H))

if __name__ == "__main__":
    main()
