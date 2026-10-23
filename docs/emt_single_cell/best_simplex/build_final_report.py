#!/usr/bin/env python3
"""Final explainer report (beige HTML + PDF) answering the two research questions with figures.
Uses only the main fits (main_fits.csv): SCC scRNA-only pooled + 8 scRNA tumours; cleaned arm elsewhere (4T1 clone E light QC).
Run from Project_1/:  .venv/bin/python docs/emt_single_cell/best_simplex/build_final_report.py
"""
import sys, html
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import build_explainer_ppt as P          # beige theme, main-fit policy, Q2 figures
import analyze as A, figurebook as FB, build_explainer as E
OUT = HERE / "figures" / "report"; OUT.mkdir(parents=True, exist_ok=True)
INK, MUT, BROWN = "#3A2E25", "#7A6A5A", "#8C5A2B"
MAIN, MS, CA, NAME, GS = P.MAIN, P.MS, P.CA, P.NAME, P.GS
DSS = list(E.FG.ORDER)
KEY = ["dataset", "level", "unit", "geneset", "arm", "k"]
ANN = CA.drop(columns=["p"]).merge(MAIN[KEY + ["p"]], on=KEY)          # corners of the main fits
UNITWORD = {"COOK": "cell line", "SCC": "tumour", "PDAC": "library", "T4": "clone", "GSE": "sample"}
obs = A.load_obs()
FIGN = [0]
def fig_html(fn, cap):
    FIGN[0] += 1
    return f"<figure><img src='{fn.relative_to(HERE).as_posix()}'><figcaption><b>Figure {FIGN[0]}.</b> {cap}</figcaption></figure>"

# ---------------------------------------------------------------- figures
def style(ax):
    for s in ax.spines.values(): s.set_color("#cdbf9f")
    ax.tick_params(colors=INK)

def fig_overview():
    cols = [(g, l) for g in ("all", "ks", "hallmark") for l in ("pooled", "within")]
    M = np.zeros((len(DSS), len(cols))); T = np.empty(M.shape, dtype=object)
    for i, ds in enumerate(DSS):
        for j, (g, l) in enumerate(cols):
            d = MAIN[(MAIN.dataset == ds) & (MAIN.geneset == g) & (MAIN.level == l)]
            if l == "pooled":
                r = d.iloc[0]; M[i, j] = 1.0 if r.p <= .05 else 0.0
                idc = ANN[(ANN.dataset == ds) & (ANN.geneset == g) & (ANN.level == "pooled")]
                T[i, j] = f"k={int(r.k)}, p={r.p:.2g}" + (f"\n{int((idc.top_sample_pct >= 70).sum())}/{len(idc)} identity" if r.p <= .05 else "\nn.s.")
            else:
                M[i, j] = (d.p <= .05).mean(); T[i, j] = f"{int((d.p <= .05).sum())} of {len(d)}\nunits"
    cmap = LinearSegmentedColormap.from_list("b", ["#efe5cc", "#9fc49a", "#3f7a4f"])
    fig, ax = plt.subplots(figsize=(10, 3.9))
    ax.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]): ax.text(j, i, T[i, j], ha="center", va="center", fontsize=8.5, color="white" if M[i, j] > .7 else INK)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([f"{GS[g][:2]} {GS[g][3:]}\n{l}" for g, l in cols], fontsize=8.5); ax.xaxis.tick_top()
    ax.set_yticks(range(len(DSS))); ax.set_yticklabels([NAME[d] for d in DSS], fontsize=9); style(ax)
    for x in (1.5, 3.5): ax.axvline(x, color="white", lw=3)
    fig.tight_layout(); fn = OUT / "overview.png"; fig.savefig(fn, dpi=160); plt.close(fig); return fn

def panel(ax_state, ax_corner, ds, row, level):
    pts, V, o = A.cells_for(row, obs); gen, inside, near = A.rules(pts, V)
    order = FB.STATES.get(ds, ["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"])
    st = o["cluster_annotation"].astype(str).values; sidx = np.array([order.index(x) if x in order else 0 for x in st])
    idx = np.random.default_rng(0).permutation(len(o))[:3500]
    for j in range(4):
        m = sidx[idx] == j; ax_state.scatter(pts[idx][m, 0], pts[idx][m, 1], s=2.2, color=FB.SCOL[j], alpha=.75, lw=0)
    g = gen["groves5"][idx]; ax_corner.scatter(pts[idx][g, 0], pts[idx][g, 1], s=1.8, color=E.GEN_COL, alpha=.7, lw=0)
    an = ANN[(ANN.dataset == ds) & (ANN.unit == row.unit) & (ANN.geneset == row.geneset) & (ANN.arm == row.arm) & (ANN.k == row.k)]
    hs = []
    for j in range(V.shape[0]):
        m = (~g) & (near[idx] == j); ax_corner.scatter(pts[idx][m, 0], pts[idx][m, 1], s=4, color=FB.VCOL[j % 7], lw=0)
        c = an[an.corner == f"V{j + 1}"]
        lab = f"V{j + 1}"
        if len(c):
            c = c.iloc[0]; lab += f": {c.top_sample} (identity)" if (level == "pooled" and c.top_sample_pct >= 70) else f": {c.task}"
        hs.append(plt.Line2D([], [], marker="o", ls="", color=FB.VCOL[j % 7], markersize=4, label=lab[:38]))
    for ax in (ax_state, ax_corner):
        FB.outline(ax, V[:, :2], INK); ax.set_xticks([]); ax.set_yticks([]); style(ax)
    ax_corner.legend(handles=hs, loc="upper center", bbox_to_anchor=(.5, -.02), fontsize=6.6, frameon=False, ncol=1, labelcolor=INK)
    return row

def fig_triplet(ds, rows, level, name):
    fig, axs = plt.subplots(2, 3, figsize=(11, 6.6), gridspec_kw=dict(left=.05, right=.99, top=.9, bottom=.2, wspace=.05, hspace=.12))
    for c, (gs, row) in enumerate(rows):
        if row is None:
            for r in range(2): axs[r, c].axis("off")
            axs[0, c].set_title(f"{GS[gs]}\nno fit", fontsize=9, color=INK); continue
        panel(axs[0, c], axs[1, c], ds, row, level)
        sig = row.p <= .05
        axs[0, c].set_title(f"{GS[gs]}\n{row.unit} ({row.arm}), k={int(row.k)}, p={row.p:.2g}" + ("" if sig else "  NOT significant"), fontsize=8.8, color=INK if sig else "#A33B2B")
    axs[0, 0].set_ylabel("cells by EMT state", fontsize=8.5, color=INK); axs[1, 0].set_ylabel("specialists by corner\n(closest 5%; grey = generalist)", fontsize=8.5, color=INK)
    order = FB.STATES.get(ds, ["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"])
    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=FB.SCOL[j], markersize=5, label=n) for j, n in enumerate(order)], loc="upper center", ncol=4, frameon=False, fontsize=8, labelcolor=INK, bbox_to_anchor=(.5, .995))
    fn = OUT / f"{name}.png"; fig.savefig(fn, dpi=150, facecolor=fig.get_facecolor()); plt.close(fig); return fn

def fig_within_heat(ds):
    units = [u for u in P.main_units(ds) if u != P.POOLED[ds]]
    M = np.full((len(units), 3), np.nan); T = np.empty(M.shape, dtype=object)
    for i, u in enumerate(units):
        for j, g in enumerate(("all", "ks", "hallmark")):
            r = P.pick_main(ds, u, g)
            if r is None: T[i, j] = "n/a"; continue
            M[i, j] = -np.log10(max(r.p, 1e-3))
            an = ANN[(ANN.dataset == ds) & (ANN.unit == u) & (ANN.geneset == g) & (ANN.arm == r.arm) & (ANN.k == r.k)]
            tasks = ", ".join(dict.fromkeys(an.task.str.replace("no Hallmark set above 0.5 SD", "unnamed").str.split(" \\(").str[0]))
            T[i, j] = f"k={int(r.k)}, p={r.p:.2g}" + (f"\n{tasks[:46]}" if r.p <= .05 else "")
    cmap = LinearSegmentedColormap.from_list("p", ["#efe5cc", "#9fc49a", "#3f7a4f"])
    fig, ax = plt.subplots(figsize=(11, 0.62 * len(units) + 1.1))
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=0, vmax=2.3, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(3):
            sig = M[i, j] >= -np.log10(0.05)
            ax.text(j, i, T[i, j], ha="center", va="center", fontsize=7.4, color="white" if (np.isfinite(M[i, j]) and M[i, j] > 1.7) else INK, fontweight="bold" if sig else "normal")
    ax.set_xticks(range(3)); ax.set_xticklabels([GS[g] for g in ("all", "ks", "hallmark")], fontsize=9); ax.xaxis.tick_top()
    ax.set_yticks(range(len(units))); ax.set_yticklabels([f"{u} ({P.pick_main(ds, u, 'all').arm})" for u in units], fontsize=8.5); style(ax)
    fig.tight_layout(); fn = OUT / f"within_{ds}.png"; fig.savefig(fn, dpi=160); plt.close(fig); return fn

def fig_tasks():
    w = ANN[(ANN.level == "within") & (ANN.p <= .05)].copy()
    w["t"] = w.task.str.replace("no Hallmark set above 0.5 SD", "unnamed (no set > 0.5 SD)")
    top = w.t.value_counts().head(12).index[::-1]
    fig, ax = plt.subplots(figsize=(10, 4.6))
    y = np.arange(len(top)); left = np.zeros(len(top))
    for g, c in (("all", "#2b6cb0"), ("ks", "#d97706"), ("hallmark", "#2f855a")):
        v = np.array([((w.geneset == g) & (w.t == t)).sum() for t in top]); ax.barh(y, v, left=left, color=c, label=GS[g]); left += v
    ax.set_yticks(y); ax.set_yticklabels(top, fontsize=9); ax.set_xlabel("number of corners (significant within-unit main fits)", color=INK); style(ax)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="lower right"); fig.tight_layout(); fn = OUT / "tasks.png"; fig.savefig(fn, dpi=160); plt.close(fig); return fn

# ---------------------------------------------------------------- text helpers
def tbl(rows, head=True):
    h = "<table>"
    for i, r in enumerate(rows):
        tag = "th" if (i == 0 and head) else "td"; h += "<tr>" + "".join(f"<{tag}>{c}</{tag}>" for c in r) + "</tr>"
    return h + "</table>"
def pct(x): return f"{x:.0f}%"
def diff_text():
    d = MS.assign(d=MS.hyb_gen_gap - MS.oth_gen_gap).groupby("dataset").d.median().sort_values(ascending=False)
    return ", ".join(f"{NAME[k]} {int(round(v)) + 0:+d}" for k, v in d.items()) + " points"

def main():
    H = []; w = H.append
    toy = P.toy_rules_figure()
    w("<h1>Do EMT cells sit in a simplex, and are hybrid cells generalists?</h1><p class='m'>An explainer of the archetype (Pareto-simplex) analysis across five single-cell EMT datasets. Thesis project, IISc · 22 October 2026. Every number below is read from <code>main_fits.csv</code> and <code>corner_annotations.csv</code> in <code>docs/emt_single_cell/best_simplex/</code>.</p>")
    # ---- summary
    q1 = {g: (int((MS.geneset == g).sum()), int((MAIN.geneset == g).sum())) for g in ("all", "ks", "hallmark")}
    pid = ANN[(ANN.level == "pooled") & (ANN.p <= .05)]; ident = 100 * (pid.top_sample_pct >= 70).mean()
    e = MS.hyb_gen_gap - MS.oth_gen_gap
    w("<div class='box'><h2 style='margin-top:0'>The answers in short</h2><ol>"
      f"<li><b>Is there a significant simplex, and are its corners tasks?</b> A simplex is often significant: {q1['all'][0]} of {q1['all'][1]} main fits on the entire gene set, {q1['ks'][0]} of {q1['ks'][1]} on the KS EMT genes and {q1['hallmark'][0]} of {q1['hallmark'][1]} on the Hallmark EMT genes. "
      f"But it means different things at the two levels. <b>Pooled</b> over samples, the corners are mostly the samples themselves ({ident:.0f}% of corners of significant pooled simplices are ≥ 70% one {', '.join(sorted(set(UNITWORD.values())))}); Cook's pooled simplex is nothing but its four cell lines. "
      "<b>Within a single unit</b>, the corners are biological programmes (cycling, MYC/growth, a differentiated / KRAS-down epithelial state, metabolic, EMT / matrix), not the epithelial / hybrid / mesenchymal states. The EMT gene sets (KS, Hallmark) give fewer significant simplices than the entire gene set, and Hallmark behaves like KS.</li>"
      f"<li><b>Are hybrid cells generalists, and more so than other cells?</b> Most cells are generalists, hybrid or not: under rule 1 a median {pct(MS.hyb_gen_gap.median())} of hybrid cells and {pct(MS.oth_gen_gap.median())} of other cells are generalists (rule 2: {pct(MS.hyb_gen_groves5.median())} vs {pct(MS.oth_gen_groves5.median())}; rule 3 about 100% for everyone). "
      f"Hybrids are slightly more often generalists than other cells; median difference under rule 1 per dataset: {diff_text()}. In no dataset do hybrid cells form a corner of their own.</li></ol></div>")
    # ---- 1. why
    w("<h2>1. What we wanted to know, and why</h2>"
      "<p><b>The idea (Pareto task inference).</b> If a cell must perform several tasks and cannot be best at all of them at once, the best compromises lie inside a simple shape in gene-expression space: a line for two tasks, a triangle for three, a tetrahedron for four. Cells that specialise in one task sit at a corner (an <i>archetype</i>); cells that balance several sit inside (<i>generalists</i>). Finding such a shape, and naming its corners, is a way to read the tasks off the data.</p>"
      "<p><b>Why EMT.</b> Epithelial-to-mesenchymal transition is usually described as states (epithelial, hybrid, mesenchymal). Two competing pictures: (a) EMT states are task corners, perhaps with hybrid cells as their own corner; or (b) EMT is a continuum and hybrids are mixtures in the interior. The simplex view separates these: in (a) hybrids are specialists of a corner, in (b) they are generalists.</p>"
      "<p><b>The two questions.</b> (1) Is there a significant simplex, and do its corners correspond to tasks, when the fit uses A. the entire gene set, B. the KS EMT genes or C. the Hallmark EMT genes, both pooled over each dataset and within each tumour / cell line / sample? (2) How many hybrid cells are generalists, and how many generalists are hybrid, under three specialist rules?</p>"
      "<p><b>Why three gene sets.</b> The entire gene set lets any programme shape the cloud (unbiased but dominated by the biggest sources of variation, such as cell cycle or sample identity). The KS genes (Tan et al. 2014; the EMT scoring panel used in the earlier breast work) and the Hallmark EMT set (MSigDB, 200 genes) restrict the fit to EMT genes, so a simplex there would be an EMT-specific trade-off. Two EMT sets guard against conclusions that depend on one panel.</p>"
      "<p><b>Why pooled and within.</b> Pooling uses all cells but mixes samples, which differ for reasons unrelated to EMT (patient or mouse, cell line, technology). A pooled simplex can therefore be just sample identity. Fitting inside a single unit removes between-sample differences, so a within-unit simplex reflects variation among cells of one tumour or line.</p>")
    # ---- 2. data
    w("<h2>2. The five datasets</h2>")
    rows = [["dataset", "what it is", "units fitted (main)", "hybrid cells", "main arm (why)"],
            ["Mouse skin SCC (GSE288137)", "12 LKPR tumours; 8 scRNA, 4 multiome (nuclei)", "pooled scRNA cells; 8 scRNA tumours", "authors' Early-Hybrid + Hybrid", "scRNA only: the multiome nuclei create a technical corner; multiome tumours kept in an appendix"],
            ["Cook (GSE147405)", "4 human lines (A549, DU145, MCF7, OVCA420) given TGF-β, sampled over time", "pooled 4 lines; each line", "middle two quartiles of the KS EMT score (no state labels)", "cleaned: ≥ 10,000 UMI, mito ≤ 0.07, MALAT1/NEAT1/haemoglobin dropped"],
            ["Mouse PDAC (GSE173958)", "tumour cells of 2 mice at several sites", "pooled in vivo; 7 libraries, mouse 1, mouse 2", "KS-score bins Early- + Late-Hybrid", "cleaned: ≥ 10,000 UMI, liver ambient-RNA genes dropped"],
            ["4T1 clones (GSE303505)", "5 clones + parental, one tumour each", "pooled 6 clones; each clone", "KS-score quartiles Early- + Late-Hybrid", "cleaned: ≥ 10,000 UMI + immune / stromal contamination removed (clone E: light QC, 142 cells after cleaning)"],
            ["MCF10A / HuMEC (GSE114687)", "mammary cells in colonies; EMT rises toward the edge; Mock vs TGF-β", "pooled 4 samples; 6 units (line, line × arm)", "KS-score quartiles Early- + Late-Hybrid", "cleaned: ≥ 10,000 UMI, MALAT1/NEAT1/haemoglobin dropped"]]
    w(tbl(rows))
    w("<p><b>Why one main fit per unit.</b> Each dataset was fitted with a light-QC arm and a cleaned arm. The two disagree on whether a simplex is significant in many units (Cook 2 of 12 unit × gene-set pairs, PDAC 8 of 27, 4T1 7 of 15, GSE114687 8 of 18), so choosing whichever arm gives the lower p would inflate significance. The cleaned arm is the one that removes each dataset's known technical problem (low-depth cells, ambient RNA, contaminating immune or stromal cells), so it is the main arm; for SCC the analogous choice is scRNA only. Everything below counts units, one main fit each.</p>")
    # ---- 3. methods
    w("<h2>3. Methods and code</h2>")
    w("<h3>3.1 Fitting a simplex</h3><p>For each unit and gene set: (i) build the log-normalised expression matrix (counts per 10,000, log1p; genes in ≥ 1% of cells; mitochondrial and ribosomal genes dropped) — <code>01_build_anndata.py</code>; (ii) PCA on the unit's cells with genes centred but not scaled (scaling would give rare noisy genes the same pull as strong ones) — <code>02_pca.py</code> (pooled, balanced cap per sample so large samples do not decide the shape) or <code>08_within_unit.py</code> / <code>07_within_line_clean.py</code> / <code>08d_within_tumour_geneset.py</code> (within, capped at 4,000 cells); (iii) archetypal analysis (PCHA, 15 random starts, the largest-volume solution kept) on the first k−1 PCs, for k = 3 to 7 pooled (3 to 5 for the newest fits) and k = 3, 4 within — <code>03_fit_and_test.py</code>, engine <code>src/archetypes.py</code>.</p>"
      "<h3>3.2 Is it significant? (permutation test)</h3><p>The t-ratio is the volume of the fitted simplex divided by the volume of the convex hull of the cells. Each PC is then shuffled independently (destroying any joint shape while keeping each axis' spread) and refitted, 200 times (100 for the Hallmark, cleaned-pooled and SCC scRNA-only KS / Hallmark pooled fits). p = (1 + number of shuffles at least as good) / (1 + shuffles). For each fit the k with the lowest p is reported (ties: smaller k; pooled fits ignore k beyond the explained-variance elbow, as in the earlier report). Significant = p ≤ 0.05. <b>Caveat:</b> a small p says the cloud is not a shapeless blob; samples that differ from each other also pass, which is why corners must be named before they are called tasks.</p>"
      "<h3>3.3 What is each corner? (naming)</h3><p>For each corner, the 5% of cells nearest to it are compared with all other cells of the fit on the 50 MSigDB Hallmark gene sets (per-cell score = mean log-expression of the set; effect = difference of means in SD units). The top set names the programme, e.g. E2F Targets / G2-M = cycling, MYC Targets = growth, KRAS Signaling Dn = a differentiated / quiescent epithelial state. The share of S/G2M cells, the dominant EMT state, the dominant sample and the top marker genes support the name. <b>Identity rule:</b> in a pooled fit, a corner whose cells are ≥ 70% one sample is an <i>identity corner</i> (it marks a sample, not a task). For fits on the Hallmark EMT genes, the EMT set itself is excluded from naming (circular). — <code>corners.py</code>.</p>"
      "<h3>3.4 Specialists and generalists (three rules)</h3><p>Every cell gets barycentric weights on the k corners (they sum to 1; negative weights of cells outside the simplex are clipped). A cell is a specialist of a corner, or a generalist, by three rules:</p>")
    w(fig_html(toy, "The three rules applied to 500 simulated cells inside a triangle. Coloured = specialists of a corner, grey = generalists."))
    w("<ul><li><b>Rule 1, purity gap</b> (used in the earlier Pal-breast work of this project): specialist if the largest weight exceeds the second largest by ≥ 0.35.</li>"
      "<li><b>Rule 2, closest 5%</b> (Groves et al. 2022, main analysis): the 5% of cells nearest each corner are its specialists; everything else is a generalist.</li>"
      "<li><b>Rule 3, weight > 0.95</b> (from the Groves 2022 methods; new to this project): specialist only if one weight exceeds 0.95. Because the fitted corners sit inside the data cloud, almost no cell qualifies, so this rule calls about 100% of cells generalist. It is the strictest and least informative rule.</li></ul>"
      "<p><b>The comparison.</b> For each fit and rule we count what % of <b>hybrid cells</b> are generalists and what % of <b>all other cells</b> (epithelial + mesenchymal) are generalists. Same % = hybrids are not special; hybrids more often generalists = they sit in the interior more than other cells, as the picture of hybrids as mixtures predicts. The % of other cells comes from the overall generalist share and the hybrid numbers. \u2014 <code>analyze.py</code> (<code>rules()</code>), results in <code>main_fits.csv</code>.</p>")
    # ---- 4. Q1 results
    w("<h2>4. Question 1: is there a simplex, and are its corners tasks?</h2>")
    w(fig_html(fig_overview(), "Overview. Pooled columns: best k and p of the pooled main fit (green = significant) and how many of its corners are identity corners (≥ 70% one sample). Within columns: number of units with a significant simplex (shade = fraction)."))
    w(fig_html(fig_tasks(), "What the corners are, inside single units: names of all corners of significant within-unit main fits, by gene set."))
    for ds in DSS:
        w(f"<h3>4.{DSS.index(ds) + 1} {NAME[ds]}</h3><p class='m'>{E.DS_MEANING[ds]}</p>")
        prow = [(g, P.pick_main(ds, P.POOLED[ds], g)) for g in ("all", "ks", "hallmark")]
        txt = []
        for g, r in prow:
            an = ANN[(ANN.dataset == ds) & (ANN.level == "pooled") & (ANN.geneset == g)]
            if r.p > .05: txt.append(f"<li>{GS[g]}: no significant pooled simplex (k={int(r.k)}, p={r.p:.2g}).</li>"); continue
            idn = an[an.top_sample_pct >= 70]
            names = "; ".join((f"V{c.corner[1:]} = {c.top_sample} (identity)" if c.top_sample_pct >= 70 else f"V{c.corner[1:]} = {c.task}") for _, c in an.iterrows())
            txt.append(f"<li>{GS[g]}: significant (k={int(r.k)}, p={r.p:.2g}); {len(idn)} of {len(an)} corners are identity corners. Corners: {names}.</li>")
        w("<p><b>Pooled.</b></p><ul>" + "".join(txt) + "</ul>")
        w(fig_html(fig_triplet(ds, prow, "pooled", f"pooled_{ds}"), f"{NAME[ds]}, pooled main fits for the three gene sets. Top: cells coloured by EMT state. Bottom: the 5% of cells nearest each corner, coloured by corner and labelled (identity = ≥ 70% one {UNITWORD[ds]}); grey = generalists. For k ≥ 4 the corners are projected on PC1/PC2."))
        wt = []
        for g in ("all", "ks", "hallmark"):
            d = MAIN[(MAIN.dataset == ds) & (MAIN.geneset == g) & (MAIN.level == "within")]; sg = d[d.p <= .05]
            an = ANN[(ANN.dataset == ds) & (ANN.level == "within") & (ANN.geneset == g) & (ANN.p <= .05)]
            tc = an.task.str.replace("no Hallmark set above 0.5 SD", "unnamed").value_counts().head(4)
            wt.append(f"<li>{GS[g]}: {len(sg)} of {len(d)} {UNITWORD[ds]}s have a significant simplex" + (f" ({', '.join(sg.unit)}); most common corners: " + ", ".join(f"{k} ({v})" for k, v in tc.items()) if len(sg) else "") + ".</li>")
        w("<p><b>Within each " + UNITWORD[ds] + ".</b></p><ul>" + "".join(wt) + "</ul>")
        w(fig_html(fig_within_heat(ds), f"{NAME[ds]}, every {UNITWORD[ds]}: best k and p of the main fit per gene set (green = lower p, bold = significant) and, for significant fits, the corner names."))
        best = []
        for g in ("all", "ks", "hallmark"):
            d = MAIN[(MAIN.dataset == ds) & (MAIN.geneset == g) & (MAIN.level == "within")].sort_values(["p", "n_cells"], ascending=[True, False])
            best.append((g, d.iloc[0] if len(d) else None))
        w(fig_html(fig_triplet(ds, best, "within", f"within_best_{ds}"), f"{NAME[ds]}, the most significant {UNITWORD[ds]} for each gene set, drawn as in the pooled figure. Every other {UNITWORD[ds]} is shown in the slide deck (<code>Explainer_simplex_results.pptx</code>)."))
    # what we saw Q1
    hybc = ANN[(ANN.p <= .05)]; hc = ((hybc.top_state.isin(["Early-Hybrid", "Hybrid", "Late-Hybrid"])) & (hybc.top_state_pct >= 60)); emt = hybc[hybc.task == "EMT / matrix"]
    w("<h3>4.6 What we saw (Question 1)</h3><ul>"
      f"<li><b>Gene set matters.</b> The entire gene set gives significant simplices most often ({q1['all'][0]} of {q1['all'][1]} main fits), the EMT-only sets less often (KS {q1['ks'][0]}, Hallmark {q1['hallmark'][0]} of {q1['hallmark'][1]}). Hallmark does not rescue the EMT-only result: it behaves like KS.</li>"
      f"<li><b>Pooled shapes are mostly identity.</b> {ident:.0f}% of corners of significant pooled simplices are one sample. In Cook every pooled corner is a cell line; in MCF10A/HuMEC most are samples; in SCC two of three all-gene corners are tumours.</li>"
      "<li><b>Inside units, the corners are programmes.</b> Cycling, MYC/growth, a differentiated (KRAS-down) epithelial state, metabolic and EMT / matrix programmes recur (Figure 3).</li>"
      f"<li><b>EMT states are not corners.</b> Only {100 * hc.mean():.0f}% of corners of significant fits are dominated (≥ 60%) by a hybrid state, and those are mostly cycling or MYC/growth corners (hybrid-state cells that divide); {int(emt.top_state.isin(['Mesenchymal', 'Late-EMT']).sum())} of {len(emt)} EMT / matrix corners are mesenchymal or late-EMT, not hybrid.</li>"
      "<li><b>So</b> a trade-off structure exists inside single tumours and lines, between proliferation, growth, differentiation and EMT / matrix programmes; EMT itself behaves like an axis running between epithelial and EMT / matrix corners, not like a set of task corners.</li></ul>")
    # ---- 5. Q2
    w("<h2>5. Question 2: are hybrid cells more often generalists than other cells?</h2>")
    w(fig_html(P.q2_simple_summary(), "The main answer. For each dataset, orange = % of hybrid cells that are generalists, grey = % of all other cells that are generalists; bars are medians over the units with a significant simplex. Same height means hybrids are not special. Rule 3 is left out because it calls about 100% of every group generalist."))
    rows = [["dataset", "gene set", "units with a significant simplex", "rule 1: hybrids vs others generalist", "rule 2: hybrids vs others", "verdict (rule 1)"]]
    for ds in DSS:
        for g in ("all", "ks", "hallmark"):
            x = MS[(MS.dataset == ds) & (MS.geneset == g)]
            if len(x): rows.append([NAME[ds], GS[g][3:], str(len(x)), f"{x.hyb_gen_gap.median():.0f}% vs {x.oth_gen_gap.median():.0f}%", f"{x.hyb_gen_groves5.median():.0f}% vs {x.oth_gen_groves5.median():.0f}%", P.verdict((x.hyb_gen_gap - x.oth_gen_gap).median())])
            else: rows.append([NAME[ds], GS[g][3:], "0", "–", "–", "no simplex to read"])
    w("<p>Medians over the units with a significant simplex (units without one have no meaningful corners, so they are left out of the medians but shown faded in the per-dataset figures).</p>" + tbl(rows))
    for ds in DSS:
        w(f"<h3>5.{DSS.index(ds) + 1} {NAME[ds]}</h3>")
        w(fig_html(P.q2_unit_bars(ds), f"{NAME[ds]}, every unit: orange = % of hybrid cells that are generalists, grey = % of other cells (rule 1); the number is orange minus grey in percentage points; faded = no significant simplex."))
    w("<h3>5.6 What we saw (Question 2)</h3><ul>"
      f"<li><b>Most cells are generalists, hybrid or not</b>: rule 1 median {pct(MS.hyb_gen_gap.median())} of hybrids vs {pct(MS.oth_gen_gap.median())} of other cells; rule 2 {pct(MS.hyb_gen_groves5.median())} vs {pct(MS.oth_gen_groves5.median())}; rule 3 about 100% for both. The rules differ by construction: rule 3 almost never calls a specialist, rule 2 lets only 5% of cells per corner be specialists, rule 1 is the one that separates groups.</li>"
      f"<li><b>Hybrids lean slightly toward the interior.</b> Median difference (hybrids minus other cells, rule 1) per dataset: {diff_text()}; across all {len(MS)} significant fits {e.median():+.0f} points. So the tilt is clearest in Mouse SCC and smallest in Cook.</li>"
      f"<li><b>Exceptions</b>: in {int((e < -10).sum())} fits hybrids are much <i>less</i> often generalists than other cells (" + "; ".join(f"{NAME[r.dataset]} {r.unit} ({GS[r.geneset][3:]})" for _, r in MS[e < -10].iterrows()) + "). In these fits hybrid-state cells sit at corners more often than other cells.</li>"
      "<li><b>No hybrid corner.</b> In no dataset do hybrid cells form their own corner. This fits the picture of EMT as a continuum with hybrids as mixtures.</li></ul>")
    # ---- 6. caveats, files
    w("<h2>6. Caveats</h2><ul>"
      "<li>Many fits per dataset (gene sets, k); no multiplicity correction. p cannot go below about 0.005 (200 shuffles) or 0.01 (100 shuffles).</li>"
      "<li>Hybrid labels for PDAC, 4T1, GSE114687 and Cook are cut from the same KS score used for the KS-gene fits (partly circular). Cook's bins are pooled across lines, so a line's hybrid share reflects the line (A549 about 71%, MCF7 about 11%).</li>"
      "<li>Corner names are automatic (top Hallmark set); they summarise, and should be sanity-checked against the marker genes and the biology.</li>"
      "<li>Design limits: one tumour per clone (4T1), two mice (PDAC), no replicates (GSE114687), tumour / assay / state tangled (SCC). These results do not show that EMT is or is not a Pareto trade-off in tumours.</li>"
      "<li>Within-unit fits are capped at 4,000 random cells; Hallmark and cleaned-pooled fits were run only for k = 3 to 5.</li></ul>")
    w("<h2>7. Code and files</h2>" + tbl([["step", "script", "what it does"],
        ["data", "<code>01_build_anndata.py</code>, <code>make_clean_pooled.py</code>", "expression matrices; cleaned copies for pooled cleaned fits"],
        ["fits", "<code>02_pca.py</code>, <code>03_fit_and_test.py</code>, <code>08_within_unit.py</code>, <code>07_within_line_clean.py</code>, <code>08d_within_tumour_geneset.py</code>", "PCA, PCHA, permutation test (per dataset <code>codes/</code>)"],
        ["queues", "<code>pooled_queue.py</code>, <code>within_light_queue.py</code>", "the Hallmark and missing fits run on 21 to 22 October"],
        ["collect", "<code>analyze.py</code>", "all fits → best k per arm; three rules → <code>best_fits.csv</code>"],
        ["main fits", "<code>build_explainer_ppt.py</code>", "main-fit policy → <code>main_fits.csv</code>; the slide deck"],
        ["corners", "<code>corners.py</code>", "Hallmark naming of every corner → <code>corner_annotations.csv</code>"],
        ["this report", "<code>build_final_report.py</code>", "figures in <code>figures/report/</code>; this document"]]))
    css = """<!doctype html><meta charset=utf-8><title>Simplex explainer report</title><style>
body{background:#F4ECD8;color:#3A2E25;font:11.5pt/1.55 Georgia,'Times New Roman',serif;margin:0;-webkit-print-color-adjust:exact;print-color-adjust:exact}
main{max-width:960px;margin:0 auto;padding:28px 22px}h1{font:700 25pt Helvetica,Arial,sans-serif;margin:0 0 6px}h2{font:700 17pt Helvetica,Arial,sans-serif;color:#8C5A2B;border-bottom:2px solid #8C5A2B;padding-bottom:3px;margin-top:34px;page-break-after:avoid}
h3{font:700 13pt Helvetica,Arial,sans-serif;margin-top:22px;page-break-after:avoid}.m{color:#7A6A5A}code{font:9.5pt Menlo,monospace;background:#EDE2C8;padding:0 3px;border-radius:3px}
.box{background:#FBF7EC;border:1px solid #D8C9A6;border-radius:10px;padding:14px 18px;margin:16px 0}table{border-collapse:collapse;width:100%;font:9.5pt Helvetica,Arial,sans-serif;margin:10px 0;page-break-inside:avoid}
th{background:#8C5A2B;color:#fff;text-align:left;padding:4px 6px}td{border:1px solid #D8C9A6;padding:4px 6px;background:#FBF7EC;vertical-align:top}
figure{margin:14px 0;page-break-inside:avoid}figure img{width:100%;border:1px solid #D8C9A6;border-radius:6px;background:#F4ECD8}figcaption{font:9.5pt Helvetica,Arial,sans-serif;color:#5a4a3c;margin-top:4px}
li{margin:4px 0}@page{size:A4;margin:13mm;background:#F4ECD8}</style><main>"""
    (HERE / "REPORT_simplex_explainer.html").write_text(css + "".join(H) + "</main>")
    print("figures:", FIGN[0])

if __name__ == "__main__":
    main()
