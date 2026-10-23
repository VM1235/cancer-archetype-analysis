#!/usr/bin/env python3
"""Story deck (beige): five messages, one purpose-built figure each, then intrinsic-dimension tests and next steps.
Run from Project_1/:  .venv/bin/python docs/emt_single_cell/best_simplex/build_story_deck.py"""
import sys, io, contextlib
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
with contextlib.redirect_stdout(io.StringIO()):
    import build_explainer_ppt as P
    import build_brief as BR
import analyze as A, figurebook as FB
from families import family, FCOL, FSHORT
from pptx import Presentation
from pptx.util import Inches
OUT = HERE / "figures" / "story"; OUT.mkdir(parents=True, exist_ok=True)
BG, PANEL, INK, MUT, BROWN = "#F4ECD8", "#FBF7EC", "#3A2E25", "#7A6A5A", "#8C5A2B"
plt.rcParams.update({"figure.facecolor": BG, "savefig.facecolor": BG, "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "axes.titlecolor": INK, "axes.edgecolor": "#D8C9A6"})
MAIN, MS, ANN, KEY = P.MAIN, P.MS, BR.ANN, BR.KEY
E_, M_, Pr = "Epithelial differentiation", "Mesenchymal / matrix (EMT)", "Proliferation"
HYB = {"Early-Hybrid", "Hybrid", "Late-Hybrid"}
obs = A.load_obs()
W = ANN[(ANN.level == "within") & (ANN.p <= .05)]; fits = W.groupby(KEY).fam.apply(set)

def frame(ax):
    ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor(PANEL)
    for s in ax.spines.values(): s.set_color("#D8C9A6")

def simplex(ax, row, colour="state", level="within", title=None, highlight_edge=None):
    pts, V, o = A.cells_for(row, obs); ds = row.dataset
    idx = np.random.default_rng(0).permutation(len(o))[:4500]
    if colour == "state":
        order = FB.STATES.get(ds, ["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"]); st = o["cluster_annotation"].astype(str).values
        sidx = np.array([order.index(x) if x in order else 0 for x in st])
        for j in range(4):
            m = sidx[idx] == j; ax.scatter(pts[idx][m, 0], pts[idx][m, 1], s=3.5, color=FB.SCOL[j], alpha=.75, lw=0, label=order[j])
    elif colour == "hybrid":
        h = o["cluster_annotation"].astype(str).isin(HYB).values
        ax.scatter(pts[idx][~h[idx], 0], pts[idx][~h[idx], 1], s=3, color="#c9bea6", alpha=.7, lw=0, label="other cells")
        ax.scatter(pts[idx][h[idx], 0], pts[idx][h[idx], 1], s=4, color="#dd6b20", alpha=.8, lw=0, label="hybrid cells")
    else:
        samp = o[BR.SAMPLE[ds]].astype(str).values; us = list(pd.Series(samp).value_counts().index); pal = plt.get_cmap("tab10")(np.arange(len(us)) % 10)
        for k, u in enumerate(us):
            m = samp[idx] == u; ax.scatter(pts[idx][m, 0], pts[idx][m, 1], s=3.5, color=pal[k], alpha=.75, lw=0, label=u)
    FB.outline(ax, V[:, :2], INK)
    an = ANN[(ANN.dataset == ds) & (ANN.unit == row.unit) & (ANN.geneset == row.geneset) & (ANN.arm == row.arm) & (ANN.k == row.k)]
    if highlight_edge is not None:
        fams = {c.corner: c.fam for _, c in an.iterrows()}
        ii = [int(c[1:]) - 1 for c, f in fams.items() if f in highlight_edge]
        if len(ii) == 2: ax.plot(V[ii, 0], V[ii, 1], color=FCOL[M_], lw=4, alpha=.55, zorder=4, solid_capstyle="round")
    BR.corner_labels(ax, V, an, level)
    frame(ax)
    if title: ax.set_title(title, fontsize=10.5, color=INK)
    return an

def best_within(ds, need=2, k3=True):
    d = MAIN[(MAIN.dataset == ds) & (MAIN.level == "within") & (MAIN.p <= .05)].copy()
    sc = []
    for _, r in d.iterrows():
        f = fits.get(tuple(r[k] for k in KEY), set()); sc.append(len({E_, M_, Pr} & f))
    d["core"] = sc; d["g"] = d.geneset.map({"all": 0, "hallmark": 1, "ks": 2}); d["k3"] = (d.k == 3).astype(int)
    return d.sort_values(["core", "k3", "p", "g"], ascending=[False, False, True, True]).iloc[0]

# ---------------- figures
def fig1():
    fig = plt.figure(figsize=(13, 5.6)); a1 = fig.add_axes([.02, .1, .48, .78]); a2 = fig.add_axes([.6, .18, .37, .66])
    r = P.pick_main("COOK", P.POOLED["COOK"], "all")
    simplex(a1, r, colour="sample", level="pooled", title=f"Cook: four TGF-β cell lines pooled (all genes, k={int(r.k)})")
    a1.legend(loc="upper center", bbox_to_anchor=(.5, -.01), ncol=4, fontsize=8.5, frameon=False, markerscale=3)
    vals = []
    for ds in BR.DSS:
        x = ANN[(ANN.dataset == ds) & (ANN.level == "pooled") & (ANN.p <= .05)]
        vals.append(100 * (x.top_sample_pct >= 70).mean() if len(x) else np.nan)
    y = np.arange(5)[::-1]
    a2.barh(y, np.nan_to_num(vals), color="#8a8172")
    for yi, v in zip(y, vals): a2.text((0 if np.isnan(v) else v) + 2, yi, "no pooled simplex" if np.isnan(v) else f"{v:.0f}%", va="center", fontsize=9.5)
    a2.set_yticks(y); a2.set_yticklabels([BR.SHORTNAME[d] for d in BR.DSS]); a2.set_xlim(0, 115); a2.set_xlabel("% of pooled corners that are ≥ 70% one sample")
    a2.set_title("Across datasets: how much of the pooled shape is sample identity", fontsize=10.5)
    for s in a2.spines.values(): s.set_color("#D8C9A6")
    a2.set_facecolor(PANEL); f = OUT / "m1_identity.png"; fig.savefig(f, dpi=160); plt.close(fig); return f

def fig2():
    picks = [best_within(ds) for ds in ("GSE", "SCC", "PDAC")]
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 5.3), gridspec_kw=dict(left=.05, right=.97, top=.86, bottom=.12, wspace=.06))
    for ax, r in zip(axs, picks):
        simplex(ax, r, colour="state", title=f"{BR.SHORTNAME[r.dataset]}: {r.unit} ({BR.GSL[r.geneset]}, k={int(r.k)})")
    h = [plt.Line2D([], [], marker="o", ls="", color=FB.SCOL[j], markersize=6, label=n) for j, n in enumerate(["epithelial", "early hybrid", "hybrid / late hybrid", "mesenchymal / late EMT"])]
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, fontsize=9.5)
    f = OUT / "m2_programmes.png"; fig.savefig(f, dpi=160); plt.close(fig); return f, picks

def fig3():
    r = best_within("GSE")
    fig = plt.figure(figsize=(13, 5.6)); a1 = fig.add_axes([.02, .1, .44, .8]); a2 = fig.add_axes([.66, .2, .31, .62])
    simplex(a1, r, colour="state", title=f"{BR.SHORTNAME[r.dataset]} {r.unit}: the epithelial → mesenchymal edge (highlighted)", highlight_edge={E_, M_})
    labs, vals, cols = [], [], []
    for fam, states, lab in ((E_, {"Epithelial"}, "Epithelial-differentiation corners:\nmostly epithelial cells"), (M_, {"Mesenchymal", "Late-EMT"}, "Mesenchymal / matrix corners:\nmost often mesenchymal / late-EMT cells"), (Pr, HYB, "Proliferation corners:\nmostly hybrid-state cells")):
        x = W[W.fam == fam]; labs.append(lab); vals.append(100 * x.top_state.isin(states).mean()); cols.append(FCOL[fam])
    y = np.arange(3)[::-1]; a2.barh(y, vals, color=cols)
    for yi, v in zip(y, vals): a2.text(v + 2, yi, f"{v:.0f}%", va="center", fontsize=10, fontweight="bold")
    a2.set_yticks(y); a2.set_yticklabels(labs, fontsize=9.5); a2.set_xlim(0, 110); a2.set_xlabel("% of corners (significant within-unit fits) dominated by that state")
    a2.set_title("Which cells sit at each kind of corner (all datasets)", fontsize=10.5); a2.set_facecolor(PANEL)
    for s in a2.spines.values(): s.set_color("#D8C9A6")
    f = OUT / "m3_edge.png"; fig.savefig(f, dpi=160); plt.close(fig); return f, r

def fig4():
    r = best_within("SCC")
    fig = plt.figure(figsize=(13, 5.6)); a1 = fig.add_axes([.02, .1, .44, .8]); a2 = fig.add_axes([.53, .14, .45, .7])
    simplex(a1, r, colour="hybrid", title=f"{BR.SHORTNAME[r.dataset]} {r.unit}: hybrid cells (orange) in the interior")
    a1.legend(loc="upper center", bbox_to_anchor=(.5, -.01), ncol=2, fontsize=9, frameon=False, markerscale=3)
    for i, d in enumerate(BR.DSS):
        x = MS[MS.dataset == d]
        for off, col, c in ((-.19, "#dd6b20", "hyb_gen_gap"), (.19, "#a99f8c", "oth_gen_gap")):
            m = x[c].median(); a2.bar(i + off, m, .36, color=col); a2.text(i + off, m + 1.5, f"{m:.0f}%", ha="center", fontsize=9, fontweight="bold")
    a2.set_xticks(range(5)); a2.set_xticklabels([BR.SHORTNAME[d].replace(" ", "\n", 1) for d in BR.DSS], fontsize=9); a2.set_ylim(0, 100); a2.set_ylabel("% that are generalists (gap rule)")
    a2.set_title("Hybrid cells (orange) vs other cells (grey): share that are generalists", fontsize=10.5); a2.set_facecolor(PANEL)
    for s in a2.spines.values(): s.set_color("#D8C9A6")
    f = OUT / "m4_hybrids.png"; fig.savefig(f, dpi=160); plt.close(fig); return f, r

def fig5():
    fig = plt.figure(figsize=(13, 5.6)); a1 = fig.add_axes([.12, .16, .36, .68]); a2 = fig.add_axes([.55, .06, .36, .82])
    o = pd.read_csv(HERE.parents[2] / "Mouse_SCC/data/processed/cell_metadata.csv")
    ct = pd.crosstab(o.cluster_annotation, o.assay).reindex(["Epithelial", "Early-Hybrid", "Hybrid", "Late-EMT"]); pc = 100 * ct.div(ct.sum(1), axis=0)
    y = np.arange(4)[::-1]
    a1.barh(y, pc["multiome"], color="#8C5A2B", label="multiome (nuclei)"); a1.barh(y, pc["scRNA"], left=pc["multiome"], color="#C9B994", label="scRNA (whole cells)")
    for yi, st in zip(y, pc.index):
        a1.text(pc.loc[st, "multiome"] / 2, yi, f"{pc.loc[st, 'multiome']:.0f}%", color="white", ha="center", va="center", fontweight="bold")
        a1.text(103, yi, "ID ≈ 29" if st in ("Epithelial", "Early-Hybrid") else "ID ≈ 22–23", va="center", fontsize=10)
    a1.set_yticks(y); a1.set_yticklabels(pc.index); a1.set_xlim(0, 125); a1.set_xlabel("% of each EMT state's cells (mouse SCC)"); a1.set_facecolor(PANEL)
    a1.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(.45, -.14), ncol=2); a1.set_title("Test 1: high-ID states come mostly from multiome tumours", fontsize=10.5)
    for s in a1.spines.values(): s.set_color("#D8C9A6")
    # conceptual prediction: ID rises with mixing (1 - purity) in a toy triangle
    rng = np.random.default_rng(2); V = np.array([[0, 0], [1, 0], [.5, .87]]); wts = rng.dirichlet([.8, .8, .8], 2500); pts = wts @ V
    s_ = np.sort(wts, 1); pur = s_[:, -1] - s_[:, -2]
    sc = a2.scatter(pts[:, 0], pts[:, 1], c=1 - pur, cmap="YlOrBr", s=7, lw=0)
    FB.outline(a2, V, INK); a2.set_aspect("equal"); frame(a2)
    for v, t, c in ((V[0], "epithelial\ndifferentiation", FCOL[E_]), (V[1], "mesenchymal /\nmatrix", FCOL[M_]), (V[2], "proliferation", FCOL[Pr])):
        a2.annotate(t, v, v + (v - V.mean(0)) * .28, ha="center", va="center", fontsize=9, color="white", fontweight="bold", bbox=dict(boxstyle="round,pad=.3", fc=c, ec="none"))
    a2.set_xlim(-.3, 1.3); a2.set_ylim(-.25, 1.15)
    cb = fig.colorbar(sc, ax=a2, fraction=.035, pad=.01); cb.set_label("predicted local ID (mixing of programmes)", fontsize=9, color=INK)
    a2.set_title("Tests 2–3 (hypothesis): ID is low at corners, high where programmes mix", fontsize=10.5)
    f = OUT / "m5_id.png"; fig.savefig(f, dpi=160); plt.close(fig); return f

# ---------------- deck
P.prs = Presentation(); P.prs.slide_width, P.prs.slide_height = Inches(13.333), Inches(7.5); P.HT.clear()
para, box, title, picture, new_slide, table = P.para, P.box, P.title, P.picture, P.new_slide, P.table

def msg_slide(num, head, sub, fig, bullets, so_what):
    new_slide(); title(f"{num}. {head}", sub)
    picture(fig, 0.35, 1.35, 12.6, 4.55)
    box(0.35, 6.0, 8.6, 1.35, [para(t, 12, False, P.INK, after=4, bullet=True) for t in bullets], fill=P.PANEL, border="D8C9A6")
    box(9.1, 6.0, 3.85, 1.35, [para("So what", 11.5, True, P.BROWN, after=3), para(so_what, 11.5, False, P.INK, after=0)], fill="F1E3C6", border="D8C9A6")

new_slide()
box(0.8, 1.4, 11.7, 2.0, [para("EMT is a trade-off between programmes, not between states", 34, True, P.INK, after=8), para("Archetype (Pareto-simplex) analysis of five single-cell EMT datasets, and what it means for intrinsic dimension", 16, False, P.BROWN, after=0)], fit=False)
box(0.8, 3.9, 11.7, 2.8, [para("The story in five messages", 14, True, P.BROWN, after=6)] + [para(f"{i}. {t}", 14, False, P.INK, after=4) for i, t in enumerate([
    "Pooled simplices separate samples, not EMT states",
    "Within a tumour or cell line, cells trade off three programmes",
    "EMT is the edge between the epithelial and the mesenchymal / matrix corner",
    "Hybrid cells are interior mixtures and the proliferating cells",
    "This predicts how intrinsic dimension should vary, and gives a test"], 1)], fill=P.PANEL, border="D8C9A6")

ident = 100 * (ANN[(ANN.level == "pooled") & (ANN.p <= .05)].top_sample_pct >= 70).mean()
msg_slide(1, "Pooled simplices separate samples, not EMT states", "A pooled fit finds whatever differs most between cells; across samples, that is the sample", fig1(),
          [f"In the Cook lines every pooled corner is one cell line, in all three gene sets.", f"Across datasets, {ident:.0f}% of the corners of significant pooled simplices are ≥ 70% one sample."],
          "Questions about EMT tasks must be asked inside one tumour or line.")
f2, picks = fig2()
two = sum(len({E_, M_, Pr} & f) >= 2 for f in fits)
msg_slide(2, "Within a tumour or cell line, cells trade off three programmes", "Most significant within-unit simplex of three datasets; corners named by the Hallmark programme enriched at each", f2,
          [f"{two} of {len(fits)} significant within-unit simplices have at least two of: proliferation, epithelial differentiation, mesenchymal / matrix.", "A fourth corner depends on context: metabolism (PDAC), stress / immune signalling (SCC)."],
          "There is a shared, low-dimensional trade-off, and it is about programmes.")
f3, r3 = fig3()
epimes = sum({E_, M_} <= f for f in fits)
msg_slide(3, "EMT is the edge between epithelial and matrix corners", "Cells coloured by EMT state run from one corner to the other; hybrids sit along the edge", f3,
          [f"Epithelial and mesenchymal / matrix corners appear together in {epimes} of {len(fits)} simplices.", "Epithelial corners hold epithelial cells (87%); matrix corners most often hold mesenchymal / late-EMT cells (55%); no corner is defined by a hybrid state."],
          "EMT behaves as a continuum along one edge of a broader trade-off.")
f4, r4 = fig4()
d = MS.hyb_gen_gap - MS.oth_gen_gap
msg_slide(4, "Hybrid cells are interior mixtures and the proliferating cells", "Left: one tumour with hybrid cells highlighted. Right: hybrid vs other cells as generalists, all five datasets", f4,
          [f"{MS.hyb_gen_gap.median():.0f}% of hybrid vs {MS.oth_gen_gap.median():.0f}% of other cells are generalists; hybrids more often generalists in {(d > 0).sum()} of {len(MS)} fits (clearest in SCC, absent in Cook).",
           f"In {int(W[W.fam == Pr].top_state.isin(HYB).sum())} of {int((W.fam == Pr).sum())} proliferation corners, hybrid-state cells are the largest group."],
          "Hybrids are balanced mixtures that keep dividing, not a separate task.")
msg_slide(5, "This predicts how intrinsic dimension should vary, and gives a test", "EMT_ID: ID is higher in Epithelial / Early-Hybrid (≈ 29) than in Hybrid / Late-EMT (≈ 22–23) in mouse SCC", fig5(),
          ["Test 1 (assay): compare ID between states within scRNA tumours only (SC15, SC32, SC46 have all four states).", "Tests 2–3 (geometry, programmes): model ID against simplex position; predict lower ID near corners, higher in mixed generalists."],
          "First rule out the multiome / scRNA split; then ask whether position in the simplex explains ID.")

new_slide(); title("Next steps", "Ordered by how much they change the conclusions, for the least new work")
table(0.5, 1.5, 12.3, [["#", "Step", "Question it answers", "Needs"],
    ["1", "ID within assay and within tumour (SCC)", "Is the ID–state effect biology or the multiome / scRNA split?", "EMT_ID per-cell table"],
    ["2", "ID against simplex position (purity, corner distances)", "Does position or programme commitment explain ID?", "Saved weights + ID table"],
    ["3", "Hybrid vs generalist by all four states", "Is the generalist share highest in hybrids (inverted U)?", "Existing fits"],
    ["4", "Path vs triangle inside units (SCC, PDAC)", "Is EMT a one-dimensional path within the simplex?", "Existing script"],
    ["5", "Batch-corrected pooled fits", "Is there a shared simplex across samples once identity is removed?", "Harmony / scANVI"],
    ["6", "Time and space vs simplex position (Cook pseudotime, GSE colony position)", "Do cells move along the epithelial → matrix edge?", "Existing metadata"],
    ["7", "ID for the other four datasets", "Does the ID pattern generalise?", "Compute time"]], [0.4, 4.3, 5.4, 2.2], size=11)
box(0.5, 5.3, 12.3, 1.9, [para("Caveats", 12.5, True, P.BROWN, after=4)] + [para(t, 11, False, P.INK, after=3, bullet=True) for t in [
    "Corner names are automatic (top Hallmark programme) and grouped into families by judgement, checked against the cells at each corner.",
    "Hybrid labels in four datasets come from the KS score; few samples per condition; many fits, no multiplicity correction.",
    "EMT_ID medians quoted from earlier notes; ID not recomputed here."]], fill=P.PANEL, border="D8C9A6")

P.prs.save(HERE / "Story_EMT_programmes.pptx")
css = "<!doctype html><meta charset=utf-8><style>@page{size:13.333in 7.5in;margin:0}body{margin:0;font-family:Helvetica,Arial,sans-serif}.pg{position:relative;width:13.333in;height:7.5in;background:#%s;overflow:hidden;page-break-after:always;-webkit-print-color-adjust:exact;print-color-adjust:exact}</style>" % P.BEIGE
(HERE / "Story_EMT_programmes.html").write_text(css + "".join(P.HT) + "</section>")
print("saved", len(P.prs.slides), "slides;", [f"{r.dataset} {r.unit} {r.geneset} k{int(r.k)}" for r in picks], r3.unit, r4.unit)
