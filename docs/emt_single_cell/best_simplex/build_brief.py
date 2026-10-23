#!/usr/bin/env python3
"""Brief for the supervisor: one-page at-a-glance table + one picture per dataset (beige). Run from Project_1/:
.venv/bin/python docs/emt_single_cell/best_simplex/build_brief.py   -> figures/brief/*.png and Brief_simplex_results.pdf"""
import sys, io, contextlib, textwrap
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, FancyBboxPatch
from matplotlib.backends.backend_pdf import PdfPages
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
with contextlib.redirect_stdout(io.StringIO()): import build_explainer_ppt as P
import analyze as A, figurebook as FB
from families import FAM, FCOL, FSHORT, family
OUT = HERE / "figures" / "brief"; OUT.mkdir(parents=True, exist_ok=True)
BG, PANEL, INK, MUT, BROWN = "#F4ECD8", "#FBF7EC", "#3A2E25", "#7A6A5A", "#8C5A2B"
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.facecolor": BG, "savefig.facecolor": BG, "text.color": INK})
MAIN, MS = P.MAIN, P.MS
KEY = ["dataset", "level", "unit", "geneset", "arm", "k"]
ANN = P.CA.drop(columns=["p"]).merge(MAIN[KEY + ["p"]], on=KEY); ANN["fam"] = ANN.task.map(family)
DSS = list(P.E.FG.ORDER); NAME = P.NAME
SHORTNAME = {"SCC": "Mouse skin SCC", "COOK": "Cook TGF-β lines", "PDAC": "Mouse PDAC", "T4": "4T1 clones", "GSE": "MCF10A / HuMEC"}
UNITW = {"SCC": "tumours", "COOK": "cell lines", "PDAC": "libraries", "T4": "clones", "GSE": "samples"}
SING = {"SCC": "tumour", "COOK": "cell line", "PDAC": "library", "T4": "clone", "GSE": "sample"}
SAMPLE = {"SCC": "sample", "COOK": "cell_line", "PDAC": "sample", "T4": "sample", "GSE": "sample"}
CORE = ["Proliferation", "Epithelial differentiation", "Mesenchymal / matrix (EMT)", "Metabolism & hypoxia", "Stress & immune signalling"]
GSL = {"all": "All genes", "ks": "KS EMT", "hallmark": "Hallmark EMT"}
obs = A.load_obs()

def pooled_status(ds, gs):
    r = P.pick_main(ds, P.POOLED[ds], gs)
    if r.p > .05: return "none", r
    an = ANN[(ANN.dataset == ds) & (ANN.level == "pooled") & (ANN.geneset == gs)]
    idn = int((an.top_sample_pct >= 70).sum())
    return ("identity" if idn >= max(1, len(an) - 1) else ("mixed" if idn else "programmes")), r

def within_frac(ds, gs):
    d = MAIN[(MAIN.dataset == ds) & (MAIN.geneset == gs) & (MAIN.level == "within")]
    return int((d.p <= .05).sum()), len(d)

def delta(ds):
    x = MS[MS.dataset == ds]; return float((x.hyb_gen_gap - x.oth_gen_gap).median()), float(x.hyb_gen_gap.median()), float(x.oth_gen_gap.median())

BOTTOM = {"SCC": "Proliferation vs matrix inside tumours; hybrids most interior",
          "COOK": "Pooled shape = cell lines; within lines a simplex mostly only with Hallmark genes",
          "PDAC": "Robust within-library triangles; EMT genes alone weaker",
          "T4": "Almost no structure inside clones",
          "GSE": "Clearest trade-off: differentiation / EMT-matrix / proliferation"}

# ---------------------------------------------------------------- 1. at-a-glance table
def glance():
    fig = plt.figure(figsize=(17, 8.4)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 17); ax.set_ylim(0, 8.4); ax.axis("off")
    ax.text(0.35, 8.0, "EMT archetype analysis — findings at a glance", fontsize=20, fontweight="bold", color=INK, va="center")
    ax.text(0.35, 7.6, "Within single tumours / cell lines, cells trade off a small set of shared programmes; pooled shapes mostly reflect sample identity; hybrids are interior mixtures, not a corner.", fontsize=11.5, color=BROWN, va="center")
    # column layout
    X = dict(ds=0.35, pool=2.75, wit=5.65, fam=8.55, q2=12.05, bl=13.75)
    hy = 6.95
    for x0, x1, t in ((2.6, 8.35, "QUESTION 1 · Is there a simplex, and are its corners tasks?"), (8.45, 13.6, ""), (11.95, 13.6, "QUESTION 2 · Hybrids"), ):
        pass
    ax.add_patch(FancyBboxPatch((2.6, 7.08), 9.25, 0.32, boxstyle="round,pad=0.02", fc=BROWN, ec="none")); ax.text(7.22, 7.24, "QUESTION 1  ·  Is there a simplex, and are its corners tasks?", color="white", ha="center", va="center", fontsize=11, fontweight="bold")
    ax.add_patch(FancyBboxPatch((11.95, 7.08), 1.65, 0.32, boxstyle="round,pad=0.02", fc=BROWN, ec="none")); ax.text(12.78, 7.24, "QUESTION 2", color="white", ha="center", va="center", fontsize=11, fontweight="bold")
    heads = [(X["pool"] + 1.35, "Pooled over samples\nA all  ·  B KS  ·  C Hallmark"), (X["wit"] + 1.35, "Within each unit\n(share of units with a simplex)"), (X["fam"] + 1.65, "What the corners are\n(within units, all gene sets)"),
             (X["q2"] + 0.75, "Hybrid cells vs\nother cells as generalists"), (X["bl"] + 1.0, "Bottom line")]
    for x, t in heads: ax.text(x, 6.62, t, ha="center", va="center", fontsize=9.5, color=INK, fontweight="bold")
    for j, f in enumerate(CORE): ax.text(X["fam"] + 0.35 + j * 0.66, 6.18, FSHORT[f], ha="center", va="center", fontsize=8, color=FCOL[f], fontweight="bold", rotation=0)
    for k, g in enumerate(("all", "ks", "hallmark")):
        ax.text(X["pool"] + 0.45 + k * 0.9, 6.18, "ABC"[k], ha="center", fontsize=9, color=MUT, fontweight="bold"); ax.text(X["wit"] + 0.45 + k * 0.9, 6.18, "ABC"[k], ha="center", fontsize=9, color=MUT, fontweight="bold")
    rh = 1.08; y0 = 5.55
    W_ = ANN[(ANN.level == "within") & (ANN.p <= .05)]
    for i, ds in enumerate(DSS):
        y = y0 - i * rh
        ax.add_patch(FancyBboxPatch((0.25, y - rh / 2 + 0.06), 16.5, rh - 0.12, boxstyle="round,pad=0.02", fc=PANEL if i % 2 == 0 else "#F1E8D2", ec="#D8C9A6", lw=.6))
        ax.text(X["ds"], y + 0.12, SHORTNAME[ds], fontsize=11.5, fontweight="bold", va="center"); ax.text(X["ds"], y - 0.2, UNITW[ds], fontsize=8.5, color=MUT, va="center")
        for k, g in enumerate(("all", "ks", "hallmark")):
            st, r = pooled_status(ds, g); cx = X["pool"] + 0.45 + k * 0.9
            col, lab = {"identity": ("#8a8172", "ID"), "mixed": ("#d69e2e", "mix"), "programmes": ("#3f7a4f", "prog"), "none": (None, "–")}[st]
            if col: ax.add_patch(FancyBboxPatch((cx - 0.36, y - 0.17), 0.72, 0.34, boxstyle="round,pad=0.03", fc=col, ec="none")); ax.text(cx, y, lab, ha="center", va="center", color="white", fontsize=9, fontweight="bold")
            else: ax.text(cx, y, lab, ha="center", va="center", color=MUT, fontsize=13)
            s_, n_ = within_frac(ds, g); cx = X["wit"] + 0.45 + k * 0.9; R = 0.25
            ax.add_patch(Circle((cx, y + 0.06), R, fc="white", ec=INK, lw=1))
            if s_: ax.add_patch(Wedge((cx, y + 0.06), R, 90 - 360 * s_ / n_, 90, fc=INK, ec=INK))
            ax.text(cx, y - 0.33, f"{s_}/{n_}", ha="center", va="center", fontsize=8, color=MUT)
        wd = W_[W_.dataset == ds]; fits = wd.groupby(KEY).fam.apply(set); nf = max(len(fits), 1)
        for j, f in enumerate(CORE):
            share = sum(f in s for s in fits) / nf if len(fits) else 0; cx = X["fam"] + 0.35 + j * 0.66
            if share > 0: ax.add_patch(Circle((cx, y), 0.06 + 0.24 * share, fc=FCOL[f], ec="none", alpha=.9 if len(fits) >= 5 else .35))
            else: ax.text(cx, y, "·", ha="center", va="center", color=MUT, fontsize=12)
        d_, h_, o_ = delta(ds); cx = X["q2"] + 0.75
        arrow, txt, col = (("▲", "leans interior", "#3f7a4f") if d_ >= 8 else (("▴", "slight lean", "#6b8f5e") if d_ >= 3 else ("●", "no difference", "#8a8172")))
        ax.text(cx, y + 0.13, f"{arrow} {txt}", ha="center", va="center", fontsize=10, color=col, fontweight="bold"); ax.text(cx, y - 0.2, f"{h_:.0f}% vs {o_:.0f}%  ({int(round(d_)) + 0:+d} pts)", ha="center", va="center", fontsize=8.3, color=MUT)
        ax.text(X["bl"], y, textwrap.fill(BOTTOM[ds], 30), fontsize=9.3, va="center", color=INK)
    # legend
    ly = 0.32
    items = [("#3f7a4f", "prog", "corners are programmes"), ("#d69e2e", "mix", "some corners = one sample"), ("#8a8172", "ID", "corners = samples (identity)")]
    x = 0.35
    ax.text(x, ly + 0.28, "Pooled:", fontsize=8.5, fontweight="bold", color=INK)
    for c, l, t in items:
        ax.add_patch(FancyBboxPatch((x, ly - 0.1), 0.5, 0.24, boxstyle="round,pad=0.02", fc=c, ec="none")); ax.text(x + 0.25, ly + 0.02, l, ha="center", va="center", color="white", fontsize=7.5, fontweight="bold"); ax.text(x + 0.6, ly + 0.02, t, fontsize=8, va="center", color=INK); x += 2.6
    ax.text(x, ly + 0.02, "– no significant simplex", fontsize=8, va="center", color=INK); x += 2.1
    ax.add_patch(Circle((x + 0.15, ly + 0.02), 0.13, fc="white", ec=INK)); ax.add_patch(Wedge((x + 0.15, ly + 0.02), 0.13, 0, 90, fc=INK)); ax.text(x + 0.35, ly + 0.02, "share of units with a significant simplex", fontsize=8, va="center"); x += 3.0
    ax.text(x, ly + 0.02, "Bubble = share of significant within-unit simplices\nwith a corner of that family (faded if < 5 simplices).\nQ2: gap rule, median over units with a simplex.", fontsize=8, va="center", color=INK)
    fn = OUT / "at_a_glance.png"; fig.savefig(fn, dpi=180); plt.close(fig); return fn

# ---------------------------------------------------------------- 2. dataset pictures
HEAD = {
 "SCC": "Inside SCC tumours, cells trade off proliferation against a mesenchymal / matrix programme; hybrid cells sit between the corners",
 "COOK": "In the Cook lines the pooled simplex is just the four cell lines; inside a line a simplex appears mainly with the Hallmark EMT genes (3 of 4 lines)",
 "PDAC": "Every PDAC library forms a simplex on all genes; its corners are proliferation and mesenchymal / matrix plus a metabolic or epithelial programme",
 "T4": "4T1 clones show almost no task structure: the pooled shape is mostly clone identity and single clones rarely form a simplex",
 "GSE": "MCF10A / HuMEC show the clearest trade-off: epithelial differentiation vs mesenchymal / matrix vs proliferation in every sample"}

def corner_labels(ax, V, an, level):
    ax.margins(0.2); ax.autoscale_view()
    c0 = V[:, :2].mean(0)
    for j in range(V.shape[0]):
        a = an[an.corner == f"V{j + 1}"]
        if not len(a): continue
        a = a.iloc[0]; idn = level == "pooled" and a.top_sample_pct >= 70
        lab = f"{a.top_sample}\n(identity)" if idn else FSHORT.get(family(a.task), "?") + f"\n{a.task.split(' (')[0]}"
        col = FCOL["Identity (one sample)"] if idn else FCOL.get(family(a.task), "#a0aec0")
        v = V[j, :2]; d = v - c0; d = d / (np.linalg.norm(d) + 1e-9)
        span = np.ptp(V[:, :2], axis=0).max()
        ax.annotate(lab, v, v + d * 0.12 * span, ha="center", va="center", fontsize=8, color="white", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.25", fc=col, ec="none"), arrowprops=dict(arrowstyle="-", color=col, lw=1))
        ax.scatter(*v, s=60, color=col, zorder=6, edgecolor="white")

def card(ds):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.03, 0.95, SHORTNAME[ds], fontsize=13, color=BROWN, fontweight="bold")
    fig.text(0.03, 0.915, textwrap.fill(HEAD[ds], 105), fontsize=15.5, color=INK, fontweight="bold", va="top")
    axP = fig.add_axes([0.03, 0.17, 0.30, 0.58]); axW = fig.add_axes([0.36, 0.17, 0.34, 0.58]); axQ = fig.add_axes([0.76, 0.30, 0.21, 0.45])
    # pooled (all genes), coloured by sample
    rp = P.pick_main(ds, P.POOLED[ds], "all"); pts, V, o = A.cells_for(rp, obs)
    samp = o[SAMPLE[ds]].astype(str).values; us = list(pd.Series(samp).value_counts().index)
    pal = plt.get_cmap("tab20")(np.linspace(0, 1, max(len(us), 2)))
    idx = np.random.default_rng(0).permutation(len(o))[:5000]
    for k_, u in enumerate(us):
        m = samp[idx] == u; axP.scatter(pts[idx][m, 0], pts[idx][m, 1], s=3, color=pal[k_], alpha=.75, lw=0, label=u)
    FB.outline(axP, V[:, :2], INK)
    an = ANN[(ANN.dataset == ds) & (ANN.unit == rp.unit) & (ANN.geneset == "all") & (ANN.arm == rp.arm) & (ANN.k == rp.k)]
    if rp.p <= .05: corner_labels(axP, V, an, "pooled")
    st, _ = pooled_status(ds, "all")
    axP.set_title(f"Pooled, all genes: {'no significant simplex' if st == 'none' else {'identity': 'corners mostly = ' + UNITW[ds], 'mixed': 'partly ' + UNITW[ds] + ', partly programmes', 'programmes': 'corners = programmes'}[st]}\n(cells coloured by {SING[ds]}; k={int(rp.k)}, p={rp.p:.2g})", fontsize=10, color=INK)
    axP.legend(loc="upper center", bbox_to_anchor=(.5, -.01), ncol=min(4, len(us)), fontsize=7, frameon=False, markerscale=3)
    # within: most significant all-genes unit (else best across gene sets)
    d = MAIN[(MAIN.dataset == ds) & (MAIN.level == "within")].copy(); d["g"] = d.geneset.map({"all": 0, "hallmark": 1, "ks": 2})
    rw = d.sort_values(["p", "g", "n_cells"], ascending=[True, True, False]).iloc[0]
    pts, V, o = A.cells_for(rw, obs)
    order = FB.STATES.get(ds, ["Epithelial", "Early-Hybrid", "Late-Hybrid", "Mesenchymal"]); stt = o["cluster_annotation"].astype(str).values
    sidx = np.array([order.index(x) if x in order else 0 for x in stt]); idx = np.random.default_rng(0).permutation(len(o))[:4500]
    for j in range(4):
        m = sidx[idx] == j; axW.scatter(pts[idx][m, 0], pts[idx][m, 1], s=3.5, color=FB.SCOL[j], alpha=.75, lw=0, label=order[j])
    FB.outline(axW, V[:, :2], INK)
    an = ANN[(ANN.dataset == ds) & (ANN.unit == rw.unit) & (ANN.geneset == rw.geneset) & (ANN.arm == rw.arm) & (ANN.k == rw.k)]
    corner_labels(axW, V, an, "within")
    s_, n_ = within_frac(ds, rw.geneset)
    axW.set_title(f"Inside one {SING[ds]} ({rw.unit}, {GSL[rw.geneset]}): corners = programmes\n(cells coloured by EMT state; k={int(rw.k)}, p={rw.p:.2g}; {s_} of {n_} {UNITW[ds]} have a simplex)", fontsize=10, color=INK)
    axW.legend(loc="upper center", bbox_to_anchor=(.5, -.01), ncol=4, fontsize=8, frameon=False, markerscale=3)
    for ax in (axP, axW):
        ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor(PANEL)
        for s in ax.spines.values(): s.set_color("#D8C9A6")
    # Q2 mini bars
    x = MS[MS.dataset == ds]
    for i, g in enumerate(("all", "ks", "hallmark")):
        xx = x[x.geneset == g]
        if len(xx):
            axQ.bar(i - .19, xx.hyb_gen_gap.median(), .36, color="#dd6b20"); axQ.bar(i + .19, xx.oth_gen_gap.median(), .36, color="#a99f8c")
            axQ.text(i - .19, xx.hyb_gen_gap.median() + 2, f"{xx.hyb_gen_gap.median():.0f}", ha="center", fontsize=8.5); axQ.text(i + .19, xx.oth_gen_gap.median() + 2, f"{xx.oth_gen_gap.median():.0f}", ha="center", fontsize=8.5)
        else: axQ.text(i, 5, "no\nsimplex", ha="center", fontsize=8, color=MUT)
    axQ.set_xticks(range(3)); axQ.set_xticklabels([GSL[g] for g in ("all", "ks", "hallmark")], fontsize=8.5); axQ.set_ylim(0, 100); axQ.set_facecolor(PANEL)
    for s in axQ.spines.values(): s.set_color("#D8C9A6")
    d_, h_, o_ = delta(ds)
    axQ.set_title("Hybrid cells (orange) vs other cells (grey):\n% that are generalists (gap rule)", fontsize=10, color=INK)
    fig.text(0.865, 0.2, ("Hybrids lean toward the interior" if d_ >= 8 else ("Hybrids lean slightly toward the interior" if d_ >= 3 else "No difference between hybrids and others")) + f"\n(median {int(round(d_)) + 0:+d} points)", ha="center", fontsize=10.5, color=INK, fontweight="bold")
    # takeaway strip
    fams = ANN[(ANN.dataset == ds) & (ANN.level == "within") & (ANN.p <= .05)].groupby(KEY).fam.apply(set)
    fam_txt = ", ".join(f"{f.split(' (')[0].lower()} ({sum(f in s for s in fams)}/{len(fams)})" for f in CORE if sum(f in s for s in fams)) if len(fams) else "none"
    w = [within_frac(ds, g) for g in ("all", "ks", "hallmark")]
    fig.text(0.03, 0.035, textwrap.fill(f"Significant within-{SING[ds]} simplices: all genes {w[0][0]}/{w[0][1]}, KS {w[1][0]}/{w[1][1]}, Hallmark {w[2][0]}/{w[2][1]}.   Corner families across them: {fam_txt}.", 190), fontsize=10, color=MUT)
    fn = OUT / f"card_{ds}.png"; fig.savefig(fn, dpi=150); plt.close(fig); return fn

if __name__ == "__main__":
    figs = [glance()] + [card(ds) for ds in DSS]
    with PdfPages(HERE / "Brief_simplex_results.pdf") as pdf:
        for f in figs:
            im = plt.imread(f); h, w = im.shape[:2]; fig = plt.figure(figsize=(w / 150, h / 150)); fig.figimage(im, resize=True); pdf.savefig(fig); plt.close(fig)
    print("\n".join(str(f) for f in figs))
