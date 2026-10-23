#!/usr/bin/env python3
"""Beige-theme PowerPoint (and matching PDF) version of the explainer: every result slide = simplex figure on the left, plain-language explanation on the right.
Run from Project_1/ after analyze.py, corners.py:  .venv/bin/python docs/emt_single_cell/best_simplex/build_explainer_ppt.py
Outputs: Explainer_simplex_results.pptx, Explainer_simplex_results.html/.pdf (via Chrome), figures/explainer_beige/.
"""
import sys, math, html
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from PIL import Image
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import analyze as A, figurebook as FB, build_explainer as E

# ---------------------------------------------------------------- theme
BEIGE, PANEL, INK, MUTED, BROWN = "F4ECD8", "FBF7EC", "3A2E25", "7A6A5A", "8C5A2B"
GREEN, RED, BLUE = "3F7A4F", "A33B2B", "2B6CB0"
FB.SCOL = ["#2b6cb0", "#d69e2e", "#dd6b20", "#c53030"]
FB.VCOL = ["#2b6cb0", "#d97706", "#7c3aed", "#2f855a", "#c53030", "#be185d", "#4a5568"]
E.BG, E.FGC, E.MUT, E.GEN_COL, E.OUTL_COL, E.SPINE = "#" + BEIGE, "#" + INK, "#" + MUTED, "#c9bea6", "#" + INK, "#cdbf9f"
E.OUT = HERE / "figures" / "explainer_beige"; E.OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"figure.facecolor": "#" + BEIGE, "axes.facecolor": "#" + PANEL, "savefig.facecolor": "#" + BEIGE, "text.color": "#" + INK, "axes.edgecolor": "#cdbf9f"})
rgb = lambda h: RGBColor.from_string(h)

prs = Presentation(); prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
HT = []; CUR = None

def new_slide():
    global CUR
    s = prs.slides.add_slide(prs.slide_layouts[6]); s.background.fill.solid(); s.background.fill.fore_color.rgb = rgb(BEIGE)
    HT.append(("</section>" if HT else "") + '<section class="pg">'); CUR = s; return s

def est_height(paras, w):
    h = 0
    for p in paras:
        size = p.get("size", 12); txt = "".join(r[0] for r in p["runs"])
        cpl = max(8, (w - 0.15) * 72 / (size * 0.50)); lines = max(1, math.ceil(len(txt) / cpl))
        h += lines * size * 1.22 / 72 + p.get("after", 4) / 72
    return h

def para(text, size=12, bold=False, color=INK, after=4, bullet=False, align=None):
    runs = text if isinstance(text, list) else [(text, bold, color)]
    if bullet: runs = [("• ", False, color)] + list(runs)
    return dict(runs=runs, size=size, after=after, align=align)

def box(x, y, w, h, paras, fill=None, border=None, fit=True, valign="top"):
    if fit:   # shrink the font until the estimate fits
        for k in range(10):
            if est_height(paras, w) <= h - 0.1 or all(p["size"] <= 8.5 for p in paras): break
            for p in paras: p["size"] = max(8.5, p["size"] - 0.5)
    if fill or border:
        sh = CUR.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
        sh.adjustments[0] = 0.04
        if fill: sh.fill.solid(); sh.fill.fore_color.rgb = rgb(fill)
        else: sh.fill.background()
        if border: sh.line.color.rgb = rgb(border); sh.line.width = Pt(0.75)
        else: sh.line.fill.background()
        sh.shadow.inherit = False
    tb = CUR.shapes.add_textbox(Inches(x + 0.08), Inches(y + 0.05), Inches(w - 0.16), Inches(h - 0.1)); tf = tb.text_frame; tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    style = f"left:{x}in;top:{y}in;width:{w}in;height:{h}in;" + (f"background:#{fill};" if fill else "") + (f"border:1px solid #{border};border-radius:8px;" if border else "") + "padding:5px 7px;box-sizing:border-box;overflow:hidden;position:absolute;"
    inner = ""
    for i, p in enumerate(paras):
        pp = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); pp.space_after = Pt(p["after"])
        if p.get("align") == "c": pp.alignment = PP_ALIGN.CENTER
        h_ = ""
        for t, b, c in p["runs"]:
            r = pp.add_run(); r.text = t; r.font.size = Pt(p["size"]); r.font.bold = b; r.font.color.rgb = rgb(c); r.font.name = "Helvetica"
            h_ += f'<span style="font-weight:{700 if b else 400};color:#{c}">{html.escape(t)}</span>'
        inner += f'<p style="font-size:{p["size"]}pt;margin:0 0 {p["after"]}pt 0;line-height:1.22;{"text-align:center;" if p.get("align") == "c" else ""}">{h_}</p>'
    HT.append(f'<div style="{style}">{inner}</div>')

def title(t, sub=None):
    box(0.45, 0.25, 12.4, 0.7, [para(t, 24, True, INK, after=0)], fit=False)
    if sub: box(0.45, 0.85, 12.4, 0.4, [para(sub, 12, False, MUTED, after=0)], fit=False)
    HT.append(f'<div style="position:absolute;left:0.45in;top:1.2in;width:12.4in;height:2px;background:#{BROWN}"></div>')
    ln = CUR.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0.45), Inches(1.2), Inches(12.4), Inches(0.02)); ln.fill.solid(); ln.fill.fore_color.rgb = rgb(BROWN); ln.line.fill.background()

def picture(path, x, y, w, h):
    im = Image.open(path); ar = im.width / im.height
    ww, hh = (w, w / ar) if w / ar <= h else (h * ar, h)
    CUR.shapes.add_picture(str(path), Inches(x), Inches(y), Inches(ww), Inches(hh))
    HT.append(f'<img src="{path.relative_to(HERE).as_posix()}" style="position:absolute;left:{x}in;top:{y}in;width:{ww}in;height:{hh}in">'); return ww, hh

def table(x, y, w, rows, colw, size=10.5, head=BROWN):
    n, m = len(rows), len(rows[0]); rh = 0.3
    shp = CUR.shapes.add_table(n, m, Inches(x), Inches(y), Inches(w), Inches(rh * n)); t = shp.table
    hh = f'<table style="position:absolute;left:{x}in;top:{y}in;width:{w}in;border-collapse:collapse;font-size:{size}pt;color:#{INK}">'
    for i, r in enumerate(rows):
        hh += "<tr>"
        for j, v in enumerate(r):
            c = t.cell(i, j); c.text = str(v); c.fill.solid(); c.fill.fore_color.rgb = rgb(head if i == 0 else (PANEL if i % 2 else "F1E8D2"))
            c.margin_left = c.margin_right = Inches(0.05); c.margin_top = c.margin_bottom = Inches(0.02)
            for p in c.text_frame.paragraphs:
                for rr in p.runs: rr.font.size = Pt(size); rr.font.bold = (i == 0); rr.font.color.rgb = rgb("FFFFFF" if i == 0 else INK); rr.font.name = "Helvetica"
            hh += f'<{"th" if i == 0 else "td"} style="border:1px solid #d8c9a6;padding:2px 5px;background:#{head if i == 0 else (PANEL if i % 2 else "F1E8D2")};{"color:#fff;" if i == 0 else ""}text-align:left">{html.escape(str(v))}</{"th" if i == 0 else "td"}>'
        hh += "</tr>"
    HT.append(hh + "</table>")
    for j, cw in enumerate(colw): t.columns[j].width = Inches(cw)
    for i in range(n): t.rows[i].height = Inches(rh)

# ---------------------------------------------------------------- data
B = E.B; CA = E.CA; NAME = E.NAME; GS = E.GS
for _r in ("gap", "groves5", "w95"):   # % of NON-hybrid cells that are generalists (from the overall generalist share and the hybrid numbers)
    _g, _h, _a = B[f"gen_pct_{_r}"] / 100, B.hybrid_pct / 100, B[f"hyb_gen_{_r}"] / 100
    B[f"oth_gen_{_r}"] = 100 * (_g - _a * _h) / (1 - _h)
SIG = B[B.p <= 0.05]
# ---------------------------------------------------------------- main-fit policy (agreed 22 Oct): one fit per unit and gene set
#   SCC: pooled = scRNA-only fit; within = the 8 scRNA tumours (multiome tumours and the all-cells pooled fit go to the appendix)
#   Cook, PDAC, 4T1, GSE114687: the cleaned arm (4T1 clone E has none -> light QC); light QC is shown only as a comparison
POOLED = {ds: E.FG.ORDER[ds][0] for ds in E.FG.ORDER}
APPX_UNITS = {"SCC": ["MU3", "MU4", "MU7", "MU8"]}
def pick_arm(ds, unit, gs, arm):
    d = B[(B.dataset == ds) & (B.unit == unit) & (B.geneset == gs) & (B.arm == arm)]
    return None if d.empty else d.iloc[0]
def pick_main(ds, unit, gs):
    if ds == "SCC": return pick_arm(ds, unit, gs, "scRNA only" if unit == POOLED[ds] else "light QC")
    r = pick_arm(ds, unit, gs, "cleaned")
    return r if r is not None else pick_arm(ds, unit, gs, "light QC")
def other_arm(ds, unit, gs, row):
    if ds == "SCC": return pick_arm(ds, unit, gs, "light QC") if unit == POOLED[ds] else None
    return pick_arm(ds, unit, gs, "light QC") if row.arm == "cleaned" else None
def main_units(ds): return [u for u in E.FG.ORDER[ds] if u not in APPX_UNITS.get(ds, [])]
MAIN = pd.DataFrame([r for ds in E.FG.ORDER for gs in ("all", "ks", "hallmark") for u in main_units(ds) for r in [pick_main(ds, u, gs)] if r is not None])
MAIN.drop(columns=["filt"], errors="ignore").to_csv(HERE / "main_fits.csv", index=False)
MS = MAIN[MAIN.p <= 0.05]
DS_NOTES = {
 "SCC": ["Hybrids = Early-Hybrid + Hybrid (the authors' labels); Late-EMT is not counted as hybrid.", "The four multiome tumours (MU3, MU4, MU7, MU8) were read as nuclei, not whole cells: a corner dominated by them is a technical axis, not biology.", "Within-tumour fits were run for all 12 tumours; most tumours contain one or two EMT states only."],
 "COOK": ["Cook has no state labels, so hybrids are defined from the KS EMT score (middle two quartiles).", "The cut-offs are pooled over the four lines, so a line's hybrid share mostly reflects the line (A549 about 71%, MCF7 about 11%).", "The removal arm (TGF-β withdrawn) is projected as guests and left out of the cells analysed."],
 "PDAC": ["EMT states are ours (KS-score bins), so they are partly circular with the KS-gene fits.", "Units: seven libraries, mouse 1, mouse 2. 'Cleaned' keeps cells with at least 10,000 UMI and drops liver ambient-RNA genes.", "Pooled fit = in vivo tumour cells; mouse and site are confounded."],
 "T4": ["One tumour per clone: clone, tumour and batch are the same thing.", "Clone E has only 142 cells after cleaning, so its cleaned arm is n/a.", "Only clone E is clearly separate on the EMT score; the 'spectrum' is one epithelial clone plus five similar ones."],
 "GSE": ["Colony position (inner vs outer) is an independent EMT coordinate; hybrids = Early- + Late-Hybrid (KS-score quartiles).", "Cleaning keeps only 45% of cells (HuMEC 35 to 40%), so cleaned and light-QC fits differ.", "No replicates: one experiment per line and arm."]}

def narrative(ds, gs, level, row, dall, ann):
    sig = row.p <= 0.05; n_sig = int((dall.p <= 0.05).sum()); n = len(dall); out = []
    if level == "pooled":
        if not sig: out.append(f"No simplex: the pooled fit (k={int(row.k)}, p={row.p:.3g}) is no more triangular than shuffled data, so its corners should not be read as tasks.")
        else:
            ident = ann[ann.top_sample_pct >= 70]
            if len(ident) >= max(1, len(ann) - 1):
                out.append(f"A significant simplex (k={int(row.k)}, p={row.p:.3g}), but {len(ident)} of {len(ann)} corners are almost entirely one sample ({', '.join(ident.top_sample.astype(str).unique()[:4])}). It separates samples; it is not a shared EMT trade-off.")
            else: out.append(f"A significant simplex (k={int(row.k)}, p={row.p:.3g}) whose corners mix samples (no corner is 70% or more one sample), so it is not a pure sample-identity simplex.")
    else:
        if not sig: out.append(f"No simplex inside {row.unit}: the fit (k={int(row.k)}, p={row.p:.3g}) is no more triangular than shuffled data, so its corners should not be read as tasks. {n_sig} of {n} units of this dataset have a significant simplex for this gene set.")
        else: out.append(f"Inside {row.unit} alone, cells fill a {int(row.k)}-corner simplex (p={row.p:.3g}). {n_sig} of {n} units of this dataset have a significant simplex for this gene set.")
    ident_all = level == "pooled" and len(ann) and (ann.top_sample_pct >= 70).sum() >= max(1, len(ann) - 1)
    if sig and len(ann) and ident_all:
        out.append("Each corner is essentially one " + {"COOK": "cell line", "SCC": "tumour", "PDAC": "library", "T4": "clone", "GSE": "sample"}[ds] + ". The programme names (Hallmark) only describe what distinguishes that " + {"COOK": "line", "SCC": "tumour", "PDAC": "library", "T4": "clone", "GSE": "sample"}[ds] + ", and the EMT-state colours follow the " + {"COOK": "lines (MCF7 scores epithelial, A549 mesenchymal)", "SCC": "tumours", "PDAC": "libraries", "T4": "clones", "GSE": "samples"}[ds] + ". These corners are not tasks.")
    elif sig and len(ann):
        tasks = list(dict.fromkeys(ann.task))
        out.append("Corners are programmes: " + "; ".join(tasks[:5]) + ".")
    hyb = row
    dl = row.hyb_gen_gap - row.oth_gen_gap
    out.append(f"Hybrid vs other cells: {row.hyb_gen_gap:.0f}% of hybrid cells and {row.oth_gen_gap:.0f}% of other cells are generalists (rule 1; rule 2: {row.hyb_gen_groves5:.0f}% vs {row.oth_gen_groves5:.0f}%), so {verdict(dl)}." + ("" if sig else " (No significant simplex, so this is descriptive only.)"))
    return out

def result_slide(ds, gs, level, obs, unit=None, row=None, label=None):
    unit = unit or POOLED[ds]
    if row is None: row = pick_main(ds, unit, gs)
    if row is None: return
    dall = MAIN[(MAIN.dataset == ds) & (MAIN.geneset == gs) & (MAIN.level == "within")]
    fn, ann = E.figure(ds, gs, ("pooled" if level == "pooled" else "within") + f"_{unit}_{row.arm}", row, obs)
    sig = row.p <= 0.05
    alt = other_arm(ds, unit, gs, row) if label is None else None
    cmp_ = f"   |   comparison, {alt.arm}: k={int(alt.k)}, p={alt.p:.3g}" if alt is not None else ""
    new_slide(); title((label + ": " if label else "") + f"{NAME[ds]}: {GS[gs]}, " + ("pooled fit" if level == "pooled" else f"unit {unit}"),
                       f"Fit shown: {row.arm} \u00b7 k = {int(row.k)} \u00b7 p = {row.p:.3g} \u00b7 {int(row.n_cells):,} cells \u00b7 {200 if row.n_null >= 150 else 100} shuffles" + cmp_)
    ww, hh = picture(fn, 0.4, 1.35, 7.6, 3.55)
    box(0.4, 1.35 + hh + 0.05, 7.6, 0.45, [para("Left: cells coloured by EMT state (blue epithelial → red mesenchymal; hybrids are the two middle colours). Right: the 5% of cells nearest each corner coloured by corner, all other cells grey = generalists.", 9, False, MUTED, after=0)], fit=True)
    y0 = 1.35 + hh + 0.55
    table(0.4, y0, 7.6, [["rule", "% of hybrid cells that are generalists", "% of other cells that are generalists"],
                         ["1  gap (top weight \u2212 second \u2265 0.35)", f"{row.hyb_gen_gap:.0f}%", f"{row.oth_gen_gap:.0f}%"],
                         ["2  closest 5% of cells to each corner", f"{row.hyb_gen_groves5:.0f}%", f"{row.oth_gen_groves5:.0f}%"],
                         ["3  weight > 0.95 on one corner", f"{row.hyb_gen_w95:.0f}%", f"{row.oth_gen_w95:.0f}%"]], [3.6, 2.0, 2.0], size=10)
    box(0.4, y0 + 1.25, 7.6, 0.35, [para(f"Hybrids are {row.hybrid_pct:.0f}% of all cells here; {row.inside_pct:.0f}% of cells lie inside the simplex.", 9.5, False, MUTED, after=0)])
    # right column
    tag = ("SIGNIFICANT", GREEN) if sig else ("NOT SIGNIFICANT", RED)
    box(8.2, 1.35, 4.7, 0.4, [para([(tag[0], True, tag[1]), ("   permutation test, p ≤ 0.05" if sig else "   p > 0.05", False, MUTED)], 12, after=0)], fit=False)
    nar = narrative(ds, gs, level, row, dall, ann)
    box(8.2, 1.78, 4.7, 2.15, [para("In plain words", 11, True, BROWN, after=2)] + [para(t, 10.5, False, INK, after=4) for t in nar], fill=PANEL, border="D8C9A6")
    cps = [para("Corners and their tasks" + ("" if sig else " (shown for completeness only)"), 11, True, BROWN, after=2)]
    for _, c in ann.iterrows():
        ev = f"{c.top_set} +{c.top_effect} SD" if c.top_effect > 0.5 else "no Hallmark set above 0.5 SD"
        bits = [ev]
        if not np.isnan(c.cycling_pct): bits.append(f"{c.cycling_pct:.0f}% S/G2M")
        bits.append(f"{c.top_state_pct:.0f}% {c.top_state}")
        if c.top_sample_pct >= 60: bits.append(f"{c.top_sample_pct:.0f}% {c.top_sample}")
        if not np.isnan(c.get("multiome_pct", np.nan)) and c.multiome_pct >= 70: bits.append(f"{c.multiome_pct:.0f}% multiome = technical")
        idc = level == "pooled" and c.top_sample_pct >= 70
        cps.append(para([(f"{c.corner} identity corner: {c.top_sample}" if idc else f"{c.corner} {c.task}", True, RED if idc else INK), (f" (programme: {c.task})" if idc else "", False, INK), (" — " + "; ".join(bits) + "; " + " ".join(str(c.top_genes).split()[:3]), False, MUTED)], 9.5, after=3))
    note = {"hallmark": "Hallmark EMT set left out when naming (circular: fitted on those genes).", "ks": "KS genes overlap the EMT set: an 'EMT / matrix' label is partly expected."}.get(gs)
    if note: cps.append(para(note, 8.5, False, MUTED, after=0))
    box(8.2, 4.0, 4.7, 3.3, cps, fill=PANEL, border="D8C9A6")

# ---------------------------------------------------------------- slides
def toy_rules_figure():
    rng = np.random.default_rng(3); V = np.array([[0, 0], [1, 0], [.5, .87]])
    w = rng.dirichlet([.9, .9, .9], 500); pts = w @ V
    gen, inside, near = A.rules(pts, V)
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.6), gridspec_kw=dict(left=.01, right=.99, top=.84, bottom=.02, wspace=.05))
    names = {"gap": "Rule 1: gap", "groves5": "Rule 2: closest 5%", "w95": "Rule 3: weight > 0.95"}
    for ax, r in zip(axs, ("gap", "groves5", "w95")):
        g = gen[r]; ax.scatter(pts[g, 0], pts[g, 1], s=6, color=E.GEN_COL, lw=0)
        for j in range(3): m = (~g) & (near == j); ax.scatter(pts[m, 0], pts[m, 1], s=9, color=FB.VCOL[j], lw=0)
        FB.outline(ax, V, "#" + INK); ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal")
        ax.set_title(f"{names[r]}\n{100 * (~g).mean():.0f}% specialists, {100 * g.mean():.0f}% generalists", fontsize=10, color="#" + INK)
    fn = E.OUT / "toy_rules.png"; fig.savefig(fn, dpi=150, facecolor=fig.get_facecolor()); plt.close(fig); return fn

RULES = (("gap", "Rule 1 gap", "#2b6cb0", "o"), ("groves5", "Rule 2 closest 5%", "#d97706", "s"), ("w95", "Rule 3 weight > 0.95", "#7c3aed", "D"))
def q2_figure(ds, kind):
    units = main_units(ds)
    fig, axs = plt.subplots(1, 3, figsize=(12.3, 0.42 * len(units) + 1.4), sharey=True, gridspec_kw=dict(left=.13, right=.99, top=.86, bottom=.16 if len(units) > 6 else .22, wspace=.06))
    for ax, gs in zip(axs, ("all", "ks", "hallmark")):
        for i, u in enumerate(units):
            r = pick_main(ds, u, gs)
            if r is None: continue
            y = len(units) - 1 - i; sig = r.p <= 0.05
            if kind == "gen_hyb": ax.plot([r.hybrid_pct] * 2, [y - .32, y + .32], color="#" + INK, lw=2)
            for k, nm, col, mk in RULES:
                v = r[f"{kind}_{k}"]; ax.scatter(v, y, s=46, marker=mk, facecolor=col if sig else "none", edgecolor=col, lw=1.4, zorder=3)
        ax.set_xlim(-2, 102); ax.set_title(GS[gs] + "\n(filled = significant simplex, hollow = not)", fontsize=9.5, color="#" + INK)
        ax.grid(axis="x", color="#ddd2b8"); ax.set_axisbelow(True); ax.axvline(50, color="#cdbf9f", lw=.8)
        ax.set_xlabel("% of hybrid cells that are generalists" if kind == "hyb_gen" else "% of generalist cells that are hybrid", fontsize=9, color="#" + INK)
        for s_ in ax.spines.values(): s_.set_color("#cdbf9f")
        ax.tick_params(colors="#" + INK, labelsize=8.5)
    axs[0].set_yticks(range(len(units))); axs[0].set_yticklabels([u if u != POOLED[ds] else "pooled" for u in units][::-1], fontsize=9, color="#" + INK)
    hs = [plt.Line2D([], [], marker=mk, ls="", color=col, markersize=7, label=nm) for k, nm, col, mk in RULES]
    if kind == "gen_hyb": hs.append(plt.Line2D([], [], color="#" + INK, lw=2, label="hybrid share of all cells (no-enrichment line)"))
    fig.legend(handles=hs, loc="lower center", ncol=4, frameon=False, fontsize=9, labelcolor="#" + INK)
    fn = E.OUT / f"q2_{ds}_{kind}.png"; fig.savefig(fn, dpi=150, facecolor=fig.get_facecolor()); plt.close(fig); return fn

HYC, OTC = "#dd6b20", "#a99f8c"
SHORT = {"SCC": "Mouse SCC", "COOK": "Cook", "PDAC": "Mouse PDAC", "T4": "4T1 clones", "GSE": "MCF10A /\nHuMEC"}
def q2_simple_summary():
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=True, gridspec_kw=dict(left=.07, right=.99, top=.80, bottom=.2, wspace=.08))
    for ax, (r, t) in zip(axs, (("gap", "Rule 1: purity gap"), ("groves5", "Rule 2: closest 5%"))):
        for i, d in enumerate(E.FG.ORDER):
            x = MS[MS.dataset == d]
            for off, col, c in ((-.19, HYC, f"hyb_gen_{r}"), (.19, OTC, f"oth_gen_{r}")):
                m = x[c].median(); ax.bar(i + off, m, .36, color=col, zorder=2); ax.text(i + off, m + 1.5, f"{m:.0f}%", ha="center", fontsize=10, color="#" + INK, fontweight="bold")
        ax.set_xticks(range(5)); ax.set_xticklabels([SHORT[d] for d in E.FG.ORDER], fontsize=10, color="#" + INK); ax.set_ylim(0, 100); ax.set_title(t, fontsize=12, color="#" + INK)
        ax.grid(axis="y", color="#e2d6bb", zorder=0)
        for s_ in ax.spines.values(): s_.set_color("#cdbf9f")
    axs[0].set_ylabel("% of cells that are generalists", fontsize=11, color="#" + INK)
    fig.suptitle("Are hybrid cells more often generalists than other cells?", fontsize=15, color="#" + INK, fontweight="bold")
    fig.text(.5, .855, "Median over the tumours / cell lines / samples that have a significant simplex.  Same height = hybrids are not special.", ha="center", fontsize=10, color="#" + MUTED)
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=HYC, label="hybrid cells"), plt.Rectangle((0, 0), 1, 1, color=OTC, label="all other cells (epithelial + mesenchymal)")], loc="lower center", ncol=2, frameon=False, fontsize=10.5, labelcolor="#" + INK)
    fn = E.OUT / "q2_simple_summary.png"; fig.savefig(fn, dpi=150, facecolor=fig.get_facecolor()); plt.close(fig); return fn

def q2_unit_bars(ds):
    units = main_units(ds)
    fig, axs = plt.subplots(1, 3, figsize=(12.3, 0.5 * len(units) + 1.5), sharey=True, gridspec_kw=dict(left=.12, right=.99, top=.86, bottom=.14 if len(units) > 6 else .2, wspace=.06))
    for ax, gs in zip(axs, ("all", "ks", "hallmark")):
        for i, u in enumerate(units):
            r = pick_main(ds, u, gs)
            if r is None: continue
            y = len(units) - 1 - i; sig = r.p <= 0.05; al = 1 if sig else .3
            ax.barh(y + .19, r.hyb_gen_gap, .36, color=HYC, alpha=al); ax.barh(y - .19, r.oth_gen_gap, .36, color=OTC, alpha=al)
            ax.text(max(r.hyb_gen_gap, r.oth_gen_gap) + 2, y, f"{r.hyb_gen_gap - r.oth_gen_gap:+.0f}" + ("" if sig else "  n.s."), va="center", fontsize=8.5, color="#" + INK)
        ax.set_xlim(0, 112); ax.set_title(GS[gs], fontsize=10, color="#" + INK); ax.grid(axis="x", color="#e2d6bb"); ax.set_axisbelow(True)
        ax.set_xlabel("% that are generalists (rule 1, purity gap)", fontsize=9, color="#" + INK)
        for s_ in ax.spines.values(): s_.set_color("#cdbf9f")
        ax.tick_params(colors="#" + INK, labelsize=8.5)
    axs[0].set_yticks(range(len(units))); axs[0].set_yticklabels([u if u != POOLED[ds] else "pooled" for u in units][::-1], fontsize=9, color="#" + INK)
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=HYC, label="hybrid cells"), plt.Rectangle((0, 0), 1, 1, color=OTC, label="all other cells"), plt.Rectangle((0, 0), 1, 1, color=HYC, alpha=.3, label="faded = no significant simplex (n.s.)")], loc="lower center", ncol=3, frameon=False, fontsize=9, labelcolor="#" + INK)
    fn = E.OUT / f"q2_units_{ds}.png"; fig.savefig(fn, dpi=150, facecolor=fig.get_facecolor()); plt.close(fig); return fn

def verdict(dlt):
    return "hybrids are clearly more often generalists" if dlt >= 8 else ("hybrids are slightly more often generalists" if dlt >= 3 else ("no real difference" if dlt > -3 else "hybrids are less often generalists (they sit at corners)"))

def q2_slides(ds):
    d = MS[MS.dataset == ds]
    new_slide(); title(f"Question 2 · {NAME[ds]}: are hybrid cells more often generalists than other cells?", "Orange = hybrid cells, grey = all other cells; one pair of bars per unit; the number = orange minus grey (percentage points)")
    ww, hh = picture(q2_unit_bars(ds), 0.5, 1.35, 12.3, 4.5)
    txt = [para("Reading", 11, True, BROWN, after=2)]
    for g in ("all", "ks", "hallmark"):
        x = d[d.geneset == g]
        if len(x):
            dl = (x.hyb_gen_gap - x.oth_gen_gap).median()
            txt.append(para(f"{GS[g]} ({len(x)} significant): hybrids {x.hyb_gen_gap.median():.0f}% vs others {x.oth_gen_gap.median():.0f}% generalist under rule 1 (rule 2: {x.hyb_gen_groves5.median():.0f}% vs {x.oth_gen_groves5.median():.0f}%) → {verdict(dl)}.", 10.5, False, INK, after=3, bullet=True))
        else: txt.append(para(f"{GS[g]}: no significant simplex, nothing to compare.", 10.5, False, INK, after=3, bullet=True))
    txt.append(para("Rule 3 (weight > 0.95) is not drawn: it calls about 100% of every group generalist, so it cannot show a difference.", 10, False, MUTED, after=0))
    box(0.5, 1.4 + hh + 0.05, 12.3, 7.3 - (1.4 + hh + 0.05), txt, fill=PANEL, border="D8C9A6")

def q2_summary():
    fn = q2_simple_summary()
    new_slide(); title("What it all says: Question 2", "Are hybrid cells more often generalists than other cells? (all five datasets)")
    ww, hh = picture(fn, 0.9, 1.3, 11.5, 4.5)
    dl = MS.groupby("dataset").apply(lambda x: (x.hyb_gen_gap - x.oth_gen_gap).median(), include_groups=False)
    clear = [NAME[k] for k, v in dl.items() if v >= 8]; same = [NAME[k] for k, v in dl.items() if v < 8]
    box(0.5, 1.35 + hh + 0.05, 12.3, 7.3 - (1.35 + hh + 0.05), [
        para([("Most cells are generalists, hybrid or not. ", True, BROWN), (f"Under rule 1 a median {MS.hyb_gen_gap.median():.0f}% of hybrid cells and {MS.oth_gen_gap.median():.0f}% of other cells are generalists; under rule 2 {MS.hyb_gen_groves5.median():.0f}% vs {MS.oth_gen_groves5.median():.0f}%.", False, INK)], 12, after=5),
        para([("Hybrids lean slightly toward the interior. ", True, BROWN), ("Median difference, hybrids minus other cells (rule 1): " + ", ".join(f"{NAME[k]} {int(round(v)) + 0:+d}" for k, v in dl.sort_values(ascending=False).items()) + " points. The tilt is clearest in Mouse SCC and absent in Cook.", False, INK)], 12, after=5),
        para([("No hybrid corner. ", True, BROWN), ("In no dataset do hybrid cells pile up at a corner of their own: hybrids are mixtures in the interior about as often as, or slightly more often than, other cells.", False, INK)], 12, after=0)], fill=PANEL, border="D8C9A6")

def main():
    obs = A.load_obs()
    new_slide(); box(0.9, 1.7, 11.5, 1.6, [para("Simplex results, explained", 38, True, INK, after=6), para("Five single-cell EMT datasets · three gene sets · pooled and within each unit", 18, False, BROWN, after=0)], fit=False)
    box(0.9, 4.2, 11.5, 1.8, [para("What is checked: (1) is there a significant simplex, and do its corners correspond to tasks, when the fit uses A. the entire gene set, B. the KS EMT genes, C. the Hallmark EMT genes; (2) how many hybrid cells are generalists, and how many generalists are hybrid, under three rules.", 14, False, INK, after=8), para("Thesis project, IISc · 22 October 2026", 12, False, MUTED, after=0)], fit=False)

    new_slide(); title("Key ideas in five lines", "Read this once; every later slide uses these words")
    box(0.6, 1.5, 12.1, 5.6, [
        para([("Simplex. ", True, BROWN), ("If cells must trade off k tasks, they fill a shape (a line, triangle, tetrahedron...) whose k corners are cells specialised in one task each. We fit that shape to the leading PCA axes of the gene-expression data.", False, INK)], 14, after=9),
        para([("Significant simplex. ", True, BROWN), ("The fitted shape is compared with the same data after shuffling each PCA axis independently, 100 to 200 times. If the real shape is more triangle-like than nearly all shuffles (p ≤ 0.05), it is not just a random cloud. A small p alone does not prove a trade-off: samples that differ from each other can also give a significant shape.", False, INK)], 14, after=9),
        para([("Corner = task? ", True, BROWN), ("The 5% of cells nearest each corner are scored on the 50 MSigDB Hallmark gene sets (effect in SD versus all other cells). The top set names the programme (e.g. E2F Targets = cycling). Cell-cycle phase, EMT-state mix, dominant sample and marker genes support the name. A corner that is almost all one sample (≥ 70%) is an identity corner, not a task.", False, INK)], 14, after=9),
        para([("Hybrid cells. ", True, BROWN), ("The middle EMT states (Early-Hybrid, Hybrid, Late-Hybrid). Cook has no state labels, so hybrids are the middle two quartiles of the KS EMT score.", False, INK)], 14, after=9),
        para([("Specialist vs generalist. ", True, BROWN), ("A cell claimed by one corner is a specialist; every other cell is a generalist, an interior mixture. The key comparison: what % of hybrid cells are generalists, versus what % of all other cells are generalists. Same % = hybrids are not special.", False, INK)], 14, after=0)], fill=PANEL, border="D8C9A6")

    new_slide(); title("The three specialist rules on a toy triangle", "500 simulated cells mixing three corners; the same rules are applied to the real fits")
    fn = toy_rules_figure(); ww, hh = picture(fn, 0.5, 1.4, 12.3, 3.9)
    box(0.5, 5.35, 12.3, 1.95, [
        para([("Rule 1, gap (Pal work): ", True, BROWN), ("a cell is a specialist if its largest corner weight exceeds the second largest by at least 0.35.", False, INK)], 12, after=5),
        para([("Rule 2, closest 5% (Groves 2022): ", True, BROWN), ("for each corner, the 5% of cells nearest to it are its specialists; everything else is a generalist. Caps specialists at about 15% for a triangle.", False, INK)], 12, after=5),
        para([("Rule 3, weight > 0.95 (Groves 2022, methods): ", True, BROWN), ("specialist only if one weight exceeds 0.95. Because the fitted corners lie inside the cloud (delta = 0), almost no cell qualifies, so about 100% are generalists. It is the strictest and least informative.", False, INK)], 12, after=0)], fill=PANEL, border="D8C9A6")

    new_slide(); title("How the numbers on each slide were produced", "Everything is read from saved fits; no result is refitted for the rules")
    box(0.6, 1.5, 12.1, 5.7, [
        para([("1. Fit. ", True, BROWN), ("Each fit = PCA on the unit's cells (genes centred, not scaled) then archetypal analysis on the first k−1 PCs (15 starts); the permutation test shuffles each PC and refits (200 shuffles; 100 for the new Hallmark, cleaned-pooled and SCC scRNA-only / six-tumour pooled fits, marked in tables).", False, INK)], 13, after=7),
        para([("2. Best simplex. ", True, BROWN), ("For each dataset × unit × gene set × arm the saved k values are compared; best = lowest p, ties go to the smaller k; pooled fits ignore k beyond the ESV elbow. On each slide the arm shown is the one with the lowest p (ties: cleaned arm).", False, INK)], 13, after=7),
        para([("3. Weights and rules. ", True, BROWN), ("Each cell's barycentric weights on the k corners give the three rules (previous slide). Statistics use all cells of the fit; the plots draw at most 4,000 random cells. Within-unit fits are capped at 4,000 cells.", False, INK)], 13, after=7),
        para([("4. Corner names. ", True, BROWN), ("5% nearest each corner → Hallmark scores (mean log-expression per set; effect = difference in SD) + S/G2M share + state mix + dominant sample + top genes. For Hallmark-EMT fits the EMT set is excluded from naming (circular).", False, INK)], 13, after=7),
        para([("5. Arms. ", True, BROWN), ("Light QC = standard cell QC; cleaned = at least 10,000 UMI (plus dataset-specific technical-gene / contamination filters). SCC has no cleaned arm. For k ≥ 4 the simplex is drawn by projecting the corners on PC1/PC2 and joining their convex hull.", False, INK)], 13, after=0)], fill=PANEL, border="D8C9A6")

    # ---------------------------------------------------------------- which fit is shown
    new_slide(); title("Which fit is shown on each slide", "One main fit per unit and gene set, so nothing is picked for having the lowest p")
    box(0.6, 1.5, 12.1, 5.7, [
        para([("Mouse SCC. ", True, BROWN), ("Main = the pooled fit on scRNA cells only, and the eight scRNA tumours. The four multiome tumours (MU3, MU4, MU7, MU8) were sequenced as nuclei, a different technology: mixing them in creates a technical corner, and inside a multiome tumour the signal is nuclear RNA. The all-cells pooled fit and the multiome tumours are kept in the appendix.", False, INK)], 13, after=9),
        para([("Cook, PDAC, 4T1, GSE114687. ", True, BROWN), ("Main = the cleaned arm: cells with at least 10,000 UMI, plus each dataset's technical filter (Cook: mitochondrial fraction ≤ 0.07 and MALAT1 / NEAT1 / haemoglobin genes dropped; PDAC: liver ambient-RNA and haemoglobin genes dropped; 4T1: immune / stromal contamination removed; GSE114687: MALAT1 / NEAT1 / haemoglobin dropped). 4T1 clone E has only 142 cells after cleaning, so its light-QC fit is shown.", False, INK)], 13, after=9),
        para([("Why. ", True, BROWN), ("Light QC and cleaning disagree on whether a simplex is significant in many units (Cook 2 of 12 unit × gene-set pairs, PDAC 8 of 27, 4T1 7 of 15, GSE114687 8 of 18). Picking whichever arm has the lower p would inflate significance. The light-QC p is printed on every slide as a comparison and tabulated in the appendix.", False, INK)], 13, after=9),
        para([("Counting. ", True, BROWN), ("All summary numbers below count units (one main fit each), not fits.", False, INK)], 13, after=0)], fill=PANEL, border="D8C9A6")

    new_slide(); title("At a glance: where is a simplex significant?", "Main fits. Pooled: k and p. Within: number of units with a significant simplex")
    rows = [["dataset", "A. pooled", "A. within", "B. pooled", "B. within", "C. pooled", "C. within"]]
    for ds in E.FG.ORDER:
        r = [NAME[ds]]
        for gs in ("all", "ks", "hallmark"):
            pr = pick_main(ds, POOLED[ds], gs); dw = MAIN[(MAIN.dataset == ds) & (MAIN.geneset == gs) & (MAIN.level == "within")]
            r += [f"k={int(pr.k)}, p={pr.p:.3g}" + (" *" if pr.p <= .05 else ""), f"{int((dw.p <= .05).sum())} of {len(dw)}"]
        rows.append(r)
    table(0.5, 1.6, 12.3, rows, [3.0, 1.7, 1.2, 1.7, 1.2, 1.7, 1.2], size=12)
    box(0.5, 4.2, 12.3, 2.9, [para("* significant (p ≤ 0.05). Units: SCC 8 scRNA tumours, Cook 4 lines, PDAC 9 (7 libraries, mouse 1, mouse 2), 4T1 6 clones, GSE114687 6 samples.", 11, False, MUTED, after=6)] +
        [para(f"{GS[g]}: {int((MAIN[(MAIN.geneset == g)].p <= .05).sum())} of {int((MAIN.geneset == g).sum())} main fits significant (pooled {int(((MAIN.geneset == g) & (MAIN.level == 'pooled') & (MAIN.p <= .05)).sum())} of 5, within {int(((MAIN.geneset == g) & (MAIN.level == 'within') & (MAIN.p <= .05)).sum())} of {int(((MAIN.geneset == g) & (MAIN.level == 'within')).sum())} units).", 13, False, INK, after=4) for g in ("all", "ks", "hallmark")], fill=PANEL, border="D8C9A6")

    for ds in E.FG.ORDER:
        d = MAIN[MAIN.dataset == ds]
        new_slide(); title(NAME[ds], "What the dataset is, and what to keep in mind when reading its slides")
        box(0.5, 1.5, 6.2, 5.7, [para("What the dataset means", 13, True, BROWN, after=4), para(E.DS_MEANING[ds], 12, False, INK, after=8), para("Units in the main analysis", 13, True, BROWN, after=4), para(", ".join(u for u in main_units(ds) if u != POOLED[ds]), 11, False, INK, after=8), para("Gene sets", 13, True, BROWN, after=4),
            para("A. all genes (11,000 to 14,500); B. KS EMT genes (Tan 2014; about 160 to 190 present); C. Hallmark Epithelial-Mesenchymal-Transition set (200 genes, those present in the data).", 11, False, INK, after=0)], fill=PANEL, border="D8C9A6")
        sg = d[d.p <= .05]
        box(6.9, 1.5, 6.0, 5.7, [para("Keep in mind", 13, True, BROWN, after=4)] + [para(t, 11.5, False, INK, after=6, bullet=True) for t in DS_NOTES[ds]] +
            [para("This dataset in numbers (main fits)", 13, True, BROWN, after=4)] +
            [para(f"{GS[g][3:]}: {int((sg.geneset == g).sum())} of {int((d.geneset == g).sum())} fits significant" + (f"; hybrids {sg[sg.geneset == g].hyb_gen_gap.median():.0f}% / {sg[sg.geneset == g].hyb_gen_groves5.median():.0f}% / {sg[sg.geneset == g].hyb_gen_w95.median():.0f}% generalist (median, rules 1/2/3)." if (sg.geneset == g).any() else "."), 11, False, INK, after=4, bullet=True) for g in ("all", "ks", "hallmark")], fill=PANEL, border="D8C9A6")
        box(6.9, 6.75, 6.0, 0.4, [para("Question 1 slides come first, then the Question 2 slides for this dataset.", 9.5, False, MUTED, after=0)], fit=False)
        for gs in ("all", "ks", "hallmark"):
            result_slide(ds, gs, "pooled", obs)
            units = [u for u in main_units(ds) if u != POOLED[ds]]
            dw = MAIN[(MAIN.dataset == ds) & (MAIN.geneset == gs) & (MAIN.level == "within")]
            new_slide(); title(f"{NAME[ds]}: {GS[gs]}, all units at a glance", f"Main fit of each unit: {int((dw.p <= .05).sum())} of {len(dw)} units have a significant simplex")
            rows_ = [["unit", "fit shown", "k", "p", "light-QC p", "hybrid % of cells", "hybrids→generalist (1/2/3)", "generalists→hybrid (1/2/3)"]]
            for u in units:
                r_ = pick_main(ds, u, gs)
                if r_ is None: rows_.append([u, "not run", "", "", "", "", "", ""]); continue
                al = other_arm(ds, u, gs, r_)
                rows_.append([u, r_.arm, int(r_.k), f"{r_.p:.3g}" + (" *" if r_.p <= .05 else ""), (f"{al.p:.3g}" if al is not None else "–"), f"{r_.hybrid_pct:.0f}%", f"{r_.hyb_gen_gap:.0f} / {r_.hyb_gen_groves5:.0f} / {r_.hyb_gen_w95:.0f}%", f"{r_.gen_hyb_gap:.0f} / {r_.gen_hyb_groves5:.0f} / {r_.gen_hyb_w95:.0f}%"])
            table(0.5, 1.45, 12.3, rows_, [1.6, 1.2, 0.5, 1.1, 1.1, 1.5, 2.7, 2.6], size=10.5)
            box(0.5, min(1.45 + 0.3 * len(rows_) + 0.15, 6.3), 12.3, 0.9, [para("* p ≤ 0.05 for the fit shown. Rules: 1 = gap, 2 = closest 5%, 3 = weight > 0.95. One slide per unit follows; non-significant units are included for completeness and their corners should not be read as tasks.", 10.5, False, MUTED, after=0)])
            for u in units: result_slide(ds, gs, "within", obs, unit=u)
        q2_slides(ds)
        print(ds, "slides done", flush=True)

    q2_summary()

    Cs = CA.merge(MS[["dataset", "level", "unit", "geneset", "arm", "k"]], on=["dataset", "level", "unit", "geneset", "arm", "k"])
    ident = (Cs[Cs.level == "pooled"].top_sample_pct >= 70).mean() * 100
    wc = Cs[Cs.level == "within"].task.value_counts().head(6)
    hybc = ((Cs.top_state.isin(["Early-Hybrid", "Hybrid", "Late-Hybrid"])) & (Cs.top_state_pct >= 60)); emtc = Cs[Cs.task == "EMT / matrix"]
    new_slide(); title("What it all says: Question 1", "Is there a simplex, and are its corners tasks? (main fits)")
    box(0.6, 1.5, 12.1, 5.7, [
        para([("Significance depends on the gene set and the level. ", True, BROWN), ("Main fits significant: " + "; ".join(f"{GS[g]} {int((MS.geneset == g).sum())} of {int((MAIN.geneset == g).sum())}" for g in ("all", "ks", "hallmark")) + ". Hallmark behaves like KS rather than like the entire gene set.", False, INK)], 13, after=8),
        para([("Pooled simplices mostly separate samples. ", True, BROWN), (f"{ident:.0f}% of the corners of significant pooled simplices are at least 70% one sample, tumour, line or clone.", False, INK)], 13, after=8),
        para([("Within one unit the corners are programmes, not epithelial / hybrid / mesenchymal states. ", True, BROWN), ("Most frequent corner names in significant within-unit fits: " + "; ".join(f"{k} ({v})" for k, v in wc.items()) + f". {100 * hybc.mean():.0f}% of corners are hybrid-dominated (≥ 60% of their specialists), mainly cycling or MYC/growth corners; {int(emtc.top_state.isin(['Mesenchymal', 'Late-EMT']).sum())} of {len(emtc)} 'EMT / matrix' corners are mesenchymal or late-EMT.", False, INK)], 13, after=8),
        para([("So: ", True, BROWN), ("a trade-off structure exists inside single units (proliferation, growth, differentiation, EMT / matrix), but EMT states do not form their own corners and no dedicated hybrid corner is seen.", False, INK)], 13, after=0)], fill=PANEL, border="D8C9A6")
    new_slide(); title("Caveats", "What these results do not show")
    box(0.6, 1.5, 12.1, 5.7, [para(t, 13, False, INK, after=8, bullet=True) for t in [
        "Many fits per dataset (gene sets, k, arms); no multiplicity correction. p cannot go below about 0.005 (200 shuffles) or 0.01 (100 shuffles).",
        "Hallmark pooled, cleaned-pooled and SCC scRNA-only KS / Hallmark pooled fits used k = 3 to 5 and 100 shuffles.",
        "Hybrid labels for PDAC, 4T1, GSE114687 and Cook are cut from the same KS score used for the KS-gene fits (circular); Cook bins are pooled across lines.",
        "Within-unit fits are capped at 4,000 random cells; corner names are automatic (top Hallmark set) and should be sanity-checked against biology.",
        "Design limits: one tumour per clone (4T1), two mice (PDAC), no replicates (GSE114687), tumour / assay / state tangled (SCC). Nothing here shows that EMT is or is not a Pareto trade-off in tumours.",
        "Files: main_fits.csv (the fits shown), best_fits.csv (all arms), corner_annotations.csv; scripts analyze.py, corners.py, build_explainer_ppt.py in docs/emt_single_cell/best_simplex/."]], fill=PANEL, border="D8C9A6")

    # ---------------------------------------------------------------- appendix
    new_slide(); box(0.9, 2.4, 11.5, 2.5, [para("Appendix", 36, True, INK, after=8), para("A1. Mouse SCC: the all-cells pooled fit (with multiome nuclei) and the four multiome tumours.   A2. Light QC versus cleaned arm, per unit, for Cook, PDAC, 4T1 and GSE114687.", 15, False, BROWN, after=0)], fit=False)
    for gs in ("all", "ks", "hallmark"):
        r_ = pick_arm("SCC", "pooled", gs, "light QC")
        if r_ is not None: result_slide("SCC", gs, "pooled", obs, unit="pooled", row=r_, label="Appendix A1")
        for u in APPX_UNITS["SCC"]:
            r_ = pick_arm("SCC", u, gs, "light QC")
            if r_ is not None: result_slide("SCC", gs, "within", obs, unit=u, row=r_, label="Appendix A1")
    for ds in ("COOK", "PDAC", "T4", "GSE"):
        new_slide(); title(f"Appendix A2. {NAME[ds]}: light QC vs cleaned", "p of each arm (* = p ≤ 0.05); the cleaned arm is the one used in the main slides")
        rows_ = [["unit", "A. light QC", "A. cleaned", "B. light QC", "B. cleaned", "C. light QC", "C. cleaned"]]
        for u in E.FG.ORDER[ds]:
            rr = [u]
            for gs in ("all", "ks", "hallmark"):
                for arm in ("light QC", "cleaned"):
                    x = pick_arm(ds, u, gs, arm); rr.append("–" if x is None else f"k={int(x.k)}, {x.p:.3g}" + (" *" if x.p <= .05 else ""))
            rows_.append(rr)
        table(0.5, 1.45, 12.3, rows_, [2.1, 1.7, 1.7, 1.7, 1.7, 1.7, 1.7], size=10.5)
    prs.save(HERE / "Explainer_simplex_results.pptx")
    css = "<!doctype html><meta charset=utf-8><style>@page{size:13.333in 7.5in;margin:0}body{margin:0;font-family:Helvetica,Arial,sans-serif}.pg{position:relative;width:13.333in;height:7.5in;background:#%s;overflow:hidden;page-break-after:always;-webkit-print-color-adjust:exact;print-color-adjust:exact}</style>" % BEIGE
    (HERE / "Explainer_simplex_results.html").write_text(css + "".join(HT) + "</section>")
    print("saved", len(prs.slides), "slides")

if __name__ == "__main__":
    main()
