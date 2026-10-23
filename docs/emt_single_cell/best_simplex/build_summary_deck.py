#!/usr/bin/env python3
"""Short visual summary (7 slides, beige) of the two research questions, for a quick read. Reuses the explainer deck's helpers and the
report figures. Run from Project_1/ after build_final_report.py:  .venv/bin/python docs/emt_single_cell/best_simplex/build_summary_deck.py"""
import sys, io, contextlib
from pathlib import Path
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
with contextlib.redirect_stdout(io.StringIO()): import build_explainer_ppt as P
from pptx import Presentation
from pptx.util import Inches
P.prs = Presentation(); P.prs.slide_width, P.prs.slide_height = Inches(13.333), Inches(7.5); P.HT.clear()
para, box, title, picture, new_slide = P.para, P.box, P.title, P.picture, P.new_slide
INK, BROWN, MUTED, PANEL = P.INK, P.BROWN, P.MUTED, P.PANEL
MAIN, MS, CA = P.MAIN, P.MS, P.CA
F = HERE / "figures" / "report"
KEY = ["dataset", "level", "unit", "geneset", "arm", "k"]
ANN = CA.drop(columns=["p"]).merge(MAIN[KEY + ["p"]], on=KEY)
ident = 100 * (ANN[(ANN.level == "pooled") & (ANN.p <= .05)].top_sample_pct >= 70).mean()
n = {g: (int((MS.geneset == g).sum()), int((MAIN.geneset == g).sum())) for g in ("all", "ks", "hallmark")}
dd = MS.assign(d=MS.hyb_gen_gap - MS.oth_gen_gap).groupby("dataset").d.median().sort_values(ascending=False)
dtxt = ", ".join(f"{P.NAME[k]} {int(round(v)) + 0:+d}" for k, v in dd.items())

def side(lines, x=8.9, y=1.45, w=4.0, h=5.7):
    box(x, y, w, h, [para(t, 13, b, BROWN if b else INK, after=8) for t, b in lines], fill=PANEL, border="D8C9A6")

# 1. title + answers
new_slide(); box(0.7, 0.4, 12, 1.6, [para("EMT archetype analysis across five single-cell datasets: summary", 28, True, INK, after=6), para("Mouse SCC · Cook TGF-β lines · Mouse PDAC · 4T1 clones · MCF10A / HuMEC  —  22 October 2026", 14, False, BROWN, after=0)], fit=False)
box(0.7, 2.1, 5.9, 4.9, [para("Q1. Is there a significant simplex, and are its corners tasks?", 15, True, BROWN, after=8), para("(fit on A. all genes, B. KS EMT genes, C. Hallmark EMT genes; pooled and within each tumour / cell line)", 11, False, MUTED, after=10),
    para(f"Often significant (A {n['all'][0]}/{n['all'][1]}, B {n['ks'][0]}/{n['ks'][1]}, C {n['hallmark'][0]}/{n['hallmark'][1]} fits).", 14, False, INK, after=6),
    para("Pooled: the corners are mostly the samples themselves (Cook: the four cell lines).", 14, False, INK, after=6),
    para("Within one tumour or line: the corners are programmes (cycling, growth, differentiated epithelial, metabolic, EMT / matrix), not EMT states. No hybrid corner.", 14, False, INK, after=0)], fill=PANEL, border="D8C9A6")
box(6.9, 2.1, 5.9, 4.9, [para("Q2. Are hybrid cells generalists, and more than other cells?", 15, True, BROWN, after=8), para("(three specialist rules: purity gap, closest 5%, weight > 0.95)", 11, False, MUTED, after=10),
    para(f"Most cells are generalists, hybrid or not: {MS.hyb_gen_gap.median():.0f}% of hybrids vs {MS.oth_gen_gap.median():.0f}% of other cells (gap rule).", 14, False, INK, after=6),
    para(f"Hybrids lean slightly toward the interior; difference by dataset: {dtxt} points.", 14, False, INK, after=6),
    para("The weight > 0.95 rule calls ~100% of cells generalist, so it cannot separate groups.", 14, False, INK, after=0)], fill=PANEL, border="D8C9A6")

# 2. Q1 at a glance
new_slide(); title("Q1 at a glance: where is a simplex significant?", "Main fits; pooled = k, p and identity corners; within = units with a significant simplex")
picture(F / "overview.png", 0.4, 1.4, 8.3, 5.8)
side([("How to read", True), ("Dark green = significant. 'Identity' = corners that are ≥ 70% one sample.", False), ("What it shows", True),
      ("All genes give a simplex most often; the two EMT gene sets less often, and Hallmark behaves like KS.", False), (f"Pooled simplices are mostly sample identity ({ident:.0f}% of their corners).", False)])

# 3. pooled = identity
new_slide(); title("Pooled fits separate samples, not EMT states", "Example: Cook, the four TGF-β cell lines (all genes / KS / Hallmark)")
picture(F / "pooled_COOK.png", 0.4, 1.4, 8.3, 5.8)
side([("What it shows", True), ("Each corner is one cell line (MCF7, A549, DU145, OVCA420), in all three gene sets.", False),
      ("The EMT colours follow the lines (MCF7 epithelial, A549 mesenchymal), which makes line corners look like EMT corners.", False), ("Why it matters", True), ("Pooled shapes must be read per sample: within-unit fits are the informative ones.", False)])

# 4. within = programmes
new_slide(); title("Within one tumour or cell line, corners are programmes", "Names of all corners of significant within-unit fits (top Hallmark gene set at each corner)")
picture(F / "tasks.png", 0.4, 1.4, 8.3, 5.8)
side([("What it shows", True), ("Recurring corners: differentiated (KRAS-down) epithelial, MYC / growth, EMT / matrix, metabolic, cycling.", False),
      ("EMT / matrix corners are mostly mesenchymal / late-EMT cells; there is no hybrid corner.", False), ("Reading", True), ("EMT behaves like an axis between epithelial and EMT / matrix corners, alongside other tasks.", False)])

# 5. example within-unit simplex
new_slide(); title("Example: within-sample simplices in MCF10A / HuMEC", "Most significant sample per gene set; top = EMT state, bottom = corners (grey = generalists)")
picture(F / "within_best_GSE.png", 0.4, 1.4, 8.3, 5.8)
bst = []
for g in ("all", "ks", "hallmark"):
    d = MAIN[(MAIN.dataset == "GSE") & (MAIN.geneset == g) & (MAIN.level == "within")].sort_values(["p", "n_cells"], ascending=[True, False]).iloc[0]
    a = ANN[(ANN.dataset == "GSE") & (ANN.unit == d.unit) & (ANN.geneset == g) & (ANN.arm == d.arm) & (ANN.k == d.k)]
    bst.append(f"{P.GS[g][:2]} {d.unit}: " + ", ".join(dict.fromkeys(a.task.str.split(" \\(").str[0])))
side([("What it shows", True), ("Inside one sample, cells fill a simplex whose corners are programmes:", False)] + [(t, False) for t in bst] + [
      ("Hybrid cells (yellow / orange) lie between corners, as mixtures.", False)])

# 6. Q2
new_slide(); title("Q2: are hybrid cells more often generalists than other cells?", "Orange = hybrid cells, grey = all other cells; median over units with a significant simplex")
picture(P.E.OUT / "q2_simple_summary.png", 0.4, 1.4, 8.4, 5.6)
side([("What it shows", True), ("Most cells are generalists either way.", False), (f"Hybrids are slightly more often generalists: {dtxt} points (gap rule).", False),
      ("Closest-5% rule: same pattern, higher overall. Weight > 0.95: ~100% everywhere.", False)])

# 7. takeaways
new_slide(); title("Takeaways and caveats", "Full methods and every tumour / cell line are in the attached report and deck")
box(0.6, 1.5, 6.0, 5.7, [para("Takeaways", 15, True, BROWN, after=8)] + [para(t, 13.5, False, INK, after=8, bullet=True) for t in [
    "A simplex is common, but pooled shapes are mostly sample identity.",
    "Within units the corners are programmes (growth, differentiation, metabolic, EMT / matrix), not EMT states.",
    "Hybrids are interior mixtures, slightly more often generalists than other cells; no hybrid corner.",
    "Consistent with EMT as a continuum running between programmes, rather than discrete task states."]], fill=PANEL, border="D8C9A6")
box(6.9, 1.5, 5.9, 5.7, [para("Caveats", 15, True, BROWN, after=8)] + [para(t, 12.5, False, INK, after=7, bullet=True) for t in [
    "Main fits: SCC scRNA-only; cleaned arm elsewhere (light QC in the appendix).",
    "Newest fits (Hallmark, cleaned pooled) used 100 shuffles and k = 3 to 5.",
    "Hybrid labels in four datasets come from the KS score (partly circular for KS fits).",
    "Few samples per condition; many fits, no multiplicity correction.",
    "Corner names are automatic (top Hallmark set)."]], fill=PANEL, border="D8C9A6")

P.prs.save(HERE / "Summary_simplex_results.pptx")
css = "<!doctype html><meta charset=utf-8><style>@page{size:13.333in 7.5in;margin:0}body{margin:0;font-family:Helvetica,Arial,sans-serif}.pg{position:relative;width:13.333in;height:7.5in;background:#%s;overflow:hidden;page-break-after:always;-webkit-print-color-adjust:exact;print-color-adjust:exact}</style>" % P.BEIGE
(HERE / "Summary_simplex_results.html").write_text(css + "".join(P.HT) + "</section>")
print("saved", len(P.prs.slides), "slides")
