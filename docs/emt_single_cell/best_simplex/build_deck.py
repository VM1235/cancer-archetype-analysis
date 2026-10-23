#!/usr/bin/env python3
"""Build the dark-theme deck from best_fits.csv and figures/. Run from Project_1/: .venv/bin/python docs/emt_single_cell/best_simplex/build_deck.py"""
from pathlib import Path
import pandas as pd, numpy as np
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from PIL import Image

HERE = Path(__file__).resolve().parent; FIG = HERE / "figures"
B = pd.read_csv(HERE / "best_fits.csv")
NAME = dict(SCC="Mouse skin SCC (GSE288137)", COOK="Cook TGF-β cell lines (GSE147405)", PDAC="Mouse PDAC (GSE173958)", T4="4T1 clones (GSE303505)", GSE="MCF10A / HuMEC pseudospace (GSE114687)")
BG, FG, MUT, ACC = RGBColor(0x14, 0x17, 0x1c), RGBColor(0xe8, 0xe8, 0xe8), RGBColor(0x9a, 0xa0, 0xa6), RGBColor(0x63, 0xb3, 0xed)
HT = ["<!doctype html><meta charset=utf-8><style>@page{size:13.333in 7.5in;margin:0}body{margin:0;background:#14171c;color:#e8e8e8;font:15px/1.45 -apple-system,Helvetica,Arial,sans-serif;-webkit-print-color-adjust:exact;print-color-adjust:exact}.pg{width:13.333in;height:7.5in;box-sizing:border-box;padding:.35in .5in;page-break-after:always;overflow:hidden;position:relative}h1{font-size:26px;margin:0 0 12px}p{margin:0 0 7px;font-size:16px}p.b{color:#63b3ed;font-weight:700}img{display:block;margin:0 auto;max-width:100%;max-height:5.9in}table{border-collapse:collapse;width:100%;font-size:11px}th,td{border:1px solid #2a2f37;padding:3px 5px}th{background:#2b4a6f}tr:nth-child(even) td{background:#171b21}.m{color:#9aa0a6}</style>"]
prs = Presentation(); prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BL = prs.slide_layouts[6]

def slide(title, notes=None):
    s = prs.slides.add_slide(BL); s.background.fill.solid(); s.background.fill.fore_color.rgb = BG
    tb = s.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(12.3), Inches(0.9)); tf = tb.text_frame; tf.word_wrap = True
    p = tf.paragraphs[0]; p.text = title; p.font.size = Pt(26); p.font.bold = True; p.font.color.rgb = FG
    HT.append(f'</section><section class="pg">' + (f'<h1>{title}</h1>' if title else ''))
    return s

def text(s, lines, x=0.6, y=1.4, w=12.1, h=5.5, size=16, color=FG):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)); tf = tb.text_frame; tf.word_wrap = True
    for i, l in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        bold = l.startswith("**"); l = l.strip("*")
        p.text = l; p.font.size = Pt(size); p.font.color.rgb = ACC if bold else color; p.font.bold = bold; p.space_after = Pt(7)
        HT.append(f'<p class="{"b" if bold else ("m" if color == MUT else "")}">{l}</p>')

def pic(s, name, x=0.5, y=1.2, w=12.3, h=5.9):
    HT.append(f'<img src="figures/{name}" style="max-height:6.9in">')
    im = Image.open(FIG / name); ar = im.width / im.height
    ww, hh = (w, w / ar) if w / ar <= h else (h * ar, h)
    s.shapes.add_picture(str(FIG / name), Inches(x + (w - ww) / 2), Inches(y), Inches(ww), Inches(hh))

def table(s, df, cols, heads, y=1.3, size=9.5, widths=None):
    HT.append("<table><tr>" + "".join(f"<th>{h}</th>" for h in heads) + "</tr>" + "".join("<tr>" + "".join("<td>" + (f"{r[c]:.3g}" if c == "p" else (f"{r[c]:.0f}" if isinstance(r[c], (float, np.floating)) else str(r[c]))) + "</td>" for c in cols) + "</tr>" for _, r in df.iterrows()) + "</table>")
    n = len(df) + 1; tbl = s.shapes.add_table(n, len(cols), Inches(0.4), Inches(y), Inches(12.5), Inches(0.3 * n)).table
    for j, h in enumerate(heads):
        c = tbl.cell(0, j); c.text = h; c.fill.solid(); c.fill.fore_color.rgb = RGBColor(0x2b, 0x4a, 0x6f)
        for p in c.text_frame.paragraphs: p.font.size = Pt(size); p.font.bold = True; p.font.color.rgb = FG
    for i, (_, r) in enumerate(df.iterrows(), 1):
        for j, cname in enumerate(cols):
            v = r[cname]; v = f"{v:.3g}" if cname == "p" else (f"{v:.0f}" if isinstance(v, (float, np.floating)) else str(v))
            c = tbl.cell(i, j); c.text = v; c.fill.solid(); c.fill.fore_color.rgb = RGBColor(0x1d, 0x23, 0x2b) if i % 2 else RGBColor(0x17, 0x1b, 0x21)
            for p in c.text_frame.paragraphs: p.font.size = Pt(size); p.font.color.rgb = FG
    if widths:
        for j, w in enumerate(widths): tbl.columns[j].width = Inches(w)

sig = B[B.p <= 0.05]
med = {r: sig[f"hyb_gen_{r}"].median() for r in ("gap", "groves5", "w95")}
medg = {r: sig[f"gen_hyb_{r}"].median() for r in ("gap", "groves5", "w95")}

s = slide("Best simplex per dataset and the hybrid / generalist question")
text(s, ["Five single-cell EMT datasets · all genes and KS EMT genes · pooled and within each unit",
         "For the most significant simplex of every fit: what share of hybrid cells are generalists, and what share of generalists are hybrid, under all three rules.", "",
         "Mouse SCC · Cook · Mouse PDAC · 4T1 clones · MCF10A/HuMEC pseudospace", "", "Thesis project, IISc • 21 October 2026"], size=18)

s = slide("What was asked and how each term is defined", "Definitions are the ones in the dataset scripts (06_groves_specialists.py, 10_describe_labels.py).")
text(s, ["**Simplex fit: PCA, then PCHA on the first k−1 PCs; the permutation test (200 shuffles) gives p. Best fit = lowest p, ties → smallest k.",
         "**Hybrid cells: Early-Hybrid, Hybrid (SCC) or Late-Hybrid states. Cook has no state labels: KS-score quartile bins (middle two = hybrid).",
         "**Three specialist rules (a cell not claimed by a corner is a generalist):",
         "   gap rule: top weight minus second weight ≥ 0.35 (Pal work)",
         "   closest-5% rule: the 5% of cells nearest each corner (Groves 2022)",
         "   weight > 0.95 rule: weight above 0.95 on one corner (Groves STAR Methods)",
         "**Two numbers per fit: % of hybrids that are generalists, and % of generalists that are hybrid (compare with the hybrid share of all cells)."], size=16)

s = slide("Method flow"); pic(s, "schematic.png", y=1.6, h=4.2)
text(s, ["Nothing is refitted for the rules: saved PCA scores and vertices give each cell's barycentric weights."], y=6.3, h=0.8, size=14, color=MUT)

s = slide("How often does a simplex pass the permutation test?"); pic(s, "summary_pass.png")

for ds in ("SCC", "COOK", "PDAC", "T4", "GSE"):
    d = B[B.dataset == ds]; sg = d[d.p <= 0.05]
    s = slide(f"{NAME[ds]}: most significant simplex per unit", f"Heat map from best_fits.csv. {int(len(sg))} of {len(d)} fits have p <= 0.05.")
    pic(s, f"p_heat_{ds}.png", y=1.15, h=6.0)
    s = slide(f"{NAME[ds]}: hybrids and generalists under three rules"); pic(s, f"rules_{ds}.png", y=1.1, h=6.2)
    s = slide(f"{NAME[ds]}: significant fits in detail (p ≤ 0.05)")
    t = sg.sort_values(["p", "k"]).head(14).copy()
    t["genes"] = t.geneset.map({"all": "all", "ks": "KS"})
    cols = ["unit", "genes", "arm", "k", "p", "inside_pct", "hyb_gen_gap", "hyb_gen_groves5", "hyb_gen_w95", "gen_hyb_gap", "gen_hyb_groves5", "gen_hyb_w95", "hybrid_pct"]
    heads = ["unit", "genes", "arm", "k", "p", "inside %", "hyb→gen gap", "hyb→gen 5%", "hyb→gen .95", "gen→hyb gap", "gen→hyb 5%", "gen→hyb .95", "hybrid % of cells"]
    if len(t): table(s, t, cols, heads, size=9, widths=[1.7, .6, 1.2, .4, .8, .8, .85, .85, .85, .85, .85, .85, .9])
    else: text(s, ["No fit reaches p ≤ 0.05."])

GSN = {"all": "A. Entire gene set", "ks": "B. KS EMT genes", "hallmark": "C. Hallmark EMT genes"}
book = sorted((HERE / "figures" / "book").glob("*.png"))
s = slide("Per-analysis figures: how to read them")
text(s, ["**For every dataset and gene set (A entire, B KS, C Hallmark EMT): one page per 3 units (pooled first, then every tumour / line / library / clone).",
         "Each row is one unit's best simplex (lowest p; ties go to the cleaned arm, then smaller k), drawn on PC1 vs PC2. Red title = not significant.",
         "Column 1: cells coloured by EMT state (hybrids = the two middle states). Columns 2 to 4: specialists coloured by corner, generalists grey, under the gap, closest-5% and weight > 0.95 rules.",
         "Under each panel: % of hybrids that are generalists | % of generalists that are hybrid.",
         "For k \u2265 4 the corners are projected onto PC1/PC2 and joined by their convex hull."], size=15)
for ds in ("SCC", "COOK", "PDAC", "T4", "GSE"):
    for gs in ("all", "ks", "hallmark"):
        for f in sorted(x for x in book if x.name.startswith(f"{ds}_{gs}_p")):
            s = slide(""); pic(s, "book/" + f.name, y=0.1, h=7.3, x=0.2, w=12.9)
s = slide("All datasets: hybrids as generalists, significant fits only"); pic(s, "cross_rules.png")
s = slide("Where hybrid cells sit in the best all-genes triangle"); pic(s, "triangles.png")
s = slide("Findings")
text(s, [f"Across {len(sig)} significant fits, the median share of hybrid cells that are generalists is {med['gap']:.0f}% (gap), {med['groves5']:.0f}% (closest 5%) and {med['w95']:.0f}% (weight > 0.95).",
         f"The median share of generalists that are hybrid is {medg['gap']:.0f}%, {medg['groves5']:.0f}% and {medg['w95']:.0f}%, close to the hybrid share of the cells. Hybrids are not over-represented among generalists beyond their base rate.",
         "The weight > 0.95 rule is near 100% in every fit by construction (corners sit inside the cloud, delta = 0); read the gap and closest-5% rules.",
         "The closest-5% rule is always the most generous to generalists because only 5% per corner can be specialists.",
         "Pooled KS-gene fits mostly fail the test; where they pass (Cook) the corners are cell lines."], size=16)
s = slide("Hallmark vs KS EMT genes, pooled fits")
pp = B[B.level == "pooled"]; rows = []
for (ds_, arm_), g_ in pp.groupby(["dataset", "arm"]):
    r_ = {gg: g_[g_.geneset == gg] for gg in ("all", "ks", "hallmark")}
    if all(len(v) for v in r_.values()):
        f_ = lambda v: f"k={int(v.k.iloc[0])}  p={v.p.iloc[0]:.3g}" + ("\u2020" if v.n_null.iloc[0] < 150 else "")
        rows.append(dict(dataset=ds_, arm=arm_, allg=f_(r_["all"]), ks=f_(r_["ks"]), hm=f_(r_["hallmark"])))
table(s, pd.DataFrame(rows), ["dataset", "arm", "allg", "ks", "hm"], ["dataset", "arm", "all genes", "KS genes", "Hallmark EMT"], size=12, widths=[2.2, 2.2, 2.7, 2.7, 2.7])
text(s, ["\u2020 = 100 shuffles. No consistent large difference between Hallmark and KS: both mostly fail pooled (PDAC, SCC all cells) and pass where corners are cell-line / sample identity."], y=6.4, h=0.9, size=13, color=MUT)
s = slide("Caveats")
text(s, ["Within-unit fits are capped at 4,000 cells (random subsample), as in the original analyses.",
         "Best-of-many selection: many fits per dataset, so 'most significant' is optimistic; p has a floor near 0.005 with 200 shuffles.",
         "Hybrid labels come from the same genes (KS score) for PDAC, 4T1, GSE114687 and Cook, so KS-gene results are partly circular. Cook bins are pooled across lines.",
         "Pooled k beyond the ESV elbow is ignored, following the report's rule; the unrestricted best is in best_fits.csv (k_any, p_any).",
         "SCC within-tumour fits were extended to all 12 tumours in this run (new script 08d), with the same method."], size=15)
prs.save(HERE / "Best_simplex_hybrid_generalists.pptx"); print("saved", len(prs.slides), "slides")
(HERE / "Best_simplex_hybrid_generalists.html").write_text("".join(HT).replace("</section>", "", 1) + "</section>")
