#!/usr/bin/env python3
"""Interpretation, presentation plan, link to the intrinsic-dimension (EMT_ID) result, and next steps. Beige HTML + PDF.
Run from Project_1/:  .venv/bin/python docs/emt_single_cell/best_simplex/build_interpretation_doc.py"""
import sys, io, contextlib
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
with contextlib.redirect_stdout(io.StringIO()): import build_explainer_ppt as P
from families import family
ROOT = HERE.parents[2]; OUT = HERE / "figures" / "interpretation"; OUT.mkdir(parents=True, exist_ok=True)
BG, INK, MUT = "#F4ECD8", "#3A2E25", "#7A6A5A"
KEY = ["dataset", "level", "unit", "geneset", "arm", "k"]
MAIN, MS = P.MAIN, P.MS
A = P.CA.drop(columns=["p"]).merge(MAIN[KEY + ["p"]], on=KEY); A["fam"] = A.task.map(family)
W = A[(A.level == "within") & (A.p <= .05)]; fits = W.groupby(KEY).fam.apply(set)
E_, M_, Pr = "Epithelial differentiation", "Mesenchymal / matrix (EMT)", "Proliferation"
HYB = {"Early-Hybrid", "Hybrid", "Late-Hybrid"}
nfit = len(fits); two = sum(len({E_, M_, Pr} & f) >= 2 for f in fits); trio = sum({E_, M_, Pr} <= f for f in fits)
epimes = sum({E_, M_} <= f for f in fits)
epi_c = W[W.fam == E_]; mes_c = W[W.fam == M_]; pro_c = W[W.fam == Pr]
epi_epi = int((epi_c.top_state == "Epithelial").sum()); mes_mes = int(mes_c.top_state.isin(["Mesenchymal", "Late-EMT"]).sum()); pro_hyb = int(pro_c.top_state.isin(HYB).sum())
d = MS.hyb_gen_gap - MS.oth_gen_gap
from scipy.stats import wilcoxon
pw = wilcoxon(d).pvalue
dd = MS.assign(d=d).groupby("dataset").d.median().sort_values(ascending=False)
dtxt = ", ".join(f"{P.NAME[k]} {int(round(v)) + 0:+d}" for k, v in dd.items())
q1 = {g: (int((MS.geneset == g).sum()), int((MAIN.geneset == g).sum())) for g in ("all", "ks", "hallmark")}
ident = 100 * (A[(A.level == "pooled") & (A.p <= .05)].top_sample_pct >= 70).mean()
o = pd.read_csv(ROOT / "Mouse_SCC/data/processed/cell_metadata.csv")
ct = pd.crosstab(o.cluster_annotation, o.assay).reindex(["Epithelial", "Early-Hybrid", "Hybrid", "Late-EMT"])
pct = (100 * ct.div(ct.sum(1), axis=0)).round(0)

# figure: SCC state composition by assay, next to the EMT_ID medians
fig, ax = plt.subplots(figsize=(8.5, 3.6)); fig.patch.set_facecolor(BG); ax.set_facecolor("#FBF7EC")
y = np.arange(4)[::-1]
ax.barh(y, pct["multiome"], color="#8C5A2B", label="multiome (nuclei)"); ax.barh(y, pct["scRNA"], left=pct["multiome"], color="#C9B994", label="scRNA (whole cells)")
for yi, st in zip(y, pct.index):
    ax.text(pct.loc[st, "multiome"] / 2, yi, f"{pct.loc[st, 'multiome']:.0f}%", color="white", ha="center", va="center", fontsize=9, fontweight="bold")
    ax.text(103, yi, {"Epithelial": "ID ≈ 29", "Early-Hybrid": "ID ≈ 29", "Hybrid": "ID ≈ 22–23", "Late-EMT": "ID ≈ 22–23"}[st], va="center", fontsize=9.5, color=INK)
ax.set_yticks(y); ax.set_yticklabels(pct.index, color=INK); ax.set_xlim(0, 122); ax.set_xlabel("% of cells of each EMT state", color=INK)
for s in ax.spines.values(): s.set_color("#D8C9A6")
ax.legend(frameon=False, fontsize=8.5, loc="lower center", bbox_to_anchor=(.45, -.42), ncol=2); ax.set_title("Mouse SCC: the high-ID states come mostly from multiome tumours", fontsize=10.5, color=INK)
fig.tight_layout(); f1 = OUT / "scc_state_by_assay.png"; fig.savefig(f1, dpi=170, facecolor=BG); plt.close(fig)

def tbl(rows):
    h = "<table>"
    for i, r in enumerate(rows):
        t = "th" if i == 0 else "td"; h += "<tr>" + "".join(f"<{t}>{c}</{t}>" for c in r) + "</tr>"
    return h + "</table>"

H = []; w = H.append
w("<h1>Interpreting the archetype analysis, linking it to intrinsic dimension, and next steps</h1><p class='m'>Thesis project, IISc · 22 October 2026 · builds on <code>REPORT_simplex_explainer.pdf</code> and the EMT_ID headline result. Numbers are from <code>main_fits.csv</code> and <code>corner_annotations.csv</code>; the EMT_ID medians are quoted from the earlier EMT_ID notes (not recomputed here).</p>")
w("<div class='box'><b>In one paragraph.</b> Inside single tumours and cell lines, cells are spread across a small number of competing programmes: proliferation, epithelial differentiation and a mesenchymal / matrix programme, with a metabolic or stress programme as a frequent fourth. EMT runs along the edge between the epithelial-differentiation and mesenchymal / matrix corners rather than forming corners of its own; hybrid-state cells lie between corners and populate the proliferation corner. Pooled over samples, these shapes are dominated by sample identity. For intrinsic dimension this has two consequences: the EMT_ID state difference in SCC is entangled with the multiome / scRNA split, and position in the simplex (corner vs interior) gives a concrete, testable explanation of why ID might differ between states.</div>")

w("<h2>1. What we found</h2><ol>"
  f"<li><b>Simplices are common but pooled shapes are mostly identity.</b> Significant main fits: all genes {q1['all'][0]}/{q1['all'][1]}, KS genes {q1['ks'][0]}/{q1['ks'][1]}, Hallmark {q1['hallmark'][0]}/{q1['hallmark'][1]}. {ident:.0f}% of the corners of significant pooled simplices are ≥ 70% one sample (Cook: the four cell lines).</li>"
  f"<li><b>Within units, corners are programmes from a few families.</b> {two} of {nfit} significant within-unit simplices contain at least two of proliferation, epithelial differentiation and mesenchymal / matrix; {trio} contain all three. Metabolism (PDAC) and stress / immune signalling (SCC) are the common other corners.</li>"
  f"<li><b>EMT is an edge, not a set of corners.</b> Epithelial-differentiation and mesenchymal / matrix corners occur together in {epimes} of {nfit} simplices. The epithelial corners are epithelial-dominated ({epi_epi} of {len(epi_c)}), the mesenchymal / matrix corners mesenchymal- or late-EMT-dominated ({mes_mes} of {len(mes_c)}).</li>"
  f"<li><b>Hybrid-state cells populate the proliferation corner.</b> In {pro_hyb} of {len(pro_c)} proliferation corners a hybrid state is the largest group of cells. The corner is defined by the programme (cycling / MYC), not by the state: hybrid cells are the proliferating cells, which matches the SCC paper's cycling hybrid population.</li>"
  f"<li><b>Hybrid cells are interior mixtures, slightly more than other cells.</b> 62% of hybrids vs 56% of other cells are generalists (purity-gap rule); hybrids are more often generalists in 40 of {len(MS)} significant fits (Wilcoxon p = {pw:.3f}; fits overlap across gene sets, so indicative). By dataset: {dtxt} points.</li>"
  "<li><b>Gene set matters.</b> The entire gene set gives a simplex most often; the EMT-only sets (KS, Hallmark) give fewer and behave alike: EMT genes alone mostly order cells along one axis.</li>"
  "<li><b>Clearest in clean, continuously varying systems.</b> MCF10A / HuMEC (6 of 6 samples) and PDAC libraries (9 of 9) show the clearest within-unit simplices; 4T1 clones (1 of 6), whose EMT range is narrow, show almost none; Cook lines mainly with Hallmark genes.</li></ol>")

w("<h2>2. How to interpret it</h2>"
  "<p><b>A three-programme trade-off with EMT as one axis.</b> The picture that fits all five datasets is a low-dimensional trade-off between dividing (proliferation), staying differentiated (epithelial programme) and remodelling / migrating (mesenchymal / matrix programme). EMT is the move along the epithelial → mesenchymal / matrix edge; proliferation is a third direction. This reconciles the 'EMT is a continuum' view (cells spread continuously along the edge) with the Pareto view (there are trade-offs, but between programmes rather than EMT states).</p>"
  "<p><b>Hybrid cells as balanced, proliferating mixtures.</b> Hybrid-state cells are the ones in between the epithelial and mesenchymal corners and the ones occupying the proliferation corner. A natural reading: the hybrid state is where cells keep proliferative capacity while partly activating the mesenchymal programme, consistent with reports that late EMT cells divide less.</p>"
  "<p><b>Context-dependent fourth task.</b> The third or fourth corner changes with the system: metabolic programmes in PDAC (tissue site, hypoxia), stress / immune signalling in SCC tumours, TGF-β response in induced EMT. This is expected if the trade-off is set by the environment, not by EMT alone.</p>"
  "<p><b>Why pooled shapes mislead.</b> Samples differ more from each other (patient / mouse, cell line, technology) than cells within a sample differ along EMT, so a pooled simplex finds the samples. Any claim about EMT tasks must be made within a unit, or after explicit batch correction.</p>"
  "<p><b>What not to claim.</b> A significant simplex is not proof of a Pareto trade-off (any non-blob cloud passes); corner names are automatic; most datasets have one sample per condition.</p>")

w("<h2>3. How I would present it</h2><p>Five messages, one figure each, in this order (all figures exist in <code>summary/</code> and the report):</p>")
w(tbl([["#", "Message (slide title)", "Figure", "Why this order"],
       ["1", "Pooled simplices separate samples, not EMT states", "Cook pooled picture (corners = four cell lines)", "Removes the obvious objection first"],
       ["2", "Within a tumour or line, cells trade off three programmes", "MCF10A / HuMEC and PDAC within-unit pictures", "The positive result"],
       ["3", "EMT is the edge between epithelial and mesenchymal / matrix corners", "Within-unit simplex coloured by EMT state", "Connects to the EMT-continuum literature"],
       ["4", "Hybrid cells are interior mixtures and the proliferating cells", "Orange vs grey bar chart + one simplex with hybrids highlighted", "Answers question 2"],
       ["5", "This predicts how intrinsic dimension should vary", "SCC state × assay figure (below) + the test design in Section 4", "Bridges to EMT_ID"]]))
w("<p>Supplement: the overview table, the per-unit grids, the light-QC vs cleaned comparison, the three-rule details. Lead with words and pictures; keep p-values and counts in captions.</p>")

w("<h2>4. Linking to the intrinsic dimension – EMT state correlation</h2>")
w("<p><b>The EMT_ID result.</b> In mouse SCC (GSE288137), local intrinsic dimension (IDEAS <code>local_2nn</code>, k = 600, raw counts) is higher in Epithelial and Early-Hybrid cells (median ≈ 29) than in Hybrid and Late-EMT cells (≈ 22–23).</p>")
w(f"<figure><img src='{f1.relative_to(HERE).as_posix()}'><figcaption>Share of each SCC EMT state coming from multiome (nuclei) vs scRNA (whole-cell) tumours. The two high-ID states are {pct.loc['Epithelial', 'multiome']:.0f}% and {pct.loc['Early-Hybrid', 'multiome']:.0f}% multiome; the two low-ID states are only {pct.loc['Hybrid', 'multiome']:.0f}% and {pct.loc['Late-EMT', 'multiome']:.0f}%.</figcaption></figure>")
w("<h3>4.1 The first thing to rule out: assay</h3>"
  "<p>The high-ID / low-ID split lines up almost exactly with the multiome / scRNA split, and the archetype analysis showed that multiome nuclei form their own technical corner (different RNA content, unspliced long transcripts, lower depth). Local ID computed on raw counts is sensitive to exactly these differences (sparsity, depth). So part or all of the ID–state correlation could be an assay effect.</p>"
  "<p><b>Test.</b> Compare ID between states <i>within the scRNA tumours only</i>, and best within the three mixed scRNA tumours (SC15, SC32, SC46) that contain all four states. Fit ID ~ state + log(depth) + (1 | tumour), separately per assay. If the Epithelial/Early-Hybrid > Hybrid/Late difference survives within scRNA tumours, it is biology; if it collapses, it was the assay.</p>")
w("<h3>4.2 Geometry: corners vs interior</h3>"
  "<p>Local ID estimators depend on where a cell sits in the cloud: a cell at a corner or edge has a one-sided neighbourhood, an interior cell a two-sided one. The archetype fit gives each cell its position (barycentric weights, purity, distance to each corner), so ID can be modelled against position.</p>"
  "<p><b>Test.</b> In each tumour's own simplex: ID ~ purity + state + (1 | tumour); Spearman correlation of ID with purity within tumour, combined across tumours with a sign test; repeat with states permuted within tumour. <b>Reading:</b> if state stops predicting ID once purity is in the model, the state difference is a geometric effect of where states sit; if it persists, state carries information beyond position.</p>")
w("<h3>4.3 Programmes: a mechanistic hypothesis</h3>"
  "<p>The archetype result suggests what ID might be measuring. Cells committed to one programme (the mesenchymal / matrix corner, mostly Late-EMT cells) use a narrow set of genes, so their local neighbourhood should be lower-dimensional; cells mixing programmes (interior generalists, proliferating cells) should have higher local ID. This predicts that ID falls with purity and is lowest near the mesenchymal / matrix corner, which is in line with low ID in Late-EMT. It also predicts that within the same EMT state, generalists have higher ID than specialists.</p>"
  "<p><b>Tests.</b> (i) ID vs distance to the mesenchymal / matrix corner and to the proliferation corner; (ii) ID of generalists vs specialists within each state; (iii) ID vs cell-cycle score, since proliferation corners are hybrid-occupied.</p>")
w("<h3>4.4 Scale: two different 'dimensions'</h3>"
  "<p>The simplex lives in 2 to 3 dimensions (k − 1 principal components, the large-scale shape); local ID (about 22 to 29) measures the small-scale complexity of each neighbourhood in raw gene space, which includes noise. They are complementary, not competing: the simplex says how cells are organised globally, ID how rich their immediate neighbourhood is. Report both; do not expect the numbers to match.</p>")
w("<h3>4.5 Beyond SCC</h3>"
  "<p>The four other datasets have no ID values yet. Computing ID the same way for them would test whether the state pattern generalises, and the archetype results give predictions: the clearest relation should appear in MCF10A / HuMEC and PDAC libraries (clear within-unit simplices), little in 4T1 clones (narrow EMT range); Cook is a test with real time ordering.</p>")

w("<h2>5. Next steps, in priority order</h2>")
w(tbl([["#", "Step", "Question it answers", "Effort", "Value"],
       ["1", "ID within assay and within tumour (SCC)", "Is the EMT_ID state effect biology or the multiome / scRNA split?", "Low (join EMT_ID cell table to metadata)", "Very high"],
       ["2", "ID vs simplex position (purity, corner distances), mixed model", "Does geometry or programme commitment explain ID differences?", "Low to medium (weights already saved)", "Very high"],
       ["3", "Hybrid vs generalist by four states (epithelial → mesenchymal)", "Is the generalist share an inverted U across states, as the mixture picture predicts?", "Low (existing fits)", "High"],
       ["4", "Path vs triangle test within units (SCC, PDAC)", "Is EMT a one-dimensional path inside the simplex?", "Low (script exists for 3 datasets)", "High"],
       ["5", "Batch-corrected pooled fits (Harmony / scANVI)", "Is there a shared simplex across samples once identity is removed?", "Medium", "High"],
       ["6", "Time and space against simplex position (Cook pseudotime, GSE114687 colony position)", "Do cells move from the epithelial to the mesenchymal / matrix corner over time / toward the colony edge?", "Medium", "High"],
       ["7", "ID for the other four datasets", "Does the ID–state pattern generalise?", "Medium (compute time)", "High"],
       ["8", "Corner naming with marker genes and an independent gene-set collection", "Are the corner families robust to the naming method?", "Low", "Medium"],
       ["9", "Planted-triangle power check", "How often would the test find a real but weak simplex?", "Medium", "Medium"],
       ["10", "Multiome replication and knockout projection (from the EMT_ID ideas)", "Independent replication; directional tests of the SCC transcription factors", "Medium to high", "Medium"]]))
w("<p><b>Suggested first milestone:</b> steps 1 to 3 together. They need no new fitting, answer the most consequential question for the ID paper (assay), and give a direct test linking simplex position to ID.</p>")
w("<h2>6. Caveats to carry into any write-up</h2><ul><li>Hybrid labels in four datasets come from the KS score (circular for KS fits); SCC labels are the authors'.</li><li>Many fits, no multiplicity correction; the newest fits used 100 shuffles and k = 3 to 5.</li><li>Corner names are automatic; families are a grouping judgement checked against cell composition (<code>families.py</code>).</li><li>The EMT_ID medians are quoted from the earlier notes; ID was not recomputed here.</li></ul>")

css = """<!doctype html><meta charset=utf-8><title>Interpretation and next steps</title><style>
body{background:#F4ECD8;color:#3A2E25;font:11.5pt/1.55 Georgia,serif;margin:0;-webkit-print-color-adjust:exact;print-color-adjust:exact}
main{max-width:900px;margin:0 auto;padding:28px 22px}h1{font:700 22pt Helvetica,Arial,sans-serif}h2{font:700 16pt Helvetica,Arial,sans-serif;color:#8C5A2B;border-bottom:2px solid #8C5A2B;padding-bottom:3px;margin-top:30px;page-break-after:avoid}
h3{font:700 12.5pt Helvetica,Arial,sans-serif;margin-top:18px;page-break-after:avoid}.m{color:#7A6A5A}code{font:9.5pt Menlo,monospace;background:#EDE2C8;padding:0 3px;border-radius:3px}
.box{background:#FBF7EC;border:1px solid #D8C9A6;border-radius:10px;padding:12px 16px;margin:14px 0}table{border-collapse:collapse;width:100%;font:9.5pt Helvetica,Arial,sans-serif;margin:10px 0;page-break-inside:avoid}
th{background:#8C5A2B;color:#fff;text-align:left;padding:4px 6px}td{border:1px solid #D8C9A6;padding:4px 6px;background:#FBF7EC;vertical-align:top}figure{margin:12px 0;page-break-inside:avoid}figure img{width:100%;border:1px solid #D8C9A6;border-radius:6px}figcaption{font:9.5pt Helvetica,Arial,sans-serif;color:#5a4a3c}
li{margin:5px 0}@page{size:A4;margin:13mm;background:#F4ECD8}</style><main>"""
(HERE / "INTERPRETATION_and_next_steps.html").write_text(css + "".join(H) + "</main>")
print("ok", pro_hyb, len(pro_c), epimes, nfit)
