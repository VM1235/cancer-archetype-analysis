#!/usr/bin/env python3
"""Plain at-a-glance table (words, thin rules, beige) -> summary/00_overview_table.png. Run from Project_1/ after build_brief.py."""
import sys, io, contextlib, textwrap
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
with contextlib.redirect_stdout(io.StringIO()): import build_brief as BR
BG, INK, MUT, RULE, BROWN = "#F4ECD8", "#3A2E25", "#7A6A5A", "#C9B994", "#8C5A2B"
CORE = BR.CORE
FAMWORD = {"Proliferation": "proliferation", "Epithelial differentiation": "epithelial differentiation", "Mesenchymal / matrix (EMT)": "mesenchymal / matrix", "Metabolism & hypoxia": "metabolism", "Stress & immune signalling": "stress / immune"}
POOLW = {"none": "No simplex", "identity": "Yes: corners are samples", "mixed": "Yes: samples + programmes", "programmes": "Yes: programmes"}

def corners_text(ds):
    W = BR.ANN[(BR.ANN.dataset == ds) & (BR.ANN.level == "within") & (BR.ANN.p <= .05)]
    fits = W.groupby(BR.KEY).fam.apply(set)
    if not len(fits): return "–"
    sh = sorted(((sum(f in s for s in fits) / len(fits), f) for f in CORE), reverse=True)
    words = [FAMWORD[f] for v, f in sh if v >= 0.4][:3]
    return ", ".join(words).capitalize() + (f"\n(only {len(fits)} simplices)" if len(fits) < 5 else "")

def q2_text(ds):
    d, h, o = BR.delta(ds)
    v = "More often generalists" if d >= 8 else ("Slightly more often generalists" if d >= 3 else "No difference")
    return f"{v}\n{h:.0f}% vs {o:.0f}% of other cells"

def build():
    cols = ["Dataset", "Pooled over samples\nA all genes", "B KS genes", "C Hallmark", "Units with a simplex\n(A \u00b7 B \u00b7 C)", "What the corners are\n(within units)", "Hybrid vs other cells\nas generalists", "Bottom line"]
    widths = [1.55, 1.75, 1.75, 1.75, 1.85, 2.25, 1.95, 2.6]
    rows = []
    for ds in BR.DSS:
        pooled = [POOLW[BR.pooled_status(ds, g)[0]] for g in ("all", "ks", "hallmark")]
        w = [BR.within_frac(ds, g) for g in ("all", "ks", "hallmark")]
        rows.append([BR.SHORTNAME[ds]] + pooled + [f"{w[0][0]}/{w[0][1]}  \u00b7  {w[1][0]}/{w[1][1]}  \u00b7  {w[2][0]}/{w[2][1]}", corners_text(ds), q2_text(ds), BR.BOTTOM[ds]])
    W = sum(widths) + 0.6; H = 2.35 + 1.0 * len(rows)
    fig = plt.figure(figsize=(W, H)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis("off")
    ax.text(0.3, H - 0.38, "EMT archetype analysis across five datasets: findings at a glance", fontsize=16, fontweight="bold", color=INK, va="center")
    ax.text(0.3, H - 0.78, textwrap.fill("Within single tumours or cell lines, cells trade off a few shared programmes (proliferation, epithelial differentiation, mesenchymal / matrix); pooled shapes mostly reflect sample identity; hybrid cells are interior mixtures, not a corner.", 190), fontsize=10, color=MUT, va="center")
    xs = [0.3]
    for w_ in widths: xs.append(xs[-1] + w_)
    top = H - 1.35; hh = 0.62
    ax.plot([0.3, xs[-1]], [top, top], color=INK, lw=1.2)
    for i, c in enumerate(cols): ax.text(xs[i] + 0.06, top - hh / 2, c, fontsize=9, fontweight="bold", color=INK, va="center")
    ax.plot([0.3, xs[-1]], [top - hh, top - hh], color=INK, lw=0.8)
    ax.text((xs[1] + xs[4]) / 2, top + 0.13, "Question 1: is there a simplex, and are its corners tasks?", ha="center", fontsize=9, color=BROWN, fontweight="bold")
    ax.text((xs[6] + xs[7]) / 2, top + 0.13, "Question 2", ha="center", fontsize=9, color=BROWN, fontweight="bold")
    rh = 1.0
    for r, row in enumerate(rows):
        y = top - hh - r * rh
        for i, v in enumerate(row):
            wrap = int(widths[i] * 13)
            txt = "\n".join(textwrap.fill(part, wrap) for part in str(v).split("\n"))
            if i == 0: ax.text(xs[i] + 0.06, y - rh / 2 - 0.2, BR.UNITW[BR.DSS[r]], fontsize=8.5, color=MUT, va="center")
            ax.text(xs[i] + 0.06, y - rh / 2 + (0.07 if i == 0 else 0), txt, fontsize=9.3 if i else 10.5, fontweight="bold" if i == 0 else "normal", color=INK if (i == 0 or not str(v).startswith("No simplex")) else MUT, va="center", linespacing=1.3)
        ax.plot([0.3, xs[-1]], [y - rh, y - rh], color=RULE, lw=0.6)
    ax.plot([0.3, xs[-1]], [top - hh - len(rows) * rh] * 2, color=INK, lw=1.2)
    ax.text(0.3, 0.2, "Main fits: Mouse SCC scRNA cells only; cleaned arm elsewhere. Significant = permutation p ≤ 0.05. 'Corners are samples' = most corners are ≥ 70% one sample. Corner families: Hallmark gene sets at each corner, grouped. Q2: purity-gap rule, median over units with a simplex.", fontsize=7.8, color=MUT, va="center")
    fn = HERE / "summary" / "00_overview_table.png"; fig.savefig(fn, dpi=180, facecolor=BG); plt.close(fig); return fn

if __name__ == "__main__":
    print(build())
