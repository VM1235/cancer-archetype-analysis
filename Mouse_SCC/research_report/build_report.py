#!/usr/bin/env python3
"""Build Research_Report.pdf from Research_Report.md.

Research_Report.md is the editable source. Steps
  1. Markdown -> HTML (python-markdown + pymdown-extensions arithmatex), with CSS
     and MathJax 3 (loaded from cdn.jsdelivr.net) for equations.
  2. HTML -> PDF with headless Google Chrome.
  3. Page numbers stamped with PyMuPDF (fitz).

Requirements: `markdown` and `pymdown-extensions` importable (set PYTHONPATH if they
are installed in a separate folder), Google Chrome, PyMuPDF, internet access for MathJax.

Usage (from Project_1/):
  PYTHONPATH=<dir with markdown> .venv/bin/python Mouse_SCC/research_report/build_report.py
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MD = HERE / "Research_Report.md"
HTML = HERE / "Research_Report.html"
PDF = HERE / "Research_Report.pdf"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

CSS = r"""
@page { size: A4; margin: 18mm 17mm 20mm 17mm; }
body { font-family: "Charter", "Georgia", "Times New Roman", serif; font-size: 10.6pt;
       line-height: 1.45; color: #1a1a1a; background: #ffffff; max-width: 100%; }
h1 { font-size: 20pt; color: #1f3a5f; border-bottom: 2px solid #1f3a5f; padding-bottom: 4px;
     margin-top: 0; page-break-before: always; font-family: "Helvetica Neue", Arial, sans-serif; }
h2 { font-size: 14pt; color: #1f3a5f; margin-top: 18px; font-family: "Helvetica Neue", Arial, sans-serif; }
h3 { font-size: 12pt; color: #2d5582; font-family: "Helvetica Neue", Arial, sans-serif; }
p, li { text-align: left; }
code { font-family: "SF Mono", Menlo, Consolas, monospace; font-size: 8.8pt; background: #f3f5f8;
       padding: 0 2px; border-radius: 2px; }
pre { background: #f3f5f8; border: 1px solid #dde3ea; border-radius: 4px; padding: 8px 10px;
      font-size: 8.3pt; line-height: 1.35; white-space: pre-wrap; page-break-inside: avoid; }
pre code { background: none; padding: 0; font-size: 8.3pt; }
table { border-collapse: collapse; margin: 10px 0 14px 0; font-size: 8.9pt; width: 100%;
        page-break-inside: auto; }
th, td { border: 1px solid #c9d1db; padding: 3px 6px; vertical-align: top; }
th { background: #e8eef5; text-align: left; }
tr { page-break-inside: avoid; }
img { max-width: 100%; display: block; margin: 10px auto 2px auto; page-break-inside: avoid; }
p.caption { font-size: 8.9pt; color: #333; margin: 2px 6px 16px 6px; line-height: 1.35; }
blockquote { border-left: 4px solid #d98e04; background: #fff8ea; margin: 10px 0; padding: 6px 12px; }
.checkpoint { border: 1px solid #7aa874; background: #f1f8ef; border-radius: 5px; padding: 6px 14px;
              margin: 18px 0; page-break-inside: avoid; font-size: 9.8pt; }
.pagebreak { page-break-after: always; }
.titlepage { padding-top: 60mm; }
.tp-kicker { font-family: "Helvetica Neue", Arial, sans-serif; letter-spacing: 2px; text-transform: uppercase;
             color: #2d5582; font-size: 10pt; }
.tp-title { border: none; font-size: 28pt; page-break-before: avoid; margin: 8px 0 14px 0; }
.tp-sub { font-size: 13pt; color: #333; margin-bottom: 26px; line-height: 1.4; }
.tp-meta { font-size: 10.5pt; color: #333; margin-bottom: 30px; }
.tp-note { font-size: 9.5pt; color: #444; border-top: 1px solid #c9d1db; padding-top: 10px; }
.toc { font-size: 10pt; }
.toc ul { list-style: none; padding-left: 14px; }
.toc > ul { padding-left: 0; }
.toc a { color: #1a1a1a; text-decoration: none; }
a { color: #1f4f8a; text-decoration: none; }
mjx-container[display="true"] { margin: 8px 0 !important; }
"""

MATHJAX = """
<script>
window.MathJax = { tex: { inlineMath: [['\\\\(', '\\\\)']], displayMath: [['\\\\[', '\\\\]']] },
                   svg: { fontCache: 'global' } };
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>
"""


def to_html():
    import markdown
    text = MD.read_text().replace('<div class="checkpoint">', '<div class="checkpoint" markdown="1">')
    body = markdown.markdown(
        text,
        extensions=["tables", "fenced_code", "toc", "attr_list", "md_in_html", "pymdownx.arithmatex"],
        extension_configs={"pymdownx.arithmatex": {"generic": True},
                           "toc": {"toc_depth": "1-2", "title": "Contents"}},
    )
    HTML.write_text(f"<!doctype html><html><head><meta charset='utf-8'><title>Mouse SCC archetype report"
                    f"</title><style>{CSS}</style>{MATHJAX}</head><body>{body}</body></html>")
    print("wrote", HTML.name)


def to_pdf():
    tmp = HERE / "_tmp.pdf"
    cmd = [CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
           "--virtual-time-budget=60000", "--run-all-compositor-stages-before-draw",
           f"--print-to-pdf={tmp}", HTML.as_uri()]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    import fitz
    doc = fitz.open(tmp)
    n = len(doc)
    for i, page in enumerate(doc):
        if i == 0:
            continue
        r = page.rect
        page.insert_text((r.width / 2 - 20, r.height - 22), f"{i + 1} / {n}", fontsize=8,
                         color=(0.35, 0.35, 0.35))
        page.insert_text((48, r.height - 22), "Mouse_SCC archetype analysis: educational report",
                         fontsize=7, color=(0.5, 0.5, 0.5))
    doc.save(PDF)
    doc.close()
    tmp.unlink()
    print("wrote", PDF.name, n, "pages")


if __name__ == "__main__":
    to_html()
    to_pdf()
