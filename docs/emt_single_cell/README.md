# Single-cell EMT archetype analysis: write-up

Archetype analysis (Pareto task inference, PCHA) applied to five single-cell EMT datasets. Two questions are asked of every dataset: whether corners (task trade-offs) exist, and whether hybrid cells sit at corners (specialists) or in the interior (generalists).

| File | Content |
|---|---|
| `EMT_archetype_analysis_FINAL_REPORT_v2.md` | **Current full report (17 Oct):** methods, code flow, every figure explained, results per dataset, trends, the cleaning audit (Section 9), gaps |
| `EMT_archetype_analysis_FINAL_REPORT.md` / `.pdf` | Version 1 (16 Oct), kept unchanged |
| `EMT_archetype_analysis_FINAL_DECK.pptx` | Slide version of the report (not tracked in git) |
| `Mouse_SCC_Archetype_Pipeline_Guide.pdf` | Short guide to the SCC pipeline |

## Datasets

| Dataset | Accession | Folder |
|---|---|---|
| Mouse skin SCC | GSE288137 | [`Mouse_SCC/`](../../Mouse_SCC/) |
| Cook cell lines, TGF-β time course | GSE147405 | [`Mouse_SCC/Cook/`](../../Mouse_SCC/Cook/) |
| Mouse PDAC | GSE173958 | [`Mouse_SCC/Mouse_PDAC_GSE173958/`](../../Mouse_SCC/Mouse_PDAC_GSE173958/) |
| 4T1 clones | GSE303505 | [`Breast Cancer/Mouse_4T1_GSE303505/`](../../Breast%20Cancer/Mouse_4T1_GSE303505/) |
| MCF10A and HuMEC pseudospace | GSE114687 | [`Breast Cancer/GSE114687_pseudospace/`](../../Breast%20Cancer/GSE114687_pseudospace/) |

## Main findings

- EMT marker genes do not form a simplex beyond sample identity in SCC, Cook, PDAC or 4T1.
- Corners that appear are programmes (proliferation, TGF-β or matrix response, hypoxia, differentiation) and sample or technical axes.
- GSE114687 is the exception: inside single samples, after cleaning, 18 of 24 fits pass with three recurring corners (cycling, differentiated epithelial, EMT / matrix) that follow colony position. The result survives depth regression and depth matching.
- Hybrid cells are mostly generalists in every dataset.
- Results depend on cell cleaning (at least 10,000 UMI), which adds triangles in PDAC and GSE114687 and removes them in 4T1.

The cleaning audit is Section 9 of the v2 report; the gaps and next steps are sections 10 and 11 (sections 9 and 10 of version 1).

## Best simplex, Hallmark genes and hybrid / generalist analysis (`best_simplex/`, 21 to 23 Oct)

Answers two questions across the five datasets: is there a significant simplex (all genes, KS genes, Hallmark EMT genes; pooled and within each unit) and are its corners tasks; and are hybrid cells generalists under three specialist rules. Main fits: SCC scRNA-only pooled + 8 scRNA tumours, cleaned arm elsewhere.

| File | Content |
|---|---|
| `summary/` | One-page overview table + one picture per dataset (sent to supervisor 22 Oct) |
| `REPORT_simplex_explainer.pdf` | Full explainer report: methods, results per dataset, figures |
| `Explainer_simplex_results.pptx` / `.pdf` | Every unit, every gene set (168 slides) |
| `Story_EMT_programmes.pptx` / `.pdf` | Seven-slide story deck |
| `Summary_simplex_results.pptx` / `.pdf`, `Brief_simplex_results.pdf` | Short summaries |
| `INTERPRETATION_and_next_steps.pdf` | Interpretation, link to intrinsic dimension (EMT_ID), next steps |
| `main_fits.csv`, `best_fits.csv`, `all_fits.csv`, `corner_annotations.csv` | Numbers behind every figure |
| `analyze.py`, `corners.py`, `families.py`, `build_*.py` | Scripts (run from `Project_1/`; `analyze.py` first) |
