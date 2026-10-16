# Single-cell EMT archetype analysis: write-up

Archetype analysis (Pareto task inference, PCHA) applied to five single-cell EMT datasets. Two questions are asked of every dataset: whether corners (task trade-offs) exist, and whether hybrid cells sit at corners (specialists) or in the interior (generalists).

| File | Content |
|---|---|
| `EMT_archetype_analysis_FINAL_REPORT.md` / `.pdf` | Full report: methods, code flow, every figure explained, results per dataset, trends, gaps |
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

The gaps and next steps are in sections 9 and 10 of the report.
