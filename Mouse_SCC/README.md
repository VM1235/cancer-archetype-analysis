# Mouse skin SCC (GSE288137)

Archetype analysis of 23,959 tumour cells from 12 LKPR mouse skin squamous cell carcinomas (Pérez-González et al., Nat Commun 2026). Eight tumours are scRNA-seq and four are multiome (nuclei). The cells carry four EMT states from the paper (Epithelial, Early-Hybrid, Hybrid, Late-EMT).

This folder also holds two further EMT datasets that reuse the same pipeline: [Cook](Cook/) (human cell lines, TGF-β time course) and [Mouse_PDAC_GSE173958](Mouse_PDAC_GSE173958/). Two more are in `Breast Cancer/`: `Mouse_4T1_GSE303505` and `GSE114687_pseudospace`. The cross-dataset write-up is in [`docs/emt_single_cell/`](../docs/emt_single_cell/).

## Layout

| Path | Content |
|---|---|
| `codes/` | Pipeline scripts 01 to 09 and `run_fits.sh` (see below) |
| `data/` | GEO download and processed matrices (ignored by git; see `data/raw/DOWNLOAD_LOG.md`) |
| `results/` | One folder per fit: `all_main`, `all_scrna`, `all_scrna6`, `ks_main`, `within_tumour_raw`, `cross_projection_SC32_raw`, others |
| `research_report/` | Long-form report, figures and presentation script |
| `hallmark_check/` | MSigDB Hallmark scoring of corner gene signatures |
| `WALKTHROUGH.md` | Concepts, decisions and results from first principles |

## Pipeline

Run from the repository root with `.venv/bin/python`.

| Step | Script | Output |
|---|---|---|
| 1 | `01_build_anndata.py` | log1p CP10k matrix, 14,502 genes |
| 2 | `02_pca.py --geneset all\|ks` | PCA on balanced fit cells (at most 800 per tumour) |
| 3 | `run_fits.sh all main 200` | PCHA for k = 3 to 7 with 200-shuffle permutation test |
| 4 | `04_describe_simplex.py --k 3` | cell weights, label tables, marker genes, figures |
| 5 | `05_robustness.py` | leave-one-tumour-out, new subsample, cell cycle removed |
| 6 | `06_groves_specialists.py` | closest-5% specialist rule |
| 8, 9 | `08_within_tumour.py --raw`, `09_cross_projection.py --raw` | per-tumour fits; projection into one tumour's triangle |

`07_regress_v1.py` removes the technical axis by regression. It also removes EMT signal, so its output is not used.

## Results in brief

- All genes: a triangle at k = 3 (p = 0.005) with corners for late EMT, proliferating epithelium and a technical axis (high in multiome samples and two scRNA tumours).
- KS EMT genes: no simplex (p = 0.83).
- Six clean scRNA tumours without the technical axis: no simplex (k = 3, p = 1.00).
- Single tumours: triangles in SC15 (p = 0.005) and SC32 (p = 0.015), not SC46 (p = 0.244).
- Hybrid cells are 98.8% generalists under the closest-5% rule.

Full numbers, figures and caveats: [`../docs/emt_single_cell/EMT_archetype_analysis_FINAL_REPORT.md`](../docs/emt_single_cell/EMT_archetype_analysis_FINAL_REPORT.md), section 3.
