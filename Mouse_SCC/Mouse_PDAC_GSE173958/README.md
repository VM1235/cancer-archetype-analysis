# Mouse PDAC across sites (GSE173958)

scRNA-seq of mouse pancreatic cancer cells from two mice at several sites (primary tumour, metastases, liver, lung, circulating cells) and two cultured populations: 12 libraries. EMT labels are not provided by GEO. They are derived here from the mouse KS EMT score (Epithelial lowest 40%, Early-Hybrid 40 to 75%, Late-Hybrid 75 to 95%, Mesenchymal top 5%), so a state's share at a corner is partly built in.

## Pipeline

Run from the repository root.

| Step | Script | Purpose |
|---|---|---|
| 1 | `codes/01_build_anndata.py` | read 12 libraries, QC (1,000 to 7,000 genes, mitochondrial fraction at most 15%) |
| 1b | `codes/01b_select_tumour_cells.py` | keep tumour cells by marker clusters (57,770 cells) |
| 1c | `codes/01c_emt_labels.py` | KS score, state bins, k-means clusters |
| 2 | `codes/02_pca.py --assay in_vivo` | PCA, at most 1,000 fit cells per library |
| 3 | `codes/run_fits.sh all\|ks in_vivo 200` | PCHA k = 3 to 7, 200 shuffles |
| 4 to 6 | `04_describe_simplex`, `05_robustness`, `06_groves_specialists` | description, leave-out refits, specialist rules |
| 8 | `codes/08_within_unit.py [--clean]` | fits inside one library, mouse or all in-vivo cells (`run_within*.sh`) |
| 8b, 9, 10 | `08b_describe_corners`, `09_cross_projection`, `10_describe_labels` | corner genes, projection into M1-SS, state tables |

Raw GEO matrices go in `data/raw/` (not tracked).

## Results in brief

- Pooled in-vivo fit: no k with p < 0.01 on either gene set.
- Within libraries, cleaned (at least 10,000 UMI, technical genes dropped): triangles in 6 of 7 libraries at k = 3.
- Corners are an epithelial/ductal programme, proliferation, and either EMT signalling (mouse 1) or hypoxia (mouse 2). Hallmark gene sets agree with these names in 19 of 19 corners.
- EMT runs along an edge between the epithelial and EMT corners. Early-hybrid cells are 64% generalists (gap rule) and 88% (closest-5% rule).

Details: [`../../docs/emt_single_cell/EMT_archetype_analysis_FINAL_REPORT.md`](../../docs/emt_single_cell/EMT_archetype_analysis_FINAL_REPORT.md), section 5.
