# 4T1 EMT clones (GSE303505)

Five single-cell-derived clones of the mouse 4T1 breast cancer line (E, EM1, EM2, EM3, M) and the parental line, one tumour each, from the Pattabiraman lab. Only the six "Total" libraries are used (12 files, about 180 MB; see `data/raw/DOWNLOAD_LOG.md`). Tumour cells are the cells the authors annotate as `Tumor` (19,930 cells). Clone, tumour and batch are the same variable.

## Pipeline

Run from the repository root.

| Step | Script | Purpose |
|---|---|---|
| 1 | `codes/01_build_anndata.py` | read counts, keep tumour cells, QC, CP10k + log1p |
| 1c | `codes/01c_emt_labels.py` | KS EMT score and state quartiles |
| 1d | `codes/01d_clean_cells.py` | cleaned copy: at least 10,000 UMI, immune and stromal cells removed |
| 2 | `codes/02_pca.py` (tag `main` or `clean`) | PCA, at most 790 fit cells per clone |
| 3 | `codes/run_fits.sh` | PCHA k = 3 to 7, 200 shuffles |
| 4 to 6 | `04_describe_simplex`, `05_robustness`, `06_groves_specialists` | description, refits, specialist rules |
| 8 | `codes/08_within_unit.py [--clean --decontam]` | fits inside one clone |
| 8b, 8c | `08b_describe_corners`, `08c_path_vs_triangle` | corner genes; triangle versus an EMT-score path |
| 13 | `codes/13_mail_figures.py` | slide figures in `results/mail_figures/` |

Cleaned pooled fits need `EXPR_NAME=expr_tumour_clean.h5ad` in the environment.

## Results

See [`RESULTS.md`](RESULTS.md). In brief: no simplex on KS genes; on all genes a weak shape that depends on cleaning (k = 4 before, k = 3 after); corners are cycling cells, epithelial-leaning cells and the M clone; within single clones the shapes disappear after cleaning (1 of 20 fits); hybrid states are 83 to 87% generalists.
