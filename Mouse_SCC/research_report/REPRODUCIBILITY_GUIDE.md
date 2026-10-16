# Reproducibility guide: Mouse_SCC archetype analysis

How to reproduce (1) the analysis results and (2) this report's figures and PDF.
All commands run from `Project_1/` unless stated. Expected numbers let you check a rerun.

## 0. Requirements

| Item | Version used | Notes |
|---|---|---|
| Machine | macOS, Apple M2, 8 cores, 16 GB RAM | Version A PCA needs ~4 GB peak |
| Python | 3.9, in `Project_1/.venv` | |
| Packages | numpy 1.26.4, scipy 1.13.1, scikit-learn 1.6.1, pandas 2.3.3, matplotlib 3.9.4, anndata 0.10.9, scanpy 1.10.3, py-pcha 0.1.3 | `py_pcha` needs numpy < 2 (uses `np.mat`) |
| Repository modules | `src/archetypes.py`, `src/enrichment.py` | imported via `sys.path` |
| Gene list | `Breast Cancer/data/processed/ks_cellline_signature_tan2014.csv` | version C only |
| Report build only | `markdown`, `pymdown-extensions`, Google Chrome, PyMuPDF (`fitz`), internet (MathJax CDN) | see §3 |

## 1. Data

1. Download GEO GSE288137 supplementary file `GSE288137_SCO.rds.gz` (sha256 `f36350196cf107505999aabb2e3e3b82ef234cda6fd04b6e66ac63957d50b3d7`; see `data/raw/DOWNLOAD_LOG.md`). It is gzip-compressed twice; unpack to `GSE288137_SCO_unpacked.rds`.
2. **Export from R (script not saved in the project; reconstruct as below).** In R with Seurat ≥ 5: load the object; write `GetAssayData(obj, assay = "RNA", layer = "counts")` as MatrixMarket to `data/processed/counts.mtx` (genes × cells); gene names to `genes.txt`; cell barcodes to `barcodes.txt`; `obj@meta.data` columns `sample, assay, cluster_annotation, clusterU, phase, nCount_RNA, nFeature_RNA, pct_counts_mt, pct_counts_ribo, LKPR_*_andrea` with a `barcode` column to `cell_metadata.csv` (same order as `barcodes.txt`). Expected: 20,358 genes × 23,959 cells, integer counts.
3. Build the analysis matrix:

```bash
.venv/bin/python Mouse_SCC/codes/01_build_anndata.py
```

Expected: `expr_log_filtered.h5ad` with 23,959 cells × 14,502 genes; printed tumour × state table equal to Table in Report §4.4.

## 2. Analysis

```bash
# Step 2: PCA (A, B, C)
.venv/bin/python -u Mouse_SCC/codes/02_pca.py --geneset all                 # A -> results/all_main
.venv/bin/python -u Mouse_SCC/codes/02_pca.py --geneset all --assay scRNA   # B -> results/all_scrna
.venv/bin/python -u Mouse_SCC/codes/02_pca.py --geneset ks                  # C -> results/ks_main
.venv/bin/python -u Mouse_SCC/codes/02b_tumour_within_state.py --geneset all --tag main
.venv/bin/python -u Mouse_SCC/codes/02b_tumour_within_state.py --geneset all --tag scrna
.venv/bin/python -u Mouse_SCC/codes/02b_tumour_within_state.py --geneset ks  --tag main

# Step 3: fits + 200-shuffle permutation tests (2-3 h each; resumable)
bash Mouse_SCC/codes/run_fits.sh all main 200
bash Mouse_SCC/codes/run_fits.sh all scrna 200
bash Mouse_SCC/codes/run_fits.sh ks main 200

# Step 4: describe simplices
.venv/bin/python -u Mouse_SCC/codes/04_describe_simplex.py --geneset all --tag main  --k 3
.venv/bin/python -u Mouse_SCC/codes/04_describe_simplex.py --geneset all --tag main  --k 4
.venv/bin/python -u Mouse_SCC/codes/04_describe_simplex.py --geneset all --tag scrna --k 3
.venv/bin/python -u Mouse_SCC/codes/04_describe_simplex.py --geneset all --tag scrna --k 4
.venv/bin/python -u Mouse_SCC/codes/04_describe_simplex.py --geneset ks  --tag main  --k 3

# Step 5: robustness (sequential; ~1-2 min per refit)
.venv/bin/python -u Mouse_SCC/codes/05_robustness.py --geneset all --tag scrna --k 3 --variant loto
.venv/bin/python -u Mouse_SCC/codes/05_robustness.py --geneset all --tag scrna --k 3 --variant nophase
.venv/bin/python -u Mouse_SCC/codes/05_robustness.py --geneset all --tag scrna --k 3 --variant seed
.venv/bin/python -u Mouse_SCC/codes/05_robustness.py --geneset all --tag main  --k 3 --variant loto

# Step 6: Groves specialists
.venv/bin/python -u Mouse_SCC/codes/06_groves_specialists.py --geneset all --tag scrna --k 3
```

**Note on version A's PCA.** The saved `all_main/pca.npz` was computed with `--pca-method svd` (the default at the time). The current default is `gram`; both give identical PCs to ~10⁻¹³. To reproduce the file bit-for-bit, add `--pca-method svd` (~13 min).

**Note on resuming.** `05_robustness.py` appends to `robustness/k3_<variant>.csv` and skips variants already present. Delete the CSV to rerun from scratch. `03_fit_and_test.py` resumes from `null_k*.npy`; delete those files to restart the null.

### Expected key numbers

| Output | Expected |
|---|---|
| `all_main/pca_info.json` | 9,600 fit cells; 0.1% rule 49 PCs |
| `all_scrna/pca_info.json` | 6,400 fit cells; 56 PCs |
| `ks_main/pca_info.json` | 176 genes; 108 PCs |
| `all_main/fit_summary.csv` | k=3 t 0.5284, p 0.005; k=4 p 0.1841 |
| `all_scrna/fit_summary.csv` | k=3 t 0.5290, p 0.005 |
| `ks_main/fit_summary.csv` | k=3 t 0.4806, p 0.8308 |
| ESV elbow | 4, 4, 5 |
| `all_scrna/robustness/k3_loto.csv` | agreement 89.7–98.9% |
| `all_scrna/groves5_k3/summary.md` | 703 specialists per corner; 0 overlaps |

Randomness is seeded (`SEED = 0`; `np.random.seed` before each PCHA call; `default_rng(SEED + 1000 + k)` for shuffles), so reruns on the same software should reproduce these numbers exactly. Different BLAS libraries can change floating-point details and therefore, rarely, a PCHA start's outcome.

## 2b. Follow-up analyses (12 Oct 2026; Report Addendum)

```bash
# 1. remove the V1 score, redo version B (PCA ~2 min; fits + 200 shuffles ~1.5 h)
.venv/bin/python -u Mouse_SCC/codes/07_regress_v1.py
bash Mouse_SCC/codes/run_fits.sh all scrna_noV1 200
.venv/bin/python -u Mouse_SCC/codes/04_describe_simplex.py --geneset all --tag scrna_noV1 --k 3
.venv/bin/python -u Mouse_SCC/codes/04_describe_simplex.py --geneset all --tag scrna_noV1 --k 4

# 2. within-tumour fits (each loads the full matrix; run <= 3 at a time on 16 GB)
for T in SC15 SC32 SC46; do .venv/bin/python -u Mouse_SCC/codes/08_within_tumour.py --fit --tumour $T --k 3; done
for T in SC15 SC32 SC46; do .venv/bin/python -u Mouse_SCC/codes/08_within_tumour.py --fit --tumour $T --k 4; done
.venv/bin/python -u Mouse_SCC/codes/08_within_tumour.py --compare

# 3. cross-projection into SC32 (needs step 2's --compare output)
.venv/bin/python -u Mouse_SCC/codes/09_cross_projection.py --host SC32

# slide figures (reads results only)
.venv/bin/python Mouse_SCC/research_report/make_slide_figures.py
```

Expected: `all_scrna_noV1/fit_summary.csv` k=3 t 0.4739, p 0.0398; within-tumour k=3 p: SC15 0.8209, SC32 0.7960, SC46 0.0050.

## 3. Report figures and PDF

```bash
# Figures (reads results; writes only into Mouse_SCC/research_report/figures/)
.venv/bin/python Mouse_SCC/research_report/make_report_figures.py

# Markdown tools, installed outside the project environment
.venv/bin/python -m pip install --target /tmp/mdlib markdown pymdown-extensions

# HTML + PDF (needs Google Chrome and internet for MathJax)
PYTHONPATH=/tmp/mdlib .venv/bin/python Mouse_SCC/research_report/build_report.py
```

`make_report_figures.py` loads the full h5ad for one figure (`emp_v1_longgene.png`; ~1–2 min). `build_report.py` converts `Research_Report.md` → `Research_Report.html` → `Research_Report.pdf` and stamps page numbers. If Chrome is elsewhere, edit `CHROME` in `build_report.py`.

## 4. Not reproducible from the project alone

- The R export commands (reconstructed in §1, not the original).
- Report-time checks run interactively (recorded in `REPORT_VALIDATION.md` with their results); the most important (V1 score) is reproduced by `make_report_figures.py::emp_v1_longgene`.
