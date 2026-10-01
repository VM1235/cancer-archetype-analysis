# Hausser Fig. 1D Reproduction — Pan-cancer per-type archetypes

Reproduction attempt of **Figure 1d** from:

> Hausser, J., Szekely, P., Bar, N., Zimmer, A., Sheftel, H., Caldas, C. & Alon, U.
> "Tumor diversity and the trade-off between universal cancer tasks."
> *Nature Communications* **10**, 5423 (2019).
> https://www.nature.com/articles/s41467-019-13195-1

Fig. 1d shows tumor transcriptomes of 8 cancer types, each falling on its own
low-dimensional polyhedron (archetype simplex) in gene-expression PC space.

## Status

- **Numerics validated on synthetic data.** The full pipeline (prep → PCA →
  PCHA fit at k=3,4,5 → t-ratio permutation test → figure export) runs
  end-to-end and was smoke-tested with a synthetic 4-archetype dataset,
  which correctly returned k=3,4 significant and k=5 not significant.
- **No real TCGA/Xena/METABRIC data has been fitted yet** — this environment
  has no network access to xenabrowser.net or cbioportal.org. You'll need to
  download the data yourself (instructions below) and run the pipeline
  locally.
- **Methodology cross-checked against the paper's Methods section** (not
  just the main text) on 2026-08-24. Two real bugs were found and fixed in
  this pass — see "Corrections made" below — and one assumption remains
  genuinely unresolved (LUAD vs LUSC).

## Pipeline

```
codes/00_cancer_types.py            registry of the 8 Fig. 1d cancer types
codes/_registry.py                  helper to import 00_ by file path
codes/01_prepare_tcga_tumor.py      Xena HiSeqV2 + clinicalMatrix -> primary-tumor CSV (per TCGA type)
codes/01b_prepare_metabric_full.py  cBioPortal METABRIC export -> primary-tumor CSV (breast, full gene set)
codes/02_run_archetypes_per_cancer.py  PCA -> PCHA (k=3,4,5) -> t-ratio permutation test, one cancer type
codes/03_run_all_cancer_types.py    loop 02_ over all 8 types, aggregate summary table
codes/04_export_figure_1d.py        render the 8-panel figure
```

Run from the repository root (i.e. one level above this folder), e.g.:

```bash
.venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/01_prepare_tcga_tumor.py" --cancer-type THCA
.venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/01b_prepare_metabric_full.py"
.venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/03_run_all_cancer_types.py" --n-perm 1000
.venv/bin/python -u "Hausser Fig1D Reproduction - Pan-cancer per-type archetypes/codes/04_export_figure_1d.py"
```

## Data download instructions

### TCGA types (via UCSC Xena "TCGA Hub")

For each of the 6 non-breast cancer types below: go to
https://xenabrowser.net/datapages/ → **TCGA Hub** → select the cohort →
download two files:

1. **`IlluminaHiSeq RNASeqV2 (unc.edu, gene RSEM log2 -- HiSeqV2)`** → save as
   `data/<CODE>/HiSeqV2`
2. **`Phenotypes -> Curated survival data`, or the cohort's
   `<cohort>_clinicalMatrix`** (must contain a `sampleID` column and, ideally,
   a `sample_type` column with values like `Primary Tumor` / `Solid Tissue
   Normal`) → save anywhere under `data/<CODE>/` — `01_prepare_tcga_tumor.py`
   recursively searches for a filename containing `clinicalMatrix`.

| Code | Xena cohort (confirm exact name on the portal) | Fig. 1d label |
|---|---|---|
| THCA | TCGA Thyroid Cancer (THCA) | Thyroid |
| BLCA | TCGA Bladder Cancer (BLCA) | Bladder |
| LIHC | TCGA Liver Cancer (LIHC) | Liver |
| COAD | TCGA Colon Cancer (COAD) | Colon |
| LGG | TCGA Lower Grade Glioma (LGG) | Glioma |
| LUAD (or LUSC — see below) | TCGA Lung Adenocarcinoma (LUAD) | Lung |
| HNSC | TCGA Head and Neck Cancer (HNSC) | Head & Neck |

Xena periodically renames cohorts; if a name above 404s, search the datapages
site for the TCGA disease code instead (e.g. "THCA").

### Breast (METABRIC, not TCGA)

Fig. 1d's caption is explicit: *"TCGA, breast cancer from Metabric."* Download
the modern cBioPortal export (`brca_metabric`) from
https://www.cbioportal.org/study/summary?id=brca_metabric → **Download** tab,
or via the cBioPortal API/GitHub data repo, and place under
`Breast Cancer/brca_metabric/` (reusing the same folder other analyses in
this repo already use):

- `data_mrna_illumina_microarray.txt`
- `data_clinical_sample.txt`
- `data_clinical_patient.txt`

`01b_prepare_metabric_full.py` reads these directly — no further conversion
needed. If you've already followed `Breast Cancer/brca_metabric/DOWNLOAD.md`
for another analysis in this repo, nothing new needs downloading.

Note: the original paper (2019) used the legacy monolithic cBioPortal export
(`data_expression.txt`, `data_CNA.txt`, `data_clinical.txt` from
`brca_metabric.tar.gz`). cBioPortal has since split the clinical file into
`data_clinical_sample.txt` / `data_clinical_patient.txt`; the prep script
already re-joins them, so this is a format change, not a data change.

## Corrections made in this pass (2026-08-24)

The pipeline's numerics were previously validated only on synthetic data, so
several assumptions in `02_run_archetypes_per_cancer.py` had not yet been
checked against the paper's actual Methods section. Checked now:

1. **Gene filtering before PCA — removed.** The paper's Methods describe PCA
   run on the full mean-centered gene matrix (20,530 TCGA genes, or the full
   METABRIC gene set) with **no variance-based gene filter**. The script
   previously filtered to the top 5,000 variance genes by default; that was
   an invented step, not something the paper does. `--n-genes` now defaults
   to `None` (= all genes) and only applies a filter if you explicitly pass
   a value.
2. **k-selection threshold — fixed from p<0.05 to p<0.01.** Methods, verbatim:
   *"We chose the smallest number of archetypes that produced a statistically
   significant polyhedron (p < 0.01)."* The script's `best_k` selection
   previously used p<0.05.
3. **Primary-tumor filter — now prefers the literal field.** Methods: primary
   tumors are those with *"field 'sample_type' set to 'Primary Tumor' in the
   TCGA clinical annotation."* `01_prepare_tcga_tumor.py` now checks for a
   `sample_type` column first and only falls back to the TCGA-barcode-suffix
   heuristic (`01`) if that column isn't present in the downloaded
   clinicalMatrix.
4. **k range (3, 4, 5), no k=6+ — confirmed unchanged.** Methods state this
   explicitly ("we did not attempt to find six or more archetypes because of
   the limited number of tumor samples"), matching the script already.
5. **Permutation count — left at 100 by default, but flagged.** Methods use
   1000 shuffles; the script's 100-shuffle default is fine for a fast dev
   pass but `03_run_all_cancer_types.py` now prints a warning if `--n-perm`
   is below 1000, and the run commands above use `--n-perm 1000`.

### Still open

- **LUAD vs LUSC for "Lung."** The paper's main-text Methods list LUAD and
  LUSC as two separate entries among the 15 cancer types tested (not a
  pooled "Lung" category), which confirms Fig. 1d's "Lung" panel is one
  specific one of the two — but the main text never says which. This can
  only be resolved from the Supplementary Methods PDF, which isn't
  fetchable from this environment. `00_cancer_types.py` defaults to LUAD;
  override with `--cancer-type LUSC` if the SI says otherwise once you can
  check it.
- **PCA scaling in `src/pca.py`.** The paper explicitly does *not* scale by
  standard deviation before PCA (only mean-centers). This review didn't have
  access to `src/pca.py`'s implementation — confirm `fit_pca` doesn't
  standardize by default before trusting final numbers.

## Interpreting your output vs. the paper

`t_ratio_summary.csv` / the figure panels report both your fitted p-value and
Hausser's reported p-value per type, for direct comparison:

| Cancer type | Hausser's reported p |
|---|---|
| Thyroid | p < 0.001 |
| Bladder | p = 0.001 |
| Liver | p = 0.002 |
| Colon | p = 0.009 |
| Glioma | p < 0.001 |
| Breast | p = 0.001 |
| Lung | borderline, p = 0.013 |
| Head & Neck | borderline, p = 0.022 |

All six non-borderline types are significant at FDR < 10%; lung and head &
neck are explicitly called out in the paper as borderline. Don't expect exact
p-value matches — PCHA has random restarts and the permutation test is
stochastic — but the same types should come out significant/borderline in
the same pattern.
