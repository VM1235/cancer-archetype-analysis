# Archetype analysis across tumor types

This repository applies Pareto Task Inference (archetypal analysis / PCHA) to bulk RNA-seq, following:

- Groves et al., *Cell Systems* **13**, 690–710 (2022). Figure 1A–C on small-cell lung cancer.  
  [doi:10.1016/j.cels.2022.07.006](https://doi.org/10.1016/j.cels.2022.07.006)
- Hausser et al., *Nature Communications* **10**, 5423 (2019). The underlying method.  
  [doi:10.1038/s41467-019-13195-1](https://doi.org/10.1038/s41467-019-13195-1)

PDFs of both papers are in `papers/`.



Each disease has its own folder with **data**, **codes**, **figures**, and **results**. Shared Python engines live in `src/` at the repo root (PCA, PCHA, enrichment, I/O). Run all commands from this root unless a disease README says otherwise.

| Folder | Status |
|---|---|
| [SCLC Reproduction - Groves Cell Systems 2022](SCLC%20Reproduction%20-%20Groves%20Cell%20Systems%202022/) | Figure 1A–C reproduction |
| [Breast Cancer](Breast%20Cancer/) | DepMap lines + TCGA-BRCA; KS gene-list + METABRIC projection |
| [Glioblastoma](Glioblastoma/) | DepMap GB lines + TCGA-GBM Panel A/B/C |
| [supplementary/emt_hybrid_analysis](supplementary/emt_hybrid_analysis/) | Optional EMT / PN-MES vs archetype weights (not Fig 1A–C) |

## Setup

Python 3.9+ (NumPy must stay **&lt; 2** because `py_pcha` still uses `np.mat`):

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Breast Panel C ComBat also needs R with Bioconductor `sva` (a local copy can live in `Breast Cancer/rlib/`). PAM50 labels used `genefu` + `org.Hs.eg.db`.

## What is and is not in git

**In git (or intended to be):** code, READMEs, figures, small/medium tables, processed SCLC matrices, processed breast Panel A matrix, Model.csv, TCGA clinical annotations.

**Kept locally, gitignored (too large for GitHub):** DepMap `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (~500 MB), UCSC Xena `HiSeqV2` (~172 MB), METABRIC cBioPortal study (`Breast Cancer/brca_metabric/`, ~1.3 GB), Panel C ComBat intermediates, the authors’ cloned `reference/` repo, `.venv/`, and `rlib/`.

Place those files in the paths listed in each disease README before rebuilding from raw data. Existing processed matrices are enough to rerun Panels A/B for SCLC and breast.

**KS + METABRIC (breast):** see [`Breast Cancer/docs/KS_METABRIC_PIPELINE.md`](Breast%20Cancer/docs/KS_METABRIC_PIPELINE.md). Download METABRIC via [`Breast Cancer/brca_metabric/DOWNLOAD.md`](Breast%20Cancer/brca_metabric/DOWNLOAD.md).

## Shared code (`src/`)

| Module | Role |
|---|---|
| `src/archetypes.py` | PCHA (ParTI-style multi-start), ESV, t-ratio |
| `src/pca.py` | PCA, sign alignment, gene-space inverse |
| `src/enrichment.py` | Distance bins, hypergeometric + BH |
| `src/io.py` | Load SCLC Groves matrices and labels |
| `src/preprocess.py` | Breast DepMap filters |
| `src/combat.py` | Python ComBat fallback |
| `src/paths.py` | Folder locations |
