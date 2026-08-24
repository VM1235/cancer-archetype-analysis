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
| [Breast Cancer](Breast%20Cancer/) | Same method on DepMap cell lines + TCGA-BRCA tumors |
| [Glioblastoma](Glioblastoma/) | Placeholder (not started) |

## Setup

Python 3.9+ (NumPy must stay **&lt; 2** because `py_pcha` still uses `np.mat`):

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Breast Panel C ComBat also needs R with Bioconductor `sva` (a local copy can live in `Breast Cancer/rlib/`). PAM50 labels used `genefu` + `org.Hs.eg.db`.

## What is and is not in git

**In git (or intended to be):** code, READMEs, figures, small/medium tables, processed SCLC matrices, processed breast Panel A matrix, Model.csv, TCGA clinical annotations.

**Kept locally, gitignored (too large for GitHub):** DepMap `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (~500 MB), UCSC Xena `HiSeqV2` (~172 MB), Panel C ComBat intermediates (~100–270 MB), the authors’ cloned `reference/` repo, `.venv/`, and `rlib/`.

Place those files in the paths listed in each disease README before rebuilding from raw data. Existing processed matrices are enough to rerun Panels A/B for SCLC and breast.

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
