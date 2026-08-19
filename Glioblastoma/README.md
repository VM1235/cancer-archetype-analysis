# Glioblastoma — Groves method on CCLE + TCGA-GBM

Same **method** as Groves et al. 2022 Figure 1A–C and the breast analysis, **not** SCLC or breast numbers.

- **Panel A/B:** DepMap glioblastoma cell lines (OncoTree `GB`).
- **Panel C:** TCGA-GBM HiSeq tumors, colored by Verhaak `GeneExp_Subtype`.

PCHA settings match the SCLC paper / ParTI `algNum=5`: `numIter=50`, **150** observed inits, **50** inits per shuffle, `delta=0`, fit in \(k-1\) PCs. Permutation uses **500** shuffles (Groves SCLC used 1000).

## Folder map

| Path | Contents |
|---|---|
| `data/` | DepMap `Model.csv` + expression CSV; Xena `HiSeqV2` + clinical matrix |
| `data/processed/` | Panel A genes × cell-line matrix |
| `codes/` | Build, Verhaak labels, Panels A/B/C |
| `figures/` | Figure 1A–C analogues |
| `results/` | `panel_a`, `panel_b`, `panel_c_tcga` |

Shared engines: repo-root `src/`.

## Data the codes expect

### Cell lines (A and B)

In `data/`:

| File | Notes |
|---|---|
| `Model.csv` | DepMap metadata |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | Already log2(TPM+1); gitignored |

`01_build_input_panelA.py` keeps OncoTree GB **cell lines** with RNA-seq, drops genes with max log2(TPM+1) &lt; 1.

### Verhaak labels (B only)

Independent of PCHA. `02_assign_verhaak.py` scores each DepMap line on Classical / Mesenchymal / Proneural marker genes (gene-wise z-scores on the cell-line matrix). Neural is not used. The polytope is fit on all lines.

### Tumors (C)

In `data/`:

| File | Notes |
|---|---|
| `HiSeqV2` | UCSC Xena TCGA-GBM, already log2 |
| `TCGA.GBM.sampleMap-GBM_clinicalMatrix` | `GeneExp_Subtype`, G-CIMP |

Keep primary (`-01` / `sample_type_id=1`). ComBat batch = cell line vs tumor, `ref.batch=cell_line`.

## How to run

From the **repository root**:

```bash
.venv/bin/python -u "Glioblastoma/codes/01_build_input_panelA.py"
.venv/bin/python -u "Glioblastoma/codes/02_assign_verhaak.py"
.venv/bin/python -u "Glioblastoma/codes/run_panelA.py"
.venv/bin/python -u "Glioblastoma/codes/run_panelB.py"
.venv/bin/python -u "Glioblastoma/codes/prepare_tcga_gbm.py"
.venv/bin/python -u "Glioblastoma/codes/run_panelC_tcga.py" --n-shuffle 20
```

Panel A is slow (150 observed inits and 50 inits × 500 shuffles for each k). It checkpoints nulls under `results/panel_a/` and can be restarted.

## Outputs

- `figures/Figure_1A_gbm.png`, `Figure_1B_gbm.png`, `Figure_1C_gbm.png`
- `results/panel_a/` — PCHA, ESV, t-ratio
- `results/panel_b/` — Verhaak labels and enrichment
- `results/panel_c_tcga/` — combined PCA; large CSVs gitignored

## Hypothesis-driven k=2 (PN–MES axis)

Separate from the Groves multi-archetype run above. Motivated by Jolly lab
iScience 2024: Verhaak’s four subtypes are not four exclusive states; the
antagonistic pair is Proneural–Mesenchymal. **Do not** overwrite
`results/panel_a/`, `panel_b/`, `panel_c_tcga/`, or `Figure_1A_gbm.png` / `1B` / `1C`.

Uses the **same** 12-PC cell-line scores and (by default) the existing ComBat
matrix. k=2 is fit on the first 1 PC as the two PC1 endpoints (`delta=0`
max-volume 1-simplex). **py_pcha is not used for k=2** (it can hang in 1-D);
k=3–7 in `panel_a/` still use PCHA. This is **not** Groves “smallest
significant k”.

```bash
.venv/bin/python -u "Glioblastoma/codes/run_panelA_k2.py"
.venv/bin/python -u "Glioblastoma/codes/run_panelB_k2.py"
.venv/bin/python -u "Glioblastoma/codes/run_panelC_k2.py"
```

- `figures/Figure_1A_gbm_k2.png`, `Figure_1B_gbm_k2.png`, `Figure_1C_gbm_k2.png`
- `results/panel_a_k2/`, `panel_b_k2/`, `panel_c_k2/`

### k=2 in the full 12-PC space (non-degenerate t-ratio)

The (k−1)-PC ParTI rule at k=2 is a 1-D fit; t-ratio is then identically 1.
`run_panelA_k2_12.py` instead fits two archetypes in the **same 12 PCs** as
Panel A, and sets t = pair length / data extent along that axis. Same 150/50
inits and 500 shuffles. This is a documented relaxation, not Groves k.

```bash
.venv/bin/python -u "Glioblastoma/codes/run_panelA_k2_12.py"
```

- `figures/Figure_1A_gbm_k2_12.png`
- `results/panel_a_k2_12/`
