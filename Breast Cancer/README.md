# Breast cancer — Groves method on CCLE + TCGA-BRCA

Same **method** as Groves et al. 2022 Figure 1A–C, **not** SCLC numbers. Cell-line archetypes are fit on DepMap invasive breast carcinoma lines; tumors are TCGA-BRCA. Panel C colors are **IHC ER/HER2**, not PAM50.

Write-up: `docs/how_we_implemented_figures.md`, `docs/Breast_CA_analysis_figures.md`. Parameter logs: `docs/PARAMS_vs_Groves.txt`, `docs/PARAMS_panelC.txt`.

Earlier experimental fits are in `archive/run1/` and are not the figures in `figures/`.

## Folder map

| Path | Contents |
|---|---|
| `data/raw/` | DepMap `Model.csv` + expression CSV |
| `data/processed/` | Panel A matrix (16,500 genes × 63 lines) and PAM50 Entrez table |
| `data/tumors/` | UCSC Xena TCGA-BRCA `HiSeqV2` + clinical matrix |
| `codes/` | Build matrix, PAM50, Panels A/B/C |
| `figures/` | Figure 1A, 1B, 1C analogues |
| `results/` | `panel_a`, `panel_b`, `panel_c_tcga` |
| `docs/` | Notes and parameter tables |
| `archive/run1/` | First pass (not official) |
| `rlib/` | Local R packages for `sva` (gitignored) |

Shared engines: repo-root `src/`.

## Data the codes expect

### Cell lines (Panels A and B)

Place in `data/raw/`:

| File | Source | Notes |
|---|---|---|
| `Model.csv` | DepMap / download.depmap.org | Model metadata |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | DepMap | **Already** log2(TPM+1). ~500 MB; gitignored |

`01_build_input_panelA.py` keeps OncoTree invasive breast / BRCA **cell lines** with RNA-seq, drops genes with max log2(TPM+1) &lt; 1, and writes:

- `data/processed/input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv` (genes × 63 lines)
- `data/processed/input_panelA_models_used.csv`

No ComBat on Panel A (single study).

### PAM50 (Panel B only)

Independent of PCHA. `03_map_and_pam50.R` needs Bioconductor `genefu` and `org.Hs.eg.db`. Labels used by Panel B:

- `results/panel_b/pam50_labels_panelA.csv`

Enrichment **drops LumA and Normal**; the polytope is still fit on all 63 lines.

### Tumors (Panel C)

Place in `data/tumors/`:

| File | Source | Notes |
|---|---|---|
| `HiSeqV2` | UCSC Xena TCGA-BRCA | Already log2; ~172 MB; gitignored |
| `TCGA.BRCA.sampleMap-BRCA_clinicalMatrix` | UCSC Xena | IHC ER/HER2 columns |

Keep primary tumors (`-01` / `sample_type_id=1`). No extra log. ComBat batch = cell line vs tumor, `mod=~1`, `ref.batch=cell_line` (Groves `bc2`).

## How to run

From the **repository root**:

```bash
.venv/bin/python -u "Breast Cancer/codes/01_build_input_panelA.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelA.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelB.py"
.venv/bin/python -u "Breast Cancer/codes/prepare_tcga_brca.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelC_tcga.py" --n-shuffle 20
```

Rebuild PAM50 only if labels are missing:

```bash
Rscript "Breast Cancer/codes/03_map_and_pam50.R" "Breast Cancer"
.venv/bin/python -u "Breast Cancer/codes/04_match_pam50_to_panelA.py"
```

`run_panelC_tcga.py --force-combat` rebuilds the ComBat matrix. Panel C prefers R `sva` in `Breast Cancer/rlib/`; otherwise it uses the Python ComBat fallback in `src/combat.py`.

If the processed Panel A CSV is already present, you can skip `01_build` (and you do not need the 500 MB DepMap expression file). If `results/panel_c_tcga/CCLE_TCGA_COMBAT.csv` is already present, Panel C will reuse it unless `--force-combat` is set.

## Outputs

- `figures/Figure_1A_breast.png`, `Figure_1B_breast.png`, `Figure_1C_tcga.png`
- `results/panel_a/` — PCHA, ESV, t-ratio nulls, PC scores
- `results/panel_b/` — PAM50 labels and enrichment
- `results/panel_c_tcga/` — combined PCA, IHC tables; large CSVs are gitignored
