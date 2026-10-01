# Breast cancer — Groves method on CCLE + TCGA-BRCA (+ METABRIC KS)

Same **method** as Groves et al. 2022 Figure 1A–C, **not** SCLC numbers. Cell-line archetypes are fit on DepMap invasive breast carcinoma lines; tumors are TCGA-BRCA (IHC ER/HER2) and/or METABRIC (Claudin subtype).

| Analysis track | Panel A genes | Tumor cohort | Panel C colors |
|---|---|---|---|
| **Default** | ~16,500 (DepMap filtered) | TCGA-BRCA | IHC ER/HER2 |
| **KS genelist** | 206 (Tan 2014 KS signature) | TCGA or METABRIC | IHC or Claudin subtype |

Write-up: `docs/how_we_implemented_figures.md`, `docs/Breast_CA_analysis_figures.md`.  
KS → METABRIC: **`docs/KS_METABRIC_PIPELINE.md`**.  
Script list: **`docs/SCRIPTS.md`**. Results map: **`results/README.md`**.

Earlier experimental fits are in `archive/run1/` and are not the figures in `figures/`.

## Folder map

| Path | Contents |
|---|---|
| `data/raw/` | DepMap `Model.csv` + expression CSV |
| `data/processed/` | Panel A matrix; KS gene lists + `input_panelA_ks_genelist.csv` |
| `data/tumors/` | UCSC Xena TCGA-BRCA `HiSeqV2` + clinical matrix |
| `data/gse161529/` | Pal et al. GEO/Figshare (**raw gitignored**) |
| `brca_metabric/` | cBioPortal METABRIC export (**gitignored**; see `DOWNLOAD.md`) |
| `codes/` | Panels A/B/C, sc projection, Pal reverse fit |
| `figures/` | Figure 1A–C analogues + Fig. 4-style sc / Pal reverse |
| `results/` | See `results/README.md` |
| `docs/` | Notes, KS/METABRIC pipeline, `SCRIPTS.md` |
| `archive/run1/` | First pass (not official) |
| `rlib/` | Local R packages for `sva` (gitignored) |

Shared engines: repo-root `src/`.

## Default pipeline (16k genes, TCGA)

From the **repository root**:

```bash
.venv/bin/python -u "Breast Cancer/codes/01_build_input_panelA.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelA.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelB.py"
.venv/bin/python -u "Breast Cancer/codes/prepare_tcga_brca.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelC_tcga.py" --n-shuffle 20
```

Figures: `Figure_1A_breast.png`, `Figure_1B_breast.png`, `Figure_1C_tcga.png`.

## KS genelist + METABRIC pipeline

See **`docs/KS_METABRIC_PIPELINE.md`** for the full run order. Short version:

```bash
.venv/bin/python -u "Breast Cancer/codes/00_build_ks_genelist_input.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelA_ks_genelist.py"              # k=3
.venv/bin/python -u "Breast Cancer/codes/run_panelA_ks_genelist_extendedk.py --only-k 4"  # k=4
# Download METABRIC → brca_metabric/ (see brca_metabric/DOWNLOAD.md)
.venv/bin/python -u "Breast Cancer/codes/prepare_metabric_brca.py"
.venv/bin/python -u "Breast Cancer/codes/run_panelC_metabric_ks_genelist.py" --n-shuffle 20
.venv/bin/python -u "Breast Cancer/codes/run_panelC_metabric_ks_genelist_k4.py" --n-shuffle 20
```

Figures: `Figure_1C_metabric_ks_genelist.png`, `Figure_1C_metabric_ks_genelist_k4.png`.

## Data the codes expect

### Cell lines (Panels A and B)

Place in `data/raw/`:

| File | Source | Notes |
|---|---|---|
| `Model.csv` | DepMap | Model metadata |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | DepMap | **Already** log2(TPM+1). ~500 MB; gitignored |

`01_build_input_panelA.py` → `data/processed/input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv` (genes × 63 lines).

KS track uses `data/processed/input_panelA_ks_genelist.csv` (206 genes × 63 lines; in git).

### PAM50 (Panel B only)

`03_map_and_pam50.R` + `04_match_pam50_to_panelA.py` → `results/panel_b/pam50_labels_panelA.csv`.  
Enrichment **drops LumA and Normal**; the polytope is still fit on all 63 lines.

### Tumors — TCGA (Panel C)

Place in `data/tumors/`:

| File | Source | Notes |
|---|---|---|
| `HiSeqV2` | UCSC Xena TCGA-BRCA | Already log2; ~172 MB; gitignored |
| `TCGA.BRCA.sampleMap-BRCA_clinicalMatrix` | UCSC Xena | IHC ER/HER2 columns |

### Tumors — METABRIC (KS Panel C)

Download cBioPortal study into `brca_metabric/` — **`brca_metabric/DOWNLOAD.md`**.  
`prepare_metabric_brca.py` writes gene-matched tables under `results/panel_c_metabric_ks_genelist/`.

ComBat: batch = cell line vs tumor, `mod=~1`, `ref.batch=cell_line` (Groves `bc2`). Large ComBat CSVs are gitignored; Panel C reuses an existing local file or rebuilds with `--force-combat`.

## Outputs (official figures)

| Track | Panel A | Panel B | Panel C |
|---|---|---|---|
| Default | `Figure_1A_breast.png` | `Figure_1B_breast.png` | `Figure_1C_tcga.png` |
| KS / TCGA | `Figure_1A_breast_ks_genelist.png` | `Figure_1B_breast_ks_genelist.png` | `Figure_1C_tcga_ks_genelist.png` |
| KS / METABRIC k=3 | extendedk Fig 1A optional | — | `Figure_1C_metabric_ks_genelist.png` |
| KS / METABRIC k=4 | `Figure_1A_breast_ks_genelist_extendedk.png` | `Figure_1B_breast_ks_genelist_k4.png` | `Figure_1C_metabric_ks_genelist_k4.png` |

`run_panelC_tcga.py --force-combat` rebuilds the TCGA ComBat matrix. Panel C prefers R `sva` in `Breast Cancer/rlib/`.

## scRNA → DepMap KS triangle, then Pal reverse

Wu and GSE173634 counts: repo-root `Sahoo/` (gitignored dumps). Pal counts: `data/gse161529/`. Full command list: `docs/SCRIPTS.md`.

```bash
# Pal → DepMap k=3 (after 00_extract + 00_build)
.venv/bin/python -u "Breast Cancer/codes/1run_panelA_ks_project_gse161529_sc.py"

# Fit Pal simplex + 50-shuffle t-ratio
.venv/bin/python -u "Breast Cancer/codes/17.9_run_panelA_ks_genelist_gse161529.py"
.venv/bin/python -u "Breast Cancer/codes/20_t_ratio_gse161529.py"

# Project Wu, GSE173634, DepMap lines into Pal k=3
.venv/bin/python -u "Breast Cancer/codes/18_project_gse176078_onto_gse161529.py"
.venv/bin/python -u "Breast Cancer/codes/19_project_gse173634_onto_gse161529.py"
.venv/bin/python -u "Breast Cancer/codes/22_project_depmap_onto_gse161529.py"
.venv/bin/python -u "Breast Cancer/codes/21_panelB_gse161529.py"
```
