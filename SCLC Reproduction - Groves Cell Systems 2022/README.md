# SCLC reproduction — Groves et al., *Cell Systems* 2022

Reproduction of **Figure 1A–C** from:

> Groves, S.M. et al. Archetype tasks link intratumoral heterogeneity to plasticity and cancer hallmarks in small cell lung cancer. *Cell Systems* 13, 690–710 (2022).

We did **not** re-derive ComBat from raw FASTQs. Minna Lab RNA-seq is dbGaP-controlled (`phs001823`). Input matrices are the authors’ published, batch-corrected tables from [QuLab-VU/Groves-CellSys2022](https://github.com/QuLab-VU/Groves-CellSys2022).

Narrative notes: `docs/01_what_the_figures_mean.md`, `docs/02_how_we_implemented_this.md`. Pipeline spec: `docs/PROJECT_CONTEXT.md`.

## Folder map

| Path | Contents |
|---|---|
| `data/` | Expression matrices and subtype labels the codes load |
| `codes/` | Panel A/B/C scripts and figure exporters |
| `figures/` | Official Figure 1A, 1B, 1C PNGs |
| `results/` | Intermediate arrays, t-ratio table, enrichment tables |
| `docs/` | Methods write-up, slides |

Shared engines: repo-root `src/`.

## Data the codes expect

All under `data/` (genes × samples unless noted):

| File | Used for | Shape (approx.) |
|---|---|---|
| `SCLC_combined_Minna_CCLE_batch_corrected_wo_lowgenes.csv` | Panels A and B | 15,950 genes × 120 cell lines, log, ComBat (Minna vs CCLE) |
| `combined_clusters_2020-05-27-MC copy.csv` | Panel B labels | 120 lines; column `NEW_10_2020` = A / A2 / N / P / Y |
| `CCLE_Minna_Thomas_COMBAT.csv` | Panel C | 14,546 genes × 201 samples (120 lines + 81 tumors), second ComBat |
| `Metadata_CCLE_Minna_Thomas_COMBAT.csv` | Panel C sample metadata | optional |

If a file is missing here, `src/io.py` falls back to `reference/data/bulk-rna-seq/` (local clone of the Groves repo; not in git).

**Do not** put dbGaP Minna FASTQs in this folder.

## How to run

From the **repository root**, with `.venv` installed (`pip install -r requirements.txt`):

```bash
# Panel A (ParTI-matched PCHA; full 1000-shuffle test is slow)
.venv/bin/python -u "SCLC Reproduction - Groves Cell Systems 2022/codes/run_panel_a_parti_full.py"

# Faster sanity check (100 shuffles)
.venv/bin/python -u "SCLC Reproduction - Groves Cell Systems 2022/codes/run_panel_a_sanity_parti.py"

# Panel B (needs Panel A outputs in results/panel_a/)
.venv/bin/python -u "SCLC Reproduction - Groves Cell Systems 2022/codes/run_panel_b.py"

# Panel C (needs Panel A archetypes)
.venv/bin/python -u "SCLC Reproduction - Groves Cell Systems 2022/codes/run_panel_c.py"

# Rebuild the three official figures
.venv/bin/python -u "SCLC Reproduction - Groves Cell Systems 2022/codes/export_sclc_reproduction_figures.py"
```

Older exploratory scripts (`run_panel_a_firstpass.py`, `run_panel_a_permutation.py`) are **not** the official result. Official calling convention: fit in \(k-1\) PCs, `delta=0`, many random starts, keep the maximum-volume simplex. See `docs/02_how_we_implemented_this.md`.

## Outputs

- Figures: `figures/Figure_1A.png`, `Figure_1B.png`, `Figure_1C.png`
- t-ratios vs paper: `results/t_ratio_official.csv`
- Fits and nulls: `results/panel_a/`, `results/panel_b/`, `results/panel_c/`
