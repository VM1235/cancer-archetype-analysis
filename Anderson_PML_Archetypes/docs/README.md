# Anderson et al. (Mol Cancer Res 2026) — Fig 3A / 3B rebuild

Rebuild of PML archetype PCA panels from:

> Anderson et al. Archetype analysis of lung adenocarcinoma premalignancy…  
> *Mol Cancer Res* 2026;24:737–53. Code: https://github.com/kelley27/PML_Archetypes

## What we reproduced

| Panel | Content |
|---|---|
| **3A** | PC1–PC2 of 9-module GSVA scores; color = signed (−) distance to each of 4 archetype vertices |
| **3B** | Same PC space; primary = closest **12.5%** to each vertex (unique only); else secondary; shape = histology |

## Data (all 4 cohorts, n = 137 PMLs)

| Study | Source | PMLs |
|---|---|---|
| Lung PCA | GEO GSE319666 counts | 68 AAH + 12 AIS/MIA |
| Sivakumar | GSE102511 TPM | 17 AAH |
| Yoo | GSE166720 FPKM | 17 AIS/MIA |
| Altorki | GSE193725 FPKM | 23 AIS/MIA |

Author **gene modules** from GitHub `outputData/pml_modules.rds` (not re-derived).

## Caveats (why this is a rebuild, not pixel-identical)

- Authors’ residual RDS matrices are request-only; we rebuilt expression (log-CPM / log2 TPM+1 / log2 FPKM+1) + limma residuals `~ Study + TopGenes`.
- Public GEO matrices are TPM/FPKM, not raw counts (paper used counts via GEOquery).
- Archetypes fit with project `src.archetypes` PCHA (`delta=0`, 150 inits), not ParetoTI bootstrap (`n=5000`).

Despite that, primary labels agree with authors on **~88%** of samples that both call primary.

## Run

From repo root (OpenMP workaround needed on this Mac):

```bash
export KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 MPLBACKEND=Agg
cd Anderson_PML_Archetypes
Rscript codes/01_build_pml_matrix.R
Rscript codes/02_compute_residuals.R
Rscript codes/02b_gsva.R
../.venv/bin/python -u codes/03_fit_archetypes_plot_fig3.py
```

Outputs: `results/fig3/fig3A_distance_panels.{png,pdf}`, `fig3B_primary_assignments.{png,pdf}`.
