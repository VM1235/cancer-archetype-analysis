# KS-restricted breast archetypes → METABRIC projection

Tan 2014 KS epithelial–mesenchymal gene list (218 genes) restricts the DepMap breast cell-line matrix. Archetypes are fit on lines only (Panel A); METABRIC tumors are projected without refitting PCHA (Panel C).

## Prerequisites

- Repo setup: root `README.md` (`.venv`, `requirements.txt`, `src/`).
- `Breast Cancer/data/processed/input_panelA_ks_genelist.csv` (206 genes × 63 lines) — in git.
- METABRIC raw files — see `brca_metabric/DOWNLOAD.md` (~1.3 GB, local only).
- R `sva` in `Breast Cancer/rlib/` for ComBat (or Python fallback in `src/combat.py`).

## Run order

From the **repository root**:

### 1. Build KS cell-line matrix (once)

```bash
.venv/bin/python -u "Breast Cancer/codes/00_build_ks_genelist_input.py"
```

### 2. Panel A — default k=3

```bash
.venv/bin/python -u "Breast Cancer/codes/run_panelA_ks_genelist.py"
```

Outputs: `results/panel_a_ks_genelist/` (`n_pcs=2`, `suggested_k=3`).

### 3. Panel A — extended k (for k=4..7)

```bash
# Full sweep (slow; ~500 shuffles per k):
.venv/bin/python -u "Breast Cancer/codes/run_panelA_ks_genelist_extendedk.py"

# Or k=4 only (~2 min):
.venv/bin/python -u "Breast Cancer/codes/run_panelA_ks_genelist_extendedk.py" --only-k 4
```

Outputs: `results/panel_a_ks_genelist_extendedk/` (`n_pcs=6`, `archetypes_k4_parti.npy`, etc.).

### 4. Panel B (optional) — PAM50 enrichment on KS archetypes

```bash
.venv/bin/python -u "Breast Cancer/codes/run_panelB_ks_genelist.py"      # k=3
.venv/bin/python -u "Breast Cancer/codes/run_panelB_ks_genelist_k4.py"  # forced k=4
```

### 5. Prepare METABRIC

```bash
.venv/bin/python -u "Breast Cancer/codes/prepare_metabric_brca.py"
```

Gene matching report → `results/panel_c_metabric_ks_genelist/prepare_report.txt`  
Prepared expression → `metabric_log_shared_genes.csv` (203 genes × 1,980 tumors).

### 6. Project METABRIC into archetype space

```bash
# k=3 (2-D simplex)
.venv/bin/python -u "Breast Cancer/codes/run_panelC_metabric_ks_genelist.py" --n-shuffle 20

# k=4 (3-D simplex; reuses ComBat from k=3 run if present)
.venv/bin/python -u "Breast Cancer/codes/run_panelC_metabric_ks_genelist_k4.py" --n-shuffle 20
```

`--force-combat` rebuilds the ComBat matrix instead of reusing `CCLE_METABRIC_COMBAT.csv`.

## Key outputs

| k | Results folder | Figure |
|---|----------------|--------|
| 3 | `results/panel_c_metabric_ks_genelist/` | `figures/Figure_1C_metabric_ks_genelist.png` |
| 4 | `results/panel_c_metabric_ks_genelist_k4/` | `figures/Figure_1C_metabric_ks_genelist_k4.png` |

Subtype breakdown tables: `claudin_subtype_by_nearest_archetype_counts.csv`, `ihc_by_nearest_archetype_counts.csv`, `metabric_sample_subtype_archetype.csv`.

## TCGA comparison (same KS genes)

```bash
.venv/bin/python -u "Breast Cancer/codes/prepare_tcga_brca.py"   # if not done
.venv/bin/python -u "Breast Cancer/codes/run_panelC_ks_genelist.py" --n-shuffle 20
```

Figure: `figures/Figure_1C_tcga_ks_genelist.png`.

## Notes

- **PCHA is never refit on tumors.** Archetypes come from cell lines; tumors enter via ComBat + combined PCA projection (Groves Panel C).
- k=3 Panel A uses 2 PCs; k=4 requires extended-k Panel A (3 PCs for the 4-vertex simplex).
- Large ComBat / merged matrices are gitignored; rerun Panel C or copy from a prior local run.
