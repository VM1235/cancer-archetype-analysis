# Project context for a new Cursor chat

Read this fully before editing code. Workspace: `/Users/apple/Desktop/thesis_iisc/Project_1`.  
GitHub (private): https://github.com/VM1235/cancer-archetype-analysis (`main`).  
Python: `.venv` at repo root. **NumPy must stay &lt; 2** (`py_pcha` still uses `np.mat`). `requirements.txt`.

Run all scripts from the **repository root** unless a README says otherwise. Disease folders have spaces or long names; quote paths.

```bash
.venv/bin/python -u "Glioblastoma/codes/run_panelA.py"
```

This document is for navigation, file meaning, caveats, and next steps. Per-disease how-to-run lives in each folder’s README.

---

## 1. What the project is

Reproduce and reuse **Figure 1A–C** from Groves et al., *Cell Systems* 13:690–710 (2022): Pareto Task Inference / PCHA on bulk RNA-seq.

- **A:** Fit archetypes on **cell lines**. Choose \(k\) from ESV elbow + t-ratio permutation (smallest significant \(k\), else DimensionFinder).
- **B:** Independent subtype labels vs distance-to-vertex bins (hypergeometric, BH, call only if peak is **bin 0**).
- **C:** Project **tumors** into that cell-line geometry after ComBat (cell line vs tumor). Do **not** refit PCHA on tumors. Scatter + variance vs tumor-only PCA.

PDFs: `papers/Groves_Cell_Systems_2022.pdf`, `papers/Hausser_NatCommun_2019.pdf`.  
SCLC science write-up: `SCLC Reproduction - Groves Cell Systems 2022/docs/` (`01_what_the_figures_mean.md`, `02_how_we_implemented_this.md`, `PROJECT_CONTEXT.md`).  
Authors’ MATLAB/data clone (gitignored, local only): `reference/` = [QuLab-VU/Groves-CellSys2022](https://github.com/QuLab-VU/Groves-CellSys2022).

**Do not copy SCLC \(k\), PC count, or subtype names onto breast or GBM.** Each matrix gets its own dim (first n PCs with ≥50% variance) and its own \(k\).

---

## 2. How to walk the tree

```
Project_1/
  README.md                          # high-level map
  PROJECT_CONTEXT.md                 # this file
  requirements.txt
  src/                               # SHARED engines (import from repo root)
  papers/
  reference/                         # gitignored Groves clone
  .venv/
  SCLC Reproduction - Groves Cell Systems 2022/
  Breast Cancer/
  Glioblastoma/
```

| Folder | Role | Status |
|---|---|---|
| SCLC | Reproduce Groves Fig 1A–C | Official figures exist; t-ratios match paper |
| Breast | Same **method** on DepMap invasive breast + TCGA-BRCA | Official = current `codes/` + `figures/` (not `archive/run1/`) |
| GBM | Same **method** on DepMap GB + TCGA-GBM | First full A/B/C run done 18 Aug 2026; **not** a Groves-like win |

`src/paths.py`: `ROOT`, `SCLC`, `BREAST`, `GBM`, `PAPERS`, `REFERENCE`.

---

## 3. Shared library (`src/`)

Disease scripts do `sys.path.insert(0, str(ROOT))` then `from src....`.

| Module | What it does | Used by |
|---|---|---|
| `archetypes.py` | `fit_pcha`, `fit_pcha_best` (multi-start, keep **max volume**), `simplex_volume` (ParTI det formula), `esv_curve`, `_shuffle_columns` | All Panel A |
| `pca.py` | `fit_pca`, `cumulative_variance`, `inverse_transform_scores`, `align_pca_signs` | A and C |
| `enrichment.py` | `distance_bins` (equal-count), `hypergeometric_enrichment` (BH, `sig_peak_at_bin0`) | Panel B |
| `io.py` | Load **SCLC** Groves CSVs (SCLC `data/` first, else `reference/`). `load_expression_csv(path)` also loads any genes×samples CSV | SCLC defaults; all diseases pass an explicit path |
| `preprocess.py` | DepMap: `pick_breast_models`, `pick_gbm_models`, `expression_to_genes_by_samples`, `drop_low_genes` | Breast/GBM `01_build` |
| `combat.py` | Python ComBat fallback if R `sva` missing | Panel C |
| `paths.py` | Folder constants | GBM Panel C rlib fallback to breast `rlib/` |

**ParTI / Groves PCHA convention (the official one):**

- Fit on first **\(k-1\)** PCs (a \(k\)-simplex is \((k-1)\)-D).
- `delta=0` (archetypes are convex combinations of real samples).
- Observed: `3 * numIter` random inits, keep max-volume simplex. Paper `numIter=50` → **150** inits.
- Each shuffle: independently permute each PC across samples; `numIter` inits (paper → **50**); keep max t-ratio.
- \( t = V_{\mathrm{simplex}} / V_{\mathrm{convex\ hull}} \) in \((k-1)\)-D.
- \( p = \) fraction of finite null \(t\) ≥ observed \(t\).
- Groves SCLC used **1000** shuffles. Breast Panel A used a **reduced** `numIter=5` (15/5 inits) and 500 shuffles. **GBM used paper inits (150/50) and 500 shuffles.**

Do not use SCLC `run_panel_a_firstpass.py` / `run_panel_a_permutation.py` as the scientific result (`delta=0.1`, 12-D fit). Official SCLC is `run_panel_a_parti_full.py` + `export_sclc_reproduction_figures.py`.

---

## 4. SCLC folder (reproduction)

Path: `SCLC Reproduction - Groves Cell Systems 2022/`

| Path | Role |
|---|---|
| `data/` | Groves processed matrices (copied from `reference/`) |
| `codes/` | Drivers |
| `results/panel_{a,b,c}/` | Intermediates |
| `figures/Figure_1A.png` … `1C` | Official panels |
| `docs/` | Methods notes, old `PROJECT_CONTEXT.md`, slides |

**Data (do not re-log, do not rebuild ComBat):**

- `data/SCLC_combined_Minna_CCLE_batch_corrected_wo_lowgenes.csv` — A/B, ~15950 genes × 120 lines
- `data/combined_clusters_2020-05-27-MC copy.csv` — B labels, column `NEW_10_2020`
- `data/CCLE_Minna_Thomas_COMBAT.csv` — C, lines + 81 tumors
- `data/Metadata_CCLE_Minna_Thomas_COMBAT.csv`

Minna RNA is dbGaP `phs001823`. We never downloaded FASTQs.

**Codes:** `run_panel_a_parti_full.py` (1000 shuffles, checkpoints), `run_panel_b.py`, `run_panel_c.py`, exporters. Older first-pass scripts are not official.

**Key intermediates:** `results/panel_a/pc_scores_12.npy`, `archetypes_k{k}_parti.npy`, `S_k{k}_parti.npy`, `null_t_ratios_k{k}_parti_n1000.npy`, `t_ratio_official.csv`. B/C read Panel A archetypes; C does **not** refit PCHA.

---

## 5. Breast folder

Path: `Breast Cancer/` (space in the name).

Official analysis is **not** `archive/run1/`. Use `codes/`, `figures/`, `results/`.

| Path | Role |
|---|---|
| `data/raw/` | `Model.csv`; huge `OmicsExpression…csv` (gitignored) |
| `data/processed/` | Panel A matrix 16500 genes × **63** lines |
| `data/tumors/` | Xena BRCA `HiSeqV2` (gitignored), clinical matrix |
| `codes/` | Build → PAM50 → A → B → prepare TCGA → C |
| `rlib/` | Local Bioconductor `sva` (gitignored); GBM C can reuse this |
| `docs/` | `PARAMS_vs_Groves.txt`, `PARAMS_panelC.txt`, figure notes |

**Run order:**

1. `01_build_input_panelA.py` — OncoTree invasive BRCA cell lines, drop all-low genes  
2. `03_map_and_pam50.R` + `04_match_pam50_to_panelA.py` — genefu PAM50 (only if labels missing)  
3. `run_panelA.py` — **numIter=5**, 500 shuffles (lighter than paper)  
4. `run_panelB.py` — PAM50; **drops LumA and Normal** from enrichment only; polytope still all 63  
5. `prepare_tcga_brca.py`  
6. `run_panelC_tcga.py` — IHC ER/HER2 colors, not PAM50  

Suggested \(k=4\) from DimensionFinder; **no k with p&lt;0.05** (like GBM, unlike SCLC).

---

## 6. Glioblastoma — files, codes, how they chain

This is the latest work. Data sit **directly** under `Glioblastoma/data/` (not `data/raw/` / `data/tumors/` like breast).

### 6.1 Raw / input data

| File | Source | Git | Role |
|---|---|---|---|
| `data/Model.csv` | DepMap | yes | Metadata |
| `data/OmicsExpressionProteinCodingGenesTPMLogp1.csv` | DepMap, **already log2(TPM+1)** | **ignored** (~500 MB) | Cell-line RNA |
| `data/HiSeqV2` | UCSC Xena TCGA-GBM, **already log2** | **ignored** | Tumor RNA, 20530 genes × **172** samples (RNA-seq subset; clinical has 629 rows) |
| `data/TCGA.GBM.sampleMap-GBM_clinicalMatrix` | Xena | yes | `GeneExp_Subtype` (Verhaak), `G_CIMP_STATUS`, `sample_type_id` |

### 6.2 Codes (run from repo root, this order)

| Script | Writes | Reads |
|---|---|---|
| `codes/01_build_input_panelA.py` | `data/processed/input_panelA_glioblastoma_ccle_logtpm_filtered.csv`, `input_panelA_models_used.csv`, `input_panelA_build_report.txt` | DepMap Model + Omics |
| `codes/02_assign_verhaak.py` | `results/panel_b/verhaak_labels_panelA.csv` | Processed matrix + (markers only; no TCGA) |
| `codes/run_panelA.py` | `results/panel_a/*`, `figures/Figure_1A_gbm.png` | Processed matrix |
| `codes/run_panelB.py` | `results/panel_b/enrichment_*`, `figures/Figure_1B_gbm.png` | Matrix, `pc_scores.npy`, `archetypes_k{k}_parti.npy`, Verhaak labels, `suggested_k.txt` |
| `codes/prepare_tcga_gbm.py` | `results/panel_c_tcga/tcga_primary_log_shared_genes.csv` (gitignored), `tcga_primary_metadata.csv`, `prepare_report.txt` | HiSeq + clinical + Panel A gene list |
| `codes/combat_cellline_tumor.R` | Called by C | merged + batch CSVs |
| `codes/run_panelC_tcga.py` | ComBat matrix (gitignored), PCA tables, `figures/Figure_1C_gbm.png` | Panel A archetypes + prepared tumors |

`pick_gbm_models` in `src/preprocess.py`: OncoTree **GB/GBM** or disease/subtype containing “glioblastoma”. **Not** all CNS/Brain (that would add astrocytoma, medulloblastoma, etc.). Then `ModelType` must be cell line.

**Build outcome:** 67 GB models in metadata, **54 with RNA**, 19205 → **15833** genes after max log2(TPM+1) &lt; 1 filter.

### 6.3 Panel A intermediates (`results/panel_a/`)

| File | Meaning |
|---|---|
| `pc_scores.npy` | 54 × n_pcs scores. **n_pcs=12** (51.4% variance), same 50% rule as Groves |
| `n_pcs.txt` | `12` |
| `pca_variance.csv`, `pca_variance_full.csv` | Scree |
| `esv_curve.csv` | ESV vs k=2..13 in 12-D |
| `archetypes_k{k}_parti.npy` | k × (k−1) archetype coordinates in PC space |
| `S_k{k}_parti.npy` | 54 × k mixture weights |
| `null_t_ratios_k{k}_parti_n500.npy` | 500 null t-ratios (checkpointed every 25) |
| `t_ratio_parti_500_k{k}.csv` | One-row result per k (parallel workers) |
| `t_ratio_parti_500.csv` | Merged table |
| `suggested_k.txt` | First line integer k; second line reason |

`run_panelA.py` runs k=3..7 **in parallel** (`ProcessPoolExecutor`). Log lines like `100/500` **do not print which k**; look at `null_t_ratios_k*_n500.npy` sizes or `t_ratio_parti_500_k*.csv`. Restart is safe: complete k’s are skipped; partial nulls are resumed.

**Settings in this GBM run:** `NUM_ITER=50`, `N_INIT_OBS=150`, `N_INIT_NULL=50`, `N_PERM=500`, `DELTA=0`. ~20 min on this machine for all k.

**Result of the 18 Aug 2026 run:**

| k | t-ratio | p (500 shuffles) |
|---|---|---|
| 3 | 0.580 | 0.224 |
| 4 | 0.297 | 0.238 |
| 5 | 0.102 | 0.510 |
| 6 | 0.034 | 0.694 |
| 7 | 0.027 | 0.068 |

**No k with p&lt;0.05.** `suggested_k.txt` = **7** (DimensionFinder elbow). Figure 1A still draws the k=7 polytope; that is a fallback, not a Groves-style significant simplex.

### 6.4 Panel B intermediates (`results/panel_b/`)

| File | Meaning |
|---|---|
| `verhaak_labels_panelA.csv` | One row per ACH- ID: `verhaak_subtype`, `confidence_score`, marker scores |
| `verhaak_labels_all.csv` | Same aligned to Panel A column order |
| `verhaak_labels_panelB_input.csv` | After optional drops (none in this run) |
| `distances_verhaak.csv`, `bins_verhaak.csv` | Distance / bin 0–4 to each of 7 arcs, in **(k−1)=6 D** |
| `enrichment_verhaak.csv` | Hypergeometric + BH; `sig_peak_at_bin0` is the Groves call |

**Labels are not TCGA Verhaak.** First attempt (Spearman to TCGA centroids) assigned **all 54 lines Mesenchymal** (platform batch). Replaced by `02_assign_verhaak.py`: gene-wise z-score on the **DepMap matrix**, mean of marker sets:

- Classical: EGFR, NES, NOTCH3, SMO, GAS1, GLI2, NFKBIA  
- Mesenchymal: CD44, CHI3L1, RELB, TRADD, TNFRSF1A, NFKB1, STAT3  
- Proneural: OLIG2, PDGFRA, DLL3, NKX2-2, SOX10, ASCL1, NKX2-1  

This run: Proneural 21, Mesenchymal 20, Classical 13. Neural not used (Wang 2017 three-type scheme). Enrichment **n_bins=5**. **No significant bin-0 peaks.**

### 6.5 Panel C intermediates (`results/panel_c_tcga/`)

`prepare_tcga_gbm.py`: primary tumors only, intersect genes with Panel A → **13677 genes × 154 tumors**. Clinical kept: sample type, histology, `GeneExp_Subtype`, G-CIMP, etc.

`run_panelC_tcga.py` chain:

1. Rebuild cell-line PCA (12 PCs), align signs to saved `pc_scores.npy`.
2. Map k=7 archetypes from (k−1)-D → 12-D (pad zeros) → **gene space** (`archetypes_gene_space_shared.csv`).
3. Merge lines + tumors on shared genes → `merged_uncorrected_shared_genes.csv` + `combat_batch.csv`.
4. ComBat R `sva` in `Breast Cancer/rlib/` (or Python fallback): batch = cell_line vs tumor, `mod=~1`, **`ref.batch=cell_line`**. Output `CCLE_TCGA_COMBAT.csv` (**gitignored**).
5. **New** PCA on ComBat matrix (up to 20 PCs). Transform gene-space archetypes into that PCA (`archetypes_in_combined_pca.csv`). **No PCHA on tumors.**
6. Inside-simplex: barycentric weights in **6-D** (k−1). This run: **0/154** inside.
7. Colors: TCGA `GeneExp_Subtype` (Classical 39, Mesenchymal 49, Proneural 38, Neural 26, 2 missing).
8. Variance: tumor-only PCA ceiling vs combined-PCA scores of tumors vs gene-row shuffles (`--n-shuffle 20`).

**Do not over-read “7 combined PCs = 97% of tumor-only ceiling.”** Combined PCA is fit on 54 lines + **154 tumors**, so tumors dominate the axes. Groves’ claim is tumors **inside the cell-line simplex**, which failed here.

`--force-combat` rebuilds the huge ComBat CSV. Without it, an existing `CCLE_TCGA_COMBAT.csv` is reused.

R wrapper: `codes/combat_cellline_tumor.R` (copy of breast). Needs `Rscript` + `sva` in rlib.

### 6.6 Figures

- `figures/Figure_1A_gbm.png`  
- `figures/Figure_1B_gbm.png`  
- `figures/Figure_1C_gbm.png`  

---

## 7. How A → B → C share objects (any disease)

```
raw expression
    → 01_build → genes × cell lines  (no extra log if DepMap/Xena already log)
         → Panel A: PCA(n_pcs) → PCHA on (k-1) PCs
              saves: pc_scores.npy, archetypes_k*_parti.npy, S_k*, suggested_k.txt
         → Panel B: independent labels + distance bins in (k-1)-D using those archetypes
         → Panel C: inverse-PCA archetypes to gene space
                    merge with tumors → ComBat
                    new combined PCA; project archetypes; do not refit PCHA
```

If you change the Panel A matrix, **rerun A then B then C**. Stale `archetypes_k7_parti.npy` with a new matrix is a silent bug. `suggested_k.txt` is what B and C use for which k file to load (GBM C reads that; do not hardcode k=4 like an older breast C snippet).

---

## 8. Caveats (especially GBM)

**Correctness of the run**

- Right DepMap GB cell lines (not all CNS). Right TCGA-GBM HiSeq + clinical. No double log. Groves ComBat model. Groves PCHA inits for GBM. Permutation 500/500 succeeded, n_fail=0.

**Why this is not “Groves on GBM”**

- Groves needs a **significant** simplex, **1-to-1** vertex–subtype map, tumors **inside** that simplex. GBM: no p&lt;0.05, no bin-0 subtype peaks, 0/154 tumors inside.
- n=54 is small vs SCLC 120; k=7 is a 6-D simplex — easy to overfit visually in PC1–PC2.
- k=7 t-ratio is **0.027**: the simplex is a small object in the hull even if p=0.068 is the smallest p.
- Cell-line “Verhaak” labels are **markers on CCLE**, not published classifications. TCGA centroid transfer failed (all Mesenchymal).
- Panel B with 7 arcs × 3 subtypes × 54 lines is underpowered; “no peak” is consistent with A but not a high-powered negative.
- Variance plot on C is **not** evidence the polytope generalizes when tumors outnumber lines.
- Breast also lacked p&lt;0.05 (and used fewer PCHA inits). SCLC is the only clean Fig 1A–C match. Do not describe GBM as a failed script; describe it as a **failed geometric hypothesis on this matrix**.

**Engineering**

- Huge files gitignored; clones without Omics/HiSeq/ComBat CSVs cannot rebuild C from raw without re-download. Processed Panel A CSV **is** in the tree (~small enough).
- `reference/` is 3.8 GB and gitignored.
- Parallel Panel A: set `OMP_NUM_THREADS=1` inside workers (already in script).
- macOS multiprocessing: `if __name__ == "__main__"` guard is required (present).

---

## 9. What a new chat might do next

Do not silently rerun 500×50 PCHA unless asked (slow). Prefer reading saved `results/`.

**Interpretation / thesis**

- Write GBM as a **contrast case** vs SCLC: same protocol, no Pareto simplex in 54 GB lines.
- Compare t-ratio tables SCLC vs breast vs GBM side by side.
- Note culture bias: CCLE GBM often mesenchymal-leaning vs TCGA Verhaak diversity.

**If the PI wants another GBM attempt (science, not plumbing)**

1. **Tumor-only Panel A** on the 154 HiSeq primaries (Groves found tumor polytopes may need **linear** not log space). Then ask whether **cell lines** fall inside the **tumor** simplex (invert C).  
2. **k=3** Wang types only, even if A is not significant, as a descriptive figure — label as exploratory.  
3. Better independent line labels: ssGSEA with Verhaak 840-gene signatures, or SVM trained on TCGA **after** ComBat (labels still not fully independent of C’s batch correction).  
4. Broader CNS filter vs GB-only: more n, dirtier biology.  
5. Linear-space PCHA on lines.  
6. Match breast inits (numIter=5) only for speed; do not call that the paper protocol.  
7. G-CIMP / IDH as Panel C colors instead of (or in addition to) Verhaak.

**If the PI wants code hygiene**

- Print `k=` on every shuffle checkpoint (GBM A logs are interleaved).  
- Unify breast/GBM `data/raw` vs `data/` layouts.  
- Stop hardcoding k=4 in any leftover breast C paths (GBM C already uses `suggested_k.txt`).  
- Git-add GBM figures if they should live on GitHub; processed matrix yes; Omics/HiSeq/ComBat no.

**Do not**

- Copy SCLC k=5 or 12 PCs into GBM.  
- Treat DimensionFinder k=7 as “the number of GBM archetypes” in a results sentence without the p-values.  
- Re-fit PCHA on tumors and call it Panel C.  
- Log-transform DepMap or Xena HiSeqV2 again.

---

## 10. Quick commands

```bash
# env
.venv/bin/python -c "import numpy; print(numpy.__version__)"   # must be 1.x

# GBM rebuild (A is slow)
.venv/bin/python -u "Glioblastoma/codes/01_build_input_panelA.py"
.venv/bin/python -u "Glioblastoma/codes/02_assign_verhaak.py"
.venv/bin/python -u "Glioblastoma/codes/run_panelA.py"
.venv/bin/python -u "Glioblastoma/codes/run_panelB.py"
.venv/bin/python -u "Glioblastoma/codes/prepare_tcga_gbm.py"
.venv/bin/python -u "Glioblastoma/codes/run_panelC_tcga.py" --n-shuffle 20
```

SCLC and breast commands: their READMEs.

---

## 11. Git / what is not on GitHub

Private repo `VM1235/cancer-archetype-analysis`. Ignored: `.venv`, `reference/`, `rlib/`, DepMap Omics CSVs, Xena HiSeqV2, Panel C ComBat/merged/tumor-expression CSVs.  
Present if committed: code, READMEs, figures, processed Panel A matrices, labels, t-ratio tables, small metadata.

---

## 12. One-sentence status

SCLC Fig 1A–C is a successful Groves reproduction; breast is the same method with mixed/weak significance; GBM (54 GB lines, paper PCHA inits, 500 shuffles, 154 TCGA tumors) **does not form a significant cell-line simplex, does not put Verhaak-style markers at vertices, and does not place tumors inside that simplex** — with the protocol and files above.
