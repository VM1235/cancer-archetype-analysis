# Pipeline audit: PCHA / ParTI across SCLC, breast, GBM

**Scope.** Trace executed code paths and saved artifacts. SCLC numerical agreement with Groves et al. 2022 is taken as already established; this audit asks whether that agreement is from the *correct* protocol, and whether breast and GBM apply the *same* protocol.

**What was not done.** No PCHA/ComBat/PAM50/PCA refits, no `numIter` changes, no figure regeneration. Counts below come from existing CSVs, npy headers, and metadata (`Model.csv`), not from rebuilding expression matrices from Omics.

**Git snapshot.** `HEAD` = `62a37eeb82ede59d28118094b42d2f83aca4b6fa` (2026-08-18 11:19 +0530, “Organize the thesis repo…”). Working tree: GBM codes/results **untracked**; `src/preprocess.py`, READMEs, `.gitignore` modified.

---

## CRITICAL ISSUES

Findings that would change a reported number, figure, or scientific claim if fixed. Recommended fixes are listed; **corrected numbers were not recomputed** (no reruns).

1. **SCLC Panel B and Panel C load the discarded 12-PC / `delta=0.1` fit, not the ParTI-matched fit.**
   - Current loads: `run_panel_b.py:104` and `run_panel_c.py:54` → `results/panel_a/archetypes_k5.npy`.
   - That file is **(5, 12)** (npy header). It is written by `run_panel_a_firstpass.py:89` from a 12-D, `delta=0.1`, single-start fit (`DELTA=0.1` at line 32; `permutation_t_ratio` → `fit_pcha` default path).
   - The ParTI-matched file is `archetypes_k5_parti.npy` **(5, 4)**, written by `run_panel_a_sanity_parti.py:61` / reused by `run_panel_a_parti_full.py:59`.
   - Official Figure 1A uses the *parti* file (`export_sclc_reproduction_figures.py:65`). Official Figure 1B replots `enrichment_author.csv` produced by `run_panel_b.py`. Official Figure 1C replots `panel_c/` CSVs produced by `run_panel_c.py`.
   - **Implication.** A/B/C are not the same simplex. Panel A t-ratios can match Groves while B/C describe a different polytope. Matching Groves’ 1-to-1 subtype map on B does **not** prove the ParTI vertices were used.
   - **Fix.** Point B and C at `archetypes_k5_parti.npy`; in B/C restrict scores to the first `k-1=4` columns before distances / inverse-PCA (pad zeros to 12-D for gene-space inverse, as breast/GBM already do). Then regenerate B, C, and Figures 1B–1C. Do not treat a coincidental match as validation.

2. **SCLC Panel B distances are not in the PCHA fitting space even under the official convention.**
   - `run_panel_b.py:103–104,76`: `distance_bins(scores, archetypes)` with `pc_scores_12.npy` **(120, 12)** and no `[:, :k-1]` slice.
   - Breast (`run_panelB.py:106–129`) and GBM (`run_panelB.py:103–125`) *do* slice to `n_vol = arcs.shape[1]`.
   - SCLC docs (`docs/02_how_we_implemented_this.md` Panel B) even *state* 12-PC distances. That is a protocol divergence from Groves/ParTI (fit in `k-1`) and from breast/GBM, independent of which npy is loaded.
   - **Fix.** Same as (1): distances in 4-D for k=5.

3. **SCLC Panel C never computes containment.**
   - `run_panel_c.py` has no barycentric / `inside_simplex` logic (grep empty). “Mostly contained” is a scatter-plot claim only.
   - Breast reports **0/1097** and GBM **0/154** from the same barycentric test. SCLC has no comparable fraction, so cross-cancer containment statements are not on equal footing.
   - **Fix.** Port `barycentric_weights` from breast/GBM Panel C; report `n_inside / n_tumor` in the same `(k-1)`-D combined-PCA space.

4. **`numIter` is not the same method across cancers, and the write-ups disagree with ParTI.m and with each other.**
   | Source | `numIter` | observed inits (`3*numIter`) | null inits |
   |---|---|---|---|
   | `reference/ParTI/ParTI.m:110–113` (`algNum==5`) | **5** | 15 | 5 |
   | SCLC `run_panel_a_parti_full.py:44–46` | **50** | 150 | 50 |
   | GBM `run_panelA.py:45–47` | **50** | 150 | 50 |
   | Breast `run_panelA.py:47–49` | **5** (hardcoded) | 15 | 5 |
   | `py_pcha.PCHA` | no `numIter` | n/a | n/a |

   Breast `NUM_ITER=5` is **not** a default of the library being called (`py_pcha`). It is a manual constant, copied from MATLAB `ParTI.m`’s PCHA branch. That MATLAB default **is real** (verified in this repo’s `reference/ParTI/ParTI.m` and AlonLabWIS/ParTI). Groves `params_lognorm.txt` records `algNum=5`, so a stock ParTI run would have used 15/5, not 150/50.

   SCLC nonetheless reproduced paper t-ratios with **150/50** (`t_ratio_parti_1000.csv`). GBM README claims `numIter=50` “matches the SCLC paper / ParTI `algNum=5`” — the second half of that sentence is **false**.

   **Implication.** Breast is closer to stock ParTI.m; SCLC/GBM are closer to each other and to Sisal’s `numIter=50`. They are not one protocol. Breast p-values (none &lt; 0.05) are not strictly comparable to SCLC/GBM p-values.
   - **Fix (science).** Pick one spec and apply it to all three: either ParTI.m PCHA (5) or the SCLC Python run that matched Groves’ table (50). Re-run only the disease(s) that change. This audit did not re-run.

5. **GBM hypothesis-driven k=2: degenerate 1-PC result (`t=1`, `p=1`) is what Panel B/C k=2 consume; the 12-PC correction is isolated.**
   - Saved separately: `results/panel_a_k2/archetypes_k2_parti.npy` **(2, 1)** vs `results/panel_a_k2_12/archetypes_k2_12pc.npy` **(2, 12)**. They do not overwrite `panel_a/` or each other. **PASS** on file isolation.
   - `run_panelB_k2.py:28` and `run_panelC_k2.py:39` load **only** `panel_a_k2/` (1-D). `panel_c_k2/subtype_inside_simplex_counts.csv` is **152/152 inside** (all True) — expected if vertices are 1-D hull endpoints.
   - Official GBM Figure 1B/C (`run_panelB.py`, `run_panelC_tcga.py`) use `suggested_k.txt` **k=7**, not k=2. Main A/B/C are not silently the degenerate fit.
   - **Risk.** Any prose that cites k=2 tumor containment or enrichment as “the corrected k=2 analysis” is wrong unless it names `panel_*_k2` vs `panel_a_k2_12`.

6. **GBM cell-line “Verhaak” labels are not Verhaak classifications.**
   - `02_assign_verhaak.py:26–46`: 7-gene sets, gene-wise z-score **across the 54 DepMap lines**, mean of available markers, **argmax**. Neural omitted. `ASCL1` is **absent** from the Panel A matrix (not in the CSV; `genes_Proneural` lists 6 genes, no ASCL1).
   - This is weaker than breast’s `genefu` PAM50 and is not independent published metadata. Do not describe Panel B as testing Verhaak subtypes on cell lines.

---

## CONFIRMED CORRECT

Do not re-litigate these unless the inputs above change.

- Shared engines exist and are what breast/GBM Panel A/B actually call: `src/archetypes.py` (`fit_pcha_best` slices to `k-1`, `delta` passed in, max **volume**), `simplex_volume` (ParTI det formula), `_shuffle_columns` (independent PC permutation), `src/pca.py` `fit_pca` (sklearn PCA, mean-centered, **not** unit-variance scaled, `whiten` unset), `src/enrichment.py` (`distance_bins` equal-count; hypergeometric `sf(k-1)`; BH on the whole table; `sig_peak_at_bin0`).
- Official SCLC Panel A permutation (`run_panel_a_parti_full.py`) uses `delta=0`, `k-1` PCs, 150/50 inits, 1000 shuffles, `p = mean(finite null ≥ observed)`, `n_fail=0` for k=3–6. Saved t-ratios match Groves to ~3 decimals (`t_ratio_parti_1000.csv` vs `params_lognorm.txt`).
- Breast and GBM official Panel A also use `fit_pcha_best` / `hull_t_ratio` = `simplex_volume` / `ConvexHull` of data (same formula as SCLC official A, not `src.archetypes.t_ratio` which uses hull-of-vertices).
- `algNum` is not a Python argument; the algorithm is PCHA via `py_pcha` (ParTI `algNum=5` equivalent). `py_pcha.PCHA` default `delta=0`; wrappers that omit `delta` on **`fit_pcha`** still default to **0.1** (`src/archetypes.py:17`) — official A scripts pass `delta=0` explicitly.
- ComBat for breast/GBM Panel C: genes × samples; `ref.batch=cell_line` is an actual argument (`combat_cellline_tumor.R:14,43–44`; Python callers pass `"cell_line"`). Batch vectors are `cell_line`×n_lines then `tumor`×n_tumors (`run_panelC_tcga.py` breast:176, GBM:155). Groves tumor ComBat used `ref.batch='m'` (Minna) in `reference/notebooks/bulk/Cell-line-tumor-batch-correction-and-clustering.Rmd:205`.
- Archetypes are not refit on tumors in breast/GBM C (or SCLC C): inverse-PCA to genes, new combined PCA, `transform` vertices.
- Breast/GBM Panel B: distances in `k-1` columns; equal-count bins; BH within the panel table; bin-0-and-peak rule from `hypergeometric_enrichment`.
- SCLC Panel B **author** labels are `NEW_10_2020` from the Groves cluster CSV (`src/io.py:79–84`), not the TF-clustering branch (that branch is extra).
- Saved PAM50 counts on 63 lines: Basal 27, LumB 17, Her2 14, LumA 4, Normal 1 (`pam50_labels_panelA.csv`).
- PCA ≥50% rule on **saved** scree: breast 8 PCs = 50.01%; GBM 12 PCs = 51.38%. SCLC 12 PCs = **47.28%** (hardcoded 12, Groves `dim=12`).
- Permutation n_fail = 0 in all official t-ratio CSVs inspected (SCLC 1000; breast/GBM 500; GBM k=2 both variants 500).
- No TCGA patient double-counting in prepared primary tables: 1097 BRCA barcodes / 1097 patients, all `-01`; 154 GBM / 154 patients, all `-01`.
- No ACH–TCGA ID collision. SCLC cell-line vs tumor split is name prefix `m.`/`c.` vs `t.` (`src/io.py:53–61`); combined matrix 70 Minna + 50 CCLE + 81 Tumor.
- GBM k=2 1-PC vs 12-PC outputs live in different folders; official k=7 A/B/C do not read them.
- Low-expression filter for DepMap: `drop_low_genes` on **log2(TPM+1)** with `max >= 1` (`src/preprocess.py:92–99`), applied **after** subsetting samples, **no extra log**. Breast/GBM build scripts do not log again. Xena HiSeqV2 is not logged again (`prepare_tcga_*.py`).

---

## Top-line summary

| Phase | Cancer | Status | Key issues |
|---|---|---|---|
| 0 Inventory | SCLC | PASS (with stale-figure caveat) | Official A is `parti_full`; B/C scripts still wired to first-pass npy; figures 14 Aug vs code mtime 18 Aug (repo organize). |
| 0 Inventory | Breast | PASS | Separate drivers + shared `src/`; `archive/run1/` is not official. |
| 0 Inventory | GBM | PASS | Separate drivers + shared `src/`; entire GBM tree untracked at `HEAD`. |
| 1 Provenance | SCLC | CANNOT VERIFY (filters) / PASS (final shapes) | Author matrices 15950×120 and 14546×201; mito/NA/low-gene steps not in this repo. |
| 1 Provenance | Breast | PASS | 77 invasive BRCA cell lines in `Model.csv`; 63 in processed matrix; no 63-vs-65 split in this DepMap table. |
| 1 Provenance | GBM | PASS | 67 GB cell lines in metadata; 54 in processed matrix; 154 primaries / 152 with `GeneExp_Subtype`. |
| 2 Preprocess/PCA | SCLC | PASS (PCA) / CANNOT VERIFY (ComBat rebuild) | Uses author ComBat; 12 PCs = 47.3% not ≥50%. |
| 2 Preprocess/PCA | Breast | PASS | Log-scale gene filter; ComBat orientation and `ref.batch` correct; 8 PCs. |
| 2 Preprocess/PCA | GBM | PASS | Same as breast; 12 PCs. |
| 3 PCHA params | SCLC | PASS (A) / FAIL (B/C pointers) | A: 150/50/1000, δ=0, k−1. B/C: discarded 12-D file. |
| 3 PCHA params | Breast | FAIL (cross-cancer parity) | Hardcoded `NUM_ITER=5`; matches ParTI.m, **not** SCLC/GBM. |
| 3 PCHA params | GBM | PASS (k=3–7) / PASS isolation (k=2) | 150/50/500; k=2 1-PC vs 12-PC separated; B/C k=2 use 1-PC. |
| 4 Permutation | All | PASS (shuffle + formula) | Independent PC shuffles; same volume ratio; `n_fail` excluded from Python p; MATLAB uses `>` and `/maxRuns`. |
| 5 Panel B | SCLC | FAIL | Wrong npy + 12-D distances; author labels OK; bins=10 equal-count. |
| 5 Panel B | Breast | PASS | 3-D distances, 5 bins, PAM50 counts match save; LumA/Normal dropped from tests only. |
| 5 Panel B | GBM | PASS (mechanics) / FAIL (label meaning) | 6-D distances, 5 bins; labels are marker z-score argmax, not Verhaak. |
| 6 Panel C | SCLC | FAIL (containment) / PASS (no tumor PCHA) | No numeric containment; variance recipe present; 120:81 ≈ 60:40. |
| 6 Panel C | Breast | PASS | 0/1097 inside from saved flags; 63:1097 ≈ 5:95; k=4 path hardcoded. |
| 6 Panel C | GBM | PASS | 0/154 inside from saved flags (computed here from CSV, was missing as a quoted fraction); 54:154 ≈ 26:74. |

---

## Phase 0 — Inventory

### Shared vs separate implementations

**Finding: hybrid. Shared mathematical engines; separate per-cancer drivers.**

| Layer | Shared? | Evidence |
|---|---|---|
| PCA, PCHA multi-start, shuffle, enrichment, DepMap helpers, Python ComBat fallback | Yes | `src/{pca,archetypes,enrichment,preprocess,combat,io,paths}.py`; every disease driver does `sys.path.insert(0, ROOT)` and imports these. |
| Panel A/B/C orchestration, paths, `NUM_ITER`, k grid, tumor IDs, figures | No | Copied scripts: SCLC `codes/run_panel_*.py`; Breast `codes/run_panelA.py` etc.; GBM `codes/run_panelA.py` etc. |
| ComBat R wrapper | Duplicated | Breast and GBM `combat_cellline_tumor.R` are the same logic (GBM copies breast). |
| PAM50 / Verhaak | Disease-specific | `03_map_and_pam50.R` vs `02_assign_verhaak.py`. |

A bug fix in `src/archetypes.py` hits all three. A wrong npy path or `NUM_ITER` in a driver hits only that cancer. **Breast and GBM must be checked independently against the spec, not assumed identical because they import `src/`.**

### SCLC pipeline map (executed order)

```
Groves repo matrices (not rebuilt)
  data/SCLC_combined_Minna_CCLE_batch_corrected_wo_lowgenes.csv     # A/B
  data/combined_clusters_2020-05-27-MC copy.csv                    # B labels
  data/CCLE_Minna_Thomas_COMBAT.csv                                 # C

  run_panel_a_firstpass.py     → pc_scores_12.npy, archetypes_k{4,5,6}.npy  [NOT official]
  run_panel_a_sanity_parti.py  → archetypes_k*_parti.npy, S_*_parti.npy
  run_panel_a_parti_full.py    → null_t_ratios_*_n1000.npy, t_ratio_parti_1000.csv
  run_panel_b.py               → results/panel_b/*   [LOADS archetypes_k5.npy]
  run_panel_c.py               → results/panel_c/*   [LOADS archetypes_k5.npy]
  export_sclc_reproduction_figures.py → figures/Figure_1A.png (parti), 1B/1C (B/C CSVs)
  export_figure_1a_paper_style.py     → 1A only, parti file
```

Also present, not official: `run_panel_a_permutation.py` (loads `archetypes_k{k}.npy`).

### Breast pipeline map

```
data/raw/Model.csv + OmicsExpression…Logp1.csv
  → 01_build_input_panelA.py → data/processed/input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv (16500×63)
  → 03_map_and_pam50.R → pam50_genefu_raw.csv
  → 04_match_pam50_to_panelA.py → pam50_labels_panelA.csv
  → run_panelA.py → results/panel_a/*, figures/Figure_1A_breast.png
  → run_panelB.py → enrichment_pam50.csv, Figure_1B_breast.png
data/tumors HiSeqV2 + clinical
  → prepare_tcga_brca.py → results/panel_c_tcga/tcga_primary_*
  → run_panelC_tcga.py → ComBat, combined PCA, Figure_1C_tcga.png
```

`archive/run1/` is explicitly unofficial.

### GBM pipeline map

```
data/Model.csv + OmicsExpression… + HiSeqV2 + clinicalMatrix
  → 01_build_input_panelA.py → 15833×54 matrix
  → 02_assign_verhaak.py → verhaak_labels_panelA.csv
  → run_panelA.py → panel_a/*, Figure_1A_gbm.png  (k=3..7 parallel)
  → run_panelB.py → Figure_1B_gbm.png
  → prepare_tcga_gbm.py → 13677×154 tumors
  → run_panelC_tcga.py → Figure_1C_gbm.png

Side chain (not official Groves k):
  run_panelA_k2.py → panel_a_k2/ (2×1)
  run_panelB_k2.py / run_panelC_k2.py → panel_{b,c}_k2/
  run_panelA_k2_12.py → panel_a_k2_12/ (2×12)  # no B/C
```

### Timestamps vs reported numbers

| Artifact | mtime (local) | Notes |
|---|---|---|
| SCLC `t_ratio_parti_1000.csv` | 2026-08-14 21:25 | Reported A table |
| SCLC `figures/Figure_1A.png` | 2026-08-14 15:06 | |
| SCLC `Figure_1B.png`, `1C.png` | 2026-08-14 12:15 | **Before** the t-ratio CSV; B/C from first-pass geometry |
| SCLC `run_panel_{a_parti_full,b,c}.py` | 2026-08-18 11:11 | Organize commit; load paths still `k5.npy` for B/C |
| Breast `t_ratio_parti_500.csv`, Fig 1A | 2026-08-17 15:37 | Matches `NUM_ITER=5` columns in the CSV |
| Breast `run_panelA.py` | 2026-08-18 11:12 | Organize; `NUM_ITER=5` still in file |
| GBM `t_ratio_parti_500.csv`, Fig 1A | 2026-08-18 14:27 | After `HEAD`; **not in git** |
| `t_ratio_official.csv` | **missing** at path named in `run_panel_a_parti_full.py:36` (`SCLC/results/t_ratio_official.csv`); copy target exists only if last full run copied it. Table lives at `results/panel_a/t_ratio_parti_1000.csv`. |

**Stale-code vs numbers.** SCLC B/C **code still points at the same discarded file that produced the Aug 14 B/C figures** — not a later silent edit of the pointer. Breast numbers match current `NUM_ITER=5`. GBM results are reproducible from current untracked scripts only if those files stay as audited.

---

## Phase 1 — Raw data provenance and integrity

### SCLC

| Check | Status | Evidence |
|---|---|---|
| 120-line RNA-seq = CCLE + Minna, Groves processed | PASS (file-level) | `data/SCLC_combined_Minna_CCLE_batch_corrected_wo_lowgenes.csv` is **15950 genes × 120 samples**; columns start `m.DMS153`, `m.NCIH60`, … (`src/io.py` Minna=`m.`, CCLE=`c.`). README + `docs/02_how_we_implemented_this.md` cite Groves `data/bulk-rna-seq/` / dbGaP `phs001823` (not downloaded). Combined C matrix: **70 Minna + 50 CCLE + 81 Tumor**. |
| 81 George/Thomas tumors | PASS | `CCLE_Minna_Thomas_COMBAT.csv` **14546 × 201**; `sample_sources.csv` Tumor=81. |
| NA / mito / low-expression filter steps | CANNOT VERIFY | Not implemented here. Authors’ `*_wo_lowgenes` matrix is used as-is. Cannot confirm intermediate counts. |
| Final A matrix shape | PASS | 15950 × 120 as stated. |
| Gene IDs | PASS | HGNC-style symbols (`A1BG`, …), no Ensembl/Entrez in the A matrix. |
| Line vs tumor leakage | PASS | Distinct prefixes; tumors not in the 120-line A matrix. Training = cell lines only. |

### Breast

| Check | Status | Evidence |
|---|---|---|
| Filter logic | PASS | `pick_breast_models` (`src/preprocess.py:18–51`): OncoTree code `BRCA` **or** disease string containing “invasive breast” / “breast invasive”; `01_build_input_panelA.py:54–56` additionally requires `ModelType` matching `"cell line"`. |
| 63 lines | PASS vs **this** `Model.csv` | Metadata: **82** invasive models (77 cell line + 5 organoid). Invasive **cell lines = 77**. Processed `input_panelA_models_used.csv` = **63**. The 14 invasive CLs not in the used list: `ACH-001065, 001249, 001358, 001397, 001514, 002163, 002179, 002208, 002320, 002322, 002324, 002326, 002331, 002400`. Build report: `n_models_in_metadata: 77`, `n_with_rnaseq: 63`. **No 65** appears. A 63-vs-65 discrepancy is **not** present in the current DepMap table; it would require another release or a lineage-wide (not invasive) filter (Breast lineage CLs = 91). RNA absence of those 14 was **not** re-checked against Omics in this audit. |
| Genes after filter | PASS (artifact) | Build report 19205 → 16500; processed CSV **16500 × 63**. |
| Gene IDs / PAM50 mapping | PASS (code) | DepMap headers stripped `SYMBOL (Entrez)` (`preprocess.py:12–15,81`). `03_map_and_pam50.R:55–90` maps SYMBOL→ENTREZ via `org.Hs.eg.db`; unmapped logged (`n_unmapped 66` in `step0bc_report.txt`); duplicate Entrez collapsed by variance. Drops are Entrez-mapping, not silent format mismatch. Fresh genefu was **not** re-run; saved labels used. |
| Line vs tumor leakage | PASS | ACH- IDs vs TCGA barcodes. |
| TCGA double-count | PASS | 1097 unique `-01` primaries, 1097 patients (`prepare_tcga_brca.py:111–116` barcode `01` and `sample_type_id==1`). |

### GBM

| Check | Status | Evidence |
|---|---|---|
| Filter | PASS | `pick_gbm_models` (`preprocess.py:54–74`): OncoTree `GB`/`GBM` or disease/subtype containing “glioblastoma”; plus cell line (`01_build_input_panelA.py:48–52`). Not all CNS. |
| 54 lines | PASS | Metadata **67** GB cell lines; used **54**. 13 GB CLs not in used list (no RNA per build report): `ACH-001118, 001214, 001606, 002223, 002225, 002227–002231, 002259, 002268, 002349`. Report: 67 / 54 / 19205 → 15833. Processed CSV **15833 × 54**. Omics not re-scanned. |
| 154 tumors / 152 Verhaak | PASS | `prepare_report.txt`: 172 HiSeq samples → 154 primary; `GeneExp_Subtype` non-null **152**. Per-sample table: Mesenchymal 49, Classical 39, Proneural 38, Neural 26, nan 2. |
| Gene IDs | PASS | Symbols after Entrez strip, same DepMap convention. Verhaak markers matched by symbol; **ASCL1 missing** from matrix (see Phase 5). |
| Leakage / double-count | PASS | ACH vs TCGA; 154 unique patients, all `-01`. |

---

## Phase 2 — Preprocessing and normalization

### Order of log vs gene filter

| Cancer | Status | Evidence |
|---|---|---|
| SCLC | CANNOT VERIFY (author pipeline) / PASS (this repo) | No log or gene filter in Python. Docs: Groves dropped genes with all values &lt; 2 **then** log2. Not re-derived. |
| Breast | PASS | DepMap already log2(TPM+1). `drop_low_genes` uses `max(axis=1) >= 1` on that log matrix (`preprocess.py:92–99`; `01_build_input_panelA.py:74`). Threshold is log-space, equivalent to TPM&lt;1 in every line if the file is truly log2(TPM+1). |
| GBM | PASS | Same function and threshold (`01_build_input_panelA.py:70`). |

**Cross-cancer:** SCLC’s author threshold (linear &lt; 2 before log) is **not** the same as DepMap max log2(TPM+1)&lt;1. Silent gene-set difference vs Groves, but breast and GBM match each other.

### Batch correction

| Check | SCLC | Breast | GBM |
|---|---|---|---|
| `ref.batch` = cell-line batch | CANNOT VERIFY in Python (author `CCLE_Minna_Thomas_COMBAT.csv`). Groves Rmd uses `ref.batch='m'` (Minna), not a generic `cell_line` label. | PASS — `run_panelC_tcga.py:118,125`; R `ref.batch` arg default/`cell_line`. | PASS — `run_panelC_tcga.py:96,103`. |
| Batch labels | PASS prefixes in saved C matrix | PASS 63 `cell_line` + 1097 `tumor` (`combat_batch.csv`) | PASS 54 + 154 |
| Matrix orientation | Author matrix genes × samples | PASS `ComBat(dat=merged)` genes × samples (`combat_cellline_tumor.R:23,43`); Python `combat()` documents genes × samples (`src/combat.py:55–68`) and checks `dat.shape[1]==len(batch)` | PASS same R file |

### PCA fitting / scaling / ≥50% PCs

All three call `fit_pca` → `sklearn.decomposition.PCA(..., svd_solver="full")` (`src/pca.py:6–16`). sklearn **centers** features; **does not** divide by standard deviation (`whiten` default False). **PASS** same convention. Difference vs a scaled PCA would be intentional only if documented; it is not documented, but it **is** consistent.

| Cancer | n_pcs | Cumulative variance (saved `pca_variance.csv`) | How chosen | Status |
|---|---|---|---|---|
| SCLC | 12 | **47.28%** | Hardcoded `N_PCS=12` in `run_panel_a_firstpass.py:28`; `parti_full` **reuses** `pc_scores_12.npy` | PASS as Groves `dim=12`; **FAIL** if the spec is “first n with ≥50%” (12 never reaches 50%) |
| Breast | 8 | **50.01%** | `choose_n_pcs` in `run_panelA.py:56–60,128–131` | PASS |
| GBM | 12 | **51.38%** | same helper `run_panelA.py:54–58,328–331` | PASS |

PCA is fit on the Panel A expression matrix in the same script that writes `pc_scores*.npy` (not a stale other variable), except SCLC official permutation **reloads** the first-pass 12-PC scores — same matrix, first-pass PCA. **PASS** with that dependency noted.

SCLC Figure 1A **ESV curve** is `esv_curve.csv` from first-pass (`delta=0.1`, PCHA in **all 12 PCs**, `esv_curve` in `archetypes.py:85–101`). `parti_full` does not rewrite it. **FAIL** as a description of the official fit; the elbow panel is the discarded protocol.

---

## Phase 3 — PCHA / archetype fitting parameters

Python never sets MATLAB `algNum`; PCHA is selected by calling `py_pcha.PCHA`.

### Actual arguments (official Panel A)

| Cancer | k tested | alg | delta | PCs in fit | n_init obs | n_init null | n_perm | seed |
|---|---|---|---|---|---|---|---|---|
| SCLC | 3,4,5,6 | PCHA | 0.0 `run_panel_a_parti_full.py:48` | `k-1` via `fit_pcha_best` line 64 | 150 | 50 | 1000 | shuffle: `default_rng([0,k,i])` line 117; **PCHA inits unseeded** |
| Breast | 3–7 | PCHA | 0.0 line 52 | `k-1` | **15** | **5** | 500 | shuffle: `default_rng(SEED+1000+k)` line 197; PCHA unseeded |
| GBM | 3–7 | PCHA | 0.0 line 50 | `k-1` (npy shapes 3×2, 4×3, … 7×6) | 150 | 50 | 500 | shuffle: `default_rng([0,k,i])` line 155; PCHA unseeded |

`k-1` is **enforced in code** (`fit_pcha_best`: `scores = pc_scores[:, :k-1]`), not only in comments — except GBM `run_panelA_k2_12.py:82–91`, which **intentionally** does not slice (documented deviation).

### Open flag (a) — breast `numIter`

**Finding:** Hardcoded in `Breast Cancer/codes/run_panelA.py:47` (`NUM_ITER = 5`), docstring lines 6–9, figure title line 307. **Not** read from `py_pcha`. **Is** the MATLAB ParTI.m default when `algNum==5` (`reference/ParTI/ParTI.m:112–113`; `findMinSimplex.m` observed loop `3*numIter`; `CalculateSimplexTratiosPCHA.m` null loop `numIter`).

Breast write-up (`docs/PARAMS_vs_Groves.txt:13–15`) is **correct about ParTI.m** and **incorrect if it claims that is what the SCLC Python reproduction used**.

Per user instruction this audit did **not** change `NUM_ITER` or re-run Panel A. Whether t/p would move is **unknown** (PCHA is described as stable in ParTI.m, which is why they lowered `numIter`).

### Open flag (b) — SCLC file pointer

**FAIL, current code.** See Critical #1. Shapes: `archetypes_k5.npy` (5, 12) vs `archetypes_k5_parti.npy` (5, 4). First-pass t for k=5 was 0.092 vs paper 0.107 (`t_ratio_firstpass.csv`); parti t is 0.1075 (`t_ratio_parti_1000.csv`). B/C used the 0.092-class geometry.

`run_panel_c.py:56–64` inverse-transforms archetypes with a **12-component** PCA. A (5, 4) array would **error** in `inverse_transform` unless padded. Current C **requires** the 12-D file. That is independent proof C is not using the parti fit.

### Open flag (c) — GBM k=2 degeneracy

**PASS** on separate saves (Critical #5). Degenerate: `panel_a_k2/t_ratio_k2.csv` t=1.0, p=1.0, `n_pcs_fit=1`. Corrected: `panel_a_k2_12/t_ratio_k2_12pc.csv` t=0.704, p=0.49, `n_pcs_fit=12`. Official k=3–7 table does not include k=2. Downstream k=2 B/C **do** use the degenerate fit **by filename**, labeled “not Groves k” in titles.

### Observed-init success counts

| Run | Saved? | Status |
|---|---|---|
| SCLC parti | Printed only if refit (`observed_fit`); files reused if present | CANNOT VERIFY exact `n_ok/150` from artifacts |
| Breast | Printed `inits_ok=n_ok/15`; not in CSV | CANNOT VERIFY from files; perm `n_fail=0` |
| GBM k=3–7 | Same reuse pattern | CANNOT VERIFY `n_ok/150`; perm `n_success=500`, `n_fail=0` |
| Swallow | `fit_pcha_best` `continue` on `RuntimeError` (`archetypes.py:78–79`) | Failures drop that init; best-of-N is among successes. If all fail, raise. |

---

## Phase 4 — Permutation / significance test

### Shuffle

**PASS all official A scripts.** `_shuffle_columns` (`archetypes.py:119–123`) permutes **each PC column independently** across samples. Matches `CalculateSimplexTratiosPCHA.m:32–35` (`randperm` per axis). Input to shuffle is `scores[:, :k-1]` (SCLC line 117, breast 200, GBM 155).

Not used for official p-values: gene-row shuffle (that is Panel C variance **null PCA**, different question).

### t-ratio formula

Official A: `hull_t_ratio` = `simplex_volume(archetypes) / ConvexHull(data).volume` with data in `k-1` D (SCLC `run_panel_a_parti_full.py:52–55`; breast 84–87; GBM 76–79). `simplex_volume` (`archetypes.py:54–59`) matches ParTI `abs(det)/ (k-1)!`.

**Not used for official tables:** `src.archetypes.t_ratio` (ConvexHull of vertices / hull of data) — first-pass only.

GBM k=2 1-D: segment length / (max−min) (`run_panelA_k2.py:92–104`). k=2 12-PC: pair length / axis extent (`run_panelA_k2_12.py:75–79`) — **different statistic**, documented.

### p-value

Python official: `p = mean(finite_null >= observed)` (`parti_full.py:144`; breast 225; GBM 182). Failed shuffles stored as NaN and **excluded from the denominator** (`finite = null[isfinite]`). With `n_fail=0`, denominator = `n_perm`.

MATLAB ParTI (`findArchetypes.m:214–215`): `sum(tRatioRand > tRatioReal) / maxRuns` — **strict `>`**, denominator **always** `maxRuns` (failed/NaN comparisons are false).

**Status:** PASS as internally consistent in Python across three cancers. CANNOT VERIFY bit-identity with MATLAB p due to `>=` vs `>` (immaterial if no exact ties) and unseeded PCHA.

**Reproducibility:** shuffle RNGs are seeded; **PCHA random inits are not**. Exact p-values cannot be bit-reproduced. Flag for audit credibility, not a logic bug.

---

## Phase 5 — Subtype enrichment (Panel B)

### Distance space

| Cancer | Fit space | Distance space in B | Status |
|---|---|---|---|
| SCLC k=5 official A | 4 D | **12 D** to **12-D discarded** vertices | FAIL |
| Breast k=4 | 3 D | `scores[:, :n_vol]` with `n_vol=3` (`run_panelB.py:106–129`) | PASS |
| GBM k=7 | 6 D | `scores[:, :n_vol]` (`run_panelB.py:103–125`) | PASS |
| GBM k=2 B | 1 D | first 1 PC (`run_panelB_k2.py:106–129`) | PASS relative to 1-PC fit; not the 12-PC k=2 fit |

### Bins

`distance_bins` (`enrichment.py:93–117`): equal-**count** via argsort ranks and cumulative sizes (`base, rem = divmod`; leftover to first bins). **PASS** implementation.

| Cancer | n_bins | Evidence |
|---|---|---|
| SCLC | 10 | `run_panel_b.py:76`; `bins_author.csv` values 0–9 |
| Breast | 5 | `N_BINS=5` (`run_panelB.py:33`); bins 0–4 |
| GBM official | 5 | `N_BINS=5` (`run_panelB.py:33`); 7 arc columns |

### Hypergeometric + BH

`hypergeometric_enrichment` (`enrichment.py:120–170`): one table, `q_value = BH(all p in table)`, `significant = q < fdr` (0.1), `sig_peak_at_bin0` requires significant **and** `bin==0` **and** `peak_bin==0` (peak = max fold among that subtype×archetype).

FDR is **within that function call**, not across cancers. Breast passes `subtype_levels=KEEP` after dropping LumA/Normal — 3×4×5 = 60 tests (`enrichment_pam50.csv`). GBM 3×7×5 = 105. SCLC uses default `SUBTYPES` A/A2/N/P/Y × 5 arcs × 10 bins = 250. **PASS** mechanics. Breast’s drop changes the test universe (intentional; documented `run_panelB.py:1–6`).

Match rule is in shared code, not prose-only. **PASS** for breast/GBM; SCLC uses the same function on the **wrong distances**.

### Subtype label provenance

**SCLC — PASS.** `load_author_subtypes(..., column="NEW_10_2020")`. Figure 1B uses `enrichment_author.csv`. Extra TF clustering (`spearman_average_clusters`) is computed but is not the official B plot. Labels are Wooten-style published names in the Groves CSV, not a new clustering for the reported panel.

**Breast — PASS on saved counts; CANNOT VERIFY a fresh genefu run.** Pipeline: symbols → Entrez → `genefu::molecular.subtyping(sbt.model="pam50")` (`03_map_and_pam50.R:122–128`); match to ACH IDs (`04_match_pam50_to_panelA.py`). Saved: Basal 27, LumB 17, Her2 14, LumA 4, Normal 1 (63/63). `step0bc_report.txt`: 0 failed calls. This audit did not re-execute R.

**GBM — method (actual code), not Verhaak.**

1. Marker lists (`02_assign_verhaak.py:26–30`):  
   Classical: EGFR, NES, NOTCH3, SMO, GAS1, GLI2, NFKBIA  
   Mesenchymal: CD44, CHI3L1, RELB, TRADD, TNFRSF1A, NFKB1, STAT3  
   Proneural: OLIG2, PDGFRA, DLL3, NKX2-2, SOX10, ASCL1, NKX2-1  
   Compact “Wang 2017 three-type” style; **Neural not scored**.
2. For each gene, z-score **across the 54 cell lines** (not tumors): `(x - mean_line) / sd_line` (`line 36`).
3. Subtype score = **mean z of markers present** in the matrix (`lines 42–45`). Missing markers are skipped, not imputed.
4. Assignment = **argmax** of the three scores (`line 46`). No threshold, no ssGSEA, no TCGA centroid, no SVM.
5. **ASCL1 is not in the Panel A CSV**; Proneural uses 6 genes. That is a silent panel change.

Saved assignments: Proneural 21, Mesenchymal 20, Classical 13 (n=54). TCGA Panel C uses clinical `GeneExp_Subtype` (four classes including Neural) — **different label source than Panel B**.

---

## Phase 6 — Tumor projection / containment (Panel C)

### No PCHA on tumors

| Cancer | Status | Evidence |
|---|---|---|
| SCLC | PASS | Inverse 12-PC cell-line PCA (`run_panel_c.py:56–67`); new PCA on author combined matrix; `pca_all.transform(gene_arcs)` line 86. No `fit_pcha` in this file. **But gene-space vertices come from discarded 12-D archetypes.** |
| Breast | PASS | Pad k−1 → n_pcs, inverse, ComBat, new PCA, transform (`run_panelC_tcga.py:148–191`). Prints “PCHA is NOT refit”. Hardcodes `archetypes_k4_parti.npy` (line 42) rather than `suggested_k.txt`; currently k=4 so consistent, **fragile**. |
| GBM | PASS | Reads `suggested_k.txt` (line 119); same inverse/ComBat/new PCA path. |

### Containment (barycentric w ≥ 0)

Formula (breast/GBM): affine `lstsq` in combined-PCA **first `n_vol` PCs**, `inside = (weights >= -1e-6).all(axis=1)`.

| Cancer | Reported | This audit (saved `inside_simplex`, no refit) | Status |
|---|---|---|---|
| SCLC | “mostly contained” (visual) | **Not computed in code** | FAIL as a number |
| Breast | 0/1097 | **0 True / 1097** in `tcga_sample_subtype_archetype.csv` | PASS vs saved output |
| GBM official k=7 | “0/154” in PROJECT_CONTEXT; not a figure annotation | **0 True / 154** in `tcga_sample_subtype_archetype.csv` | PASS; **fraction to quote: 0/154 = 0%** |
| GBM k=2 C | — | **152/152** labeled tumors inside (`panel_c_k2/subtype_inside_simplex_counts.csv`) | Degenerate 1-D; do not cite as Groves containment |

### Combined vs tumor-only variance ceiling

Same recipe: tumor-only PCA on ComBat tumor columns vs variance of those tumors on combined PCA axes; ratio of **cumulative** % (`ratio_combined_over_tumor`). Gene-row shuffles for a null curve (SCLC/breast/GBM C).

| Cancer | n_lines : n_tumors | Ratio | Combined/tumor-only at k PCs (saved CSV) |
|---|---|---|---|
| SCLC | 120 : 81 (**59.7% : 40.3%**) | ~60:40 | 5 PCs: **80.1%** (`variance_explained.csv` n_components=5) |
| Breast | 63 : 1097 (**5.4% : 94.6%**) | ~5:95 | 4 PCs: **99.92%** |
| GBM | 54 : 154 (**26.0% : 74.0%**) | ~26:74 | 7 PCs: **97.2%** |

**These ceiling percentages are not comparable across cancers.** Tumors dominate the combined PCA wherever they outnumber lines (breast almost entirely; GBM majority; SCLC more balanced). High breast/GBM ratios do **not** mean tumors lie in the cell-line simplex (containment is 0%).

---

## Recommended next actions (no work done here)

1. Retarget SCLC `run_panel_b.py` / `run_panel_c.py` to `archetypes_k5_parti.npy`, slice/pad PCs correctly, regenerate B/C and Figures 1B–1C; treat current 1-to-1 B map as **untrusted** until then.
2. Add SCLC barycentric containment; report the fraction.
3. Decide a single `numIter` (ParTI.m PCHA=5 vs SCLC-matched 50) and apply it to all three; only then compare p-values.
4. In GBM prose, never call cell-line labels “Verhaak”; state marker z-score argmax and missing ASCL1.
5. Commit or freeze GBM artifacts so the 18 Aug numbers stay attached to a git SHA.
6. Stop hardcoding breast `archetypes_k4_parti.npy`; read `suggested_k.txt` like GBM.

---

*Audit method: line-level reads of drivers and `src/`, npy headers, CSV counts, `Model.csv` filters, `reference/ParTI/*.m`, Groves ComBat Rmd. No model refits.*
