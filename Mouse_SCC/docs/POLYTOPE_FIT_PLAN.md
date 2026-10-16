# Mouse skin SCC polytope plan (GSE288137)

**Status:** protocol only. No fit has been run.  
**Question:** On the Pérez-González et al. 2026 scRNA-seq atlas (GEO **GSE288137**), does gene-expression geometry form a statistically supported simplex, and if so, where do the paper’s four EMT tumor states sit on that simplex?

**Paper:** Pérez-González, Windels, et al., *Nat. Commun.* (2026) 17:8561, doi:10.1038/s41467-026-75521-8 (`papers/s41467-026-75521-8.pdf`).

This document is the working protocol. It follows **one** calling convention: ParTI PCHA as used for Groves 2022 and Hausser 2019 in this repository (`src/archetypes.py`, `src/pca.py`, breast/SCLC Panel A). EMT labels from the paper are **not** used to choose genes, PCs, or \(k\).

---

## 0. What we are doing (and what we are not)

### We are doing two stages, in this order

**Stage A — Fit.**  
Build a cells × genes matrix from GSE288137. Run PCA. For several \(k\), run PCHA, record explained sample variance (ESV), t-ratio, and a permutation \(p\)-value. Choose a single \(k^\star\). Freeze the PCA, the \(k^\star\) vertices, and the barycentric weights.

**Stage B — Label.**  
Take those **same cells** and the **already computed** weights. Attach the paper’s `emt_state` and `tumor` IDs. Make tables and PC1–PC2 plots. Do not refit.

Stage B is **not** a guest projection into a breast Pal host. The mouse cells *are* the host. “Projecting the dataset into it” only means: the cells already have PC scores and weights from Stage A; we look at those numbers with the paper’s metadata.

### We are not doing in this plan

- Fitting on GSE288220 (multiome) or GSE288196 (bulk).
- MAGIC imputation.
- Harmony / Seurat integration **as the PCHA input** (UMAP in the paper used integration; the simplex is fit on expression).
- Tan 2014 KS 206-gene list, human Pal triangle, or mouse→human ortholog mapping.
- Blais 2024 SDVMM / raw TPM protocol (different fitter, bulk lymphoma). Mentioned only as something we are not copying.
- Choosing \(k=4\) because the paper named four EMT states.

### Scientific decision rule (written down before looking at labels)

\[
k^\star = \min\{ k \in \{3,4,5,6,7\} : p(k) < 0.01 \}
\]

If no \(k\) has \(p < 0.01\), report the smallest \(k\) with Groves’ looser \(p < 0.05\) as a **secondary** line, and say so. If several \(k\) are significant, take the **smallest**. If the ESV elbow (DimensionFinder) disagrees with \(k^\star\), report both; **\(k^\star\) from \(p\) is the official call.**

---

## 1. Folder layout (create before any download)

All new work lives under `Mouse_SCC/` at the repo root, parallel to `Breast Cancer/` and `SCLC Reproduction - Groves Cell Systems 2022/`.

```
Mouse_SCC/
  docs/
    POLYTOPE_FIT_PLAN.md     ← this file
  data/
    raw/                     ← GEO downloads (gitignored if large)
    processed/               ← frozen analysis matrix + metadata
  codes/                     ← scripts to be written from this plan
  results/
    panel_a/                 ← PCA, ESV, t-ratio, archetypes, weights
    labels/                  ← EMT × vertex tables
  figures/
```

Add `Mouse_SCC/data/raw/` to `.gitignore` if count matrices are large. Keep small CSVs (metadata, t-ratio table) tracked.

**Engine (do not fork a second PCHA):**

| Piece | Path |
|--------|------|
| PCA | `src/pca.py` → `fit_pca`, `cumulative_variance` |
| PCHA, volume, ESV | `src/archetypes.py` → `fit_pcha`, `fit_pcha_best`, `simplex_volume`, `esv_curve` |
| t-ratio (report) | `simplex_volume(arcs) / ConvexHull(data[:, :k-1]).volume` as in breast `run_panelA.py` `hull_t_ratio` (ParTI MATLAB). Do **not** use `src.archetypes.t_ratio` for the official number: that function uses `ConvexHull` for both simplex and data, which is a slightly different volume formula. |
| PC shuffle | `src.archetypes._shuffle_columns` |
| Barycentric weights | same least-squares as `Breast Cancer/codes/host_utils.py` `barycentric` |

**Environment**

- Python from repo `.venv`.
- **NumPy &lt; 2** (`py_pcha` still uses `np.mat`). If NumPy 2 is installed, shim `np.mat = np.asmatrix` as in `host_utils.py`, but pinning is better.
- Packages: `numpy`, `pandas`, `scipy`, `scikit-learn`, `py_pcha`, `matplotlib`, `scanpy`, `anndata`. Optional: `geoparse` / `requests` for GEO.
- Seed: **`SEED = 0`** everywhere (PCA `random_state`, `np.random.seed` before PCHA, permutation RNG).
- Before PCHA: `np.random.seed(SEED)` because `py_pcha` uses the **global** NumPy RNG (`host_utils.fit_host_matrix` already does this).

---

## 2. Stage A0 — Get GSE288137 and the paper labels

### 2.1 What GEO contains

Accession: **GSE288137** (10x scRNA-seq of primary mouse skin SCC).  
The paper: 12 primary tumors from `Lgr5-CreER; KrasG12D; Trp53cKO; Rosa-YFP`; tumors with 1–65% EPCAM⁺; **54,081 YFP⁺ tumor cells** after QC; unsupervised clustering gave four EMT states.

GEO may ship:

- per-sample 10x `matrix.mtx` / `barcodes.tsv` / `features.tsv`, or
- HDF5, or
- a processed Seurat/h5ad in supplementary files.

**Labels (`emt_state`, `tumor`) are not guaranteed in the count matrix.** Without them Stage B is impossible. Search, in this order:

1. GEO supplementary: cell metadata, cluster IDs, RDS/h5ad.
2. Paper “Source data” / Nature Communications supplementary tables (cell barcodes × state).
3. Author-processed object (email if needed).
4. Last resort: re-derive states from **the paper’s marker rules** on a reprocessed atlas (Section 2.4). Treat that as a weaker Stage B.

### 2.2 Download (minute steps)

1. Open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288137.
2. Record: Series matrix, supplementary file list, number of GSM samples, genome (mm10).
3. Download into `Mouse_SCC/data/raw/GSE288137/`.
4. Write `Mouse_SCC/data/raw/DOWNLOAD_LOG.md`: date, URLs, file names, md5 if available.
5. If files are 10x per tumor: one folder per GSM, each with matrix/barcodes/features.

### 2.3 Identify tumor cells (if starting from counts, not an author object)

The paper distinguished tumor cells from untransformed keratinocytes using:

- YFP (Rosa reporter),
- elevated CNV (inferCNV),
- *Trp53* exon 2–10 deletion,
- *Kras*G12D.

**Preferred:** use the author’s tumor-cell barcode list so the polytope is fit on the **same 54,081 cells**.

**If rebuilding QC from 10x (only if author barcodes are missing):** follow their scRNA methods as closely as possible (Scanpy, SoupX ambient RNA, scDblFinder doublets, then the genetic/transcriptional tumor filter). Document every threshold. Do not invent a new EPCAM-only mask.

Do **not** throw away EPCAM⁻ tumor cells: late EMT is EPCAM⁻.

### 2.4 Attach `emt_state` and `tumor`

Create `Mouse_SCC/data/processed/cell_metadata.csv`:

| column | meaning |
|--------|---------|
| `barcode` | unique cell ID (prefix with tumor/GSM if barcodes collide) |
| `tumor` | one of the 12 primary tumors |
| `gsm` | GEO sample ID if known |
| `emt_state` | `epithelial` / `early_hybrid` / `hybrid` / `late` |
| `early_hybrid_sub` | `I` / `II` / `NA` if available (multiome had I vs II; RNA atlas may not) |
| `label_source` | `author` or `marker_rule` |

**Author states (Fig. 1 / immunostaining logic):**

| State | RNA (paper) |
|--------|-------------|
| Epithelial | Epcam⁺, Cdh1⁺, Krt14⁺, Esrp1⁺, S100a14⁺; Vim⁻ |
| Early hybrid | Epcam/Cdh1 lost; Krt14, Krt5, Itga6 retained; Vim⁺ |
| Hybrid EMT | epidermal keratins lost; Krt8⁺; mesenchymal genes higher |
| Late EMT | Col3a1, Pdgfrb, Aspn, Postn, fibroblast-like |

If labels are Leiden clusters in an h5ad, map cluster → state using the paper’s marker UMAPs / supplementary, and save the mapping table. Do not rename clusters after seeing PCHA vertices.

**Freeze metadata now.** Stage A scripts must load expression only; they may store barcode order but must not read `emt_state` for filtering or \(k\).

### 2.5 Sanity on metadata (before PCA)

Print and save:

- \(n\) cells, \(n\) tumors (expect 12),
- cells per tumor (watch one tumor dominating, Pal TN-0135 lesson),
- cells per `emt_state`,
- fraction of barcodes in the expression matrix that match metadata (must be 100% of cells used in PCA).

If one tumor is &gt;15–20% of cells, **do not change the primary fit yet**; record it. A later sensitivity (Section 10) can cap per tumor. The official first fit is **all QC tumor cells**, to match the paper’s atlas.

---

## 3. Stage A1 — Expression matrix for PCA/PCHA

### 3.1 Load

1. Read counts into AnnData: cells × genes, gene symbols as `var_names` (mouse: *Epcam*, not *EPCAM*).
2. Restrict to barcodes in `cell_metadata.csv`.
3. Same gene world for every cell (inner join of genes across samples if samples were separate 10x runs).
4. Save raw counts: `Mouse_SCC/data/processed/counts_tumor.h5ad`.

### 3.2 Normalize (matches the paper’s RNA atlas language: log-normalized)

For each cell:

1. Library-size normalize to **10,000** counts per cell (scanpy `normalize_total(target_sum=1e4)`).
2. \(\log(x+1)\) natural log (`scanpy.pp.log1p`).

Do **not** z-score genes (`scale`). Hausser 2019: mean-center only; “We did not scale log2 fold changes by the standard deviation prior to PCA.” `sklearn.PCA` mean-centers features; that is enough.

Do **not** MAGIC.  
Do **not** regress cell cycle / mito for the **primary** matrix (the paper noted cycling is similar except late EMT; residualizing would mix biology into the residual). Optional sensitivity later.

Do **not** Harmony-correct the matrix that goes into PCHA.

### 3.3 Gene filter (Hausser used all genes; scRNA has too many zeros)

Primary filter (pre-registered):

1. Remove genes detected (count &gt; 0) in **&lt; 1% of cells**.
2. Remove mitochondrial genes (`mt-` / `Mt-`) and ribosomal (`Rpl`, `Rps`, `Mrpl`, `Mrps`) **only if** they would otherwise dominate PC1. Record the list. Default: **drop mito + ribo** (standard scRNA; Hausser bulk did not have this problem).
3. Keep all remaining genes. **No** “top 2000 HVG” for the official fit (that would be a different analysis). Hausser ran PCA on the full mean-centered matrix.

Save:

- `Mouse_SCC/data/processed/expr_log_filtered.h5ad` (or genes × cells CSV/npz),
- `Mouse_SCC/data/processed/genes_kept.txt`,
- a one-line report: \(n_{\text{cells}}\), \(n_{\text{genes}}\), % zeros.

### 3.4 Export the numeric matrix PCHA will see

Let \(X\) be **cells × genes**, `float64`, no NaNs.

- Row order = barcode order in metadata (stable).
- Column order = gene list (stable).
- Save `Mouse_SCC/data/processed/X_cells_by_genes.npy` and `barcodes.txt`, `genes.txt`.

This file is the analysis. Changing it later means a new run, not a continuation.

---

## 4. Stage A2 — PCA

### 4.1 Fit

```text
X  : (n_cells, n_genes)
pca, scores = fit_pca(X, n_components=n_pcs_store, random_state=0)
```

`src/pca.py`: `sklearn.decomposition.PCA(svd_solver="full", random_state=0)`.

This mean-centers genes. It does **not** divide by gene SD.

### 4.2 How many PCs to **store**

1. Fit a generous PCA (e.g. 50 components) once, compute cumulative variance `cumulative_variance(pca)`.
2. Let \(n_{50}\) = smallest \(n\) with cumulative variance **≥ 0.50**.  
   (Groves SCLC used 12 PCs because that was ~47–50%; breast Panel A uses this same 50% rule on each matrix.)
3. Set  
   \(n_{\text{store}} = \max(n_{50}, 6)\).  
   Six PCs is enough to try \(k\) up to 7 (needs 6 dimensions).
4. If \(n_{50}\) is huge (e.g. &gt; 30), still store \(n_{50}\) for the ESV curve, but **PCHA/t-ratio never use more than \(k-1\) PCs**.

Save:

- `results/panel_a/pca_variance.csv` (PC index, % var, cumulative),
- `results/panel_a/n_pcs.txt` (the integer \(n_{\text{store}}\)),
- `results/panel_a/pc_scores.npy` shape `(n_cells, n_store)`,
- `results/panel_a/pca_mean.npy`, `pca_components.npy` (to freeze the map),
- printed line: `var_explained_first_2`, `var_explained_first_(k-1)` for each k you will test.

### 4.3 Sign convention (reproducibility)

Flip each PC so the loading with largest absolute value is **positive** (`host_utils.fit_host_matrix`). Save scores **after** this flip. Plots then survive reruns.

### 4.4 What PCA is not

PC1–PC2 is only a picture. A \(k=4\) tetrahedron is 3-D. Fitting in 2 PCs for every \(k\) is wrong. Fitting in all 12–50 PCs for t-ratio is also wrong (that was the discarded SCLC first pass).

---

## 5. Stage A3 — ESV elbow (DimensionFinder)

### 5.1 What ESV is

For each \(k\), PCHA returns `varexpl`: fraction of variance in the **fitting space** explained by reconstructing points as convex combinations of \(k\) archetypes. Groves/ParTI plot this vs \(k\) and run DimensionFinder.

### 5.2 Curve

1. Take `pc_scores` with \(n_{\text{store}}\) columns (ESV in Groves is computed in the stored PCA space in `findArchetypes.m`; breast `esv_curve` uses the scores you pass in).
2. For \(k = 3, 4, \ldots, 15\) (or up to \(\min(15, n_{\text{store}}+1)\)):
   - `fit_pcha(pc_scores, k, delta=0.0)` **or** `fit_pcha_best` on `pc_scores[:, :k-1]` with 15 inits.
   - Official for the **elbow plot**: one PCHA in the stored PC space is what ParTI’s ESV loop did; for consistency with t-ratio, **prefer `fit_pcha_best` on the first \(k-1\) PCs, 15 inits, `delta=0`**, record `varexpl`.
3. Skip \(k=2\) (`py_pcha` can hang in 1-D; t-ratio undefined).
4. Write `results/panel_a/esv_curve.csv`: `k, esv, n_ok, volume`.

### 5.3 DimensionFinder (copy breast `run_panelA.dimension_finder`)

ParTI `DimensionFinder.m`: farthest point from the chord of the ESV curve.

If `esv[0]` is \(k=3\) (we skipped \(k=2\)), **do not** blindly add `+2` as in the MATLAB code that assumed `TotESV1(1)=k=2`. Implement:

- x-axis = the actual \(k\) values,
- chord from first to last point,
- \(k_{\text{elbow}}\) = \(k\) with maximum perpendicular distance to that chord.

Save `suggested_k_elbow` in `results/panel_a/best_k.json` as a **hint only**.

Plot: ESV vs \(k\), vertical line at \(k_{\text{elbow}}\).

---

## 6. Stage A4 — Observed PCHA for each \(k\) (the real fit)

### 6.1 Protocol table (do not change between \(k\))

| Parameter | Value | Source |
|-----------|--------|--------|
| Algorithm | PCHA (`py_pcha` via `fit_pcha_best`) | Groves `algNum=5` |
| Fitting coordinates | `scores[:, : k-1]` | Groves, Hausser, `fit_pcha_best` |
| `delta` | **0.0** | vertices = convex combinations of observed cells |
| `conv_crit` | `1e-6` | `py_pcha` / ParTI default |
| `maxiter` | 500 | same |
| Observed inits | **15** | ParTI.m PCHA: `numIter=5`, observed loop `3*numIter` |
| Keep | simplex with **maximum** `simplex_volume` | ParTI `findMinSimplex.m` |
| Volume formula | \(\lvert \det(A_{:,1:k-1} - A_{:,k})\rvert / (k-1)!\) | `simplex_volume` |
| \(k\) grid | **3, 4, 5, 6, 7** | Groves/breast Panel A; Hausser only 3–5 on bulk |

Breast `PARAMS_vs_Groves.txt` is the reference: Groves notebooks sometimes used a larger `numIter`, but ParTI’s PCHA branch sets `numIter=5`. We follow that PCHA branch (15 / 5), which is also breast Panel A.

### 6.2 Per-\(k\) loop (exact)

For \(k\) in `(3, 4, 5, 6, 7)`:

1. `np.random.seed(SEED)` (or `SEED + k` if you must avoid identical inits; then **record which**). Prefer `SEED` plus documenting that `py_pcha` draws internally.
2. `arcs, weights, varexpl, vol, n_ok = fit_pcha_best(scores, k, n_init=15, delta=0.0)`.
3. Abort that \(k\) if `n_ok < 15` only if **zero** inits worked; if some failed, still keep max volume among successes and record `n_ok`.
4. `data = scores[:, :k-1]`.
5. `t_obs = vol / ConvexHull(data).volume`.
6. Save:
   - `archetypes_k{k}.npy` shape `(k, k-1)`,
   - `weights_k{k}.csv` cells × \(k\) (rows sum ≈ 1),
   - print `varexpl`, `t_obs`, `n_ok`.

`fit_pcha_best` already slices to \(k-1\) PCs. Passing full `scores` is correct.

### 6.3 Mixture weights from PCHA vs barycentric

PCHA returns `S` (mixture weights). For plotting and “inside the simplex”, also compute **barycentric** least-squares in the \(k-1\) plane (`host_utils.barycentric`), which can go slightly negative for points outside.

Save both:

- `w_pcha_*` from PCHA,
- `w_bary_*` from `barycentric`.

**Official “inside”:** all barycentric weights \(\ge -10^{-6}\) (`host_utils.describe_scores`).  
**Official nearest vertex:** \(\arg\max\) of clipped-renormalized barycentric weights (`composition` then `argmax`).

---

## 7. Stage A5 — Permutation t-ratio and \(p\)

### 7.1 Null

For a \(k\)-simplex in \((k-1)\)-D, independently permute each of the \(k-1\) PC columns across cells (`_shuffle_columns`). That destroys correlation among PCs while keeping each PC’s marginal distribution. Refit PCHA on the shuffled cloud. If the real data are a filled simplex, real t-ratio should exceed most shuffled t-ratios.

### 7.2 Settings

| Parameter | Value |
|-----------|--------|
| Shuffles (debug) | 100 |
| Shuffles (reported) | **1000** (Hausser Methods; Groves) |
| Inits per shuffle | **5**, keep **max t-ratio** among those 5 | ParTI PCHA null |
| `delta` | 0 |
| \(p\) | \(\#\{\, t_{\text{null}} \ge t_{\text{obs}} \,\} / n_{\text{success}}\) |

ParTI MATLAB PCHA used **strict** `>`; we use **`≥`** as in `permutation_t_ratio` and breast/SCLC Python. With 1000 shuffles the difference is at most one count. State this in the methods sentence.

Failed shuffles (Qhull / PCHA init) are **dropped from the denominator** (repo convention). Record `n_fail`.

### 7.3 Do not call `permutation_t_ratio` as-is for the official \(p\)

`src.archetypes.permutation_t_ratio` currently does **one** `fit_pcha` per shuffle. That is the discarded SCLC first-pass null (too small \(p\)).

Official null, for each shuffle \(i=1\ldots 1000\):

1. `shuffled = _shuffle_columns(scores[:, :k-1], rng)`.
2. `arcs_n, _, _, vol_n, _ = fit_pcha_best(shuffled, k, n_init=5, delta=0.0)` — pass already-sliced array; `fit_pcha_best` will slice `[:, :k-1]` again, which is a no-op if width is already \(k-1\).
3. `t_n = vol_n / ConvexHull(shuffled).volume`.
4. Append `t_n`.

RNG: `np.random.default_rng([SEED, k, i])` per shuffle (SCLC `run_panel_a_parti_full` style) so reruns match.

### 7.4 Compute cost (54k cells)

1000 × 5 PCHA fits × 5 values of \(k\) on 54k points is heavy.

**Official observed fit:** all cells.  
**Official permutation:** same cells if it finishes in reasonable time.

If a 10-shuffle pilot is already too slow:

1. Draw a **tumor-balanced** subset once (e.g. at most 1,500 cells per tumor, `SEED=0`): ~12 × 1500 = 18,000 cells, or all cells from small tumors.
2. Recompute PCA **on that subset** (new \(n_{50}\), new scores).  
   **Do not** permute the full-data PCs after subsetting rows; PCs would be wrong.
3. Run the full A2–A5 pipeline on the subset and treat it as a **sensitivity**, or, if full-data permutation is impossible, declare the subset the official matrix **before** looking at EMT labels.

Never mix: full-data \(t_{\text{obs}}\) with subset-data null.

Save `results/panel_a/null_t_ratios_k{k}.npy` and a histogram figure per \(k\).

### 7.5 t-ratio summary table

`results/panel_a/t_ratio_summary.csv`:

`k, n_cells, n_genes, n_pcs_stored, n_pcs_fit, esv, t_ratio, p_value, n_perm, n_success, n_fail, n_ok_obs`

This is the table that chooses \(k^\star\).

---

## 8. Stage A6 — Choose \(k^\star\) and freeze

### 8.1 Rules (apply in order)

1. From `t_ratio_summary.csv`, list \(k\) with `p_value < 0.01`.
2. \(k^\star\) = minimum of that list.
3. If the list is empty: list \(k\) with `p_value < 0.05`; if nonempty, \(k_{\text{secondary}}\) = min of those; write that Hausser’s 0.01 rule failed.
4. If still empty: **no significant polyhedron**. Still freeze \(k=3\) only for exploratory plots, and say the geometry claim did not pass.
5. Write `k_elbow` from DimensionFinder next to \(k^\star\). Do not override (4) with the elbow.

`results/panel_a/best_k.json`:

```json
{
  "k_star": 4,
  "rule": "smallest k with p < 0.01",
  "k_elbow": 3,
  "alpha_hausser": 0.01,
  "alpha_groves": 0.05,
  "n_pcs_stored": 8,
  "n_cells": 54081,
  "n_genes": null,
  "delta": 0.0,
  "n_init_obs": 15,
  "n_init_null": 5,
  "n_perm": 1000,
  "seed": 0
}
```

### 8.2 Freeze the host (needed even though we are not projecting a guest)

Write `results/panel_a/host_k{k_star}.npz` containing:

- `genes`, `barcodes`,
- `pca_mean`, `components`, `n_pcs`,
- `archetypes` \((k^\star, k^\star-1)\),
- `train_scores`,
- `weights_pcha`, `weights_bary`,
- `t_ratio`, `p_value`, `esv`,
- `seed`, `delta`, `n_init`.

After this file exists, Stage B only **reads** it.

---

## 9. Stage B — EMT and tumor labels on the frozen simplex

### 9.1 Join

Inner-join `weights` / barycentric table to `cell_metadata.csv` on `barcode`. Every analysis cell must have a state. Drop cells without labels and record the count (should be 0).

### 9.2 Per-cell columns to compute

From barycentric weights \(w\) (clip at 0, renormalize = `composition`):

| column | definition |
|--------|------------|
| `w1` … `wk` | raw barycentric (can be slightly negative) |
| `inside` | all raw \(w_j \ge -10^{-6}\) |
| `nearest` | \(1\ldots k\), `argmax` of clipped \(w\) |
| `purity` | largest clipped weight − second largest |
| `specialist` | `purity >= 0.35` (Pal convention; state the cut) |
| `emt_state` | paper label |
| `tumor` | paper / GEO tumor ID |

Optional: Shannon evenness of clipped weights, so the 0.35 cut is not the only interior definition.

### 9.3 Vertices do not get EMT names until after the contingency table

Do **not** rename archetype 1 “epithelial” because you hope so. First make:

**Table B1** — `emt_state` × `nearest` (counts and row %).  
**Table B2** — mean \(w_1,\ldots,w_k\) per `emt_state`.  
**Table B3** — % inside and % specialist per `emt_state`.  
**Table B4** — `tumor` × `nearest` (row %). If one tumor is 100% nearest one vertex, write that down (TN-0135-type warning).  
**Table B5** — `emt_state` × `tumor` counts (does one tumor own one state?).

Only after B1: if vertex \(j\) is overwhelmingly epithelial cells, you may **display**-name it epithelial in figures, with the contingency as the evidence.

### 9.4 Figures (required)

1. **PC1–PC2**, gray all cells, vertices as large markers, colored by **nothing from EMT** (geometry only). Caption: fit was in \(k^\star-1\) D.
2. Same axes, points colored by **`emt_state`** (four colors). Vertices overlaid.
3. Same axes, colored by **`tumor`** (12 colors or small multiples).
4. Stacked bar: for each EMT state, fraction nearest each vertex.
5. Optional ternary if \(k^\star=3\); tetrahedron sketch if \(k^\star=4\) (PC1–PC3).

No UMAP in the official geometry figures.

### 9.5 How to read Stage B (pre-registered interpretations)

Write these in the notebook **before** plotting:

- **Support for four archetypes:** \(k^\star=4\) and the four states enrich four different vertices (B1 off-diagonal small).
- **Support for a lower-\(k\) continuum:** \(k^\star=3\) (or 4 with two states sharing a vertex); states ordered along an edge (epithelial → early hybrid → hybrid → late) in mean weights (B2).
- **Interior hybrids:** hybrid (paper) has low specialist % and weights spread across ≥2 vertices.
- **Confounded by tumor:** B4 shows a vertex that disappears if that tumor is removed (then Stage B is not a clean EMT geometry).

### 9.6 Optional Groves Panel B-style enrichment

For each vertex \(j\), Euclidean distance in the **\((k^\star-1)\)-D fit space** from cell to vertex. 10 equal-count distance bins. Hypergeometric enrichment of each `emt_state` in each bin; BH FDR 0.1; a “hit” if the significant peak is **bin 0** (closest). This is `src/enrichment.py` if present in the SCLC Panel B code. Optional, not required to answer the main question.

---

## 10. Sensitivities (after the official run, not instead of it)

Run only if Stage A+B finished. Each is a **separate** `best_k` / B1 table. Do not pick the prettiest.

| ID | Change | Why |
|----|--------|-----|
| S1 | Cap ≤1500 cells per tumor, then full A2–B | Pal-style dominance |
| S2 | HVG ~2000 genes instead of all genes after 1% filter | scRNA convention vs Hausser |
| S3 | MAGIC on the same genes | Pal pipeline; expect higher t-ratio |
| S4 | Leave-one-tumor-out: refit without each tumor, match vertices, % agreement | script 31 idea; 12 refits |
| S5 | No mito/ribo drop | PC1 housekeeping |

S1 and S4 are the important ones. S3 is optional and easy to over-interpret.

---

## 11. Script order (to implement)

Names are proposals; keep them linear.

| Script | Does |
|--------|------|
| `00_download_gse288137.py` | fetch GEO, write DOWNLOAD_LOG |
| `01_build_anndata.py` | counts h5ad + join metadata |
| `02_normalize_filter.py` | log-normalize, gene filter, write `X_cells_by_genes.npy` |
| `03_run_pca.py` | 50% variance rule, save scores + loadings |
| `04_esv_elbow.py` | ESV curve + DimensionFinder |
| `05_pcha_observed.py` | \(k=3\ldots7\), 15 inits, save archetypes + weights |
| `06_t_ratio_perm.py` | 100 then 1000 shuffles, 5 inits, write summary + `best_k.json` |
| `07_freeze_host.py` | npz for \(k^\star\) |
| `08_label_overlay.py` | join EMT/tumor, tables B1–B5, figures |

Run from **repo root**, `.venv/bin/python -u`.

Pilot: `06` with `--n-perm 10` on a 5,000-cell random subset **only** to test code. Discard those \(p\)-values.

---

## 12. Methods paragraph (draft, to paste later)

> We fit a Pareto simplex on the Pérez-González et al. GSE288137 YFP⁺ tumor-cell atlas (cells × genes, library-size normalized to 10⁴ UMI and log1p-transformed; genes detected in ≥1% of cells; mitochondrial and ribosomal genes removed). PCA was mean-centered without gene-wise scaling. The number of stored principal components was the smallest \(n\) explaining ≥50% of gene-space variance (at least 6). For \(k=3,\ldots,7\) we ran PCHA on the first \(k-1\) PCs with \(\delta=0\), 15 random initializations, retaining the maximum-volume simplex (ParTI PCHA convention). Explained sample variance versus \(k\) was summarized with ParTI DimensionFinder. Significance used the t-ratio (simplex volume / convex-hull volume in \(k-1\) dimensions) against 1000 independent PC-column shuffles with 5 PCHA initializations per shuffle. We selected the smallest \(k\) with permutation \(p<0.01\) (Hausser et al. 2019). EMT state and tumor identifiers from the source paper were joined after this choice. Mixture weights were barycentric coordinates in the fitting plane.

---

## 13. Checklist before claiming a result

- [ ] 54k-scale (or documented subset) matrix frozen; barcodes match metadata 100%.
- [ ] `emt_state` was not used in gene filter, PCA, PCHA, or \(k^\star\).
- [ ] Fit in \(k-1\) PCs, `delta=0`, 15 / 5 inits, NumPy &lt; 2.
- [ ] Official \(p\) from 1000 shuffles (or written exception).
- [ ] \(k^\star\) = smallest \(p<0.01\), not “four states.”
- [ ] B1–B4 saved; tumor dominance inspected.
- [ ] Geometry figures are PC plots, not UMAP.
- [ ] Breast Pal projection not mixed into this result.

---

## 14. First actions (human + code, in order)

1. Open GSE288137 on GEO; list supplementary files; see whether **cell-level EMT labels** exist.
2. Create `Mouse_SCC/data/raw/` and download.
3. Build `cell_metadata.csv` (`author` source if at all possible).
4. Implement `02` → `03` on the full matrix; print \(n_{50}\) and cells per tumor.
5. Implement `05`–`06` with `--n-perm 10` on a subset to verify PCHA runs.
6. Full observed fit + 1000-shuffle \(p\) → `best_k.json`.
7. Only then run `08_label_overlay.py`.

That is the entire project as specified: **polytope on mouse SCC expression, then the paper’s EMT and tumor labels on that polytope.**
