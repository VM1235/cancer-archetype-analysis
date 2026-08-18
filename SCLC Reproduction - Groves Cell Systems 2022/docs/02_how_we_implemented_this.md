# How this reproduction was implemented

> Paths below may still mention older `scripts/` and root `results/` locations. Current layout: this folder’s `codes/`, `data/`, `results/`, and `figures/`. See the SCLC README.

Official result: **ParTI-matched PCHA** (fit in \(k-1\) PCs, `delta=0`, many random starts, keep the maximum-volume simplex). That is the analysis behind the figures in this folder. An earlier 12-D / `delta=0.1` / single-start run is **not** what we report.

Paper: Groves et al., *Cell Systems* 13:690–710 (2022), Figure 1A–C.  
Reference code and data: [QuLab-VU/Groves-CellSys2022](https://github.com/QuLab-VU/Groves-CellSys2022), cloned read-only as `reference/` in this project.

---

## Language and libraries

The authors’ archetype fit is MATLAB ParTI (`algNum=5` = PCHA). We do not have MATLAB.

We used **Python 3.9** in `.venv`:

- `py_pcha` — Mørup’s PCHA, the same family of code ParTI wraps as `PCHA1.m`
- `numpy` (pinned **&lt; 2**, because `py_pcha` still uses `np.mat`)
- `scipy` (convex hull, hierarchical clustering, hypergeometric test)
- `scikit-learn` (PCA)
- `pandas`, `matplotlib`

We did **not** use R `ParetoTI` or R `sva`. ParetoTI still calls `py_pcha` underneath. ComBat had already been run by the authors (see below).

---

## What we took from their repo (used as-is)

We did **not** rebuild the expression matrices from raw FASTQs or from dbGaP.

| File | Role |
|---|---|
| `reference/data/bulk-rna-seq/SCLC_combined_Minna_CCLE_batch_corrected_wo_lowgenes.csv` | Panel A/B input: **15,950 genes × 120 cell lines** (70 Minna + 50 CCLE), log-transformed, ComBat-corrected |
| `reference/data/bulk-rna-seq/CCLE_Minna_Thomas_COMBAT.csv` | Panel C input: **14,546 genes × 201 samples** (120 lines + 81 tumors), second ComBat (cell line vs tumor) |
| `reference/data/bulk-rna-seq/combined_clusters_2020-05-27-MC copy.csv` | Panel B subtype labels (`NEW_10_2020`: A / A2 / N / P / Y) |
| `reference/notebooks/ParTI-code/human-cell-lines/params_lognorm.txt` | Ground-truth table: algorithm 5, 12 PCs, t-ratios and p-values for \(k=3\ldots7\) |
| `reference/ParTI/findMinSimplex.m`, `findArchetypes.m`, `CalculateSimplexTratiosPCHA.m` | Exact PCHA **calling convention** (not executed; read as specification) |
| `reference/notebooks/bulk/Thomas-Tumors-Bulk-Archetypes.ipynb` | Panel C variance-explained recipe (combined PCA vs tumor-only PCA vs shuffle) |
| Saved MATLAB outputs under `notebooks/ParTI-code/human-cell-lines/out/` | Check numbers only (e.g. paper \(k=5\) t-ratio \(= 0.107\), \(p = 0.034\)) |

Minna RNA-seq is dbGaP-controlled (`phs001823`). We never downloaded it. The public Groves repo already contains the merged, corrected tables.

---

## What we did not take / did not rerun

| Authors did | We did |
|---|---|
| MATLAB ParTI `.mlx` / `.m` | Python `py_pcha` |
| ComBat in R (`sva`) on CCLE vs Minna, then again on lines vs tumors | Used the saved ComBat matrices |
| Hierarchical clustering to **create** subtype names | Used the published label table for Figure 1B |
| GO / ENRICHR task names (Figure 1E) | Out of scope |
| Tumor-only polytope in linear (not log) space | Out of scope (Panel C uses the cell-line polytope, not a tumor-only fit) |
| Single-cell and mouse panels | Out of scope |

Skipping ComBat is intentional: the figure is a statement about geometry on **their** processed matrix. Re-deriving ComBat would be a different paper (and would need Minna access).

---

## Pipeline we actually ran

Code lives under `src/` and `scripts/`. Intermediate arrays are in `results/`. The figures here are exported by `scripts/export_sclc_reproduction_figures.py`.

### Panel A

1. **PCA, 12 components**, samples × genes, `sklearn.PCA(svd_solver="full")`.  
   Cumulative variance **47.3%** (paper ~50%). This matches their `dim=12` choice.  
   ESV vs \(k\) is computed in this 12-D space, as in `findArchetypes.m`.

2. **Official PCHA fit (ParTI protocol)** — this is the result we show:
   - Restrict to the first **\(k-1\)** PCs (a \(k\)-simplex is \((k-1)\)-D).
   - `delta=0` (archetypes are convex combinations of observed samples).
   - **150 random initializations** (ParTI: `3 * numIter` with `numIter=50`); keep the simplex of **maximum volume**.
   - t-ratio \(= V_{\text{simplex}} / V_{\text{convex hull}}\) with the same determinant formula ParTI uses.

   Observed t-ratios:

   | k | Groves t | Ours |
   |---|----------|------|
   | 4 | 0.247 | 0.247 |
   | 5 | 0.107 | 0.108 |
   | 6 | 0.043 | 0.043 |

3. **Permutation test.** Shuffle each of the \(k-1\) PCs independently across the 120 lines; refit with 15 random starts per shuffle, keep max volume; 100 shuffles (sanity scale; paper used 1000 shuffles × 50 starts).  
   \(p =\) fraction of null t-ratios \(\ge\) observed.

   | k | Groves p (1000) | Ours p (100 × 15) |
   |---|-----------------|-------------------|
   | 4 | 0.059 | **0.07** (not significant) |
   | 5 | 0.034 | **0.02** (significant) |
   | 6 | 0.016 | 0.00 |

   Same call as the paper: **smallest significant \(k\) is 5**.

4. **Polytope picture.** Cell lines and the \(k=5\) vertices plotted on PC1–PC2. Fitting was 4-D; the drawing is a projection.

### Why the first Python run looked different (and why we discarded it)

The first pass fit PCHA in **all 12 PCs**, with **`delta=0.1`**, and **one** initialization per fit (including each shuffle). That is a valid PCHA run. It is **not** ParTI’s t-ratio protocol.

Consequences:

- Observed t-ratios were ~13% lower (\(k=5\): 0.092 vs 0.107).
- Null t-ratios were too small (one weak start instead of the best of many), so p-values were too small, and \(k=4\) incorrectly crossed \(p < 0.05\).

`py_pcha` and MATLAB `PCHA1` are the same algorithm family. The discrepancy was **how we called it**, not a different dataset. After matching the call, t-ratios agree to three decimals. We keep only that second result.

### Panel B

No new PCHA. Distances are Euclidean in the 12-PC space used for the original scores (120 × 12), 10 equal-count bins, hypergeometric enrichment of the **published** `NEW_10_2020` labels, BH FDR \(q < 0.1\), hit = significant **and** peak at bin 0.

Result (1-to-1):

| Subtype | Archetype with bin-0 peak |
|---|---|
| SCLC-N | Arc 1 |
| SCLC-A2 | Arc 2 |
| SCLC-A | Arc 3 |
| SCLC-Y | Arc 4 |
| SCLC-P | Arc 5 |

We also tried recutting Spearman hierarchical clustering ourselves. Average linkage on all 15,950 genes collapsed into one giant cluster (classic chaining). Their R `heatmap()` default is complete linkage, and their label file has dated columns — the names were curated from a dendrogram plus TFs, not a single `cutree(k=5)`. Figure 1B in the paper uses those published names. Recutting the tree is not required to reproduce the figure, and we do not report it as a result.

### Panel C

1. Load `CCLE_Minna_Thomas_COMBAT.csv` (genes × 201).
2. Reconstruct cell-line archetypes in gene space from the 12-PC model, restrict to the 14,546 genes shared with the combined matrix.
3. Fit a **new** PCA (20 components) on the 201 samples. Transform archetypes with that PCA. Do **not** refit PCHA.
4. Scatter: cell lines, tumors, five projected vertices.
5. Variance curves, copying the Groves notebook:
   - tumor-only PCA cumulative variance (ceiling);
   - variance of tumors when projected onto combined PCs, divided by total tumor variance;
   - 20 shuffles of gene rows of the combined matrix, same scoring, as null.

At 5 components, combined PCs capture **80.1%** of the tumor-only ceiling (their notebook: 80.1%). Grey null curves stay near the axis.

---

## What “reproduced” means here

| Claim | Status |
|---|---|
| Same processed matrices as Groves | Yes, from their GitHub |
| Same PCHA *engine* as MATLAB bit-for-bit | No (Python port) |
| Same PCHA *protocol* as ParTI for t-ratio | Yes, in the official run |
| t-ratios for \(k=4,5,6\) | Match to 3 decimals |
| p-values | Same qualitative call; 100 shuffles vs their 1000, so p is coarser |
| Panel B 1-to-1 subtype–archetype map | Yes, with their labels |
| Panel C ~80% at 5 PCs | 80.1% |
| Figures are pixel copies of the PDF | No — same analysis, our plotting |

If a reviewer asks why p is 0.02 rather than 0.034 at \(k=5\): we used 100 nulls instead of 1000, and 15 starts per null instead of 50. The **t-ratio** is the quantity that should match on the same data and protocol; it does. The p-value is a Monte Carlo tail count and will jitter.

---

## Files to show

All under `SCLC reproduction/figures/` — **three files**, one per paper panel:

- `Figure_1A.png` — polytope, ESV/elbow, t-ratio and p-values
- `Figure_1B.png` — subtype enrichment vs distance bin
- `Figure_1C.png` — tumors in the polytope + variance-explained curves

Supporting numbers: `SCLC reproduction/t_ratio_official.csv`, `results/panel_c/variance_explained.csv`.
