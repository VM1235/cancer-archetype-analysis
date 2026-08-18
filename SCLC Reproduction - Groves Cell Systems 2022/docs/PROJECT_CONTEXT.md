# Project Context: Reproducing Figure 1A–C (SCLC Archetype Analysis)

This document is for an AI coding assistant (Cursor) working in this repo. Read this
fully before writing or editing any code. It explains the science, the exact pipeline,
the data sources, and what "done" looks like for each panel.

---

## 1. What we're doing

We are reproducing **Figure 1, panels A, B, and C** from:

> Groves, S.M. et al. "Archetype tasks link intratumoral heterogeneity to plasticity
> and cancer hallmarks in small cell lung cancer." *Cell Systems* 13, 690–710 (2022).
> https://doi.org/10.1016/j.cels.2022.07.006

The method (Archetypal Analysis / Pareto Task Inference) is the same one introduced in:

> Hausser, J. et al. "Tumor diversity and the trade-off between universal cancer
> tasks." *Nat Commun* 10, 5423 (2019). https://doi.org/10.1038/s41467-019-13195-1

Both PDFs are available in this project if needed for methods detail — `p1.pdf`
(Hausser et al., the general method paper) and `p2.pdf` (Groves et al., the SCLC
paper whose Figure 1 we're reproducing).

**After finishing Fig 1A–C, the next phase (separate task, not covered by this doc)
is applying the same pipeline to breast cancer cell lines from CCLE.** Keep the code
modular/reusable with that in mind — don't hardcode SCLC-specific logic where it's
avoidable.

---

## 2. The core idea (for context, not to re-derive)

Given a matrix of cell lines (rows) × gene expression (columns), Archetypal Analysis
(via the PCHA algorithm) finds the smallest number of "extreme" points (archetypes)
such that every data point can be approximated as a weighted mixture of those
archetypes, and the archetypes themselves are weighted mixtures of real data points.
Geometrically, the archetypes are the vertices of a low-dimensional polytope
enclosing the data cloud.

- **Panel A**: determines the right number of archetypes (k) for 120 human SCLC cell
  lines via (1) an explained-variance elbow curve across k=2–15, and (2) a
  permutation-based t-ratio significance test. Result: **k=5** is the smallest
  statistically significant number (p=0.034 for k=5; p=0.059 for k=4 — not significant).
- **Panel B**: shows the 5 geometrically-derived archetypes are enriched for the 5
  already-known SCLC subtypes (SCLC-A, -A2, -N, -P, -Y), independently labeled via
  hierarchical clustering. This is a validation step, not a new archetype-fitting run.
- **Panel C**: shows the archetype space (derived from cell lines only) generalizes to
  81 real human tumor samples — tumors project inside the same polytope, and 5
  components from a combined cell-line+tumor PCA explain ~80% of the variance that a
  tumor-only PCA explains.

---

## 3. Data sources

### 3.1 Cell line bulk RNA-seq (feeds Panel A & B)

Two source datasets, merged:
- **CCLE**: 50 human SCLC cell lines (publicly available via CCLE/DepMap).
- **Minna Lab / cBioPortal**: 70 human SCLC cell lines. **This is controlled-access**
  via dbGaP (accession `phs001823.v1.p1`) — do NOT assume this is trivially
  downloadable. Check first whether we already have local access or a pre-merged file
  (see §4) before attempting a fresh download.

Combined dataset after preprocessing: **120 cell lines × ~15,950 genes**,
log-transformed, batch-corrected (ComBat, via R `sva` package).

### 3.2 Human tumor bulk RNA-seq (feeds Panel C only)

- **George et al. 2015** dataset: 81 human SCLC tumor samples. Available via the
  paper's supplementary data / GEO.

### 3.3 Pre-processed data availability — CHECK THIS FIRST

The original authors' code repo may already contain fully preprocessed, batch-corrected
matrices, which would let us skip most of §3.1–3.2 entirely:

**Repo**: https://github.com/QuLab-VU/Groves-CellSys2022

Known relevant paths inside it:
```
data/bulk-rna-seq/parti-input/CCLE_Minna_Thomas_COMBAT_no_names.csv
    → likely the merged, batch-corrected CCLE + Minna + tumor ("Thomas") matrix

notebooks/ParTI-code/human-cell-lines/Human-cell-line-RNA-seq.mlx
    → MATLAB script that runs PCHA on the cell-line-only data (Panel A's core code)

notebooks/ParTI-code/human-cell-lines/params.txt
    → exact parameters used for each ParTI run (dimensions, algorithm number, etc.)
    → USE THIS to match settings if re-implementing in Python/R

notebooks/ParTI-code/human-cell-lines/out/
    → saved original archetype results (5 archetypes + other k values tested)
    → USE THIS to validate our own reproduction's numbers against ground truth

notebooks/ParTI-code/thomas-tumors/
    → ParTI results for the 81-tumor dataset alone (relevant to Panel C)

notebooks/ParTI-code/combined-data/
    → ParTI results for the combined cell-line + tumor dataset (relevant to Panel C)

notebooks/bulk/SCLC RNA-seq batch correction-CCLE-Minna.Rmd
    → the actual batch correction code for §3.1 (Panel A/B input)

notebooks/bulk/Cell-line-tumor-batch-correction-and-clustering.Rmd
    → batch correction for cell lines + tumors combined (Panel C input)

notebooks/bulk/Compare_cell_lines_tumors_archetypes.Rmd
    → comparing cell-line archetypes to combined-dataset archetypes (Panel C logic)

notebooks/bulk/Thomas-Tumors-Bulk-Archetypes.ipynb
    → Python notebook comparing 81 tumors to the 5 cell-line archetypes (Panel C)

environment.yml
    → conda environment spec used by the original authors — check this for exact
      package versions before assuming our own environment matches
```

**First task for Cursor**: clone/fetch this repo into a `reference/` subfolder (read-only,
don't edit), inspect what's actually inside `data/bulk-rna-seq/`, `ParTI-code/out/`, and
`environment.yml`, and report back what's available before writing new preprocessing code.
Do not duplicate work the authors already did if their output files are usable directly.

---

## 4. Tooling decision: MATLAB vs. Python/R port

The original PCHA/ParTI fitting was done in **MATLAB** (`ParTI` package, Hart et al. 2015).
We need to decide one of:

- **Option A (preferred if available)**: Run their actual `.mlx` MATLAB files directly.
  Fastest path to guaranteed-correct results. Check if MATLAB is installed/licensed
  in this environment before ruling this out.
- **Option B**: Port to Python using `py_pcha` (`pip install py_pcha`), matching exact
  parameters from `params.txt`.
- **Option C**: Port to R using the `ParetoTI` package (the authors themselves used this
  for single-cell panels later in the paper, so it's confirmed compatible with this method).
  Install via `remotes::install_github("vitkl/ParetoTI")`.

**Cursor: check what's available in this environment (MATLAB license, Python packages,
R packages) and propose which option to use before starting Panel A implementation.**

---

## 5. Pipeline — Panel A (run this first; everything else depends on it)

### Step 1: Preprocess raw data → clean matrix
1. Load CCLE and Minna cell line expression tables.
2. Remove all-NA genes/samples and mitochondrial genes.
3. Normalize to TPM; merge on overlapping gene set between the two sources.
4. Log-transform: `log(TPM + 1)`.
5. Filter low-expression genes: drop where `log(TPM) < 1` across all samples.
6. **Batch correct** (CCLE vs. Minna) using ComBat (R `sva` package), preserving
   biological signal tied to the four marker TFs: `ASCL1`, `NEUROD1`, `YAP1`, `POU2F3`.
7. Output: **120 cell lines × ~15,950 genes** matrix. Save this as an intermediate
   artifact (e.g., `data/processed/cell_line_matrix_combat.csv`) so later steps and
   Panels B/C can reuse it without recomputation.

> **Skip this step and go straight to Step 2 if a usable pre-batch-corrected matrix is
> found in the reference repo (see §3.3).**

### Step 2: Dimensionality reduction
1. Run PCA on the 120 × 15,950 matrix.
2. Keep top **12 principal components** (this specific number matches the paper's
   choice, based on an elbow in their own PCA scree plot explaining ~50% of variance —
   don't just hardcode 12 blindly, verify with our own elbow check but expect ~12).
3. Output: **120 × 12** matrix — this is PCHA's actual input.

### Step 3: Run PCHA across candidate k values
For k = 2 through 15 (minimum: k = 4, 5, 6 to save time on first pass):
1. Fit PCHA to the 120×12 matrix with `noc=k`.
2. Record: archetype coordinates (k × 12), the sample-to-archetype weight matrix `S`
   (120 × k — needed later for Panel B distance binning), and the explained sample
   variance (ESV).

### Step 4: ESV elbow curve (top-left plot of Panel A)
1. Plot ESV(k) − ESV(k−1) vs. k (or plot ESV(k) directly).
2. Identify the elbow — expect candidates around k=4, 5, or 6.

### Step 5: t-ratio permutation test (top-right plot of Panel A)
For each candidate k (4, 5, 6):
1. Compute convex hull volume of the 120 points (12-D) — `scipy.spatial.ConvexHull`.
2. Compute polytope volume from the k archetype coordinates — same `ConvexHull` on
   just the k points.
3. `t_ratio = polytope_volume / convex_hull_volume`.
4. Shuffle the 120×12 matrix column-wise (permute each PC's values independently
   across the 120 samples) — this destroys cross-PC correlation while preserving each
   PC's marginal distribution.
5. Refit PCHA on the shuffled data with the same k, recompute t-ratio.
6. Repeat shuffle+refit **1000 times** (use 100 for a fast first-pass sanity check,
   then scale to 1000 for final numbers — this step is slow).
7. `p_value = fraction of shuffled t-ratios >= real t-ratio`.
8. Expected result: k=4 → p≈0.059 (not significant); k=5 → p≈0.034 (significant).
   **Report our own p-values and flag any large deviation from these reference values.**

### Step 6: Select k=5, visualize polytope (bottom plot of Panel A)
1. Take k=5 as final (smallest significant k).
2. Fit a **separate, standalone 2-component PCA** (PC1, PC2 only) to the original
   120×15,950 matrix — this is just for visualization, distinct from the 12-PC space
   used for fitting.
3. Project the 5 archetype coordinates into this 2D space.
4. Plot: 120 cell line dots + 5 archetype points + all pairwise connecting edges
   (complete graph, 10 edges) between archetypes.

**Panel A deliverables to save for reuse in B/C:**
- The 5 archetype coordinates (in 12-PC space)
- The sample-to-archetype weight/mixture matrix `S` (120 × 5)
- The batch-corrected 120×15,950 matrix (for Panel B's fresh clustering + Panel C's merge)

---

## 6. Pipeline — Panel B (reuses Panel A's output; no new PCHA fit)

1. **Independently label cell lines**: hierarchical clustering (Spearman correlation
   distance, average linkage) on the same 120×15,950 matrix → assign each of the 120
   cell lines one subtype label from {SCLC-A, -A2, -N, -P, -Y}.
2. **Bin by distance to archetype**: for each of the 5 Panel-A archetypes, compute
   each cell line's Euclidean distance to that archetype (in the 12-PC space), then
   sort into 10 equal-sized bins (bin 0 = closest 12 cell lines).
3. **Hypergeometric enrichment test**: for each archetype × bin combination, test
   whether a given subtype label is over-represented in that bin vs. the rest of the
   dataset. Apply Benjamini-Hochberg FDR correction (q < 0.1). Significant only if the
   enrichment peak is at bin 0.
4. **Plot**: 5 subplots (one per subtype), x-axis = distance bin (0–9), y-axis =
   enrichment (fold-enrichment or similar), 5 colored lines per subplot (one per
   archetype). Expect exactly one line to spike sharply at x=0 in each subplot.

---

## 7. Pipeline — Panel C (new data + new batch correction; reuses Panel A's archetypes)

1. **Load tumor data**: 81 human SCLC tumor bulk RNA-seq samples (George et al. 2015).
2. **Batch correct combined data**: merge the 81 tumors with the 120 cell lines
   (on overlapping genes) and run ComBat again — this is a *separate* batch correction
   from Panel A's CCLE-vs-Minna correction; here the batch variable is cell-line-vs-tumor.
   Output: **201 samples × genes** matrix.
   (Check reference repo: `CCLE_Minna_Thomas_COMBAT_no_names.csv` may already be this file.)
3. **Fresh PCA**: fit a new PCA to this 201-sample combined matrix (independent of
   Panel A's PCAs).
4. **Project archetypes**: take the 5 archetype coordinates from Panel A (found using
   cell lines only) and project them into this new combined-PCA space as fixed
   reference points — do NOT refit archetypes here.
5. **Plot left panel**: scatter of cell lines (light green) + tumors (dark green) + the
   5 projected archetype points, in the new combined-PCA space (PC1 vs PC2 or similar).
6. **Variance-explained comparison (right panel)**:
   - Fit a separate PCA to tumors only (81 samples) → cumulative variance explained
     per component = the "ceiling" curve.
   - Using the combined-data PCA components (step 3), compute how much of the
     tumor-only variance they explain → this is the "actual result" curve.
   - Shuffle the combined data, refit PCA → null baseline curve.
   - Plot all three cumulative-variance-explained curves vs. number of components
     (up to ~15–20). Expect the combined-data curve to reach ~80% of the tumor-only
     ceiling by 5 components.

> **Methodological note to remember**: if fitting a *fresh* polytope directly to the
> tumor-only data fails to reach significance, try switching from log-space to
> linear expression space before assuming something is broken — the original authors
> found tumor mixtures only form a significant polyhedron in linear space, not log
> space, because physical cell-type mixing is a linear-space phenomenon.

---

## 8. Definition of done for this phase

- [ ] Reference repo cloned/inspected; report on what's reusable vs. what needs rebuilding
- [ ] Tooling decision made (MATLAB / py_pcha / ParetoTI) and justified
- [ ] Panel A: ESV curve, t-ratio/p-value table, and polytope visualization reproduced,
      with our own k=5 p-value compared against the paper's reported p=0.034
- [ ] Panel B: 5-subplot enrichment figure reproduced, confirming 1-to-1 archetype-subtype
      correspondence
- [ ] Panel C: combined PCA scatter + variance-explained comparison plot reproduced,
      confirming ~80% variance-explained-by-5-components result
- [ ] All intermediate matrices saved to disk with clear filenames so Panels B/C don't
      require recomputing Panel A's output
- [ ] Code structured modularly (data loading / preprocessing / PCHA fitting /
      significance testing / plotting as separable functions) since the same pipeline
      will be reused on breast cancer CCLE data in the next phase

---

## 9. Known pitfalls (from prior analysis of the paper)

- Minna Lab data requires dbGaP controlled access — don't assume it's a simple download.
- ComBat batch correction is easy to get subtly wrong; validate against the four marker
  TFs (ASCL1, NEUROD1, YAP1, POU2F3) retaining expected biological signal post-correction.
- The bottom polytope plot in Panel A uses a *different* 2-component PCA than the
  12-component PCA used for actual PCHA fitting — don't conflate the two.
- t-ratio permutation testing (1000 shuffles × PCHA refit) is slow — prototype with
  ~100 shuffles first, scale up only once the pipeline is verified correct.
- Archetype labels ("Proliferation", "Signaling & Secretion", etc.) are NOT produced by
  PCHA itself — they come from downstream GO-term/gene enrichment analysis (Panel E,
  out of scope for this phase) and subtype-matching (Panel B). Don't expect PCHA output
  to be pre-labeled.
- Tumor-only polyhedron fitting may require linear (not log) expression space to reach
  significance (see §7 note).
