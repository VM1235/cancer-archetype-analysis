# Addendum — Author Code Repository Is Now Public (2026-09-09)

**Companion to:** `SESSION_2026-09-09_FIG3.md`
**Trigger:** The paper's preprint states "All code will be made available upon publication." I checked, and the repo referenced in your own session note — `kelley27/PML_Archetypes` on GitHub — is now live. I cloned it and read the actual scripts the authors used to produce Figure 3. This addendum records what that changes.

---

## 1. Executive summary

| Question | Answer |
|---|---|
| Is the author repo really public now? | **Yes.** Confirmed via direct clone (`Step0`–`Step5` R scripts + `Libraries.R` + `outputData/`). |
| Was your Path B (v2) residual design correct? | **Yes, exactly.** `~Study + TopGenes` design with `duplicateCorrelation(block=Patient)` on voom/TMM matches the author's `Step0_compute_residuals.R` line for line. |
| Is the raw residual/count data public now too? | **No.** The repo's own README still lists all residual `.rds` files as "request from author," including `PML_resid_duplicateCor.rds`. Code ≠ data here. |
| Is there a shortcut that doesn't need the email? | **Possibly.** The repo's `outputData/PML_arcA4_pml.rds` is the *actual fitted archetype object* — and your session note says you already have a copy of this file locally. You may be able to project your existing GSVA scores onto the real published vertices without re-fitting your own polytope. |
| Does this explain the layout mismatch? | **Largely, yes.** `voom()` requires true counts; running the paper's exact formula on log(TPM/FPKM+1) substitutes for 42% of samples isn't an approximation of their pipeline, it's a structurally different one. This was already your root-cause #1 — now confirmed directly from their code rather than inferred. |

---

## 2. What I did

1. Pulled the bioRxiv preprint full text (`doi:10.1101/2025.10.04.680402`) to check the Methods section wording and the Data/Code availability statements.
2. Checked whether `github.com/kelley27/PML_Archetypes` resolves (it does — HTTP 200, not a 404).
3. Downloaded the repo as a tarball (`codeload.github.com`, no auth needed) and inspected every script.
4. Cross-referenced the repo's `outputData/` file list against the "author assets" your own session note says you already have on disk.

---

## 3. Repo contents

```
kelley27/PML_Archetypes/
  README.md
  Libraries.R
  Step0_compute_residuals.R      # limma/voom residuals
  Step1_compute_ADJ_mat.R        # gene adjacency matrix (module discovery)
  Step2_build_modules.R          # igraph clustering -> 9 modules
  Step3_optimizeK.R              # elbow + t-ratio permutation test for k=4
  Step4_compute_archetypes.R     # GSVA + ParetoTI bootstrap fit + PC projection + labeling
  Step5_assign_archetypes_tumor.R# projection of archetypes into new/independent datasets
  outputData/
    pml_modules.rds                          # <- you already have this
    PML_arcA4_pml.rds                        # <- you already have this (the fitted polytope!)
    PML_Archetypes.rds                       # <- csv version already in your data/author/
    ENS_SYM_ID.rds, ENS_ENT_ID.rds           # <- you already have these
    combo_adj_pml4.rds
    PML_GSVA_permutation_archetpyes3to6_5000perm.rds
    LUAD_PCA_All_Archetypes.rds, LUAD_PCA_Tumor_Archetypes.rds
    GSE166720_Tumor_Archetypes.rds, TRACERx_Archetypes.rds
    GSE68465_Archetypes.rds, GSE31210_Archetypes.rds, GSE30219_Archetypes.rds
```

The README's own "Input data (request from author)" list names nine residual/count-matrix `.rds` files — including `PML_resid_duplicateCor.rds`, the exact file your session note flags as needed for Path A. **That list hasn't changed just because the code is public.** The data gate is still there.

Notably, you already hold the small set of files that *aren't* gated (modules, ID maps, and — importantly — the fitted archetype object itself). What's gated is specifically the sample-level expression/residual data, not the archetype geometry.

---

## 4. Confirmed: your methodology was already right

`Step0_compute_residuals.R` does, in order:

- `DGEList()` + `filterByExpr(group = Histology)` + `calcNormFactors(method="TMM")`
- Computes each sample's "percent of counts from top 1000 genes" (`TopGenes`) — the same QC covariate your `sample_info_v2.csv` already carries
- `design <- model.matrix(~Study + TopGenes, ...)`
- `voom()` → `duplicateCorrelation(block=Patient)` → **second** `voom()` pass with that correlation folded in → second `duplicateCorrelation` → `lmFit(..., block=Patient, correlation=...)` → `eBayes()` → `residuals(fit, v)`

This is a two-pass voom/duplicateCorrelation refinement, and it is *exactly* what `07_pathB_resid_gsva_archetypes.R` was going for. There's no methodological correction needed here — your instinct to add patient blocking was right, not a deviation from the paper.

`Step4_compute_archetypes.R` confirms the rest of your Path B assumptions too:
- GSVA call: `kcdf="Gaussian", mx.diff=1, abs.ranking=FALSE` — matches what you used.
- Archetype fit: `fit_pch_bootstrap(GSVAScores, n=5000, sample_prop=0.8, noc=4, delta=0, seed=42792, conv_crit=1e-6)` — your 500-bootstrap Python fallback was a scaled-down but faithful stand-in for this exact call.
- Octile assignment: closest 12.5% (`bin_prop` sequence up to 1, first non-zero bin wins) — matches your `assign_octile`.

**Bottom line:** none of the discrepancy is explained by "you built the wrong model." It's explained by what that correct model was fed.

---

## 5. What's newly clarified (not previously knowable without the code)

### 5.1 `voom()` needs real counts — this is now a hard, not soft, constraint

The residual step operates on a raw count `DGEList`, TMM-normalized, then voom-transformed. Voom's precision weights are derived from the mean-variance trend of **counts**, not of already-log-transformed TPM/FPKM. Feeding it log(TPM+1) as if it were log-CPM doesn't just add noise — it breaks the weighting model voom is built on. This confirms (with the actual code in hand, not just inference) that root-cause #1 in your original note — TPM/FPKM standing in for counts on 3 of 4 cohorts — is the dominant, structural reason the fitted polytope doesn't look like theirs, more so than any bootstrap-count or ParetoTI-vs-PCHA difference.

### 5.2 The vertex→archetype-name mapping is a hardcoded relabel, not an algorithm

The final step of `Step4` is literally:

```r
arch_PCs_bin$Combined <- gsub("Archetype 1", "A1: Normal-like", ...)
arch_PCs_bin$Combined <- gsub("Archetype 2", "A4: Proliferation", ...)
arch_PCs_bin$Combined <- gsub("Archetype 3", "A2: Inflammation", ...)
arch_PCs_bin$Combined <- gsub("Archetype 4", "A3: Cell Adhesion", ...)
```

In other words: ParetoTI returned four unordered vertices ("Archetype 1–4" in whatever order it converged to for their seed/data), and the authors manually inspected the biology and relabeled them once. There is no generalizable "vertex → biological name" function to recover — it's an artifact of one specific fit. Your data-driven approach (label a fitted vertex by majority author-primary-label overlap among its nearest neighbors) is the right kind of workaround for this, given you can't reproduce their exact seed and input matrix.

### 5.3 PC projection has no imposed sign convention

`project_to_pcs(..., pc_method="svd", log2=FALSE, zscore=FALSE)` is a plain SVD-based PCA, same in spirit as your `sklearn.decomposition.PCA`. Neither the paper's code nor yours pins down a canonical sign for each axis — PCA/SVD sign is inherently arbitrary. This confirms your session note's finding that no single PC1/PC2 flip fully reconciles the two layouts: it's not just an axis-orientation mismatch, it's a genuinely different fitted geometry from genuinely different input data (Section 5.1).

---

## 6. The opportunity: you may already hold the real polytope

You don't need `PML_resid_duplicateCor.rds` to get the *paper's actual archetype coordinates* — you apparently already have `PML_arcA4_pml.rds`, which **is** that fitted object (the output of `average_pch_fits()` on the real 5000-bootstrap ParetoTI run). What you're missing is only the sample-level GSVA scores computed from the real residuals.

The paper's own Methods describe exactly this scenario for external datasets (their "Projection of archetypes into LUAD datasets" section, implemented in `Step5_assign_archetypes_tumor.R`): take a new sample's GSVA scores, project them against the **stored, fixed** archetype coordinates, and compute distances — no re-fitting.

You could do the same thing with your own re-derived GSVA scores (from your v1 or v2 residuals) instead of re-fitting your own tetrahedron from scratch:

1. Read `PML_arcA4_pml.rds` in R (or via `pyreadr`/`rds2py` in Python) and pull out the archetype coordinate matrix (the object's `$XC` slot, or equivalent — 9 modules × 4 archetypes).
2. Compute Euclidean distance from each of your samples' 9 GSVA scores to each of those 4 fixed vertices.
3. Run your existing octile/primary-secondary logic and PCA-for-plotting exactly as before, but color/label using distance-to-*their*-vertices instead of distance-to-*your own refit* vertices.

This does **not** require `ParetoTI`'s heaviest dependencies (`fit_pch_bootstrap` is what needs those) — steps 1–2 are just reading a serialized R object and doing linear algebra you already have code for. It also doesn't require any email or waiting.

**What it won't fix:** the underlying TPM/FPKM-substitute problem for 42% of samples. Your GSVA scores for those samples will still be somewhat off, so distances to the real vertices will be somewhat off too. But it removes the "independently-fit, arbitrarily-rotated tetrahedron" confound entirely, which is likely the single biggest visual discrepancy in your current Fig 3A/3B — and it's a cheap, fast experiment to run before committing to the email-and-wait path.

---

## 7. Updated recommended next steps

Supersedes Section 10 of the original session note:

1. **Fastest, try first:** extract `PML_arcA4_pml.rds`'s coordinates and project your existing v2 GSVA scores onto them (Section 6). Low effort, uses files already on disk, directly tests whether the layout mismatch is fixable independent of the counts problem.
2. **Cheap, parallel:** double check GEO supplementary files for GSE102511 / GSE166720 / GSE193725 for a raw-counts matrix you may have missed (some GEO RNA-seq series deposit both a normalized table *and* a raw counts table). If real counts exist there, you can rebuild the true 137-sample count matrix and re-run `Step0` verbatim — no SRA requantification, no email needed.
3. **Still the highest-ceiling fix, but gated:** email for `PML_resid_duplicateCor.rds` (and/or the four per-cohort residual files). The repo's README confirms this remains request-only even post-publication, so this path hasn't gotten any faster — but you now know exactly which files to name in the request, and that the rest of the pipeline (your Step 07/08 code) doesn't need to change to use them.
4. **Lower priority for now:** finishing the full `ParetoTI` Bioconductor install. It's still worth doing eventually (for a faithful 5000-bootstrap re-fit once you have real counts), but it doesn't help until #2 or #3 land — fitting on still-imperfect data with the "real" package won't move the layout.
5. **Deprioritized, unchanged:** SRA/Salmon requantification. Only revisit if #2 turns up nothing and #3 doesn't pan out.

---

## 8. Follow-up executed in this workspace (same day)

Saved this addendum as `docs/ADDENDUM_2026-09-09_AUTHOR_REPO.md`, mirrored author repo under `reference/author_repo/`, and ran the Section 6 projection:

| Item | Path |
|---|---|
| Script | `codes/09_project_author_archetypes.py` |
| Author XC export | `data/author/PML_arcA4_XC.csv` (`Angiogenesis` row renamed to `Regulation Of Endothelial Proliferation` to match modules) |
| Figures | `results/fig3_v2_author_arc/fig3A_*.png`, `fig3B_*.png` |

**QC vs Path B re-fit (`fig3_v2`):**

| Metric | Re-fit PCHA (v2) | Project onto author XC |
|---|---|---|
| Primary agreement (both primary) | 82.4% (n=34) | **93.5% (n=46)** |
| Signed-distance r (A1–A4) | 0.54–0.83 | **0.91–0.95** |

**Layout:** author XC projected into PCA-of-*our*-GSVA still has A1 on the **left** and A3 on the **right** (paper is the opposite). A single PC1 flip would put A1/A3 on the paper’s sides; A2/A4 still won’t match paper’s top vs bottom-left exactly because the sample cloud (and thus PCA basis) is still ours.

**Takeaway:** most of the *label/distance* discrepancy was from re-fitting our own polytope. Most of the remaining *page-layout* discrepancy is still the TPM/FPKM → non-voom residual GSVA space (and PCA sign).

### Code bug found in re-fit (`codes/08`) and patched (`codes/10`)

`08_pathB_bootstrap_pcha_plot.py` fitted PCHA on the first **3 PCs** and measured distances in that PC subspace. Author `PML_arcA4_pml.rds` / Step4 use archetypes and distances in **full 9-module GSVA** (`XC` is 9×4).

Fixed re-fit: `codes/10_refit_pcha_gsva9d.py` → `results/fig3_v2_refit9d/`

| Metric | `08` (PC subspace) | `10` (9D GSVA) | `09` (author XC, no re-fit) |
|---|---|---|---|
| Primary agreement | 82.4% | 78.6% | **93.5%** |
| Distance r | 0.54–0.83 | 0.59–0.88 | **0.91–0.95** |

So the PC-subspace bug was real and worth fixing, but **re-fitting on our GSVA still cannot recover paper layout or author XC**; locking their polytope (`09`) remains the best public-data figure for labels/distances.

---

## 9. Pointers

- Author repo: `github.com/kelley27/PML_Archetypes` (local mirror: `reference/author_repo/`)
- Preprint: `biorxiv.org/content/10.1101/2025.10.04.680402`
- Step0 / Step4 / Step5 under `reference/author_repo/`
