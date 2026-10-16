# Report validation record

What was verified while writing `Research_Report.md` (2026-10-11), what could not be verified, and the limitations that remain. "Verified" means checked against a file, the code, or a primary source in this session. It does not mean independently reproduced by another person.

## 1. Verified

### Code behaviour (read from source)
- `scc_common.py`: balanced fit cells (≤ 800/tumour, seed 0), gram/SVD PCA, projection, barycentric weights via `np.linalg.solve`, clipping, purity, `fit_simplex` (best-of-15 PCHA on k−1 PCs, t = simplex volume / hull volume).
- `src/archetypes.py`: `fit_pcha_best` slices to k−1 PCs and keeps max volume; `simplex_volume` = |det|/(k−1)!.
- Installed `py_pcha/PCHA.py`: FurthestSum initialisation from a random cell; random S; alternating projected-gradient S and C updates (10 inner steps); stop at relative ΔSSE < 1e-6, 500 iterations or varexpl ≥ 0.9999; varexpl = (SST − SSE)/SST.
- `03_fit_and_test.py`: independent shuffles of the first k−1 PCs; 5 starts per null fit; p = (1 + #null ≥ obs)/(1 + n); ESV in 10 PCs, best of 5; elbow = max distance from the chord.
- `04`–`06` logic as described in Report Ch. 6–7.

### Outputs
- All 15 p-values recomputed from `null_k*.npy` and `fit_k*.json`: identical to `fit_summary.csv`. All 3,000 null fits finite.
- Cell counts, assay counts, phase counts and tumour × state table from `cell_metadata.csv`.
- `expr_log_filtered.h5ad`: 23,959 × 14,502; `Rplp0`, `Rplp1`, `Rplp2`, `Rpsa` present (regex gap).
- QC medians by assay and tumour (multiome % mito 13.4 vs scRNA 4.6).
- PC drivers and within-state drivers read from CSVs; top PC genes from `pc_top_genes.csv`.
- Version B tables recomputed on the 14,059 scRNA cells; 5% cut-off radii (32.56, 7.92, 8.57); scRNA-only tumour spread per corner.
- Robustness CSVs (23 refits) read and plotted.

### Report-time analyses (new, not part of the pipeline)
- **5 vs 15 PCHA starts, k = 3:** six reruns with 5 starts gave t = 0.52897 (B) and 0.52842 (A), identical to the 15-start fits.
- **V1 long-gene score** per tumour, state and assay (`emp_v1_longgene.png`).

### Sources
- Pérez-González et al. 2026 (PDF text): 54,081 cells; QC thresholds; SoupX, scDblFinder, scVI/scANVI, Harmony, Leiden; state markers; Slingshot branches ("another ending at the hybrid EMT state"); cycling statement ("no cluster was dominated by cell cycle genes"); multiome nuclei protocol; reference list entries.
- GEO series matrix: LKPR genotype, tamoxifen protocol, YFP⁺ Lin⁻ sorting, nuclei isolation for multiome, sample/library names.
- Groves et al. 2022 (PDF text): specialist definitions (Fig. 3E closest 5%; Methods weight > 0.95); MAGIC use; 1,000 randomisations; bootstrap.
- Hausser 2019, Groves 2022, Blais & Jeukens 2024: titles and authors from PDFs.

### Follow-up analyses (12 Oct 2026, Report Addendum)
- `07_regress_v1.py`: fit cells asserted identical to version B; projection of fit cells asserted to reproduce the PCA scores; V1 score r² with PCs 1–10 = 0.000 after regression.
- All p-values of `all_scrna_noV1` recomputed from saved shuffles: match. Within-tumour p-values read from `fit_k*.json` (200 shuffles each, all finite).
- `09_cross_projection.py`: host cells asserted to reproduce the host PCA scores.
- Interpretation check run before stating it: the proliferating-hybrid corner is not a generic cell-cycle corner (dividing epithelial cells reach it 1%).
- **Withdrawn claim:** `Malat1` as evidence for a nuclear-RNA reading of V1 (r = 0.19 with the V1 score; high in all tumours). Corrected in §9.5, §9.11, §11.1, Q2.

## 2. Discrepancies found (also Report §18.4)

1. Pre-registration's "cycling hybrid state" branch is not in the paper.
2. Ribosomal regex misses `Rplp0/1/2`, `Rpsa`.
3. `02_pca.py` docstring: 0.1% rule "used for the ESV curve"; actually ESV uses 10 PCs.
4. `WALKTHROUGH.md`: "nuclei poorer in mitochondrial RNA" contradicts this dataset.
5. Version B tables in `describe_k3/summary.md` include projected multiome cells; `WALKTHROUGH.md` §6.4 quoted them as if fitted cells.
6. "Multiome cells land at V1" holds for the closest-5% rule only (gap rule: 79–98% interior).
7. 23,959 vs 54,081 cells; 4 multiome samples vs "three primary skin SCCs" in the paper.
8. Version A `pca.npz` computed by full SVD; others by gram (shown equivalent).
9. `match_vertices` unused.
10. First report version cited `Malat1` as evidence for V1 being nuclear RNA (withdrawn).
11. First report version's main conclusion (k = 3, two biological corners) is revised by the Addendum.

## 3. Could not be verified

- The R export step (commands not saved).
- Why the object contains 23,959 of 54,081 cells.
- How `phase`, `clusterU` and `LKPR_*_andrea` were computed by the authors.
- The mapping of `Early Hybrid EMT_1/2/3` to the paper's EMT-I/II.
- Whether V1 reflects nuclear/unspliced RNA (needs read-level data not on GEO).
- Citations 5, 6, 10–15, 19–26, 30, 31: bibliographic details from the author's knowledge, not re-checked online in this session.
- PCHA start-count asymmetry for k ≥ 4.

## 4. Remaining limitations of the analysis (summary)

Confounding of tumour, assay and state; one technical corner; null model rejects independence, not "no Pareto front"; modest effect sizes; rule-dependent specialist fractions; gene-set dependence (no simplex on KS genes); depth correlation; label circularity; no bootstrap uncertainty; MAGIC, Harmony and "arm A" versions not run.

## 5. Report build checks

- `Research_Report.pdf` opens with PyMuPDF; 47 pages (after the 12 Oct addendum); text extractable on every page.
- Spot-checked pages render title page, table of contents, equations (MathJax SVG), tables, figures and checkpoint boxes.
- All 29 figures referenced in the Markdown exist; captions are numbered 1–29 in order and match their image labels; equations are labelled Eq. 1–11 (checked with a script after assembly).
- Original files in `Mouse_SCC/codes/`, `Mouse_SCC/results/` and `Mouse_SCC/data/` were not modified while writing the report. (`WALKTHROUGH.md` and the pipeline were written earlier in the same project session; this report corrects some of their statements instead of editing them.)

## 6. Figure inventory

Schematic (toy or diagram): `schem_pipeline`, `schem_emt_states`, `schem_pareto`, `schem_pca`, `schem_barycentric`, `schem_pcha_tratio`, `schem_rules`, `schem_confounding`.
Empirical (from saved outputs): `emp_data_composition`, `emp_pc_drivers`, `emp_scree`, `emp_null_esv`, `emp_state_classes`, `emp_robustness`, `emp_v1_longgene` (report-time analysis).
Project (copied unchanged): 14 files in `figures/project/`.
