# Connecting EMT\_ID with archetype analysis of EMT

*Ideas, methods, and a ranking by ease of implementation and value to the EMT\_ID project*

Based on a read of `Sampoorn135/EMT_ID`, `VM1235/cancer-archetype-analysis`, POLYTOPE\_FIT\_PLAN.md, and Pérez-González et al., Nat. Commun. 17:8561 (2026). No expression matrix was available, so no archetype fit has been run. Everything below is a proposal, not a result. Numbers quoted from EMT\_ID were computed from its `cell_table.csv.gz`.

## 1. Where the two projects meet

**EMT\_ID** (GSE288137 plus multiome) estimates a per-cell local intrinsic dimension (IDEAS, `local_2nn`, k=600, on raw counts), signalling entropy (CCAT), KS epithelial and mesenchymal ssGSEA scores (mouse signatures, 164 and 48 genes), and an HPCS programme score. It covers 23,959 cells from 12 tumours (8 scRNA, 4 multiome) with four EMT labels: Epithelial, Early-Hybrid, Hybrid, Late-EMT. Its headline is that ID is higher in Epithelial and Early-Hybrid (medians about 29) than in Hybrid and Late-EMT (about 22 to 23).

**Your repository** fits a Pareto simplex (PCHA on the first k-1 PCs, delta=0, permutation t-ratio), projects other datasets into a frozen host simplex by barycentric weights, names vertices with KS gene sets, and classifies cells as specialists (largest weight exceeds the second by at least 0.35) or interior generalists. Your Pal result is the key precedent: the hybrid corner was a single patient (TN-0135), and nearest vertex did not always agree with gene-score class.

The two are complementary. EMT\_ID describes how complex the local neighbourhood of each cell is. Archetype analysis describes where each cell sits in a low-dimensional polytope and which "task" it leans toward. Combining them lets each project test the other's blind spot.

**Shared prerequisite for every idea below.** All six ideas need the same Stage A fit from your plan: log-normalised cells x genes matrix, PCA, PCHA for k=3 to 7, permutation t-ratio, and a frozen host with barycentric weights. Request `emt_scc_ideas_input.h5ad` (raw counts, the same 23,959 cells) from the coworker. Their repository contains only CSV score tables, and the labels join exactly on the h5ad index. Joining to raw GEO files is error-prone because about 1% of their barcodes carry a `-1-1` suffix.

Three adjustments to the plan apply to everything that follows:

- Use the 23,959 labelled cells, not the paper's 54,081. Labels for the full atlas may not be public; this could not be confirmed.
- Treat leave-one-tumour-out and an scRNA-only fit as primary analyses. Four scRNA tumours (SC2b1, SC2b3, SC35, SC43) are 89 to 98% hybrid or late, while the multiome tumours are 93% epithelial or early hybrid, so tumour, assay and state are heavily confounded.
- Add a batch-corrected arm. The paper defined its states on scANVI and Harmony space, so a raw-expression simplex alone cannot separate "late-EMT vertex" from "one tumour's vertex".

## 2. The ideas and the method for each

### Idea 1. Geometry-controlled ID and entropy

**Why.** Local ID can depend on where a cell sits in a point cloud. Cells at a corner have one-sided neighbourhoods, and interior cells have two-sided ones. If ID tracks simplex position, part of the Epithelial and Early-Hybrid versus Hybrid and Late difference is geometric. If it survives control, the EMT\_ID claim is much stronger. This addresses their central result directly.

**Method.**

1. Join `ID_score_percell_raw`, `CCAT_score`, `transcriptome_entropy` to the frozen weights by barcode.
2. Per cell, compute purity (largest minus second-largest clipped barycentric weight) and the weights w1 to wk, using your existing definitions in `Breast Cancer/codes/26_specialist_vs_interior.py`.
3. Model ID as a function of purity and weights, adding `emt_state` as a second predictor. Use a mixed model with a random intercept for tumour and a fixed effect for assay. The unit of biological replication is the tumour (n=12), not the cell, as the EMT\_ID README itself notes.
4. As a non-parametric check, compute Spearman correlation of ID with purity within each tumour, then summarise across tumours with a sign test.
5. Compare the simplex dimension k\*-1 with the ID values (about 22 to 29). They should differ by a lot, since ID on raw counts is noise-dominated. Report both rather than expecting agreement.
6. Repeat with the state labels permuted within tumour to show the effect is not an artefact of state sizes.

**Caveats.** `compute_ID` takes no seed, so the deposited ID values are one random draw. IDEAS works in gene space with 600 neighbours, while PCHA works in a few PCs, so the two see different geometry.

### Idea 2. Vertex identity and robustness (tumour, assay, state)

**Why.** Without this, no other result is interpretable. It is the TN-0135 check, and EMT\_ID already has the mouse KS scores needed.

**Method.**

1. Name each vertex from the KS\_Epi and KS\_Mes ssGSEA scores of the cells nearest it, using the Pal rule: both scores high (75th percentile or above) means hybrid. Do this after the contingency tables B1 to B5 in your plan, not before.
2. Build `emt_state` x nearest-vertex and tumour x nearest-vertex tables. Report how many tumours supply each vertex.
3. Run leave-one-tumour-out by adapting `Breast Cancer/codes/28_pal_leaveout_tn0135.py`: refit k=3 without each of the 12 tumours (15 starts, delta=0) and match vertices by distance.
4. Run an scRNA-only fit (8 tumours) and a batch-corrected fit (scANVI or Harmony PCs as the input to PCHA) and compare vertex stability.
5. Report the gap between nearest vertex and gene-score class, as in your TN-B1-0554 finding.

**Caveat.** With only 12 tumours and 4 supplying most late cells, a vertex may genuinely depend on a few tumours. Report it as such rather than trying to fix it.

### Idea 3. Gene and TF gradients across the simplex

**Why.** The paper's biology is a TF hierarchy: Klf5 and Pitx1 in epithelial and early hybrid, Zeb1, Twist1 and Prrx1 rising progressively, Nfatc1 and Creb3l1 in late EMT. Archetype weights give a continuous coordinate on which to test whether each TF is vertex-specific or edge-specific. This is the Hausser-style "which genes peak near each vertex" analysis.

**Method.**

1. For each vertex, compute the Euclidean distance from every cell in the (k\*-1)-D fit space and form 10 equal-count bins with `src/enrichment.py` (`distance_bins`).
2. For genes, test whether expression peaks in the closest bin (Spearman against distance, FDR 0.1). Run GO or hallmark enrichment on the peaking genes to name the "tasks" (for example proliferation, ECM, differentiation).
3. For the paper's TFs, regress expression (and SCENIC regulon AUC if the coworker or authors supply it) on barycentric weights. Predictions: Klf5 and Pitx1 decline along the epithelial-to-late edge; Nfatc1 and Creb3l1 peak at the late vertex; Zeb1 is broad.
4. Score cell cycle and include it as a covariate, since the paper reports cycling is reduced specifically in late EMT.

**Caveats.** Ranking genes by distance to a vertex that was fit on the same genes is partly circular; hold out a gene subset or use the multiome as replication (Idea 5).

### Idea 4. Guest projection of the knockout bulk RNA-seq

**Why.** The paper's knockouts (Pitx1, Klf5, Nfatc1, Creb3l1; GSE288196) make directional predictions that a polytope can test. Pitx1 and Klf5 loss should move cells toward hybrid and late. Nfatc1 and Creb3l1 loss should move them toward epithelial and early hybrid.

**Method.**

1. Freeze the host (PCA mean, components, vertices) from Stage A.
2. Apply the host PCA to the bulk log-normalised matrix (same mouse gene symbols, no ortholog step) and compute barycentric weights, following `Breast Cancer/codes/22_project_depmap_onto_gse161529.py`.
3. Report the knockout-minus-control shift vector for each factor, separately for EPCAM-positive and EPCAM-negative fractions. Shift vectors cancel much of the constant offset between bulk in vitro and single-cell in vivo.
4. Compare shift direction with the predictions; no significance test is possible with n=2.

**Caveats.** Bulk and single-cell have different scales, the cells were TGFβ-treated in culture, and late EMT was not observed in vitro. I did not verify that count matrices are deposited for GSE288196.

### Idea 5. Replication on the multiome and comparison with trajectories

**Why.** The multiome (GSE288220) has five states, including early hybrid I and II, and the paper reports a stepwise I-to-II progression. If the simplex is real, I should sit nearer the epithelial vertex than II. The paper's Slingshot analysis also gives two branches from the epithelial root, which predicts a triangle (epithelial, late EMT, cycling hybrid) with early hybrid along the edge.

**Method.**

1. Treat the scRNA-only fit as host and project the multiome RNA cells in as a guest, as in `18_project_gse176078_onto_gse161529.py`. Report the fraction inside the simplex, as in your Wu and Pal tables.
2. Compare mean weights of early hybrid I and II.
3. If pseudotime is available from the authors, correlate it with the coordinate along the epithelial-to-late edge.
4. Write the branched-triangle prediction down before looking at the weights.

**Caveats.** Nuclei (multiome) and whole cells (scRNA) differ in chemistry and gene detection, so inside-the-simplex rates will be lower. The EMT\_ID repository does not contain pseudotime values.

### Idea 6. Archetypes on the reprogramming and lung datasets (stretch)

**Why.** EMT\_ID also analyses mouse reprogramming (181,124 cells, with Waddington-OT fate probabilities) and lung HPCS data. Archetype weights could be compared with WOT fate probabilities and HPCS scores to ask whether potency, as measured by ID, relates to simplex position outside EMT.

**Method.** Fit a separate simplex per dataset, never a shared one, then repeat Idea 1 and correlate weights with `w_*` and `P_*` fate columns.

**Caveats.** The reprogramming object is 8.6 GB and a PCHA permutation test at that size is expensive. This is outside the EMT scope of the question.

## 3. Ranking

Scores run from 1 (hard or low value) to 5 (easy or high value). "Ease" assumes the shared Stage A fit already exists and asks how much new code is needed beyond your current pipelines. "Value" asks how much the result would strengthen or change EMT\_ID's conclusions.

| Rank | Idea | Ease | Value to EMT\_ID | Main reusable code | Why |
| --- | --- | --- | --- | --- | --- |
| 1 | 2. Vertex identity and robustness | 4 | 4 | `26_specialist_vs_interior.py`, `28_pal_leaveout_tn0135.py`, KS scores already in EMT\_ID | Nearly a re-run of existing scripts, and a precondition for trusting everything else |
| 2 | 1. Geometry-controlled ID and entropy | 4 | 5 | purity and weights from `26_`; small new regression code | Tests the coworker's main claim directly with a control they do not currently have |
| 3 | 3. Gene and TF gradients | 3 | 4 | `src/enrichment.py`, `05_emt_generalist_panelD.py` | Links the paper's TF biology to geometry; needs new gene-level enrichment code and a circularity guard |
| 4 | 5. Multiome replication and trajectory | 3 | 4 | `18_project_gse176078_onto_gse161529.py` | Independent replication and a pre-registered topology test; needs multiome counts and possibly pseudotime |
| 5 | 4. Knockout projection | 3 | 3 | `22_project_depmap_onto_gse161529.py` | Directional predictions, but n=2, in vitro, and bulk versus single-cell offsets limit what it can show |
| 6 | 6. Reprogramming and lung | 1 | 2 | none directly | Large data, separate fits, and off the EMT question |

**Suggested order.** Run Ideas 2 and 1 together as the first milestone, since they use the same fit and answer "is the geometry real, and does it explain the ID result". Then Idea 3, then Idea 5. Idea 4 is a good figure but weak evidence. Idea 6 should wait.

## 4. Practical notes on the repository

- `POLYTOPE_FIT_PLAN.md` refers to `host_utils.py` (for `barycentric`, `fit_host_matrix`, `describe_scores`). That file is not in the GitHub clone I read. Barycentric functions are defined separately inside individual scripts (for example `28_pal_leaveout_tn0135.py` and `18_project_gse176078_onto_gse161529.py`). Confirm whether `host_utils.py` is local and uncommitted, and commit it or consolidate the helper into `src/`.
- `src/archetypes.py` matches the plan's description: `fit_pcha_best` slices to k-1 PCs and keeps the maximum-volume simplex, and `permutation_t_ratio` does a single fit per shuffle, so the multi-init null in the plan is a real change. Its `t_ratio` uses a hull for both volumes, unlike the ParTI formula the plan specifies.
- A significant permutation p-value with tens of thousands of cells means the PCs are not independent. It does not on its own show a simplex. Batch structure alone can produce it, which is another reason to run Idea 2 first.
- Runtime: 1000 shuffles x 5 inits x 5 values of k is thousands of PCHA fits. Time a 10-shuffle pilot before committing.
