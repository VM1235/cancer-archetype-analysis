# Breast cancer Run_2 — how Figures 1A, 1B, 1C were made

This note is only for the **invasive breast carcinoma** analysis in `Breast Cancer/Run_2/`. It does not change the SCLC reproduction write-up.

**Why this exists.** Groves et al. 2022 (Figure 1A–C) argued that SCLC cell-line transcriptomes fill a simplex whose vertices match known subtypes, and that human tumors sit in the same space. The PI asked us to apply that **method** to breast: DepMap cell lines for A/B, then TCGA-BRCA tumors whose **ER/HER2 status is already given by histopathology** (not PAM50 LumA/LumB from Panel B).

**What “today’s” pipeline is.** Official breast outputs are **Run_2**. We reused shared engines in `src/` (`pca.py`, `archetypes.py`, `enrichment.py`, `io.py`, `preprocess.py`) and wrote breast-specific drivers under `Breast Cancer/Run_2/codes/`. We did **not** copy SCLC numbers, k, or PC count. We did **not** overwrite `results/` (SCLC) or original breast Run_1 matrices.

**Stack.** Python 3.9 `.venv`, NumPy &lt; 2 (`py_pcha` still uses `np.mat`), scikit-learn PCA, scipy (hull, hypergeometric), matplotlib. PCHA via `py_pcha` with the **ParTI `algNum=5` calling convention** (not MATLAB itself). Panel C ComBat via R `sva` 3.58 in `Breast Cancer/Run_2/rlib`. PAM50 via R `genefu` (Run_1 labels reused).

---

## Code map (what each file is for)

| File | Role |
|---|---|
| `Breast Cancer/codes/01_build_input_panelA.py` | DepMap → genes × 63 lines matrix |
| `Breast Cancer/codes/03_map_and_pam50.R` | Symbol → Entrez, genefu PAM50 (Run_1; not rerun in Run_2) |
| `Breast Cancer/Run_2/codes/run_panelA.py` | PCHA, ESV, t-ratio, Figure 1A |
| `Breast Cancer/Run_2/codes/run_panelB.py` | Distance bins + hypergeometric, Figure 1B |
| `Breast Cancer/Run_2/codes/prepare_tcga_brca.py` | Primary TCGA matrix + clinical columns |
| `Breast Cancer/Run_2/codes/combat_cellline_tumor.R` | Groves `sva::ComBat` (bc2: `mod=~1`, `ref.batch`) |
| `Breast Cancer/Run_2/codes/run_panelC_tcga.py` | Combined PCA, IHC scatter + variance, Figure 1C |
| `src/archetypes.py` | PCHA multi-start, simplex volume, ESV, column shuffle |
| `src/enrichment.py` | Equal-count distance bins, hypergeometric + BH |
| `src/pca.py` | PCA, sign alignment, inverse transform to gene space |

Parameter logs: `PARAMS_vs_Groves.txt` (A), `PARAMS_panelC.txt` (C).

---

## Data

### Cell lines (Panels A and B)

- Raw: `Breast Cancer/raw data/Model.csv` and `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (DepMap **log2(TPM+1)**).
- Filter: OncoTree invasive breast / BRCA **and** `ModelType=Cell Line` **and** present in the RNA table → **63 lines** (metadata had 77 invasive models; 63 have RNA-seq). Portal “65” vs 63 is a release/overlap issue, not a coding bug.
- Genes: start 19,205 protein-coding; drop those with **max log2(TPM+1) &lt; 1** (TPM &lt; 1 in every line) → **16,500 × 63**.
- **No ComBat** for Panel A: one study, unlike Groves’ Minna+CCLE merge.
- Matrix: `Breast Cancer/input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv`.

### PAM50 labels (Panel B only)

- Independent of PCHA. `genefu::molecular.subtyping(sbt.model="pam50")` after `org.Hs.eg.db` mapping (66/16,500 symbols unmapped).
- File: `Breast Cancer/results/panel_b/pam50_labels_panelA.csv` (ACH- IDs, 63/63).
- Counts: Basal 27, LumB 17, Her2 14, LumA 4, Normal 1.
- Run_2 **drops LumA + Normal** (5 lines) from the **enrichment test only**. The polytope is still fit on all 63.

### Tumors (Panel C)

- `Breast Cancer/Panel C/HiSeqV2` — UCSC Xena, **already log2**, 20,530 genes × 1,218 samples.
- `Breast Cancer/Panel C/TCGA.BRCA.sampleMap-BRCA_clinicalMatrix`.
- Keep **primary** (`-01` barcode and `sample_type_id=1`) → **1,097** tumors.
- Shared symbols with Panel A: **14,234**. No extra log.
- **Plot colors:** IHC ER and HER2  
  `breast_carcinoma_estrogen_receptor_status`  
  `lab_proc_her2_neu_immunohistochemistry_receptor_status`  
  Groups: ER+/HER2− (436), ER+/HER2+ (123), ER−/HER2+ (40), ER−/HER2− (126), incomplete (372).  
  **Not** `PAM50Call_RNAseq`, **not** Panel B LumA/LumB.

---

## Shared geometry (why A, B, and C are not three clusterings)

PCHA finds **vertices** of a simplex. A sample can be a mixture (weights sum to 1). Clustering would partition; we did not use clustering to define A or C.

A k-vertex simplex is **(k−1)-dimensional**. For k=4 the fit, t-ratio, and Panel B distances live in **3 PCs**. 2-D figures are cartoons of that 3-D object.

Archetype **numbers are arbitrary** (random inits). “Arc 3” is an index, not a biology name. Names would come from enrichment (B) or IHC (C).

---

# Figure 1A — `Figure_1A_breast.png`

**Question.** How many extremes do 63 invasive breast lines support, and is that simplex tighter than chance?

**Why Groves’ SCLC settings are not copied blindly.** They used 12 PCs because that was ~47–50% of **their** gene-space variance, and k=5 because it was the smallest k with permutation p&lt;0.05. We used the **same rules** on this matrix.

### Technical protocol (all three 1A plots share this fit)

1. PCA on 63 × 16,500. Smallest n with cumulative variance ≥ 50% → **8 PCs (50.0%)**. Saved `n_pcs.txt`, `pc_scores.npy` (63 × 8).
2. PCHA (`delta=0`) in the first **k−1** of those PCs. ParTI PCHA: `numIter=5` → **15** observed inits (keep **max volume**), **5** inits per shuffle (keep **max t-ratio**). 500 shuffles (paper 1000).
3. ESV curve k=2…9 in 8-D; DimensionFinder elbow (farthest point from the chord of the ESV curve, MATLAB `+1` indexing).
4. t-ratio = simplex volume / convex-hull volume of the points, in (k−1)-D. p = fraction of null t-ratios ≥ observed.
5. **k we plot:** smallest k with p&lt;0.05 if any; else DimensionFinder. Result: **k=4**, elbow, **no k with p&lt;0.05**.

t-ratios (500 shuffles, 15/5 inits):

| k | t-ratio | p |
|---|---|---|
| 3 | 0.683 | 0.44 |
| 4 | 0.354 | **0.082** |
| 5 | 0.162 | 0.12 |
| 6 | 0.058 | 0.15 |
| 7 | 0.015 | 0.32 |

**Interpretation we actually claim.** k=4 is the descriptive elbow. It is **not** a significant Pareto simplex. Breast looks more continuum-like than Groves SCLC (their k=5 was p=0.034).

**Wrong protocol we discarded.** Sisal-style 150 observed / 15 null inits. ParTI explicitly sets `numIter=5` for PCHA.

---

### Plot 1 — top left: Explained sample variance (ESV)

**What you see.** X = number of archetypes N = 2…15. Y = **% ESV on top of the N−1 model** (marginal gene-space ESV). Dashed line at N=4. Arrow: “Suggested number of archetypes by elbow.”

**What it is.** After each PCHA, ESV is the fraction of PC-space variance the simplex explains. We multiply by the 8-PC cumulative (50%) so the y-axis is a percent of **total gene-space** variance, matching ParTI’s extra ESV plot. The **increment** from N−1 to N is what is drawn. The elbow is where extra vertices stop buying much new variance (here the jump 3→4 is still visible; 4→5 flattens).

**Why this plot.** Groves (and ParTI DimensionFinder) pick k from this diminishing-return curve, not from a silhouette score.

**How it is drawn.** `esv_curve.csv`; `tot_esv = esv * cumvar[8]`; `delta_esv` from that; `dimension_finder()` on `tot_esv`.

---

### Plot 2 — top right: t-ratio vs k

**What you see.** X = N = 3…7. Y = t-ratio. Each point labeled with permutation p. Dashed line at N=4. Title notes 500 shuffles and `numIter=5`.

**What it is.** A large t-ratio means the data hull is almost the fitted simplex (little “empty space” outside the vertices). Shuffling columns of the (k−1)-D scores destroys that geometry; if the real t-ratio is still in the null bulk, the simplex is not special.

**Why this plot.** Elbow alone can overfit. Groves required p&lt;0.05 to **call** k. We show both: elbow at 4, p=0.08, not &lt;0.05.

**How it is drawn.** `t_ratio_parti_500.csv`. Nulls in `null_t_ratios_k*_parti_n500.npy`.

---

### Plot 3 — bottom: archetype space (63 lines)

**What you see.** Grey dots = cell lines. Four numbered blue vertices, grey edges of the tetrahedron, in **PC1 vs PC2**.

**What it is / is not.** This is **not** the 8-D fit and **not** the 3-D PCHA space. It is a **separate 2-PC PCA** of the 63 × 16,500 matrix. Vertices from 3-D PCHA are inverse-transformed to gene space (`inverse_transform_scores` on the 8-PC model, extra PCs zero-padded) then projected with this 2-PC model so the cartoon can be drawn.

**Why this plot.** Same visual Groves used: “here is the polytope.” Most grey points sit near the projected hull; that is expected if PCHA worked, and is **not** the significance test (that is plot 2).

**How it is drawn.** `pc_scores.npy` is **not** what is scattered here; `fit_pca(..., n_components=2)` is. Vertices: `archetypes_k4_parti.npy` (4 × 3).

---

# Figure 1B — `Figure_1B_breast.png`

**Question.** Do **independent** PAM50 classes sit at the Panel A vertices?

**Why PAM50 here and not IHC.** CCLE lines do not come with clinical ER/HER2 IHC. PAM50 is the closest established RNA subtype. The PI’s TCGA request (histopathology) is Panel C, not this figure.

**Why drop LumA and Normal.** n=4 and n=1 cannot support 5-bin enrichment. Groves had five well-populated SCLC classes. We keep LumB, Her2, Basal (n=58). The **vertices stay the 63-line k=4 fit**.

### Technical protocol (all three 1B panels share this)

1. Distances: Euclidean in the **first 3 PCs** of the 8-PC cell-line scores, to the four PCHA vertices.
2. For each archetype, sort the 58 lines by distance; cut into **5 equal-count bins** (remainder allowed). Bin 0 = closest. Typical bin size 11–12.
3. Hypergeometric: is subtype S over-represented in that bin vs the other 57? BH q across all subtype × archetype × bin tests. **q &lt; 0.1**.
4. **Hit** (open circle on the curve): significant **and** the fold-enrichment **peak for that subtype–archetype pair is at bin 0**. A significant bump at bin 1 does not count as “this vertex is that subtype.”

Fold enrichment = (k/n) / (K/N), i.e. observed subtype fraction in the bin over the global fraction.

---

### Plot 1 — left: LumB

**What you see.** Four colored lines (Arc 1–4). Y = fold enrichment vs distance bin. Horizontal dashed line at 1 (no enrichment). Open circle on Arc 3 at bin 0.

**What it means.** LumB is enriched among lines **closest to Arc 3** (fold ≈ 2.56, q &lt; 0.1, peak at bin 0). Other arcs do not pass the bin-0 hit rule for LumB.

---

### Plot 2 — middle: Her2

**What you see.** Same four arcs. No open circle at bin 0.

**What it means.** Her2 is not a vertex specialist under Groves’ rule. There is a significant enrichment of Her2 in **Arc 3, bin 1** (near but not at the vertex). That is consistent with Her2 lines sitting **inside** an edge/face, or with n=14 being too small for a clean bin-0 spike. We report **no Her2 hit**, not “Her2 = Arc 3.”

---

### Plot 3 — right: Basal

**What you see.** Open circle on **Arc 4 at bin 0**. Arc 1 is high at bin 0–1 but its peak is **not** bin 0 (peak at bin 1), so it is not a hit. Arc 3 is **depleted** near the vertex and high at bins 3–4 (Basal is far from the LumB vertex).

**What it means.** Basal ↔ Arc 4 is the second 1–1-style match. Arcs 1 and 2 have no subtype that peaks significantly at bin 0.

**Run_2 hits (0-based archetype in the CSV, 1-based on the figure):** LumB → Arc 3; Basal → Arc 4; Her2 none; empty vertices Arc 1 and 2.

Code: `run_panelB.py`. Table: `results/panel_b/enrichment_pam50.csv`.

---

# Figure 1C — `Figure_1C_tcga.png`

**Question.** Do TCGA-BRCA primaries, labeled by **IHC ER/HER2**, occupy the **same** cell-line archetype geometry?

**Why this is not Panel B transferred.** Panel B labels are PAM50 on cell lines. The PI asked for existing histopathology on tumors. We never assigned LumA/LumB to TCGA for this figure.

**Why we do not refit PCHA on tumors.** Groves’ claim is generalization of the **cell-line** polytope. A new tumor-only simplex would be a different paper (and they noted tumor-only polytopes can need linear, not log, space).

### What we tried first and threw away

Transform tumors with the **cell-line-only** 8-PC PCA (missing genes filled with the PCA mean). Almost all tumors nearest Arc 4, distances ~400–700. That is DepMap log2(TPM+1) vs Xena HiSeqV2 as a **batch axis**, not “tumors don’t fit biology.” Groves never used that projection for Figure 1C.

### Groves protocol we implemented (both 1C plots)

Same as their `Cell-line-tumor-batch-correction-and-clustering.Rmd` (they **wrote `bc2`**, ComBat, not the fsva object) and `Thomas-Tumors-Bulk-Archetypes.ipynb`, and our SCLC `scripts/run_panel_c.py`:

1. Merge 63 lines + 1,097 tumors on 14,234 genes (`merged_uncorrected_shared_genes.csv`).
2. `sva::ComBat(dat, batch, mod=~1, par.prior=TRUE, ref.batch="cell_line")`. Cell-line batch is the reference (their SCLC `ref.batch` was Minna). Output `CCLE_TCGA_COMBAT.csv`. One gene with uniform expression in a batch was left unadjusted (sva warning).
3. **New** PCA, 20 components, on the ComBat matrix (1,160 samples). This is independent of Panel A’s 8 PCs.
4. Map k=4 vertices: 8-PC cell-line PCA → gene space → restrict to 14,234 genes → `pca.transform` into the combined 20-PC space. **PCHA not run on TCGA.**
5. Distances / nearest vertex (for tables, not a separate plot): Euclidean in the first **3 combined PCs**.

---

### Plot 1 — left: tumors in cell-line archetype space

**What you see.**

- Grey: 63 cell lines.
- Light grey: 372 tumors missing ER or HER2 Positive/Negative.
- Blue / pink / orange / green: IHC groups (n=436 / 123 / 40 / 126).
- Four colored vertices + grey tetrahedron edges.
- Axes: **PC1 / PC2 of the combined ComBat PCA** (not Panel A PC1/PC2).

**What it is.** After ComBat, lines and tumors share a scale (roughly ±100), unlike the failed projection. Vertices are the **same biological corners as Figure 1A**, drawn in this new basis. IHC is the only tumor color.

**What it is not.** It is not proof that every tumor lies **inside** the tetrahedron. Many points sit outside the wireframe. Barycentric weights in 3-D: **0 / 1,097** tumors have all weights ≥ 0. Groves SCLC tumors sat inside a 5-vertex polytope; we do not force that story here. The cell-line hull is small relative to TCGA spread.

**Biology we do report (nearest vertex, IHC complete n=725):**

| IHC | Arc 1 | Arc 2 | Arc 3 | Arc 4 |
|---|---|---|---|---|
| ER+/HER2− | 25 | 68 | 187 | 156 |
| ER+/HER2+ | 4 | 20 | 54 | 45 |
| ER−/HER2+ | 12 | 9 | 5 | 14 |
| ER−/HER2− | **82** | 28 | 1 | 15 |

ER−/HER2− (TNBC-like) toward Arc 1; ER+ toward Arc 3/4. Tendency, not a 1–1 map. Arc indices still match Panel A/B (Basal was Arc 4 on lines; tumor ER−/HER2− is Arc 1 — **do not** equate those without a joint analysis; combined PCA rotated the space).

---

### Plot 2 — right: variance explained

**What you see.** X = number of components 1…20. Y = cumulative **% of total tumor variance**. Blue: PCA fit on tumors only (ceiling). Orange: tumors scored on the **combined** PCA axes. Grey: 20 shuffles of **gene rows** of the ComBat matrix, PCA refit, original tumors projected (null). Vertical dashed line at **4** (k, Groves used 5 because k=5). Title: “4 combined PCs = 100% of tumor-only ceiling.”

**What Groves meant.** Not “4 PCs explain 100% of tumors.” They meant: the first k axes of the **combined** PCA capture ~80% of **whatever a tumor-only PCA captures with the same k**. If that ratio is high and shuffles stay low, tumors live in the same linear subspace as the lines.

**How we computed it.** Tumor total variance = sum of gene-wise sample variances. Combined-on-tumors EV = variance of tumor scores on those PCs / total. Ratio of **cumulative** curves at k=4 is ~0.999. Absolute variance at 4 PCs is only ~33% of all tumor variance (the ceiling itself is ~33% at 4 PCs); the orange and blue curves overlap, so the **ratio** is ~100%.

**Why this analogue is weak for breast.** Groves: 120 lines vs 81 tumors. Ours: **63 lines vs 1,097 tumors**. Combined PCA is almost a tumor PCA, so orange ≈ blue by construction. Grey shuffles stay near the axis (gene identity scrambled vs unshuffled tumors), which still shows the ceiling is not a coding bug. **Do not quote “100% of ceiling” as a Groves-style 80% result.**

---

## How A, B, and C fit together (breast)

| Figure | Plots | Question | Answer we have |
|---|---|---|---|
| 1A | ESV elbow; t-ratio+p; 2-D polytope | How many cell-line extremes? | k=4 by elbow; **not** p&lt;0.05 |
| 1B | LumB; Her2; Basal enrichment vs distance | Do PAM50 classes occupy vertices? | LumB–Arc 3, Basal–Arc 4; Her2 no bin-0 hit |
| 1C | ComBat combined scatter (IHC); variance vs shuffle | Do IHC-labeled TCGA tumors live there? | Same protocol as Groves; same **space** after ComBat; not all **inside** the 4-simplex; IHC groups differ by nearest vertex |

Panel A is a **descriptive** 4-pole model. B and C still follow Groves’ pipeline; they do not upgrade A into a significant Pareto front.

---

## Commands (repo root)

```bash
.venv/bin/python -u "Breast Cancer/codes/01_build_input_panelA.py"
.venv/bin/python -u "Breast Cancer/Run_2/codes/run_panelA.py"
.venv/bin/python -u "Breast Cancer/Run_2/codes/run_panelB.py"
.venv/bin/python -u "Breast Cancer/Run_2/codes/prepare_tcga_brca.py"
.venv/bin/python -u "Breast Cancer/Run_2/codes/run_panelC_tcga.py" --n-shuffle 20
```

`--force-combat` rebuilds the ComBat matrix. PAM50 R step is only needed if labels are missing: `Breast Cancer/codes/03_map_and_pam50.R`.
