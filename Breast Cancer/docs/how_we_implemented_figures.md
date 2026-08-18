# How we implemented the breast Figure 1A–C analogues

This is the current breast analysis (Groves *method*, not SCLC numbers). Official codes, figures, and results now live under `Breast Cancer/{codes,figures,results}/` (formerly `Run_2/`).  
Paper: Groves et al., *Cell Systems* 13:690–710 (2022), Figure 1A–C.  
Reference: `reference/` (QuLab-VU/Groves-CellSys2022). SCLC reproduction notes live separately under `SCLC reproduction/`.

Figures:

| File | Panel |
|---|---|
| `Breast Cancer/Run_2/figures/Figure_1A_breast.png` | A — cell-line polytope |
| `Breast Cancer/Run_2/figures/Figure_1B_breast.png` | B — PAM50 vs vertices |
| `Breast Cancer/Run_2/figures/Figure_1C_tcga.png` | C — TCGA tumors in that space |

Python env: `.venv`, NumPy **&lt; 2** (`py_pcha` still uses `np.mat`). Do not overwrite `results/` (SCLC) or original breast Run_1 outputs.

---

## Shared idea (all three panels)

Groves’ claim is a **simplex**, not a clustering:

1. **A** — cell-line transcriptomes sit in a low-dimensional polytope whose vertices are archetypes (PCHA).
2. **B** — known subtypes enrich at those vertices (distance bins + hypergeometric test).
3. **C** — human tumors occupy the **same** geometry after batch-correcting lines vs tumors. Archetypes are **not** refit on tumors.

Breast differences from SCLC that we did **not** hide: one DepMap study (no Minna/CCLE ComBat for Panel A); k from this matrix (not hardcoded 5 or 12 PCs); Panel C labels are **histopathology ER/HER2**, not PAM50.

---

## Figure 1A — cell-line archetypes

### Concept

If expression is mixtures of a few extremes, samples fill a simplex whose corners are those extremes. PCHA finds the corners. Significance asks whether that simplex is tighter than in shuffled data (t-ratio).

### Data

- DepMap `Model.csv` + `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (already **log2(TPM+1)**).
- Keep **invasive breast carcinoma** cell lines with RNA-seq: **16,500 genes × 63 lines**.
- Drop genes with max log2(TPM+1) &lt; 1.
- **No ComBat** (single source). Built by `Breast Cancer/codes/01_build_input_panelA.py` → `Breast Cancer/input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv`.

### Methods

Same *protocol* as ParTI `algNum=5` (PCHA), not Sisal:

- Fit in the first **(k−1)** PCs of a gene-space PCA.
- `delta=0`; `py_pcha` maxiter=500, conv=1e-6.
- Observed fit: **15** random inits, keep **max volume**.
- Each permutation: **5** inits, keep **max t-ratio**.
- **500** column shuffles (Groves used 1000).
- **dim**: smallest n PCs with **≥50%** variance on *this* matrix → **8 PCs** (Groves used 12 on SCLC because that was ~47–50% there).
- Suggested k: ParTI DimensionFinder elbow on the ESV curve.

Code: `Breast Cancer/Run_2/codes/run_panelA.py`. Params: `PARAMS_vs_Groves.txt`.

### Technicality

Earlier breast runs used 150/15 inits (Sisal `numIter=50`). ParTI **lowers `numIter` to 5 for PCHA**, so Groves’ null is 15 observed / 5 per shuffle. We matched that.

The 2-D picture is a projection. The fit and t-ratio live in **3-D** for k=4.

### What we report

- Elbow **k=4**.
- t-ratio at k=4 ≈ **0.354**, p=**0.082** (500 shuffles). **No k in 3–7 has p&lt;0.05.**
- Groves SCLC k=4 was also NS (p=0.059); they called k=5 because p=0.034. Breast has no such k.
- So this is a **descriptive 4-pole model**, not a significant simplex. Breast looks more continuum-like than SCLC.

Vertices: `results/panel_a/archetypes_k4_parti.npy` (4 × 3). Scores: `pc_scores.npy` (63 × 8).

---

## Figure 1B — PAM50 enrichment at vertices

### Concept

If a vertex is a real subtype extreme, that subtype should pile up in the **closest** distance bin to that vertex (bin 0), not just be “somewhat associated” somewhere.

Cell lines do not have clinical ER/HER2 IHC, so Panel B uses **PAM50** (RNA classifier). That is **not** the label set for Panel C.

### Data

- Same 63-line matrix and **k=4 vertices from Panel A** (all 63 lines still define the polytope).
- PAM50 from genefu `molecular.subtyping(sbt.model="pam50")` after symbol→Entrez (`org.Hs.eg.db`). Labels reused from Run_1: `Breast Cancer/results/panel_b/pam50_labels_panelA.csv`.
- Counts on 63: Basal 27, LumB 17, Her2 14, LumA 4, Normal 1.
- Enrichment **drops LumA + Normal** (n=5); n=58 for the test.

### Methods

Groves-style:

- Euclidean distance to each vertex in the **first 3 PCs** (the PCHA space).
- **5** equal-count bins per archetype (bin 0 = closest).
- Hypergeometric enrichment per (archetype × bin × subtype); BH **q &lt; 0.1**.
- **Hit** = significant **and** the enrichment peak is at bin 0.

Code: `Breast Cancer/Run_2/codes/run_panelB.py`.

### Technicality

Archetype **index is arbitrary**. Run_1 vs Run_2 can flip Arc 1 vs Arc 4 without changing geometry. Report subtype→vertex pairs, not “Arc 4 means X” as a biological name.

Bins of 5 on n=58 are coarse; Her2 can be significant in a non-zero bin without counting as a hit.

### What we report

- **LumB → Arc 3** (bin-0 peak).
- **Basal → Arc 4** (bin-0 peak).
- **Her2**: no bin-0 peak (signal exists off bin 0).
- Arcs 1 and 2 have no significant subtype match under this rule.

---

## Figure 1C — TCGA tumors (IHC ER/HER2)

### Concept

Cell lines might be a culture artifact. Panel C asks whether **primary tumors** live in the **same** polytope. Groves’ answer is not “transform tumors with the cell-line PCA.” Assay shift (lines vs tumors) would then look like a new biological axis. They **ComBat** lines and tumors together, fit a **new** PCA on that matrix, and **overlay the old vertices**.

Labels here are **histopathology**, as requested: ER+/HER2−, ER+/HER2+, ER−/HER2+, ER−/HER2− — **not** PAM50 LumA/LumB from Panel B.

### Data

Under `Breast Cancer/Panel C/`:

- UCSC Xena **HiSeqV2** (already log2; 20,530 genes × 1,218 samples).
- `TCGA.BRCA.sampleMap-BRCA_clinicalMatrix`.

Prep (`codes/prepare_tcga_brca.py`): primary tumors only (barcode `-01` and `sample_type_id=1`) → **1,097** tumors, **14,234** genes shared with Panel A. No extra log.

IHC columns (Positive/Negative):

- `breast_carcinoma_estrogen_receptor_status`
- `lab_proc_her2_neu_immunohistochemistry_receptor_status`

Complete ER and HER2: **725 / 1,097**. Incomplete IHC is grey on the plot.

### Methods (Groves, not the first attempt)

Matches `reference/notebooks/bulk/Cell-line-tumor-batch-correction-and-clustering.Rmd` (saved matrix **bc2**, not fsva) and `Thomas-Tumors-Bulk-Archetypes.ipynb`, and our SCLC `scripts/run_panel_c.py`:

1. Merge 63 lines + 1,097 tumors on shared genes.
2. **`sva::ComBat`**, batch = `cell_line` vs `tumor`, `mod = ~1`, **`ref.batch = cell_line`** (their SCLC `ref.batch` was Minna). R 4.5, sva 3.58 in `Run_2/rlib`. Wrapper: `codes/combat_cellline_tumor.R`.
3. PCA (**20** components) on the ComBat matrix.
4. Inverse-transform Panel A vertices (8-PC cell-line PCA) to gene space → `pca.transform` into the combined PCA. **No PCHA on TCGA.**
5. Left: scatter PC1–PC2, cell lines + tumors by IHC + 4 vertices.  
   Right: tumor-only PCA cumulative variance (ceiling) vs combined PCs scored on tumors vs **20** gene-row shuffles.

Code: `Breast Cancer/Run_2/codes/run_panelC_tcga.py`. Note: `PARAMS_panelC.txt`.

### Technicality — what went wrong first

Projecting Xena tumors **directly** into the DepMap 8-PC PCA (filling missing genes with the PCA mean, no ComBat) put almost every tumor nearest Arc 4 at distances ~400–700. That is **batch**, not biology. Groves never used that projection for Figure 1C.

After ComBat, lines and tumors share an axis scale. Combined PCA is still **tumor-dominated** (1,097 vs 63), unlike Groves (81 tumors vs 120 lines). So “4 combined PCs ≈ 100% of the tumor-only ceiling” is **not** comparable to their “5 PCs ≈ 80%.” The left scatter is the Panel C figure; the right panel is a weak analogue of their variance claim.

Barycentric weights in the 3-D simplex: **0 / 1,097** tumors have all weights ≥ 0. Many points sit **outside** the 4-vertex wireframe in PC1–PC2. Groves SCLC tumors sat inside a 5-vertex polytope. We report that; we do not force tumors into the hull.

### What we report (IHC, after ComBat)

Nearest vertex in combined (k−1)-D, n=725 with IHC:

| IHC | n | Arc 1 | Arc 2 | Arc 3 | Arc 4 |
|---|---|---|---|---|---|
| ER+/HER2− | 436 | 25 | 68 | 187 | 156 |
| ER+/HER2+ | 123 | 4 | 20 | 54 | 45 |
| ER−/HER2+ | 40 | 12 | 9 | 5 | 14 |
| ER−/HER2− | 126 | **82** | 28 | 1 | 15 |

ER−/HER2− (TNBC-like) toward Arc 1; ER+ toward Arc 3/4. That is a **tendency**, not a 1–1 vertex map, and not “all tumors inside the tetrahedron.”

---

## How the three panels lock (breast)

| Panel | Question | Breast answer we actually have |
|---|---|---|
| A | How many extremes in DepMap lines? | Elbow **k=4**; t-ratio **not** significant |
| B | Do PAM50 classes sit at those corners? | LumB–Arc 3, Basal–Arc 4; Her2 unmatched; LumA/Normal dropped |
| C | Do TCGA tumors (ER/HER2 IHC) live there? | After Groves ComBat, yes **same coordinate system**; not all **inside** the 4-simplex; IHC groups separate somewhat by nearest vertex |

If A is only descriptive, B and C are still the Groves *pipeline* on breast data. They are not a claim that invasive breast carcinoma is a significant 4-task Pareto front.

---

## Commands (repo root)

```bash
.venv/bin/python -u "Breast Cancer/Run_2/codes/run_panelA.py"
.venv/bin/python -u "Breast Cancer/Run_2/codes/run_panelB.py"
.venv/bin/python -u "Breast Cancer/Run_2/codes/prepare_tcga_brca.py"
.venv/bin/python -u "Breast Cancer/Run_2/codes/run_panelC_tcga.py" --n-shuffle 20
```

`--force-combat` rebuilds `results/panel_c_tcga/CCLE_TCGA_COMBAT.csv`.
