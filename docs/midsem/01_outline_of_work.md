# Outline of work done

**Student:** Vedashree Mahajan (ID f20220216), BITS Pilani, K K Birla Goa Campus  
**Supervisor:** Dr. Mohit Kumar Jolly, Department of Bioengineering, Indian Institute of Science  
**Period covered:** 10 August 2026 – 4 October 2026  
**Repo:** `Project_1` (private GitHub: `VM1235/cancer-archetype-analysis`)

This note is the step-by-step record of what was actually done. Numbers below are taken from saved result tables, not from the earlier emails, where a few percentages were later revised after MAGIC imputation.

---

## a. Concepts and science

### Step 0. The question that started the project

On 10 August 2026 the supervisor set the motivating idea: hybrid epithelial/mesenchymal (E/M) cells may be **generalists** (adaptable, sitting in the interior of a multi-task polytope), while extreme epithelial or mesenchymal cells may be **specialists** (sitting at corners, less plastic).

Two papers define the method.

- Hausser et al., *Nature Communications* 2019. Tumors lie near low-dimensional polytopes because tasks trade off. Cells near a vertex specialize in one task. Cells in the interior are generalists.
- Groves et al., *Cell Systems* 2022. On small-cell lung cancer (SCLC) cell lines, Pareto Task Inference finds a significant 5-vertex simplex. Known SCLC subtypes sit at vertices. Patient tumors, after batch correction, occupy that same cell-line geometry.

The first assignment was to reproduce Groves Figure 1A–C, then ask the same question of breast cancer cell lines, with hybrid E/M as the biological target.

### Step 1. Learn the geometry on SCLC, where the method is known to work

**What a simplex is.** A triangle is a 3-vertex simplex and lives in 2 dimensions. A tetrahedron is 4 vertices in 3 dimensions. A fit with \(k\) archetypes is always done in the first \(k-1\) principal components.

**What an archetype is.** A corner of the expression cloud. Each sample is a mixture of corners. The mixture weights are non-negative and sum to 1 (barycentric coordinates). If any weight is negative, the sample is outside the simplex.

**How the fit is judged.**

- PCHA (principal convex hull analysis), with `delta = 0`, so each corner is itself a mixture of real samples.
- Many random starts; keep the simplex with the largest volume.
- The **t-ratio** is the volume of that simplex divided by the volume of the convex hull of the data. A value near 1 means the cloud really is simplex-shaped. A value near 0 means a small polytope inside a round blob.
- Significance comes from shuffling each principal component independently across samples, refitting, and asking how often the shuffled t-ratio is at least as large as the real one.
- The working rule, taken from Groves: the smallest \(k\) with \(p < 0.05\). If none qualifies, the elbow of the explained-sample-variance curve is reported as a fallback and is not called a significant simplex.
- Panel B does not cluster. Independent labels (subtype, histology) are tested for enrichment in distance-to-vertex bins. A subtype “matches” a vertex only if enrichment peaks in the closest bin.
- Panel C does not refit archetypes on tumors. Cell-line archetypes are carried into a joint cell-line-plus-tumor space after ComBat, and one asks what fraction of tumors fall inside.

**Result.** On the Groves SCLC matrix (about 15,950 genes × 120 lines), with the paper’s initialization settings and 1,000 shuffles:

| \(k\) | t-ratio | \(p\) (this repo) | paper \(p\) | paper t-ratio |
|---|---:|---:|---:|---:|
| 3 | 0.520 | 0.471 | 0.508 | 0.52 |
| 4 | 0.247 | 0.051 | 0.059 | 0.247 |
| 5 | 0.108 | 0.045 | 0.034 | 0.107 |
| 6 | 0.043 | 0.019 | 0.016 | 0.043 |

The t-ratios match the paper. The smallest significant \(k\) is 5, which is the published choice. An early 100-shuffle run had given \(p = 0.02\) at \(k = 5\); the 1,000-shuffle table above is the one to cite.

**What this step was for.** Proof that the Python pipeline (PCA, PCHA, t-ratio, enrichment, ComBat projection) reproduces a published positive result before any new biology is claimed.

### Step 2. Same recipe on breast cell lines, full transcriptome — no significant simplex

DepMap invasive breast carcinoma lines with RNA-seq: 63 lines, about 16,500 genes after filtering. Already log2(TPM+1); not logged again.

500 shuffles, a lighter PCHA search than Groves (15 observed starts and 5 null starts, versus 150 and 50):

| \(k\) | t-ratio | \(p\) |
|---|---:|---:|
| 3 | 0.683 | 0.438 |
| 4 | 0.354 | 0.082 |
| 5 | 0.162 | 0.118 |
| 6 | 0.058 | 0.150 |
| 7 | 0.015 | 0.318 |

The elbow prefers \(k = 4\). No \(k\) has \(p < 0.05\). An early 100-shuffle run had looked borderline (\(p \approx 0.05\)); more shuffles removed that.

Panel B, with PAM50 labels and Luminal A and Normal dropped because of tiny counts: Basal and Luminal B enrich near vertices; HER2 does not. Panel C on TCGA-BRCA (histopathology ER/HER2): tumors occupy related space, but they do not fall inside the four-corner cell-line shape the way SCLC tumors fell inside theirs.

**Meaning.** Breast is not “SCLC with a different name.” A full-transcriptome cell-line simplex is not there to host a hybrid-generalist story.

### Step 3. Glioblastoma as the supervisor’s suggested contrast — also no simplex

DepMap glioblastoma lines (54 with RNA) and TCGA-GBM. Paper-grade PCHA settings (150/50 starts) and 500 shuffles. No \(k\) from 3 to 7 has \(p < 0.05\). Closest is \(k = 7\), \(p = 0.068\), t-ratio \(0.027\) (a tiny simplex). A \(k = 2\) test aimed at the proneural–mesenchymal axis was not significant (\(p = 0.49\)). Marker-based Verhaak-style labels on the lines did not peak at vertices. Zero of 154 TCGA primary tumors fell inside the \(k = 7\) shape.

Wang 2017 gene-list restriction did not create a significant simplex either.

**Meaning.** The pipeline does not invent a simplex whenever it is run. GBM is a negative control, not a failed script. It is not the thesis.

### Step 4. The turn: restrict breast to the KS epithelial and mesenchymal genes

The supervisor corrected the gene list. The intended list was not PAM50. It was the KS epithelial and KS mesenchymal signature of Tan et al. 2014 (206 genes).

On the same 63 DepMap lines, restricted to those genes, \(k = 3\):

- t-ratio \(0.587\), \(p = 0.022\) (500 shuffles, 2 principal components).
- \(k = 4\) is also significant (\(p = 0.018\)) but in a higher-dimensional fit. The Groves rule keeps the smallest significant \(k\), so the working geometry is the triangle.
- PAM50: one vertex matches Basal, another Luminal B. HER2 still has no vertex of its own.
- A 50-gene PAM50-only matrix was not significant (\(p = 0.76\)). The signal is the EMT signature, not “any short gene list.”

Bulk tumors were then projected into this **cell-line** triangle (ComBat, no refit): METABRIC, 1,049 / 1,980 tumors inside (53%), with claudin-low enriched at one vertex; TCGA tumors also separate by ER/HER2, with a larger inside fraction than the unrestricted \(k = 4\) fit. These bulk projections were encouraging and are kept as background. They are not the result the thesis now turns on.

### Step 5. Single cells into the cell-line triangle

Hausser’s Figure 4 logic: freeze the archetypes, project a new dataset in, do not refit.

Two single-cell sets used in Sahoo et al., *iScience* 2024, after MAGIC imputation on the 206 KS genes:

| Guest | What it is | Inside the DepMap \(k = 3\) triangle |
|---|---|---:|
| GSE173634 (Gambardella) | breast **cell-line** scRNA-seq, 35,271 cells | 53.6% |
| GSE176078 (Wu et al.) | breast **tumor** cancer-epithelial cells, 24,162 cells | 45.5% |

Pre-MAGIC runs had looked more contained (about 84% and 69%). MAGIC is the number to cite, because it is the preprocessing the field uses on these sparse matrices and because it is what every later comparison uses. A \(k = 4\) tetrahedron contained far fewer cells and was dropped.

A third atlas was added because two datasets cannot tell you whether the cell-line triangle is general.

### Step 6. Pal et al. as a guest — the cell-line triangle fails

GSE161529 (Pal et al., *EMBO J* 2021): 84,602 epithelial cells from breast tumors, after MAGIC, projected into the frozen DepMap KS triangle.

**11,217 / 84,602 = 13.3% inside.**

Subtype coloring still had a pattern, but the containment claim failed. The supervisor agreed to reverse the fit: let the tumor cells define the simplex, and ask whether the other datasets fall inside that.

### Step 7. Pal as the host — the simplex that contains the others

Fit PCHA on Pal’s own MAGIC KS matrix (84,602 cells × 206 genes).

- Elbow of explained variance suggested \(k = 4\).
- t-ratio at \(k = 3\) is \(0.660\); at \(k = 4\) it is \(0.325\).
- Both beat all 50 shuffles (0/50). With only 50 shuffles the honest statement is “0 of 50 nulls,” which means \(p < 0.02\), not a precise \(p = 0\).
- The \(k = 4\) tetrahedron does not contain the other datasets (Pal self-containment about 38%, Wu about 44%, GSE173634 about 51%). The triangle does.

Fraction of each dataset **inside** each host triangle:

| Dataset | \(n\) | Inside DepMap triangle | Inside Pal triangle |
|---|---:|---:|---:|
| Pal GSE161529 (tumor sc) | 84,602 | 13.3% | 75.8% (self) |
| Wu GSE176078 (tumor sc) | 24,162 | 45.5% | 78.4% |
| GSE173634 (cell-line sc) | 35,271 | 53.6% | 96.1% |
| DepMap (63 bulk lines) | 63 | 52.4% (self) | 90.5% |

**Asymmetry.** The tumor simplex contains cell lines. The cell-line simplex does not contain this tumor atlas. For this signature, who is allowed to define the corners changes the biology. The lines look like a truncated version of the tumor geometry: they cover the epithelial side and miss the mesenchymal extreme.

### Step 8. Name the Pal corners

Two different labels, which must not be collapsed into one sentence.

1. **Geometry.** Which vertex is nearest (largest archetype weight), and whether the cell is a specialist (largest weight beats the second by at least 0.35) or interior.
2. **Gene score.** Within a dataset, hybrid if both KS epithelial and KS mesenchymal scores are at or above the 75th percentile. Otherwise epithelial or mesenchymal by the mesenchymal-minus-epithelial score, split at the tertiles. The remainder is intermediate.

On Pal, vertices from the KS genes of cells nearest each corner:

| Vertex | Clinical enrichment | Gene program at the vertex |
|---|---|---|
| Arc 2 | ER+ (and HER2+ lives here) | luminal epithelial: FOXA1, AGR2, XBP1 |
| Arc 1 | TNBC, about 2.8-fold in the closest bin | mixed; not a pure VIM program |
| Arc 3 | no single IHC peak | mesenchymal (VIM) |

HER2 has no corner of its own. TNBC is not the same thing as EMT: about two-thirds of Pal TNBC cells are nearest Arc 1 and about one-third nearest Arc 3, and that split is largely between patients (TN-0135 entirely nearest Arc 1; TN-B1-0131 entirely nearest Arc 3 and 95.5% mesenchymal by gene score).

### Step 9. The original hypothesis, tested

The 26 September update scored hybrid cells on the Pal triangle and reported that they sit at a corner, especially Arc 1, so the hybrid looks like a third specialist rather than the compromise in the middle. The supervisor’s reply on 27 September asked for the two conditional probabilities, and asked whether TNBC can be both hybrid and mesenchymal (consistent with basal tumors being more EMT-heterogeneous than luminal tumors). The priority was to score all four datasets the same way. That is the 4 October analysis (`26_specialist_vs_interior.py`).

**If a cell is interior, what is its gene score?** Hybrid is 13.4% of Pal interior cells and 14.6% of Wu interior cells. The interior is mostly epithelial, mesenchymal, or intermediate. A generalist is not a hybrid.

**If a cell is hybrid, is it interior?**

| Dataset | Hybrid cells | At a corner | Interior | Where the hybrid mass actually is |
|---|---:|---:|---:|---|
| Pal | 8,719 (10.3%) | 65.6% | 34.4% | Arc 1 specialist. 92% of hybrid Arc-1 specialists are one patient, TN-0135 |
| Wu | 2,413 (10.0%) | 6.3% | 93.7% | Arc 1–Arc 3 edge, not the centroid. Half are one patient, CID4515, and none of that patient’s cells pass the corner cut |
| GSE173634 | 65 (0.2%) | 0% | 100% | the centroid; all from the line CAL851 |
| DepMap | 0 | — | — | no line is in the top quartile of both scores |

On Pal, a hybrid cell is about as likely to be interior (34%) as an epithelial cell (31%) or a mesenchymal cell (27%). The hybrid label does not pick out the middle.

On Wu it does pick out the interior, but that interior is an edge between the hybrid pole and the mesenchymal pole (median weights about 0.51, 0.11, 0.38), not the center of the triangle.

Epithelial samples can be specialists at Arc 2 (luminal / ER+). Mesenchymal samples are interior in Wu (80%), GSE173634 (100%), and DepMap (100%). Basal and TNBC samples are the ones that leave the epithelial corner, and outside Pal they usually do not lock onto a vertex.

**Reading that survives all four datasets.** Hybrid-as-its-own-corner is a Pal result, and inside Pal it is TN-0135. It does not recur. The reproducible structure is an epithelial specialist corner, with mesenchymal and basal programs spread through the interior or along an edge, often in a patient-specific way.

### Step 10. Side threads that are not the thesis spine

- **Anderson et al. 2026**, archetype analysis of premalignant lung lesions, Figure 3. Rebuilt from public counts. Primary calls agree on 44/45 lesions where both the paper and this rebuild assign a primary archetype. Useful as a methods check. Different organ, different question. Not part of the breast story.
- **Hausser Figure 1d** pan-cancer pipeline: written, only partly run, limited by data download. It supplied the language for “freeze the host, project the guest.” It is not a result of this thesis.
- **Survival / NCI-60 generalist scoring** was discussed on 20 August and was not carried out. The geometric specialist-versus-interior score is the version of that idea that was actually finished.

### Where the science stands on 4 October

1. The method reproduces Groves on SCLC.
2. It does not find a simplex in full-transcriptome breast lines or in GBM. That negative is real.
3. On the Tan 2014 KS genes, breast tumors define a triangle that contains cell lines. The cell-line triangle does not contain the Pal atlas.
4. The corners are an epithelial program, a mesenchymal program, and a third vertex that is TNBC-enriched and, in Pal, hybrid-enriched — but that third vertex is dominated by one patient.
5. Hybrid E/M cells are not, in general, the interior generalists the project set out to find. On the one dataset where they look like a specialist corner, the corner is one patient. On the next tumor atlas they sit on an edge. In cell lines they are almost absent, and the mesenchymal state is interior.

---

## b. Tech stack — what the code is, and what each file does

### Environment

- Python 3.9, NumPy kept below 2 because `py_pcha` still uses `np.mat`.
- Virtualenv at the repo root, packages in `requirements.txt`.
- R with Bioconductor `sva` for ComBat on bulk Panel C (a local library can live in `Breast Cancer/rlib/`). PAM50 labels used `genefu`.
- Commands are run from the repository root.

### Shared engine, `src/`

Every disease script adds the repo root to the path and imports these. Disease folders do not reimplement the math.

| File | What it does |
|---|---|
| `src/pca.py` | PCA on samples × genes, cumulative variance, sign alignment, inverse transform back to genes |
| `src/archetypes.py` | PCHA multi-start, keep maximum volume, simplex volume, explained-sample-variance curve, column-shuffle null |
| `src/enrichment.py` | Equal-count distance bins, hypergeometric test, Benjamini–Hochberg, call only if the peak is the closest bin |
| `src/io.py` | Load the Groves SCLC matrices; also a generic genes × samples CSV loader |
| `src/preprocess.py` | Choose DepMap breast or GBM models, orient the matrix, drop low genes |
| `src/combat.py` | Python ComBat if R `sva` is missing |
| `src/paths.py` | Folder locations |

### SCLC reproduction

Folder: `SCLC Reproduction - Groves Cell Systems 2022/`.

| File | What it does |
|---|---|
| `codes/run_panel_a_parti_full.py` | Official Panel A: 1,000 shuffles, paper initializations |
| `codes/run_panel_b.py` | Subtype enrichment at vertices |
| `codes/run_panel_c.py` | Project the 81 tumors; do not refit |
| `codes/export_sclc_reproduction_figures.py` | Figures 1A–C |

Older scripts (`run_panel_a_firstpass.py`, a 12-dimensional fit, `delta = 0.1`) are not the result. Official numbers are `results/t_ratio_official.csv`.

### Breast, default and KS tracks

Folder: `Breast Cancer/`. Official code is `codes/`, not `archive/run1/`.

| File | What it does |
|---|---|
| `codes/01_build_input_panelA.py` | Build the 63-line full-transcriptome matrix from DepMap |
| `codes/00_build_ks_genelist_input.py` | Restrict that matrix to the 206 Tan 2014 genes |
| `codes/run_panelA.py` | Full-transcriptome PCHA and t-ratio |
| `codes/run_panelA_ks_genelist.py` | KS genes, \(k = 3\), 500 shuffles. This is the significant cell-line triangle |
| `codes/run_panelA_ks_genelist_extendedk.py` | KS genes, \(k = 4\) to 7 |
| `codes/03_map_and_pam50.R`, `04_match_pam50_to_panelA.py` | PAM50 labels for the 63 lines |
| `codes/run_panelB.py`, `run_panelB_ks_genelist.py` | Vertex enrichment |
| `codes/prepare_tcga_brca.py`, `prepare_metabric_brca.py` | Bulk tumor matrices matched to the Panel A genes |
| `codes/run_panelC_tcga.py` and the `run_panelC_*ks*` scripts | Project tumors into the **cell-line** space |

### Single cells into the DepMap triangle, then the reverse

| File | What it does |
|---|---|
| `codes/1run_panelA_ks_project_gse173634_sc.py` | MAGIC, then project Gambardella cells into the DepMap triangle |
| `codes/1run_panelA_ks_project_gse176078_sc.py` | Same for Wu cancer-epithelial cells |
| `codes/00_extract_metadata_gse161529.R` | Pal sample metadata from GEO |
| `codes/00_build_matrix_gse161529.py` | Pal count matrix. ER-0001 dropped (barcode mismatch) |
| `codes/1run_panelA_ks_project_gse161529_sc.py` | Project Pal into the DepMap triangle (the 13.3% result) |
| `Sahoo/files/magic_utils.py` | MAGIC imputation used by the `1run_*` scripts |
| `codes/17.9_run_panelA_ks_genelist_gse161529.py` | Fit PCHA on Pal itself |
| `codes/20_t_ratio_gse161529.py` | 50-shuffle t-ratio at \(k = 3\) and \(k = 4\) |
| `codes/18_project_gse176078_onto_gse161529.py` | Wu into the frozen Pal triangle |
| `codes/19_project_gse173634_onto_gse161529.py` | Gambardella into Pal |
| `codes/22_project_depmap_onto_gse161529.py` | The 63 lines into Pal |
| `codes/21_panelB_gse161529.py` | Name Pal vertices by IHC enrichment and by KS genes |
| `codes/23_host_swap_comparison_and_arc1_audit.py` | Side-by-side inside-fractions and the Arc 1 patient audit |
| `codes/24_panelB_emt_hybrid_gse161529.py` | Hybrid gene-score versus barycentric region on Pal |
| `codes/25_host_swap_triangle_panels.py` | The four-dataset, two-host figure |
| `codes/26_specialist_vs_interior.py` | The same specialist/interior rule and the same gene-score rule on all four datasets. This is the analysis the supervisor asked to prioritize |

`1run_*` means MAGIC on. Older `run_panelA_ks_project_*` scripts are the pre-MAGIC versions. Cite `1run_*`.

### Glioblastoma, Anderson, Hausser

| Location | Role |
|---|---|
| `Glioblastoma/codes/01_build_input_panelA.py` through `run_panelC_tcga.py` | Same A/B/C on GB lines and TCGA-GBM |
| `Glioblastoma/codes/02_assign_verhaak.py` | Marker z-score labels. Not published Verhaak calls |
| `Anderson_PML_Archetypes/` | Figure 3 rebuild for the lung premalignancy paper |
| `Hausser/draft 1/`, `Hausser/draft 2/` | Pan-cancer and Figure-4-style drafts. Not the breast result |

### Where the numbers live

| Claim | File |
|---|---|
| SCLC t-ratios | `SCLC Reproduction - Groves Cell Systems 2022/results/t_ratio_official.csv` |
| Breast full transcriptome | `Breast Cancer/results/panel_a/t_ratio_parti_500.csv` |
| Breast KS \(k = 3\) | `Breast Cancer/results/panel_a_ks_genelist/t_ratio_parti_500.csv` |
| Pal t-ratio | `Breast Cancer/results/panel_a_ks_genelist_gse161529/t_ratio_parti_50.csv` |
| Inside fractions, two hosts | `Breast Cancer/results/host_swap_comparison/inside_fractions_k3.csv` |
| Specialist versus interior | `Breast Cancer/results/specialist_interior/summary.json` and `Breast Cancer/docs/RESULT_compare_four_datasets.md` |
| Threshold sensitivity of the 0.35 cut | `Breast Cancer/results/specialist_interior/threshold_sensitivity.csv` |
| Per-dataset write-ups | `Breast Cancer/docs/RESULT_pal_triangle.md`, `RESULT_wu_gse176078.md`, `RESULT_gse173634.md`, `RESULT_depmap_cell_lines.md` |

### Figures that belong in a talk

Use these. They are the ones the later analyses agree on.

- `Breast Cancer/figures/Figure_1A_breast_ks_genelist.png` — cell-line KS triangle
- `Breast Cancer/figures/Figure_4_ks_panelA_gse161529_sc_magic.png` — Pal failing to sit inside it
- `Breast Cancer/figures/Figure_host_swap_triangles_k3.png` — both hosts, four datasets
- `Breast Cancer/figures/Figure_host_swap_inside_fractions_k3.png` — the same comparison as bars
- `Breast Cancer/figures/Figure_4_gse176078_on_gse161529_k3.png`, `Figure_4_gse173634_on_gse161529_k3.png`, `Figure_4_depmap_on_gse161529_k3.png` — guests inside Pal
- `Breast Cancer/figures/Figure_1B_gse161529_named_scatter.png`, `Figure_1B_gse161529_ks_vertices.png` — what the corners are
- `Breast Cancer/figures/Figure_emt_hybrid_pal_triangle.png` — hybrid cells on the Pal triangle
- `Breast Cancer/figures/Figure_specialist_vs_interior_hybrid.png`, `Figure_program_given_specialist.png` — the conditional-probability result

SCLC Figure 1A is the methods proof, one slide. GBM and full-transcriptome breast are one “what did not work” slide, not a tour of every negative panel.

### What is not in git

DepMap expression (~500 MB), Xena HiSeq, METABRIC, raw GEO dumps under `Sahoo/`, MAGIC `.npz` intermediates, the authors’ MATLAB clone in `reference/`, and ComBat matrices. Processed Panel A matrices, small tables, figures, and code are in the repo. A clone can replot the saved results. It cannot rebuild Panel C from raw data without those downloads.
