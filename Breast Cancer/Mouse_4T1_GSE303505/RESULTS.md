# 4T1 EMT clones (GSE303505, Pattabiraman lab): archetype analysis

Run 2026-10-14/15 with the Mouse SCC / PDAC / Cook pipeline. Predictions written before any fit: `results/PREREGISTERED_PREDICTIONS.md`.
Download record: `data/raw/DOWNLOAD_LOG.md`. The cleaned pooled fit, cleaned within-clone fits, path test and Groves rule were added AFTER
the first results (not pre-registered).

## Data
Six "Total" tumour libraries, one tumour per clone: E, EM1, EM2, EM3, M and parental 4T1. Tumour cells = the authors' `Annotation == Tumor`
(19,930 cells; 11,911 genes after filtering). **Clone = tumour = batch** (n = 1 per clone). CD45 libraries not used.
The EMT score (KS_EMT) separates E clearly (mean -0.83) but EM1, EM2 and parental 4T1 overlap (about -0.03), as do EM3 and M (+0.17, +0.16);
only 37% of its variance is between clones.

## 1. Pooled fits (200 shuffles each; p; "k*" = smallest k with p<0.01 at or before the ESV elbow)
| Cells / genes | k=3 | k=4 | k=5 | k=6 | k=7 | k* |
|---|---|---|---|---|---|---|
| Light QC, KS genes (174) | 0.93 | 0.31 | 0.17 | 0.83 | 0.18 | none |
| Light QC, all genes | 0.12 | **0.005** | 0.015 | 0.005 | 0.020 | 4 |
| **Cleaned**, KS genes | 0.52 | 0.03 | 0.52 | 0.36 | 0.34 | none |
| **Cleaned**, all genes | **0.010** | 0.33 | 0.94 | 0.75 | 0.90 | 3 |
Cleaned = at least 10,000 UMI, immune/stromal panel score <= 0.3, technical genes dropped (13,256 cells; E keeps only 142).
Marker genes never give a simplex. On all genes the k=4 tetrahedron of the light-QC fit disappears after cleaning and a triangle appears.
p = 0.005 is the floor with 200 shuffles; k=3 here is 2 of 201.

## 2. What the corners are
- **Light QC, k=4:** V1 clone E, V2 parental 4T1 (epithelial, partly cycling), V3 clone M (76% M), V4 cycling cells (S+G2M 68%).
  Closest-5% (Groves) specialists: V3 is 98.8% clone M. Hybrid states are 83-85% generalists (epithelial 72%, mesenchymal 80%).
- **Cleaned, k=3 (all genes):** V1 cycling (S/G2M 95%; 4T1 41%, M 31%), V2 epithelial-leaning, G1, 4T1 50% + EM2 26%, V3 clone M (89% M, 49% mesenchymal state).
  Every EMT state is 65-68% interior; no state has a majority at any corner.

## 3. Robustness
- Light-QC k=4: seeds 95-98%; cell cycle removed 80%; leaving out E removes V1 (agreement 65%), leaving out M replaces V3 (heldout 69%).
- Cleaned k=3: seeds 93-96%; leaving out EM1, EM2, EM3, 4T1 or E changes little (88-96%). **Leaving out M** drops agreement to 79% and the M cells to 57%; the mesenchymal corner is replaced by an Early-Hybrid one.
  So the mesenchymal corner is the M clone; the cycling and epithelial-leaning corners survive any one clone.
- Path test (is a triangle needed beyond a 1-D EMT path; cubic in KS_EMT rank, 200 path-null shuffles): pooled cleaned fit p = 0.32 (KS_EMT explains 1.5% of PC1, 6% of PC2), so the pooled triangle is not shown to beat an EMT path.

## 4. Within one clone
Light QC, all clones (p): E k3 0.010 / k4 0.045, M k4 0.010, 4T1 k4 0.010 (all genes); EM1 KS k4 0.005. 5 of 24 below 0.05; their corners were stress genes, immune
(Cd14, S100a8, Lyz2) or stromal (Col1a1) cells, and cycling cells (`results/within_unit/corner_description.csv`).
**Cleaned + decontaminated** (E dropped: 169 cells at >= 10,000 UMI):
| Clone | Cells | all k=3 | all k=4 | KS k=3 | KS k=4 |
|---|---|---|---|---|---|
| EM1 | 1,365 | 0.62 | 0.74 | 0.83 | 0.22 |
| EM2 | 2,701 | 0.63 | 0.98 | **0.015** | 0.10 |
| EM3 | 1,507 | 0.51 | 0.40 | 0.86 | 0.46 |
| M | 3,324 | 0.50 | 0.79 | 0.16 | 0.45 |
| 4T1 | 4,000 (capped) | 0.28 | 0.06 | 0.10 | 0.06 |
One pass in 20 tests (about one expected by chance). Path test inside clones: PCs are not along KS_EMT (R2 <= 0.06); only EM1 KS (0.025) and 4T1 KS (0.010) beat a path, and neither passes the shuffle test.

## 5. Against the predictions
| Prediction | Result |
|---|---|
| 1. k*=3 on KS genes | Fails: no simplex on KS genes at any k, with or without cleaning |
| 2. KS axis orders clones E<EM1<EM2<EM3<M | Partly: E is clearly separate; EM1/EM2/4T1 overlap, EM3/M overlap |
| 3. Epithelial, mesenchymal, and a cycling or intermediate third corner | Yes on cleaned all-genes k=3: cycling, epithelial-leaning, M clone |
| 4. Intermediate clones more interior | Yes: no intermediate clone has its own corner; hybrid states are 65-85% interior/generalist |
| 5. Corners are largely single clones | True for M (89-99% of its specialists); E only in the light-QC fit |

## 6. Comparison with the other datasets (from their own results)
- SCC: KS genes no simplex; pooled all-genes shape partly technical; hybrids generalists; within tumours 2 of 3.
- PDAC: pooled no simplex; within libraries, cleaning REVEALED triangles in 6 of 7 libraries (cycling / epithelial / mesenchymal corners).
- Cook: pooled shape = cell-line identity; within lines 2 of 8 after cleaning; corners = proliferation and TGF-b response.
- 4T1: like Cook and SCC (KS null, pooled structure tied to a clone, within-clone shapes disappear after cleaning); unlike PDAC, where cleaning created shapes.

## 7. Caveats
- One tumour per clone; the M corner cannot be separated from tumour M. The EMT score barely separates five of six clones.
- Cleaning choices (10,000 UMI, panel cut 0.3) are borrowed from Cook/PDAC and the panel cut was set by eye; E cannot be cleaned (142 cells left).
- Light-QC and cleaned pooled fits disagree on k (4 vs 3); both p-values are near the floor. Treat the all-genes simplex as weak.
- PC1 of the light-QC all-genes fit is 71% sequencing depth.
- Not done: cross-projection, state-in-simplex tables inside clones, corner programme scores, MAGIC, CD45 libraries.

## Figures (`results/mail_figures/`, made by `codes/13_mail_figures.py`)
1. `fig1_clone_emt_spectrum.png`: EMT score per clone. E is separate (mean -0.83); EM1, EM2 and parental 4T1 overlap (about -0.03), as do EM3 and M (+0.17, +0.16).
2. `fig2_pooled_permutation_p.png`: permutation p against k for four pooled fits (KS or all genes; light QC or cleaned). Only all-genes fits cross p = 0.01: k=4 and 6 with light QC, k=3 after cleaning. KS genes never do.
3. `fig3_pooled_triangle_cleaned.png`: cleaned all-genes k=3 triangle (p = 0.010), cells coloured by clone, cell-cycle phase and EMT score. V1 cycling, V2 epithelial-leaning, V3 clone M.
4. `fig4_pooled_tetrahedron_lightQC.png`: light-QC k=4 fit (p = 0.005) in PC1-PC2, by clone and phase. A flat view of a 3-D shape, so corners overlap.
5. `fig5_within_clone_p.png`: within-clone p-values, light QC next to cleaned. The five light-QC passes disappear after cleaning except EM2 on KS genes at k=3; E cannot be tested after cleaning.
6. `fig6_within_clone_triangles.png`: each clone's own k=3 triangle, coloured by phase, with p in the titles. The triangles are drawn whether or not they are significant; the cells form one dense mass, not three corners.
7. `fig7_specialists_generalists.png`: Groves closest-5% specialists and generalists for the light-QC k=4 fit (top) and the cleaned k=3 fit (bottom): position on the map, where each EMT state goes (71-86% generalist), and which clone makes up each class (V3 is 89-99% clone M).

## Files
`codes/` 01 build, 01c labels, 01d clean cells, 02 PCA, 02b, 03 fits (`run_fits.sh`), 04 describe, 05 robustness, 06 Groves, 08 within-clone (`--clean --decontam`), 08b corners, 08c path test, 13 figures.
Cleaned pooled runs need `EXPR_NAME=expr_tumour_clean.h5ad`. Results: `results/{ks,all}_{main,clean}/`, `results/within_unit/`, `results/*.csv`.
