# GSE114687: MCF10A and HuMEC "pseudospace" EMT (McFaline-Figueroa et al.), archetype analysis

Run 2026-10-15/16 with the SCC / Cook / PDAC / 4T1 pipeline. No predictions file was written; choices added after looking are listed in section 7.
Download record: `data/raw/DOWNLOAD_LOG.md`. Questions: (1) are there task trade-offs / corners, (2) where do hybrid cells lie (specialists or generalists).

## Data
Two human mammary epithelial systems grown as colonies, so EMT progresses with distance from the colony edge ("pseudospace"), each with an untreated (Mock, spontaneous EMT) and a TGF-β arm, inner and outer colony.
30,590 cells after light QC (>= 500 genes, mito <= 15%): MCF10A Mock 4,991, MCF10A TGFB 4,230, HuMEC Mock 11,212, HuMEC TGFB 10,157; 12,002 genes (the CROP-seq screen is not used).
30 of the 9,268 MCF10A metadata cells had no matching matrix barcode and were dropped. Batch unit (`sample`) = cell line + arm. Position (inner/outer) is the independent EMT axis (the role time plays in Cook).
In the 4 samples, cleaning to >= 10,000 UMI keeps 2,536 / 3,127 / 3,889 / 4,102 cells (HuMEC loses about 60%).

## 1. Pooled over the four samples (200 shuffles)
| Genes | k=3 | k=4 | k=5 | k=6 | k=7 | Rule's choice |
|---|---|---|---|---|---|---|
| KS (191) | 0.85 | 0.12 | 0.005 | 0.005 | 0.005 | none (k = 5-7 are past the elbow at 4) |
| All genes | 0.96 | 0.99 | 0.24 | 0.005 | 0.005 | k = 6 (elbow 6); only 1.2% of cells inside the shape |
Corners here are the samples themselves (KS k=4: MCF10A Mock, MCF10A TGFB 94% mesenchymal, two HuMEC corners). Cell line and arm dominate the leading PCs; inner/outer position explains <= 6% of any PC.
Pooled closest-5% rule: 67-89% of the epithelial, early-hybrid and late-hybrid cells are generalists.

## 2. Inside single samples (k=3 and 4; all genes and KS genes)
| Unit | Light QC all k3 / k4 | Light KS k3 / k4 | Cleaned all k3 / k4 | Cleaned KS k3 / k4 |
|---|---|---|---|---|
| MCF10A (both arms) | 0.005 / 0.005 | 0.76 / 0.015 | 0.005 / 0.005 | 0.73 / 0.005 |
| MCF10A Mock | 0.005 / 0.43 | 1.00 / 0.29 | 0.005 / 0.005 | **0.005 / 0.005** |
| MCF10A TGFB | 0.050 / 0.005 | 0.995 / 0.045 | 0.005 / 0.005 | 0.56 / 0.39 |
| HuMEC (both arms) | 1.00 / 0.005 | 0.005 / 0.010 | 0.93 / 0.005 | 0.075 / 0.17 |
| HuMEC Mock | 0.99 / 0.09 | 0.64 / 0.24 | 0.005 / 0.005 | **0.005 / 0.005** |
| HuMEC TGFB | 1.00 / 0.75 | 0.08 / 0.010 | 0.025 / 0.005 | **0.005 / 0.005** |
Fits below p = 0.05: light QC 11 of 24, cleaned 18 of 24. This is the first dataset where the KS EMT genes themselves give a shape inside samples (3 of 6 cleaned units at k=3).
(Cleaning here is >= 10,000 UMI only; no contamination filter, since these are cultured cells.)

## 3. The corners recur in every sample (cleaned, k=3)
| Corner | Cells | Top genes | Hallmark | Colony position |
|---|---|---|---|---|
| Cycling | 92-100% S/G2M | PTTG1, HMGB2, CCNB1, UBE2C, CDC20 | E2F / G2-M top set in 8 of 8 fits | 58-98% outer (HuMEC TGFB 58%, the rest 85-98%) |
| Differentiated epithelial | KS_EMT -0.45 to -0.9 (MCF10A TGFB -0.17); 89-100% epithelial state (MCF10A TGFB 42%) | SPRR1A/B, SBSN, PI3, KRT16, S100A8/9, KRT10 | no listed Hallmark set (KRAS-up weak) | 4-49% outer (mostly inner) |
| EMT / matrix (TGF-β arms) | KS_EMT +0.36 to +0.77, 74-99% mesenchymal state | FN1, LAMC2, TAGLN, CXCL1-3, TGM2, VIM, LAMA3 | EMT top set (MCF10A pooled, MCF10A TGFB, HuMEC TGFB) | 82-97% outer |
| Basal / partial (Mock arms) | KS_EMT about 0 | KRT15, MIR205HG, TGFBI, NDRG1 (CXCL14, CKB, COL17A1 in HuMEC) | weak (top set Hypoxia / EMT at <= +0.07) | mixed |
So the shape is a trade-off between differentiation, proliferation and EMT/matrix programmes, aligned with colony position, and it appears in two cell lines and both arms.

## 4. Specialists and generalists (cleaned k=3; gap rule >= 0.35 / closest-5% rule; % of the state that is generalist)
Early-hybrid 47-77% / 72-99%; late-hybrid 40-86% / 82-97%; epithelial 7-64% / 18-85%; mesenchymal 30-70% / 66-88%.
A hybrid state (early or late) is the most generalist state in 7 of 8 fits by the gap rule and 8 of 8 by the closest-5% rule. The exception under the gap rule is HuMEC Mock (epithelial 64%, early-hybrid 59%, late-hybrid 47%, mesenchymal 32%).
In both Mock fits the mesenchymal state is the least generalist (MCF10A Mock: 55 / 68 / 56 / 44% for epithelial / early / late / mesenchymal). The TGF-β arms are less regular (MCF10A TGFB: epithelial 7%, late-hybrid 86%).
Table: `results/within_unit/generalist_share_by_state.csv`.

## 5. Path versus triangle
Question: is a triangle needed beyond a one-dimensional EMT path? Null = a cubic path in the EMT score (KS_EMT, rank-scaled) with residuals permuted; cleaned fits, k = 3, 200 path-null shuffles.
| Unit | All genes: EMT score explains PC1 / PC2 | p vs path | KS genes: PC1 / PC2 | p vs path |
|---|---|---|---|---|
| MCF10A Mock | 0.85 / 0.16 | **0.005** | 0.88 / 0.36 | **0.005** |
| MCF10A TGFB | 0.36 / 0.35 | **0.005** | 0.76 / 0.35 | 0.050 |
| HuMEC Mock | 0.75 / 0.28 | **0.005** | 0.83 / 0.05 | **0.005** |
| HuMEC TGFB | 0.58 / 0.22 | 0.17 | 0.76 / 0.20 | **0.005** |
| MCF10A (both arms) | 0.81 / 0.24 | 0.075 | 0.90 / 0.11 | 0.30 |
| HuMEC (both arms) | 0.71 / 0.25 | 0.65 | 0.82 / 0.24 | 0.25 |
The EMT score explains 36-90% of PC1 inside a sample (it explained under 6% in 4T1), so there is a genuine EMT axis. In 6 of 12 fits the triangle still beats that path (p = 0.005), including both Mock arms with both gene sets, so a second and third dimension (cycling versus differentiation) adds structure beyond an EMT line.

## 6. Caveats
- One experiment per arm and line; no biological replicates. MCF10A and HuMEC are not tumours.
- Cleaning at >= 10,000 UMI roughly doubles the number of significant fits; it keeps only about 40% of HuMEC cells. Depth correlates with cell cycle, so a cycling corner may partly be a depth effect. Not tested.
- Position explains little of the leading PCs pooled; inner/outer is a coarse two-level axis.
- Corner names come from top genes; the differentiated-epithelial corner is not named by any Hallmark set tried.
- The pooled all-genes k=6 shape contains 1.2% of cells and its corners are samples; do not read it as biology.
- Many fits were run (48 within-unit plus 10 pooled); the significant ones exceed chance, but p = 0.005 is the floor.

## 7. Added after looking (post hoc)
Cleaned within-unit fits (>= 10,000 UMI); the Hallmark check; the path test with the EMT score as the path (not the authors' pseudospace trajectory, which is stored in an R object that was not read); the cell-line-level units (MCF10A, HuMEC).

## Files
`codes/` 01 build, 01c labels, 02 PCA, 03 fits (`run_fits.sh`), 04 describe, 06 Groves, 08 within-unit (`--clean`), 08b corners, 08c path test, 09 Hallmark, 10 label tables, 13 figures.
`results/{ks,all}_main/`, `results/within_unit/` (fit folders, `corner_description_clean.csv`, `hallmark_corners_clean.csv`, `generalist_share_by_state.csv`), `results/within_unit_summary.csv`, `results/mail_figures/` (6 figures).
