# Pre-registered predictions: archetype fit on 4T1 EMT clones (GSE303505)

Written 2026-10-14 12:39 IST, after building the matrix (step 1) and BEFORE any PCA, PCHA fit or weight.
Do not edit after the first fit; later notes go in a separate file.

## Primary analysis (fixed in advance)
- Cells: 19,930 tumour cells (authors' `Annotation == Tumor`, QC nFeature>=500, mito<=15%) from the six Total tumour libraries:
  E 794, EM1 2,913, EM2 3,381, EM3 1,986, M 4,845, parental 4T1 6,011. One tumour per clone.
- Matrix: log1p CP10k; genes detected in >=1% of cells; mito/ribo removed (11,911 genes). No imputation, no batch correction. Genes centred, not scaled.
- Fit cells: clone-balanced, at most 790 per clone (4,740 cells), seed 0; all other cells projected.
- Gene sets: `all` and `ks` (Tan 2014 KS genes, mouse symbols). Same engine, k = 3..7, 15 inits, 200 shuffles x 5 inits. Specialist = weight gap >= 0.35.
- Labels (never used in fits): clone (E < EM1 < EM2 < EM3 < M per the paper; parental 4T1 unordered), KS_EMT-derived quartile states, KMeans clusters, cell-cycle phase.
- k*: smallest k with permutation p < 0.01 that also sits at or before the ESV elbow.

## Predictions
1. k* = 3 on the KS genes; on all genes a significant p will appear for several k because clone identity is a strong batch axis.
2. On KS genes PC1 is an EMT axis ordering the clones E < EM1 < EM2 < EM3 < M; the parental 4T1 spans the range.
3. Corners: one epithelial (E / EM1 cells), one mesenchymal (M / EM3 cells), and a third that is either an intermediate-clone state (EM2) or cell-cycle driven (G2M enrichment).
4. Intermediate clones (EM1-EM3) are more often interior (purity < 0.35) than E or M.
5. On all genes, corners are largely single clones (clone = tumour, inseparable). Risk stated in advance.

## Outcomes that would count against the predictions
- No k with p < 0.01, or k* > 3 on KS genes.
- Intermediate clones sit at their own corners (then EMT "hybrid" states are clone identity, not mixtures).
- A corner is >50% one clone and disappears when that clone is left out (leave-one-clone-out).
- The third corner vanishes when cell-cycle phase is regressed out (then it is cell cycle, not EMT).
- Within a single clone no triangle (within-clone permutation p >= 0.05) means EMT heterogeneity inside a clone is not simplex-shaped.
