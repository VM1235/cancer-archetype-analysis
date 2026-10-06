# Pal triangle refit without TN-0135

Run: `Breast Cancer/codes/28_pal_leaveout_tn0135.py`.
Outputs: `Breast Cancer/results/pal_leave_tn0135/`.
Figure: `Breast Cancer/figures/Figure_pal_leaveout_tn0135.png`.

The original Pal k=3 fit is untouched.

## What was fit

Patient TN-0135 (14,322 cells) was removed from the PCA and from PCHA. The fit is k=3 only, 15 starts, delta=0, same as the full Pal fit. All 15 starts succeeded. Observed t-ratio is 0.498 (the full-Pal t-ratio was 0.660). Shuffles were not rerun, so this is not a new p-value.

New vertex numbers are a new fit. They are not the old Arc 1 / Arc 2 / Arc 3.

## Where the old corners went

Among the training cells (Pal without TN-0135):

| Old nearest vertex | Lands on new vertex |
| --- | --- |
| Old arc 1 (the TNBC-enriched corner, this patient already removed) | 99.4% on new vertex 2 |
| Old arc 3 (mesenchymal) | 99.4% on new vertex 2 |
| Old arc 2 (luminal-epithelial) | 72.3% on new vertex 3, 27.0% on new vertex 1 |

By IHC, in the training cells:

- New vertex 2 holds 97.8% of the remaining TNBC cells.
- New vertex 3 holds 71.3% of ER+.
- New vertex 1 holds 45.2% of HER2+ and 15.1% of ER+. Its closest 5% of training cells are the epithelial tip (mean epithelial score +0.40, mesenchymal −0.20). New vertex 2 is the only tip with a positive EMT score.

The third corner does not stay once TN-0135 is out of the fit. The other cells that used to sit nearest that corner join the mesenchymal vertex. The old epithelial corner splits between two new vertices, one mostly ER+ and one that also takes a large share of HER2+.

## TN-0135 projected back

Same MAGIC matrix, training PCA, no second mean/SD match. 6,206 / 14,322 cells are inside (43.3%). Every cell is nearest new vertex 2 (median weights −0.007, 0.75, 0.25). The cells outside the triangle are just past the edge between new vertices 2 and 3: the median of the most negative weight is −0.007, and only 2.6% have a weight below −0.1. They sit on the mesenchymal vertex. They do not open a corner of their own.

## Other datasets, matched to the training cells

| Dataset | Inside the new triangle | Inside the original full-Pal triangle |
| --- | --- | --- |
| Pal without TN-0135 (70,280) | 65.8% | 75.8% was the full Pal cohort, not this subset |
| TN-0135 (14,322) | 43.3% | same caveat |
| Wu GSE176078 (24,162) | 76.0% | 78.4% |
| GSE173634 (35,271) | 68.8% | 96.1% |
| DepMap (63) | 74.6% | 90.5% |

Wu is essentially unchanged. GSE173634 and the cell lines drop. The host triangle still contains most of each guest, and it contains them less completely than the triangle that was fit with TN-0135 inside it.
