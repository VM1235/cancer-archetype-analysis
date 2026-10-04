# What was weak, and what the better approach would have been

This is a critique of the work in `01_outline_of_work.md`, not a list of things to apologize for in the presentation. Several of these points are already visible in the result files. The midsem report should say them in one paragraph each. Hiding them would be worse than the gaps themselves.

---

## a. Concepts and science

### 1. The first breast and GBM fits asked the Groves question in the wrong feature space

Groves found a simplex in the SCLC transcriptome because SCLC subtypes are the dominant axes of that matrix. Breast subtypes and EMT are not the dominant axes of 16,500 genes in 63 lines. Running the identical pipeline and treating a non-significant \(k = 4\) (\(p = 0.082\)) as a puzzle to be “fixed” spent a week. The better first breast experiment, given the supervisor’s actual question, was the KS epithelial/mesenchymal list. That list was in the email’s motivating idea from day one. It was not used until 20 August, after a PAM50 gene-list detour that the supervisor had to correct.

GBM was a reasonable contrast to run once. It was not a reasonable place to keep tuning \(k\), gene lists, and Verhaak markers after the full transcriptome and the Wang list had both failed. The negative result is the result.

### 2. “Significant” was used at three different standards

| Fit | Shuffles | PCHA starts (observed / null) | How it was described |
|---|---:|---|---|
| SCLC official | 1,000 | 150 / 50 | matches the paper |
| Breast, all tracks | 500 | 15 / 5 | “\(p = 0.022\)” at KS \(k = 3\) |
| Pal | 50 | 15 / 5 | easy to say “\(p = 0\)” |

The breast KS triangle is the first positive breast result, and 500 shuffles is enough to report \(p = 0.022\). The lighter PCHA search means a different random-start budget than Groves. That should be one sentence in the methods, not a footnote discovered by a reviewer.

Pal’s \(p = 0.0\) in `t_ratio_parti_50.csv` is 0 successes out of 50. The smallest resolvable \(p\) is \(< 0.02\). Writing \(p = 0\) overstates it. Elbow-preferred \(k = 4\) and t-ratio-preferred \(k = 3\) were both real; choosing \(k = 3\) because the tetrahedron does not contain the other datasets is a good reason, and it should be stated as a containment decision, not as “the elbow was wrong.”

### 3. Early containment numbers were preprocessing-dependent, and the emails got ahead of the final pipeline

On 28 August, before MAGIC, GSE173634 was described as ~84% inside and Wu as ~69% inside. After MAGIC the same projections are 53.6% and 45.5%. The scientific conclusion (partial containment) survived. The percentages did not. Any report or slide that still quotes the August numbers is wrong.

Bulk TCGA/METABRIC “percent inside” also depends on ComBat and on which \(k\) is drawn. Those bulk numbers are supporting context. They are not interchangeable with the single-cell host-swap table.

### 4. Two labels were allowed to sound like one label

“Arc 1 is the hybrid corner” came from KS genes on cells nearest that vertex. “This cell is hybrid” came later from a within-dataset double-quartile rule. They disagree often. In Pal, TN-B1-0554 is 96.7% nearest Arc 1 and 78.4% mesenchymal by the gene score. Nearest the hybrid corner is not “this cell turned both programs on.”

The 4 October write-up keeps the two labels apart. Earlier emails did not. The report must use the later distinction.

### 5. The hybrid definition does not mean the same cells would be hybrid in another dataset

Quantiles are computed inside each dataset. A DepMap line can be extreme on a mesenchymal score relative to tumors and still fail to be “hybrid” among the 63 lines, or the reverse. That is why DepMap has zero hybrid lines while Pal has 10.3% hybrid cells, and why epithelial–mesenchymal score correlations have opposite signs (Pal \(r = +0.29\); Wu \(-0.29\); GSE173634 \(-0.77\); DepMap \(-0.69\)). Positive correlation in Pal is partly what creates a large joint-high set. A better sensitivity, not yet done, is one shared score scale (z-score against a common reference, or a fixed gene-set cutoff) so “hybrid” is a biological call rather than a rank within the file.

### 6. The specialist cut is a choice, and one alternative would reverse the Wu sentence

The rule used, because it reproduces the “66% / 53%” figures already sent on 26 September, is: after clipping negative weights at 0 and renormalizing, specialist if the largest weight beats the second by at least 0.35.

`threshold_sensitivity.csv` shows the qualitative split is stable for Pal (hybrid cells stay mostly specialists from a 0.20 cut through a 0.50 cut) and for GSE173634 (the 65 hybrid cells are never specialists). It is not stable for Wu if someone replaces purity with “largest weight ≥ 0.50”: under that rule 87% of Wu hybrid cells would be called specialists, mostly at Arc 1. Their median purity is only 0.12, so the largest weight is about 0.51 with Arc 3 close behind. They are on an edge. Calling them specialists because the leading weight exceeds one half would erase the edge, which is the actual observation.

The better approach is to lead with the continuous weights (median weight vector, purity, entropy) and treat 0.35 as one displayed cut, with the sensitivity table in the supplement. The continuous numbers are already in `summary.json`. The talk should not depend on a single threshold.

### 7. Patient structure can manufacture a vertex

TN-0135 is 14,322 of 84,602 Pal cells, all nearest Arc 1, all specialists, median weights about 0.88, 0.01, 0.11. Of the hybrid cells that are Arc 1 specialists, 92% are this patient. Wu’s lookalike, CID4515, is the same clinical story (TNBC, nearest Arc 1, high hybrid fraction) and the opposite geometry (purity too low to be a corner).

Until Pal is refit without TN-0135, “Arc 1 is a hybrid task” is an observation about one library that happens to be large. PCHA will put a corner on a tight, extreme cloud of thousands of cells from one patient. That is the single most important scientific gap. The audit in `23_host_swap_comparison_and_arc1_audit.py` and `RESULT_pal_triangle.md` already shows the dependence. It does not yet remove the patient and refit.

### 8. The generalist the papers define is not quite the generalist that was scored

Hausser’s generalist is a position in the polytope: low specialization, high task entropy, often associated with worse survival in their pan-cancer analysis. The score used here is a hard purity cut, crossed with a hybrid gene-set call. That answers the supervisor’s conditional-probability question, which was the right next measurement. It does not yet answer “are hybrid cells generalists?” in the Hausser sense, because:

- interior includes edges, where a cell is a mixture of two tasks, which is a different biological claim from “equal mixture of all three”;
- there is no survival, metastasis, or plasticity readout attached to the score;
- the NCI-60 / EMT-score analysis discussed on 20 August was never run.

The better statement, which the data already support, is narrower: under this gene list and this triangle, hybrid E/M is not the interior class. Epithelial cells can be specialists. Mesenchymal and basal cells are usually not.

### 9. Cell-line self-containment near 50% was easy to over-read

DepMap is only 52.4% inside its own triangle. A simplex fit in 2-D with `delta = 0` is not supposed to leave half the points outside if the cloud is a filled triangle; points near the hull but off the faces, and the difference between the PCHA simplex and the convex hull, both contribute. Pal’s self-containment is 75.8%, which is healthier, and the comparison that matters is the off-diagonal (Pal inside DepMap 13.3%, DepMap inside Pal 90.5%). Slides should show the off-diagonal, not “52% of lines sit in their own fit” as a success.

### 10. Side projects diluted the claim without adding a control the claim needed

Anderson Figure 3 shows the same family of methods can be rebuilt. Hausser Figure 1d was requested and only drafted. Neither tests the hybrid-generalist idea. A methods reproduction is appropriate in a thesis introduction chapter. It is not a second result. The methylome/MIDAA idea raised on 26 September was correctly deprioritized by the supervisor. Opening it before the leave-one-patient-out would be the same mistake as the August gene-list scatter.

### 11. What the work got right, so the critique stays attached to the real argument

The host swap was the correct experiment, and it was proposed before being told to do it. The four-dataset conditional probabilities are exactly the quantities asked for on 27 September. The patient audit stops the hybrid-specialist claim from being oversold. The GBM negative stops “PCHA always finds a triangle” from being an objection. Those are the parts a thesis can stand on. The gaps above are what keep it from being a paper this semester.

---

## b. Technical and engineering

### 1. The script tree records the history of the project instead of the analysis

Names like `1run_`, `17.9_`, `run_panelA_ks_genelist_extendedk.py`, and both `Hausser/draft 1` and `draft 2` are recoverable only because `Breast Cancer/docs/SCRIPTS.md` and `docs/LEARNING_WALKTHROUGH.md` were written later. A reader who starts from `PROJECT_CONTEXT.md` will miss the Pal reverse result entirely: that file still describes GBM as the latest work and breast as “mixed/weak significance.” It was not updated when the science moved.

Better approach, still worth doing before the final thesis: one `README` status block of ten lines that says which script is canonical, and a banner at the top of `PROJECT_CONTEXT.md` pointing at the midsem notes and `RESULT_compare_four_datasets.md`. Do not rename the historical scripts now. Renaming breaks every path in the docs.

### 2. Reduced PCHA settings were copied forward as if they were the method

Breast Panel A used `numIter = 5` to finish 500 shuffles in reasonable time. That setting then became the default for the KS triangle and for Pal. GBM, run with the paper’s 150/50, is the only non-SCLC fit at the published search depth. A reviewer will ask whether the breast \(p\)-values move if the search is deepened. The better approach would have been: paper settings for the one fit the claim depends on (DepMap KS \(k = 3\), and Pal \(k = 3\)), and the light settings only for scans over \(k\).

Pal’s 50 shuffles were a time compromise on 84,602 cells. That is understandable. It has to be written as 0/50, and a 500-shuffle confirmation is a methods item for the next two months, not a new biological aim.

### 3. There was no single results table until the end

Each email attached a new PNG and a new percentage. MAGIC changed the percentages. \(k = 4\) was tried and dropped. The comparison table in `inside_fractions_k3.csv` and the four `RESULT_*.md` files are the first time the claim is one claim. The better habit, from the METABRIC week onward, would have been to append a row to one CSV (`dataset, host, preprocessing, k, n, frac_inside, t, n_perm`) every time a projection finished, and to quote only that CSV.

### 4. Figure folders contain the rejected analyses next to the result

Pre-MAGIC Figure 4s, \(k = 4\) projections, PAM50-genelist panels, and GBM \(k = 6\) panels are all valid lab notebook. They are dangerous in a talk, because a figure titled like the real result will be read as the real result. The presentation prompt lists the figures that are allowed. The thesis should put the others in a supplement titled as sensitivity or as negative controls.

### 5. Large data and gitignore are correct, and they are also a reproducibility hole

A clean clone cannot rebuild ComBat or MAGIC. The saved projection CSVs and the specialist score tables can rebuild the argument. The better state for submission is: a `results/` manifest that says which files are derived and which raw inputs they need, and a one-command replot script for the host-swap and specialist figures from those CSVs. That script does not exist as a single entry point. `25` and `26` are close.

### 6. Batch correction and gene matching were done differently for every guest

Bulk Panel C uses ComBat with the cell-line batch as reference. Single-cell projection into DepMap uses MAGIC plus mean/SD matching and a joint PCA. Projection into Pal uses Pal’s PCA transform and fills the missing gene `LHFPL6` with the Pal mean. These are reasonable, and they are not the same operation. “Percent inside” is therefore not a pure biological distance between datasets. The host-swap comparison is still informative because both directions were computed, so a global “tumors always fall outside cell lines” artifact would have to explain why the reverse direction goes inside. The better sentence is: the asymmetry survives two different projection recipes. A still better check, not done, is to project with one frozen recipe in both directions (same scaling, same genes, no joint PCA refit) and show the asymmetry remains.

### 7. Tests that would have been cheap were not automated

- Sign-flip of principal components changes which corner is “left” on a plot. Naming was done carefully in the later scripts; a unit check that vertex gene ranks are stable under a sign flip was not.
- The 0.35 cut’s sensitivity was computed (good) but only after the claim had already been emailed.
- No leave-one-patient-out hook was built into `17.9` even after TN-0135 was recognized.

### 8. What not to call a technical failure

NumPy being pinned below 2, R `sva` living in a local library, and DepMap files being gitignored are constraints, not mistakes. The SCLC t-ratios matching the paper to three digits is evidence the shared `src/` engine is doing the ParTI calculation. Do not reopen that engine unless a new fit disagrees with a saved table.
