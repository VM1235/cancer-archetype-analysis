# Mid-semester report

**Who defines the corners? Archetype geometry of an epithelial–mesenchymal signature in breast cancer, and a test of the hybrid-generalist hypothesis**

Vedashree Mahajan  
ID No. f20220216  
BITS Pilani, K K Birla Goa Campus  
First Degree Thesis, Semester I, 2026–27  
Supervisor: Dr. Mohit Kumar Jolly, Associate Professor, Department of Bioengineering, Indian Institute of Science, Bengaluru  
Work period: 10 August 2026 – 4 October 2026  
Repository: private, `cancer-archetype-analysis` (local folder `Project_1`)

---

## Abstract

Pareto Task Inference represents a cohort as a simplex whose vertices are extreme expression programs. Specialists are samples near a vertex. Generalists are samples in the interior. This project asked whether that geometry, established for small-cell lung cancer, describes epithelial–mesenchymal (EMT) programs in breast cancer, and whether hybrid epithelial/mesenchymal cells are the generalists.

I reproduced the published SCLC analysis: the volume ratio at five archetypes matches Groves et al. (t-ratio 0.107), and five is the smallest significant number of vertices. The same procedure on the full transcriptome of 63 breast cancer cell lines, and on 54 glioblastoma cell lines, finds no significant simplex. Restricting the breast lines to the 206-gene KS epithelial and mesenchymal signature of Tan et al. yields a significant triangle (\(p = 0.022\)).

Single-cell tests of that cell-line triangle diverge by cohort. After MAGIC imputation, 53.6% of GSE173634 cells and 45.5% of Wu et al. tumor cells fall inside, but only 13.3% of 84,602 epithelial cells from Pal et al. do. Reversing the host — fitting the triangle on Pal and projecting the others in without refitting — puts 78.4% of Wu cells, 96.1% of GSE173634 cells, and 90.5% of the DepMap lines inside. The asymmetry is the main geometric result: for this signature, the tumor atlas contains the cell-line geometry, and the cell-line geometry does not contain the tumor atlas.

The Pal vertices are not the three clinical subtypes. One vertex is a luminal epithelial program and contains both ER+ and HER2+ cells. A second is mesenchymal. A third is enriched for TNBC. HER2 has no vertex. Scoring hybrid status from the KS genes and specialist status from the archetype weights rejects a simple hybrid-as-generalist rule. Interior cells are rarely hybrid. Pal hybrid cells are mostly specialists at the third vertex, and that specialist set is dominated by one patient (TN-0135). Wu hybrid cells lie on the edge between that vertex and the mesenchymal vertex. Cell-line datasets have almost no jointly high hybrid samples, and their mesenchymal samples are interior. The next required analysis is to refit Pal without TN-0135.

---

## 1. Introduction

### 1.1 Task trade-offs and archetypes

Tumors have to perform several tasks, and those tasks can trade off: a program that is good at one can be poor at another. Hausser et al. (2019) argued that this trade-off leaves tumors near a polytope. The vertices are specialist programs. Points inside are compromises, and those compromises were linked to generalist behavior. The fitting method is archetypal analysis, in the form used by Pareto Task Inference: each sample is a convex combination of a small number of archetypes, and each archetype is a convex combination of real samples.

Groves et al. (2022) made this concrete in small-cell lung cancer. Cell lines form a five-vertex simplex that is unlikely under a permutation null. The known SCLC subtypes sit at vertices. Patient tumors, placed into that geometry after batch correction, occupy the same simplex. The cell-line fit is the host. Tumors are the guests. That host–guest order is part of the scientific claim, not only a computational convenience.

### 1.2 Why EMT, and why hybrid cells

Epithelial–mesenchymal transition is a natural place to ask the same question. Full epithelial and full mesenchymal states are often treated as extremes. Hybrid E/M states are observed, are associated with plasticity, and have been discussed as more adaptable than either extreme. The hypothesis at the start of this project, from the supervisor, was specific: on an EMT polytope, hybrid cells would be the generalists in the interior, and the extreme states would be the specialists at the corners.

Breast cancer is a useful test. It has epithelial and mesenchymal variation, clinical subtypes that are not the same thing as EMT, public cell lines, and several single-cell atlases. It is also a test of whether the SCLC pattern is general. If every carcinoma formed a cell-line simplex that tumors then filled, the method would be a protocol. If not, the interesting object is the condition under which a simplex appears and which dataset is allowed to set the corners.

### 1.3 Aims for this semester

1. Reproduce Groves Figure 1A–C well enough that later fits can be trusted.
2. Ask whether breast, and as a contrast glioblastoma, have the same cell-line-then-tumor structure.
3. If the full transcriptome does not, ask whether the KS EMT signature does.
4. Test that signature in single cells, and test which dataset should be the host.
5. Score hybrid E/M against specialist versus interior, which is the original hypothesis.

---

## 2. Methods

### 2.1 Geometry

Samples are columns of a gene-expression matrix. Principal component analysis reduces that matrix. A simplex with \(k\) vertices is fit in the first \(k-1\) components, because that is the dimension in which \(k\) points span a simplex.

Fitting uses principal convex hull analysis (PCHA) with `delta = 0`, so archetypes stay inside the data. The search is multi-start. The start with the largest simplex volume is kept. The **t-ratio** is that volume divided by the volume of the convex hull of the samples, in the same \(k-1\) dimensions. Values near 1 mean the cloud fills a simplex. Values near 0 mean the fitted simplex is a small object inside a larger cloud.

The null destroys dependence between axes: each component is permuted across samples independently, the fit is repeated, and the t-ratio is recomputed. The \(p\)-value is the fraction of finite null t-ratios at least as large as the observed one. The decision rule follows Groves: take the smallest \(k\) with \(p < 0.05\). If none qualifies, the elbow of the explained-sample-variance curve is reported as a description, not as a discovered number of tasks.

**Barycentric coordinates** are the mixture weights. They are non-negative and sum to 1 for a point inside the simplex. A negative weight means the point is outside. “Percent inside” is the fraction of samples with all weights non-negative.

**Specialist versus interior**, used in the final analysis, is stricter than “which vertex is nearest.” Negative weights are clipped at 0 and the vector is renormalized. Purity is the largest weight minus the second-largest. A sample is a specialist at its leading vertex if purity is at least 0.35, and interior otherwise. Interior therefore includes both the center and the edges. The cut 0.35 was chosen because it is the rule behind the figures already discussed with the supervisor. A sensitivity table over other cuts is saved and is summarized in the discussion.

### 2.2 Gene-score labels

These are computed without using the archetype weights. Within each dataset, a sample is **hybrid** if its KS epithelial score and its KS mesenchymal score are both at or above the 75th percentile. Among the remaining samples, **epithelial** and **mesenchymal** are the lower and upper thirds of the mesenchymal-minus-epithelial score. The rest are **intermediate**. Quartiles and tertiles are inside the dataset. A hybrid call in Pal is not the same absolute expression state as a hybrid call in Wu.

### 2.3 What is projected, and what is not refit

Panel C and all single-cell overlays freeze the archetypes. Guests are scaled into the host gene space, placed in the host principal-component space (or a joint space, for the earlier DepMap projections), and assigned barycentric coordinates. PCHA is not run on the guest. This is the Hausser Figure 4 design: one shape, many datasets drawn on it.

Bulk tumor projections (TCGA, METABRIC) use ComBat with the cell-line batch as the reference, matching the Groves batch model. Single-cell projections use MAGIC on the 206 KS genes and mean/variance matching. Those are different recipes. The host-swap comparison is interpreted as an asymmetry that is visible under both, not as a single calibrated distance.

### 2.4 Datasets

| Dataset | Role | Size used |
|---|---|---|
| Groves SCLC cell lines and tumors | Methods reproduction | ~120 lines; 81 tumors in the published Panel C |
| DepMap invasive breast carcinoma | Cell-line host | 63 lines with RNA-seq, log2(TPM+1) |
| DepMap glioblastoma | Negative control | 54 lines |
| TCGA-BRCA, TCGA-GBM, METABRIC | Bulk guests | METABRIC: 1,980 tumors in the KS projection |
| GSE173634 | Cell-line scRNA-seq guest | 35,271 cells |
| GSE176078 (Wu et al.) | Tumor scRNA-seq guest, cancer epithelial | 24,162 cells |
| GSE161529 (Pal et al.) | Tumor scRNA-seq; guest, then host | 84,602 epithelial cells. Sample ER-0001 dropped because of a barcode mismatch |

The KS list is the epithelial plus mesenchymal signature of Tan et al. (2014), 206 genes, taken from that paper’s supplement. The full-transcriptome breast matrix is about 16,500 genes after a low-expression filter. DepMap and Xena matrices are already log-transformed and were not logged again.

### 2.5 Statistical settings, stated because they differ by fit

| Fit | Permutations | PCHA starts, observed / null |
|---|---:|---|
| SCLC, official | 1,000 | 150 / 50 (Groves settings) |
| Breast, full transcriptome and KS | 500 | 15 / 5 |
| Glioblastoma | 500 | 150 / 50 |
| Pal | 50 | 15 / 5 |

The breast and Pal searches are lighter than the SCLC reproduction. I treat the SCLC match as the check on the code, and I treat the breast \(p = 0.022\) as a 500-permutation result under the lighter search. For Pal, 0 of 50 nulls means \(p < 0.02\). I do not report \(p = 0\).

### 2.6 Implementation

Shared Python functions for PCA, PCHA, the volume ratio, and bin enrichment live in `src/`. Disease folders call that library. ComBat for bulk data uses R `sva`. PAM50 labels for the 63 lines use `genefu`. MAGIC follows the utility used for these datasets in Sahoo et al. (2024). Canonical scripts and result paths are listed in `Breast Cancer/docs/SCRIPTS.md` and in `docs/midsem/01_outline_of_work.md`.

---

## 3. Reproduction of the SCLC result

On the published, already batch-corrected SCLC matrix, the t-ratios are:

| \(k\) | t-ratio | \(p\) (1,000 permutations) | Paper t-ratio | Paper \(p\) |
|---|---:|---:|---:|---:|
| 3 | 0.520 | 0.471 | 0.52 | 0.508 |
| 4 | 0.247 | 0.051 | 0.247 | 0.059 |
| 5 | 0.108 | 0.045 | 0.107 | 0.034 |
| 6 | 0.043 | 0.019 | 0.043 | 0.016 |

The volume ratios match. The \(p\)-values are close but not identical, which is expected of a permutation test. The smallest \(k\) below 0.05 is 5, so the reproduced conclusion is the published one. Subtype enrichment and the tumor projection were run as Panels B and C. This section is a calibration. It is not a new SCLC finding.

An earlier email quoted \(p = 0.02\) at \(k = 5\) from a 100-permutation run. The table above replaces that number.

---

## 4. Transcriptome-wide fits that do not generalize

### 4.1 Breast cell lines

Sixty-three invasive breast carcinoma lines, about 16,500 genes, 500 permutations:

| \(k\) | t-ratio | \(p\) |
|---|---:|---:|
| 3 | 0.683 | 0.438 |
| 4 | 0.354 | 0.082 |
| 5 | 0.162 | 0.118 |
| 6 | 0.058 | 0.150 |
| 7 | 0.015 | 0.318 |

No \(k\) is significant. A 100-permutation pilot had looked borderline at \(k = 4\). The 500-permutation run is the one I use. PAM50 enrichment, after dropping Luminal A and Normal for lack of samples, still places Basal and Luminal B nearer vertices than HER2, but that description sits on a non-significant simplex. TCGA-BRCA tumors, projected after ComBat and colored by ER/HER2 immunohistochemistry, do not fall inside the four-vertex shape.

A further restriction to the 50 PAM50 genes made the fit worse (\(p = 0.76\) at \(k = 3\) and at \(k = 4\)). A short gene list is not automatically a simplex.

### 4.2 Glioblastoma

This was run because proneural and mesenchymal programs are a plausible two-task axis, and because the supervisor suggested TCGA-GBM as a comparison. With paper-grade PCHA settings, no \(k\) from 3 to 7 has \(p < 0.05\). The closest is \(k = 7\), t-ratio 0.027, \(p = 0.068\): a small simplex and a null that is not beaten. A dedicated \(k = 2\) fit has \(p = 0.49\). Marker scores in the style of the Verhaak subtypes do not peak in the closest distance bin. None of 154 TCGA primary tumors fell inside the \(k = 7\) shape. The Wang et al. (2017) gene list did not produce a significant fit either.

I read this as a negative biological result under this protocol, not as a broken pipeline. The SCLC reproduction is the evidence that the pipeline can return the published positive. GBM is not developed further in this thesis.

---

## 5. The KS signature in cell lines

On the supervisor’s correction, the breast matrix was restricted to the Tan et al. KS epithelial and mesenchymal genes (206 genes). At \(k = 3\), t-ratio \(0.587\), \(p = 0.022\) (500 permutations, fit in 2 components). At \(k = 4\), \(p = 0.018\) in a six-component setup. The smallest significant \(k\) is 3, so the working cell-line geometry is a triangle.

Panel B on that triangle associates vertices with Basal and Luminal B. HER2 remains without a vertex. This is the first breast fit I am willing to call a simplex rather than an elbow.

Bulk projections into this frozen triangle were more contained than the transcriptome-wide Panel C. In METABRIC, 1,049 of 1,980 tumors (53%) fall inside, and claudin-low tumors enrich at one vertex, a split that ER/HER2 labels in TCGA cannot show. I treat METABRIC and TCGA as supporting bulk context. The single-cell host swap below is the sharper test, because it uses the same gene list and can be run in both directions.

---

## 6. Single cells, and the decision to change the host

MAGIC-imputed KS genes were projected into the DepMap triangle.

| Guest | \(n\) | Inside DepMap \(k = 3\) |
|---|---:|---:|
| GSE173634 | 35,271 | 53.6% |
| Wu GSE176078 | 24,162 | 45.5% |
| Pal GSE161529 | 84,602 | 13.3% |

Before MAGIC, the first two cohorts had looked more contained (about 84% and 69%). Those numbers are superseded. The \(k = 4\) tetrahedron contained substantially fewer cells and was not pursued as the cell-line host.

Pal is the informative failure. It is a large epithelial tumor atlas, processed in the same way, and it does not live in the cell-line triangle. Subtype coloring of the few cells that do fall inside is not a substitute for containment. With the supervisor’s agreement I fit archetypes on Pal itself.

---

## 7. Pal as host

### 7.1 Choosing three vertices

On 84,602 Pal cells and the 206 MAGIC KS genes, the explained-variance elbow prefers \(k = 4\). The t-ratios prefer the triangle as the filled object: 0.660 at \(k = 3\) and 0.325 at \(k = 4\), each larger than all 50 permutation nulls. Self-containment and guest-containment decide the issue. At \(k = 4\), Pal contains only about 38% of its own cells, Wu about 44%, and GSE173634 about 51%. At \(k = 3\):

| Dataset | \(n\) | Inside DepMap triangle | Inside Pal triangle |
|---|---:|---:|---:|
| Pal GSE161529 | 84,602 | 13.3% | 75.8% |
| Wu GSE176078 | 24,162 | 45.5% | 78.4% |
| GSE173634 | 35,271 | 53.6% | 96.1% |
| DepMap, 63 lines | 63 | 52.4% | 90.5% |

DepMap’s 52% self-containment is easy to misread as success. The comparison that carries the result is the off-diagonal: Pal is outside the cell-line triangle, and the cell lines are inside the Pal triangle. Cell lines look like a subset of the tumor geometry. They do not define it. In particular they under-represent the mesenchymal extreme of the Pal fit.

### 7.2 What the vertices are

IHC enrichment and KS genes at the vertices, from the Pal fit:

| Vertex | Clinical signal | Genes at the vertex |
|---|---|---|
| Arc 2 | ER+ enriched; HER2+ occupies this vertex | luminal epithelial (FOXA1, AGR2, XBP1) |
| Arc 1 | TNBC enriched (about 2.8-fold in the closest bin) | mixed; not a VIM program |
| Arc 3 | no single IHC class | mesenchymal (VIM) |

Among Pal cells, ER+ cells are 86% nearest Arc 2. HER2+ cells are 99.7% nearest Arc 2. TNBC cells are 66% nearest Arc 1 and 33% nearest Arc 3. That TNBC split is mostly a split among patients. TN-0135 (14,322 cells) is entirely nearest Arc 1. TN-B1-0131 (5,491 cells) is entirely nearest Arc 3 and 95.5% mesenchymal by the gene score. This matches the expectation that basal/TNBC tumors are more EMT-heterogeneous than luminal tumors, which sit on the epithelial vertex. It also means “TNBC” is not the name of a task.

Nearest-vertex and gene-score are different. TN-B1-0554 is 97% nearest Arc 1 and 78% mesenchymal by the KS score. I do not call a cell hybrid because it is nearest Arc 1.

---

## 8. Are hybrid cells generalists?

This section is the test of the hypothesis in Section 1.2. It uses one rule on all four datasets, applied to samples already placed in the **Pal** triangle. The supervisor asked for two conditional probabilities: given that a cell is interior, what is its EMT state, and given its EMT state, is it interior?

### 8.1 The interior is not a hybrid compartment

Share of each gene score among interior samples, P(gene score | interior):

| Dataset | Interior \(n\) | Epithelial | Mesenchymal | Hybrid | Intermediate |
|---|---:|---:|---:|---:|---:|
| Pal | 22,402 | 35.6% | 31.1% | 13.4% | 19.9% |
| Wu | 15,473 | 20.9% | 36.7% | 14.6% | 27.9% |
| GSE173634 | 21,742 | 10.8% | 54.1% | 0.3% | 34.8% |
| DepMap | 45 | 28.9% | 46.7% | 0% | 24.4% |

A generalist, defined as interior, is not a hybrid cell. In the cell-line datasets the interior is where the mesenchymal samples are.

### 8.2 Hybrid cells are not reliably interior

P(interior | gene score):

| Dataset | Epithelial | Mesenchymal | Hybrid |
|---|---|---|---|
| Pal | 30.7% interior | 27.0% | 34.4% |
| Wu | 44.7% | 79.9% | 93.7% |
| GSE173634 | 20.0% | 100% | 100% (0 of 65 at a corner) |
| DepMap | 61.9% | 100% (0 of 21) | no hybrid line |

On Pal, hybrid cells are interior about as often as epithelial and mesenchymal cells. The hybrid label does not select the middle. On Wu, it does select the interior, and mesenchymal cells are also mostly interior. Epithelial cells are the specialist class in Wu (55% at a corner) and especially in GSE173634 (80%).

### 8.3 Where the hybrid cells actually sit

| Dataset | Hybrid \(n\) | At a corner | At Arc 1 | At Arc 2 | At Arc 3 | Interior | Median weights (Arc 1, 2, 3) |
|---|---:|---:|---:|---:|---:|---:|---|
| Pal | 8,719 | 65.6% | 53.2% | 10.3% | 2.1% | 34.4% | 0.68, 0.02, 0.11 |
| Wu | 2,413 | 6.3% | 0.5% | 5.8% | 0% | 93.7% | 0.51, 0.11, 0.38 |
| GSE173634 | 65 | 0% | 0% | 0% | 0% | 100% | 0.36, 0.28, 0.37 |
| DepMap | 0 | — | — | — | — | — | no line is jointly high |

Pal hybrid specialists are an Arc 1 phenomenon, and 92.2% of hybrid Arc-1 specialists are TN-0135. That patient’s median weights are about 0.91, 0.00, 0.08. Calling Arc 1 a hybrid task of breast cancer currently means calling it a property of this library.

Wu hybrid cells fail the corner cut because Arc 3 is close behind Arc 1 (median purity 0.12). They are not at the centroid. Half of all Wu hybrid cells are patient CID4515, a TNBC sample that looks like TN-0135 by nearest vertex and by hybrid fraction, and that has a specialist rate of zero. The same clinical picture produces a vertex in one atlas and an edge in the other.

GSE173634’s 65 hybrid cells are all from CAL851 and sit at equal weights. DepMap has no line in the top quartile of both scores. Epithelial and mesenchymal scores are positively correlated in Pal (\(r = 0.29\)) and negatively correlated in Wu (\(-0.29\)), GSE173634 (\(-0.77\)), and DepMap (\(-0.69\)). A large hybrid set is easier to obtain in Pal for that reason alone. Basal DepMap lines leave Arc 2: 26 of 27 are interior, spread between Arc 1 and Arc 3. The lines that pass the corner cut are Luminal B and HER2 at Arc 2.

### 8.4 Reading

Hybrid-as-its-own-corner does not replicate. The statement that does replicate is about who can be a specialist: epithelial samples can lock to Arc 2; mesenchymal samples, outside Pal, remain interior; basal and TNBC samples leave the epithelial corner and usually do not settle on a vertex. A subset of hybrid cells is interior, as the supervisor suggested. That subset is not large enough, and not specific enough, to identify the interior with the hybrid state.

---

## 9. Discussion

### 9.1 What I think the semester shows

The SCLC paper’s logic — fit on cell lines, project tumors in — does not transfer to the breast transcriptome, and it does not transfer to GBM. It partly transfers to a breast EMT gene list, and then a larger tumor atlas breaks it. The repair is not a new algorithm. It is a change of host. Once tumor cells set the corners, cell lines and a second tumor atlas are mostly inside, and the missing piece in the cell-line fit is the mesenchymal extreme.

The hypothesis that motivated the project is answered in the negative, with a narrower positive. Hybrid E/M is not the generalist class on this triangle. The generalist region is a mixture of epithelial, mesenchymal, and intermediate cells, plus, in Wu, hybrids that are really an edge between two poles. The specialist that is consistent across datasets is the epithelial program.

### 9.2 What would change my mind

Removing TN-0135 and refitting is the result that can change the third-vertex sentence. If a hybrid-enriched corner returns and is shared across patients, the third task is real and TN-0135 was only the most extreme member. If the corner collapses and the epithelial-versus-mesenchymal axis remains, the honest result is a two-program geometry with a patient-driven extra vertex, and the thesis should say that. If the whole triangle stops containing the other datasets, the host-swap claim itself was sensitive to one library and must be narrowed.

I have not yet done that refit. The patient audit shows the dependence. It does not remove it.

### 9.3 Limits

1. **Permutation depth.** Pal is 0/50. Breast KS used 15 PCHA starts rather than Groves’ 150. GBM, at the deeper setting, was negative, which is some reassurance that the light setting is not required to “find” a simplex. It is still a limit on the breast \(p\)-value.
2. **The 0.35 cut.** Pal’s hybrid-at-a-corner conclusion and GSE173634’s hybrid-at-the-centroid conclusion are stable across nearby cuts. Wu’s is not, if purity is replaced by “largest weight at least 0.5,” because those cells sit near 0.51 on an edge. I will report the weight vectors as the primary description.
3. **Within-dataset hybrid ranks.** The hybrid fraction is not a comparable prevalence. A common reference scale is a small follow-up, not a new project.
4. **Projection recipes differ** between ComBat bulk panels and MAGIC single-cell panels. The two-direction host swap is what makes the asymmetry hard to dismiss as one recipe’s artifact. A fully symmetric projection, same scaling both ways, is still worth doing once.
5. **One signature, one disease.** I am not claiming that every EMT list or every carcinoma behaves this way. Anderson et al. (2026) on premalignant lung, which I rebuilt as a methods exercise (primary archetype agreement 44/45 where both calls exist), is a different question and stays out of this claim.
6. **No functional readout.** Specialist here is geometric. It is not yet survival, metastasis, or a plasticity assay. That is appropriate for the midsemester. It is not yet a mechanistic paper.

### 9.4 Related computational work that did not become the thesis

A pan-cancer Hausser Figure 1d pipeline was drafted and is data-limited. It is how I learned the language of a frozen host. It is not a result. An extension to methylation was considered and, on the supervisor’s advice, not started. I think that advice was right.

---

## 10. Plan

The calendar for the rest of the semester is the mid-semester presentation now, Viva II in November, and the thesis report by the end of November.

1. Refit Pal \(k = 3\) without TN-0135. Project TN-0135 back in as a guest, and reproject Wu, GSE173634, and DepMap. Rebuild the specialist table. This is the only new biological run I intend before the thesis draft.
2. Show purity as a continuous score against the EMT score, and file the threshold table as a supplement.
3. Optionally repeat the DepMap KS \(k = 3\) t-ratio at the deeper PCHA setting, and optionally recompute hybrid calls on one shared scale. Each is one table.
4. Write the thesis in the order of this report: SCLC calibration, transcriptome-wide negatives, KS host swap, hybrid-versus-generalist test, leave-one-patient result, limits.
5. A paper in 2027 can be chapters 5 and 6 of that thesis. It should not be started as a separate analysis. Venue and scope are a discussion with the supervisor after the leave-one-patient figure exists.

I am not planning to refit Wu, to revive \(k = 4\), or to extend GBM.

---

## 11. Conclusion

The project began as a test of hybrid E/M cells as generalists on a breast-cancer simplex. The simplex that exists, on the KS signature, should be defined on tumor single cells rather than on cell lines. Hybrid cells do not occupy its interior in a way that distinguishes them from other states. Where they occupy a corner, the corner is largely one patient, and an independent atlas puts the analogous cells on an edge. The remainder of the thesis is that refit, and a write-up that keeps those limits in the same sentence as the result.

---

## Attendance

Quantitative attendance for this reporting period, as required by the course handout: working days present ____ / ____. Leave ____ days. The work was done at IISc. Supervision was by email and in person, including the meetings recorded in the August correspondence.

---

## References

1. Hausser, J., Szekely, P., Bar, N., Zimmer, A., Melton, S., Quake, S. R. & Alon, U. Tumor diversity and the trade-off between universal cancer tasks. *Nature Communications* 10, 5423 (2019).
2. Groves, S. M. et al. Archetype tasks link intratumoral heterogeneity to plasticity and cancer hallmarks in small cell lung cancer. *Cell Systems* 13, 690–710 (2022).
3. Tan, T. Z., Miow, Q. H., Miki, Y., Noda, T., Mori, S., Huang, R. Y.-J. & Thiery, J. P. Epithelial-mesenchymal transition spectrum quantification and its efficacy in deciphering survival and drug responses of cancer patients. *EMBO Molecular Medicine* 6, 1279–1293 (2014).
4. Pal, B. et al. A single-cell RNA expression atlas of normal, preneoplastic and tumorigenic states in the human breast. *EMBO Journal* 40, e107333 (2021).
5. Wu, S. Z. et al. A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics* 53, 1334–1347 (2021).
6. Sahoo, S. et al. Analysis of breast single-cell datasets with MAGIC and EMT scoring, as used for GSE173634 and GSE176078. *iScience* (2024). GEO accessions GSE173634 and GSE176078.
7. Wang, Q. et al. Tumor evolution of glioma-intrinsic gene expression subtypes associates with immunological changes in the microenvironment. *Cancer Cell* 32, 42–56 (2017). Used only as a GBM gene-list control.
8. Anderson, related preprint/article on archetype analysis of lung adenocarcinoma premalignancy, *Molecular Cancer Research* (2026), Figure 3. Methods reproduction only; not part of the breast claim.
9. The basal-versus-luminal difference in EMT heterogeneity discussed with the supervisor: PubMed 38974967.

PAM50 assignments used the `genefu` implementation of the Parker et al. classifier. DepMap expression and model metadata are the Cancer Dependency Map release used in `Breast Cancer/data/`. TCGA matrices are the UCSC Xena HiSeqV2 downloads. METABRIC is the cBioPortal study.
