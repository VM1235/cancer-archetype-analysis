# Mid-semester report (brief)

**Hybrid epithelial/mesenchymal cells on the archetype geometry of breast cancer**

Vedashree Mahajan  
ID f20220216  
BITS Pilani, K K Birla Goa Campus  
First Degree Thesis, Semester I, 2026–27  
Supervisor: Dr. Mohit Kumar Jolly, Department of Bioengineering, Indian Institute of Science  
Period: 10 August 2026 – 4 October 2026

---

## Abstract

This project asked whether hybrid epithelial/mesenchymal (E/M) cells behave as generalists in the interior of a task-trade-off simplex, with extreme epithelial and mesenchymal cells at the corners. I first reproduced the published archetype analysis of small-cell lung cancer. The same method does not find a significant simplex in the full transcriptome of breast or glioblastoma cell lines. On the 206-gene KS epithelial/mesenchymal signature, 63 breast cell lines do form a significant triangle. A large tumor single-cell atlas (Pal et al., 84,602 cells) mostly falls outside that triangle (13.3% inside). When the triangle is fit on the tumor cells instead, cell lines and a second tumor atlas fall inside it (78–96%). The corners are an epithelial program, a mesenchymal program, and a third vertex enriched for triple-negative cells. Hybrid cells are not the interior class. In the Pal atlas they concentrate at the third vertex, which is dominated by one patient. In an independent tumor atlas they lie on an edge rather than at a corner. The working conclusion is that, for this signature, cell lines are the wrong host for the simplex, and hybrid E/M is not in general the Pareto generalist.

## 1. Aim

Hausser et al. (2019) describe tumors as lying near low-dimensional polytopes set by trade-offs among tasks: specialists near vertices, generalists in the interior. Groves et al. (2022) showed that small-cell lung cancer cell lines form such a simplex and that tumors occupy it. The aim of this thesis is to ask the analogous question for epithelial–mesenchymal plasticity in breast cancer. The specific hypothesis, set at the start of the project, was that hybrid E/M cells would lie in the interior and that fully epithelial or mesenchymal cells would lie at vertices.

## 2. What was done

I implemented Pareto Task Inference (principal component analysis, principal convex hull analysis, a volume ratio tested by permutation, and projection of new samples into a frozen simplex) and checked it on the Groves SCLC data. The t-ratio at five archetypes is 0.107, matching the paper. With 1,000 permutations the smallest significant number of archetypes is 5.

I then ran the same analysis on 63 DepMap invasive breast carcinoma lines and on 54 glioblastoma lines, including projection of TCGA tumors. Neither transcriptome-wide matrix has a significant simplex. Restricting breast to the Tan et al. (2014) KS signature (206 genes) gives a significant triangle (\(p = 0.022\), 500 permutations). I projected three single-cell datasets into that cell-line triangle after MAGIC imputation, then reversed the fit: archetypes were estimated on GSE161529 (Pal et al. 2021) and the other datasets were projected in without refitting. Finally I scored every sample in all four datasets as specialist or interior from its archetype weights, and as epithelial, mesenchymal, or hybrid from the KS genes, which is the comparison requested after the 26 September update.

Glioblastoma and a separate rebuild of a lung-premalignancy archetype figure were controls and methods checks. They are not the result below.

## 3. Results

**The host of the simplex matters.** Fraction of samples inside each \(k = 3\) triangle:

| Dataset | \(n\) | Inside the cell-line triangle | Inside the Pal triangle |
|---|---:|---:|---:|
| Pal GSE161529, tumor cells | 84,602 | 13.3% | 75.8% |
| Wu GSE176078, tumor cells | 24,162 | 45.5% | 78.4% |
| GSE173634, cell-line cells | 35,271 | 53.6% | 96.1% |
| DepMap bulk lines | 63 | 52.4% | 90.5% |

The tumor-derived triangle contains the cell lines. The cell-line triangle does not contain the Pal tumors. A four-vertex fit does not contain the other datasets, so the triangle is the working geometry. On the Pal t-ratio test both \(k = 3\) (t-ratio 0.66) and \(k = 4\) (0.33) beat all 50 permutations. With 50 permutations that is \(p < 0.02\), not a more precise value.

**The corners are programs, not the three clinical bins.** Arc 2 is luminal-epithelial and holds ER+ and HER2+ cells. Arc 3 is mesenchymal. Arc 1 is enriched for TNBC and is not a single immunohistochemistry class. HER2 has no vertex of its own. TNBC splits across Arc 1 and Arc 3, largely by patient.

**Hybrid cells are not the generalists.** A sample is called hybrid when both KS scores are in the top quartile of its own dataset, and a specialist when its largest archetype weight beats the second by at least 0.35.

| Dataset | Hybrid samples | Of which at a corner | Of which interior |
|---|---:|---:|---:|
| Pal | 8,719 (10.3%) | 65.6% | 34.4% |
| Wu | 2,413 (10.0%) | 6.3% | 93.7% |
| GSE173634 | 65 (0.2%) | 0% | 100% |
| DepMap | 0 | — | — |

In Pal the hybrid specialists sit at Arc 1, and 92% of those cells are one patient, TN-0135. In Wu the hybrid cells have median weights about 0.51, 0.11, and 0.38: the edge between Arc 1 and the mesenchymal vertex, not the center. Half of them are one patient, and that patient does not pass the corner cut. Interior cells, in every dataset, are mostly epithelial, mesenchymal, or intermediate. Hybrid is only 13–15% of the interior in the tumor atlases. Epithelial samples are the ones that consistently reach a corner. Mesenchymal samples are interior in Wu, in GSE173634, and in DepMap.

## 4. Implication

The hypothesis that hybrid E/M cells are the interior generalists is not supported as a general rule on this signature. What is supported is a host effect: archetypes defined on cell lines miss the tumor geometry, in particular the mesenchymal extreme, whereas archetypes defined on a large tumor atlas contain both cell lines and a second tumor cohort. An apparent hybrid vertex should not be reported as a third task until the dominant patient is removed and the fit is repeated. That is the next experiment.

Two limits belong in the same paragraph as the result. Breast and Pal fits used fewer PCHA random starts than the SCLC reproduction, and the Pal permutation test used 50 shuffles. The hybrid call is a within-dataset rank, so the hybrid fraction is not directly comparable across cohorts. The specialist cut of 0.35 is a choice; the ordering (Pal hybrids at a corner, Wu hybrids on an edge, cell-line hybrids absent) is what the weight vectors show, and a sensitivity table is in the repository.

## 5. Plan for the rest of the thesis

Refit the Pal triangle without TN-0135 and reproject the other datasets. Report archetype-weight purity as a continuous score alongside the 0.35 cut. Then write the thesis around the host-swap and the hybrid-versus-generalist test, with SCLC as the methods check and the transcriptome-wide negatives as controls. A manuscript in 2027 is realistic only after the leave-one-patient refit, and only as this argument, not as a new method or a new data modality.

## References

1. Hausser, J. et al. Tumor diversity and the trade-off between universal cancer tasks. *Nature Communications* 10, 5423 (2019).
2. Groves, S. M. et al. Archetype tasks link intratumoral heterogeneity to plasticity and cancer hallmarks in small cell lung cancer. *Cell Systems* 13, 690–710 (2022).
3. Tan, T. Z. et al. Epithelial-mesenchymal transition spectrum quantification and its efficacy in deciphering survival and drug responses of cancer patients. *EMBO Molecular Medicine* 6, 1279–1293 (2014).
4. Pal, B. et al. A single-cell RNA expression atlas of normal, preneoplastic and tumorigenic states in the human breast. *EMBO Journal* 40, e107333 (2021).
5. Wu, S. Z. et al. A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics* 53, 1334–1347 (2021).
6. Gambardella, G. et al. A single-cell analysis of breast cancer cell lines to study tumour heterogeneity and drug response. Relevant GEO series GSE173634, as used in Sahoo et al., *iScience* (2024).

## Attendance

The handout asks for a quantitative attendance statement. Working days present in this period: ____ of ____. Days of leave: ____. I was at IISc for the thesis work and met the supervisor in person and by email through the period above.

## Note on the repository

Figures and tables cited here are in `Breast Cancer/figures/` and `Breast Cancer/results/`. The comparison write-up is `Breast Cancer/docs/RESULT_compare_four_datasets.md`.
