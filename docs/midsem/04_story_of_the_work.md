# The story this work is telling

This is the narrative for the report and the talk. It is chronological, and it is also an argument. Each step exists because the previous step was not enough. Side analyses (GBM details, Anderson, Hausser drafts, \(k = 4\)) are one sentence when they are controls, and absent when they are not.

**Title line.** Cell-line archetypes of the breast EMT signature do not contain tumors. The tumor-defined triangle does contain cell lines. Hybrid epithelial/mesenchymal cells are not, in general, the generalists in the middle.

---

## 1. The idea

Tumors are thought to face trade-offs among tasks. Archetypal analysis represents that as a simplex: specialists at the corners, generalists in the interior (Hausser et al. 2019). In small-cell lung cancer this is not only a picture. Cell lines form a significant five-vertex simplex, subtypes sit on vertices, and tumors occupy the same shape (Groves et al. 2022).

The question brought to this project was the EMT version of that picture. Extreme epithelial cells and extreme mesenchymal cells would be specialists. Hybrid E/M cells, which are plastic, would be the generalists in the interior.

## 2. The method had to be shown to work before it was aimed at that question

The Groves SCLC analysis was rebuilt in Python. The t-ratios match the paper, including 0.107 at \(k = 5\). With 1,000 shuffles the smallest significant \(k\) is 5, which is the published choice. From then on, a new dataset could be a biological result rather than a coding result.

## 3. The obvious next dataset did not contain the structure

The same pipeline on 63 breast cancer cell lines, full transcriptome, does not yield a significant simplex. The elbow says four archetypes; the best \(p\) is 0.082. TCGA tumors do not fall inside that shape. Glioblastoma, run as a second test with a cleaner statistical setting, also has no significant \(k\). The method is capable of saying no.

So the failure is not “breast has no geometry.” It is “the geometry is not the whole transcriptome of the cell lines.”

## 4. The gene list that matches the question does have a triangle

Restricted to the 206 KS epithelial and mesenchymal genes (Tan et al. 2014), the same 63 lines support a triangle: t-ratio 0.587, \(p = 0.022\). Basal and luminal lines separate toward different vertices. HER2 does not get a vertex. This is the first positive breast result, and it is still a **cell-line** result. A short PAM50 gene list does not do this. The signal is the EMT signature.

## 5. Single cells partly agree, until a large tumor atlas does not

Cell-line single cells (GSE173634) and one tumor atlas (Wu, GSE176078) fall partly inside that triangle, about half the cells after MAGIC. A third atlas, Pal et al. GSE161529, does not: 13.3% of 84,602 epithelial cells are inside. If the cell-line triangle were the tumor geometry, Pal should have behaved like Wu. It did not.

## 6. Reversing the host is the result

Fit the triangle on Pal. Do not refit it for anyone else. Wu, the cell-line single cells, and the 63 bulk lines now fall mostly inside (78%, 96%, 91%). Pal remains mostly outside the cell-line triangle. A four-vertex fit does not wrap the other datasets, so the working shape is the triangle even though an elbow preferred four.

The scientific point: for this signature, cell lines are a truncated tumor geometry, not the geometry tumors fail to reach. They cover the epithelial side and under-represent the mesenchymal extreme.

## 7. The corners are not the clinical subtypes

On the Pal triangle the epithelial vertex is ER+ and HER2+ together (luminal epithelial genes). The mesenchymal vertex is a VIM program and is not one IHC class. The third vertex is enriched for TNBC, and TNBC itself splits: some patients sit on that vertex, some on the mesenchymal vertex. HER2 never earns its own corner. Clinical labels and EMT tasks are different partitions.

## 8. The idea that started the project does not survive the test

Score each cell as specialist or interior from its archetype weights, and separately as epithelial, mesenchymal, or hybrid from the KS genes.

A cell in the interior is usually not hybrid. Hybrid is about 13–15% of interior cells in the two tumor atlases.

A hybrid cell is a specialist only in Pal, and there the specialist corner is one patient, TN-0135. In Wu, hybrid cells are interior, but they sit on the edge between that third pole and the mesenchymal pole, and half of them are one other patient who never reaches the corner. In the cell-line data, joint-high hybrid cells are essentially absent, and mesenchymal samples are the ones in the interior. Epithelial samples are the ones that can lock to a corner.

So hybrid E/M is not the generalist class. It is sometimes a patient-dominated vertex, sometimes an edge, and sometimes missing. The reproducible specialist is the epithelial program.

## 9. What the story refuses to say

- It does not say PCHA finds a simplex in every cancer. GBM and the breast transcriptome are the counterexamples.
- It does not say the Pal \(p\)-value is zero. Zero of fifty shuffles means \(p < 0.02\).
- It does not say Arc 1 is a hybrid task of TNBC. It says Arc 1 is where one large TNBC patient sits, and the next experiment is to refit without that patient.
- It does not say cell lines are useless. It says they are the wrong host for the corners of this signature.

## 10. What the story points to next

One experiment: refit Pal without TN-0135 and ask whether the third corner, the containment of the other datasets, and the hybrid-specialist count remain. That result, either way, is the last scientific slide of the thesis and the first figure a 2027 paper needs. The rest of the time is writing this story at thesis length.
