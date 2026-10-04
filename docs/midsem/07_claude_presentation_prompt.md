# Prompt to paste into Claude

Give Claude the repository as project context, then paste everything below the line. The prompt is written so Claude builds the deck from the figures on disk and does not invent a second analysis.

---

You are making a mid-semester thesis presentation for Vedashree Mahajan (BITS Pilani, K K Birla Goa Campus; supervisor Dr. Mohit Kumar Jolly, IISc Bioengineering). The audience is the supervisor, who already knows Hausser 2019 and Groves 2022 and has seen the email updates. The talk is 12 minutes plus questions. Build a PowerPoint (.pptx) with python-pptx if it is available, otherwise a Google-slides-ready PDF. Save it as `docs/midsem/midsem_presentation.pptx`.

The science lives in the repo. Read these files before you draw anything, and use only numbers that appear in them:

- `docs/midsem/04_story_of_the_work.md` (the narrative; follow it)
- `docs/midsem/05_midsem_report_brief.md` (the sentences you may say out loud)
- `Breast Cancer/docs/RESULT_compare_four_datasets.md`
- `Breast Cancer/results/host_swap_comparison/inside_fractions_k3.csv`
- `Breast Cancer/results/specialist_interior/summary.json`

Do not read `PROJECT_CONTEXT.md` for results. It is older than the Pal analysis and will mislead you.

## What the talk must feel like

One claim per slide. A person who has never fit PCHA should still follow the story. A person who has should not feel talked down to. Plain words first, the technical word in the same sentence. No slide should have more than about 40 words of body text plus a figure. No three-column method dumps. No screenshots of code. No logos of tools.

Visual style: white background, black text, one accent color (a deep blue). Sans-serif. Slide title is the claim, not a topic label. Good title: “Tumor cells contain the cell-line triangle. The reverse is not true.” Bad title: “Results — Host Swap Analysis.” Large figure, small caption with the file’s actual percentages. Do not redraw the data plots. Place the PNG. Do not add drop shadows, gradients, or icon rows.

14 slides. Not 20. If a figure is not in the allowed list, do not use it. In particular do not use any file with `k4`, `genelist` without `ks`, `pre-MAGIC` (the `Figure_4_ks_panelA_gse173634_sc.png` and `Figure_4_ks_panelA_gse176078_sc.png` that lack `_magic`), GBM, or Anderson.

## Allowed figures (repo-relative)

- `SCLC Reproduction - Groves Cell Systems 2022/figures/Figure_1A.png`
- `Breast Cancer/figures/Figure_1A_breast_ks_genelist.png`
- `Breast Cancer/figures/Figure_4_ks_panelA_gse161529_sc_magic.png`
- `Breast Cancer/figures/Figure_host_swap_triangles_k3.png`
- `Breast Cancer/figures/Figure_host_swap_inside_fractions_k3.png`
- `Breast Cancer/figures/Figure_1B_gse161529_named_scatter.png`
- `Breast Cancer/figures/Figure_emt_hybrid_pal_triangle.png`
- `Breast Cancer/figures/Figure_specialist_vs_interior_hybrid.png`
- `Breast Cancer/figures/Figure_program_given_specialist.png`

You may use `Figure_4_depmap_on_gse161529_k3.png` only if the host-swap triangles figure is too crowded to see the cell lines. Prefer the host-swap figure.

## The story, in the order the slides must follow

1. We asked whether hybrid epithelial/mesenchymal cells are generalists in the middle of a simplex, and extreme epithelial or mesenchymal cells are specialists at the corners.
2. The method is real: on small-cell lung cancer it recovers Groves’ five-archetype result (t-ratio 0.107 at k = 5).
3. It does not find that structure in the full breast transcriptome or in glioblastoma. Say this in one slide, no gallery of negative figures.
4. On the 206 KS EMT genes (Tan 2014), 63 breast cell lines do form a triangle (t-ratio 0.587, p = 0.022, 500 shuffles).
5. A large tumor atlas does not sit inside that triangle: Pal GSE161529, 13.3% of 84,602 cells inside.
6. Fit the triangle on Pal instead. Wu 78.4%, GSE173634 96.1%, DepMap lines 90.5% fall inside. Pal inside the cell-line triangle stays 13.3%. Who hosts the simplex changes the biology. Cell lines miss the mesenchymal extreme.
7. The corners are not ER / HER2 / TNBC. Arc 2 is epithelial and holds ER+ and HER2+. Arc 3 is mesenchymal. Arc 1 is TNBC-enriched. HER2 has no corner. TNBC splits by patient.
8. Hybrid cells are not the interior. Interior cells are only ~13–15% hybrid. Pal hybrids are mostly at Arc 1, and that is patient TN-0135. Wu hybrids are on the Arc 1–Arc 3 edge (median weights 0.51, 0.11, 0.38), not the center. Cell lines have no hybrid sample under this rule, and their mesenchymal samples are interior.
9. Next experiment, not yet done: refit Pal without TN-0135.
10. Stop. Three sentences of conclusion. Invite the question about the 0.35 cutoff and the 50 shuffles.

## Slide-by-slide

**Slide 1. Title.**  
Title: “Who defines the corners?”  
Subtitle: “Epithelial–mesenchymal archetypes in breast cancer, and a test of hybrid cells as generalists.”  
Name, ID f20220216, supervisor Dr. Mohit Kumar Jolly, mid-semester thesis presentation, October 2026. No figure.

**Slide 2. The question.**  
Title: “Hybrid cells were predicted to sit in the middle.”  
A simple triangle you draw (not a data figure): three corners labeled epithelial specialist, mesenchymal specialist, and a third corner left unnamed; the center labeled generalist. One sentence: Hausser’s trade-off picture, applied to EMT. One sentence: Groves showed the picture is real in small-cell lung cancer, with cell lines as the host and tumors as the guests. Speaker note: do not define PCHA on this slide.

**Slide 3. The only methods slide.**  
Title: “A simplex is fit once, then everyone else is drawn on it.”  
Four short lines, no more:  
(1) PCA, keep the first k−1 components.  
(2) PCHA finds k corners; each sample is a mixture with weights that sum to 1.  
(3) t-ratio = volume of the simplex / volume of the data cloud; shuffle the axes to get a p-value. Smallest k with p < 0.05.  
(4) New data are projected in. Archetypes are not refit.  
Put `Figure_1A.png` from the SCLC folder on the right, small, caption: “SCLC reproduction. t-ratio 0.107 at k = 5, matching Groves. This is the calibration, not the result.”

**Slide 4. The full transcriptome says no.**  
Title: “Breast and glioblastoma cell lines do not form this shape.”  
No figures. Three lines:  
Breast, 63 lines, ~16,500 genes: best is k = 4, p = 0.082. Not significant. Tumors do not fall inside.  
Glioblastoma, 54 lines, paper-grade settings: no k with p < 0.05 (closest k = 7, p = 0.068, t-ratio 0.027).  
A 50-gene PAM50 list was worse (p = 0.76). The failure is the feature space, not the code.  
Speaker note: if asked, the SCLC t-ratios match the paper to three digits, so this negative is biological.

**Slide 5. The KS triangle.**  
Title: “On 206 EMT genes, the cell lines do form a triangle.”  
Figure: `Breast Cancer/figures/Figure_1A_breast_ks_genelist.png`.  
Caption: “Tan 2014 KS epithelial + mesenchymal genes. k = 3, t-ratio 0.587, p = 0.022 (500 shuffles). k = 4 is also significant; we keep the smallest k.”  
One bullet: Basal and luminal lines separate; HER2 does not get a corner. This is still a cell-line geometry.

**Slide 6. Pal does not fit in it.**  
Title: “A large tumor atlas mostly falls outside the cell-line triangle.”  
Figure: `Breast Cancer/figures/Figure_4_ks_panelA_gse161529_sc_magic.png`.  
Caption: “Pal et al. GSE161529, MAGIC, 11,217 / 84,602 cells inside = 13.3%. Wu was 45.5% and GSE173634 was 53.6% after the same MAGIC step.”  
One line: earlier pre-MAGIC percentages were higher and are retired. Do not show them.

**Slide 7. The result.**  
Title: “Fit the triangle on the tumors, and the cell lines fall inside.”  
Figure: `Breast Cancer/figures/Figure_host_swap_triangles_k3.png`.  
Under it, a four-row table you typeset, taken from `inside_fractions_k3.csv`, not re-estimated:

| Dataset | Inside cell-line triangle | Inside Pal triangle |
|---|---:|---:|
| Pal, 84,602 tumor cells | 13.3% | 75.8% |
| Wu, 24,162 tumor cells | 45.5% | 78.4% |
| GSE173634, 35,271 cell-line cells | 53.6% | 96.1% |
| DepMap, 63 bulk lines | 52.4% | 90.5% |

One sentence under the table: “Cell lines are a truncated version of the tumor geometry. They miss the mesenchymal extreme.”  
Speaker note: k = 4 does not contain the guests (Pal self-containment ~38%). That is why the triangle is the working shape even though an elbow preferred 4. Pal’s permutation result is 0 of 50 nulls, so say p < 0.02, never p = 0. t-ratio 0.66 at k = 3, 0.33 at k = 4.

**Slide 8. Optional bar chart, include only if slide 7’s figure already shows the clouds and you still want a clean quantitative slide. If the triangles figure plus the table is enough, delete this slide and renumber.**  
Title: “The asymmetry is the result, not the diagonal.”  
Figure: `Breast Cancer/figures/Figure_host_swap_inside_fractions_k3.png`.  
One sentence: “DepMap being 52% inside its own triangle is not the success. Pal outside DepMap, and DepMap inside Pal, is the success.”

**Slide 9. Naming the corners.**  
Title: “The three corners are not the three clinical subtypes.”  
Figure: `Breast Cancer/figures/Figure_1B_gse161529_named_scatter.png`.  
Three bullets only:  
Arc 2: luminal-epithelial genes; ER+ and HER2+ live here.  
Arc 3: mesenchymal (VIM); not one immunohistochemistry class.  
Arc 1: enriched for TNBC; TNBC also uses Arc 3. HER2 has no corner of its own.  
Speaker note: nearest-to-Arc-1 is not the same as “this cell is hybrid.” Patient TN-B1-0554 is the example if asked.

**Slide 10. The hypothesis test.**  
Title: “Hybrid cells are not the cells in the middle.”  
Figure: `Breast Cancer/figures/Figure_specialist_vs_interior_hybrid.png`.  
If that figure is a methods sketch rather than the comparison, use `Figure_program_given_specialist.png` instead, and keep the hybrid triangle (`Figure_emt_hybrid_pal_triangle.png`) as a small inset or as the only figure if the other two are unreadable. Look at the PNGs and pick the one a person can read from the back of a room. Do not stack three unreadable plots.

Say these numbers, they are the talk:

- Rule: hybrid = both KS scores at or above the dataset’s 75th percentile. Specialist = largest archetype weight beats the second by ≥ 0.35.  
- Pal: 65.6% of 8,719 hybrid cells are at a corner, 53% of hybrid cells at Arc 1.  
- Wu: 93.7% of 2,413 hybrid cells are interior, on the edge (median weights 0.51, 0.11, 0.38), not the center.  
- GSE173634: 65 hybrid cells, all interior, all one line. DepMap: 0 hybrid lines.  
- Among interior cells, hybrid is only 13.4% (Pal) and 14.6% (Wu).

**Slide 11. The caveat that protects the claim.**  
Title: “The Pal hybrid corner is one patient.”  
No new analysis. Text only, large:  
TN-0135: 14,322 cells, all nearest Arc 1, all specialists.  
92% of Pal’s hybrid Arc-1 specialists are this patient.  
Wu’s similar patient (CID4515, TNBC, nearest Arc 1, mostly hybrid) sits on the edge and has a specialist rate of 0%.  
Bottom line in bold: “I am not calling Arc 1 a third task of breast cancer until the fit is repeated without TN-0135.”

**Slide 12. What this means.**  
Title: “Cell lines are the wrong host. Hybrid is not the generalist.”  
Three sentences, no figure:  
For this signature, a simplex fit on tumor cells contains cell lines and a second tumor atlas. The fit on cell lines does not contain the tumors.  
The reproducible specialist is the epithelial program. Mesenchymal and basal samples are usually interior.  
Hybrid E/M is a corner in one patient, an edge in the next atlas, and nearly absent in cell lines.

**Slide 13. Next, and the limits.**  
Title: “One experiment, then the thesis is this talk.”  
Left, “Next”: refit Pal k = 3 without TN-0135; reproject Wu, GSE173634, and DepMap; show the purity score as a continuous number, not only the 0.35 cut.  
Right, “Already true, so I will keep saying them”: breast PCHA used fewer random starts than Groves; Pal’s test is 50 shuffles; the hybrid label is a rank inside each dataset; MAGIC and the projection recipe affect percent-inside.  
Do not list methylome, GBM extensions, or a Wu refit. Those are out of scope.

**Slide 14. Close.**  
Title repeated: “Who defines the corners?”  
One line: “Tumor cells do. And the middle of their triangle is not where the hybrid cells are.”  
Thank you. Questions. No references slide unless asked; if you add one, it is Hausser 2019, Groves 2022, Tan 2014, Pal 2021, Wu 2021, and nothing else.

## Speaker notes

Write 40–80 words of speaker notes on every slide, in the student’s voice, first person. Notes should say what not to add verbally. On slide 7, the note must include “do not say p = 0.” On slide 10, the note must say that replacing the 0.35 purity cut with “largest weight ≥ 0.5” would falsely call Wu hybrids specialists, because their second weight is close behind.

## Hard bans

- Do not say the project found five breast archetypes, or that GBM has seven.
- Do not say 84% or 69% inside for the single-cell cohorts. Those are pre-MAGIC.
- Do not say hybrid cells are generalists. Do not say hybrid cells are specialists, without the patient clause.
- Do not say TNBC is the mesenchymal corner, or that HER2 has its own corner.
- Do not invent confidence intervals, hazard ratios, or survival results. There are none.
- Do not show code, file trees, or the Anderson lung figure.
- Do not use the phrase “delve,” “landscape,” “robust framework,” or “exciting.”
- If a number is not in the files listed at the top, leave it out.

## Check before you finish

Open every PNG you embedded and confirm the slide title does not contradict the picture. Confirm the deck has at most 14 slides. Confirm `Figure_host_swap_triangles_k3.png` and either the specialist figure or the hybrid Pal triangle are both in the deck. Confirm the words “TN-0135” appear on a slide, not only in notes.
