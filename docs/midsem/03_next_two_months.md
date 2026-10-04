# What to do in the next two months

**Window:** 5 October 2026 through early December 2026.  
**Hard dates from the First Degree Thesis handout (Semester I 2026–27):** mid-semester report and presentation by 9 October; Viva II by 15 November; thesis abstract and final report by 28 November; final viva in the 10–15 December window.  
**Publication:** not this semester. A submittable paper in 2027 is realistic if the two experiments below are done before the thesis is frozen, so the thesis and the paper are the same argument.

The thesis is already a complete narrative: method works on SCLC; it does not find a simplex in the breast or GBM transcriptome; on KS genes, the tumor atlas is the right host; hybrid E/M is not the generalist. What is missing is the check that stops a reviewer, and then the writing. New cancer types, methylation, and a Wu-as-host refit are how this stops being finishable.

---

## The claim worth defending

For the Tan 2014 KS genes in breast cancer:

1. A simplex fit on a large tumor single-cell atlas contains cell lines and a second tumor atlas. The simplex fit on cell lines does not contain the tumor atlas.
2. The reproducible corner is epithelial. Mesenchymal and basal programs lie toward the interior or along an edge.
3. Hybrid E/M cells are not the Pareto generalists. Where they look like a third specialist, the vertex is dominated by one patient.

Everything in the next two months either hardens this claim or writes it down.

---

## October (after the midsem): the one experiment the claim depends on

### 1. Refit the Pal triangle without TN-0135

TN-0135 is 14,322 cells, all specialists at Arc 1. Of Pal’s hybrid cells that are specialists at Arc 1, 92% are this patient. Until the fit is repeated without those cells, “hybrid is a third archetype” is a statement about one library.

Do this and only this as the new analysis in October.

- Drop TN-0135 from the MAGIC KS matrix.
- Refit \(k = 3\) with the same code path as `17.9_run_panelA_ks_genelist_gse161529.py`.
- Project the held-out TN-0135 cells into the new triangle, and project Wu, GSE173634, and DepMap in again.
- Recompute the specialist/interior table with script 26’s rules.

Three outcomes, all publishable if they are reported plainly:

| Outcome | What it means for the thesis |
|---|---|
| Arc 1 disappears or moves, and the epithelial-versus-mesenchymal axis remains | The third corner was a patient. The paper’s result is the epithelial specialist pole plus a mesenchymal interior. This is the stronger, more careful result. |
| A hybrid-enriched corner returns, now shared by several TNBC patients | The third task survives. TN-0135 was an extreme member, not the cause. |
| The triangle itself falls apart (t-ratio no longer beats the null, other datasets fall outside) | The Pal host was not robust. The thesis then stops at the host-swap asymmetry and does not claim three tasks. |

Also drop CID4515 from the Wu projection tables and recompute Wu’s hybrid median weights. That is a table, not a refit. Wu was never the host.

### 2. Report the continuous score next to the 0.35 cut

Script 26 already stores median weights, purity, and evenness. For the thesis, plot purity (largest weight minus second) against the mesenchymal-minus-epithelial score, for each dataset, colored by hybrid versus not. Put `threshold_sensitivity.csv` in the supplement. Do not adopt “max weight ≥ 0.5” as a new rule. Say why the Wu hybrids sit at purity about 0.12.

### 3. Write the midsem documents as the skeleton of the thesis, then stop exploring

The brief and detailed reports in this folder are that skeleton. After the leave-one-patient refit, add one results subsection. Do not open a new dataset in October.

---

## November: make the methods match the claim, then write

### 4. One cleaner statistical sentence for the two fits the talk shows

- DepMap KS \(k = 3\): if time allows, rerun the t-ratio at the Groves initialization depth (150 observed starts, 50 null starts) for \(k = 3\) only. If \(p\) stays below 0.05, say so. If it does not, the thesis uses the 500-shuffle, 15-start result and says the search was lighter than Groves. Do not scan \(k = 2\) to \(7\) again.
- Pal: either 200–500 shuffles at \(k = 3\) on the leave-one-patient matrix, or an explicit “0 of 50” with no fake precision. A full 1,000-shuffle Groves run on 84,000 cells is not required for the thesis.

### 5. One shared-scale hybrid check

Recompute the double-high call using Pal’s gene means and standard deviations as the reference for all four datasets, or a simple mean of the KS sets without within-dataset quartiles. The question is only whether “DepMap has no hybrid line” and “Pal hybrids are 10%” survive a common ruler. One table. If they do not survive, the quartile definition stays in the methods as a within-dataset rank, and the paper leans on the weight vectors instead.

### 6. Do not do these

- Do not refit PCHA on Wu. Wu is a guest. Fitting it would add another host with fewer patients and reopen the story.
- Do not promote \(k = 4\). The tetrahedron does not contain the guests.
- Do not expand GBM, Anderson, Hausser Figure 1d, methylome, or MIDAA. Anderson can be one paragraph in a methods chapter: an independent Figure 3 rebuild, 44/45 primary-call agreement.
- Do not attach survival unless the leave-one-patient result is already in the draft. A survival panel on a vertex that is one patient will not help.

### 7. Thesis shape (this is the November writing)

A first-degree thesis, not a paper with appendices.

1. Introduction. Task trade-offs, specialists and generalists, EMT and hybrid E/M, why breast and why not another SCLC figure.
2. Methods. PCA, PCHA, t-ratio, barycentric inside/outside, the two labels (gene score versus geometry), datasets, what was not refit.
3. Reproduction. SCLC Figure 1A–C, t-ratios versus the paper. Short.
4. Where the transcriptome-wide fit does not generalize. Breast full matrix and GBM, one section, as controls.
5. KS genes and the host swap. The table of inside-fractions. Vertex names. HER2 has no corner. TNBC is not one task.
6. Hybrid versus generalist. The conditional probabilities. The patient audit. The leave-one-patient refit.
7. Discussion. What this says about using cell lines to define EMT archetypes. What it says about hybrid cells as generalists. Limits: MAGIC, projection recipe, 50 shuffles, one signature, one disease.
8. Conclusion and a short “submitted toward a manuscript” paragraph.

Target a full draft by 20 November so the 28 November submission is a revision, not a first writing. Viva II on 15 November can be the same talk as the midsem, plus the leave-one-patient slide if it exists.

---

## December and the 2027 paper

The final viva defends the thesis. The paper is a cut of chapters 5 and 6, with chapter 4 as a supplement and chapter 3 as one methods sentence (“the engine reproduces Groves et al., t-ratio 0.107 at \(k = 5\)”).

A plausible venue, given the supervisor’s editorship and the style of the question, is a specialist systems-biology or cancer-biology journal rather than a new-method journal. The contribution is a biological correction: cell-line EMT archetypes are not the tumor geometry, and hybrid E/M is not the interior generalist on that geometry. It is not a new algorithm.

Before submission the paper needs:

- the leave-one-patient refit, as a main figure;
- the four-dataset conditional probabilities, as a main figure;
- the host-swap inside-fraction table, as a main figure;
- methods that state shuffles, starts, MAGIC, the 0.35 cut, and the within-dataset quartiles;
- a short related-work passage: Hausser 2019, Groves 2022, Tan 2014, Pal 2021, Wu 2021, Sahoo 2024, and the basal-versus-luminal EMT heterogeneity result the supervisor cited (the 2024 paper linked as PMID 38974967).

It does not need a new modality. If the leave-one-patient result is clean by December, a preprint in the first half of 2027 is a reasonable aim. If the refit overturns the triangle, the preprint is the negative-plus-asymmetry paper, which is smaller and still worth writing, and the thesis is honest either way.

---

## Order of work, as a checklist

1. Midsem talk and report from the documents in this folder. No new runs before the presentation.
2. Refit Pal \(k = 3\) without TN-0135. Project the four datasets. Rebuild the specialist table.
3. Purity-versus-EMT scatter and the sensitivity table in the supplement.
4. Optional: deeper t-ratio on the DepMap KS triangle only; common-scale hybrid call, one table.
5. Thesis draft by 20 November, submit 28 November.
6. After the viva, cut chapters 5–6 into a paper outline and discuss venue with the supervisor. Do not start that cut before the leave-one-patient figure exists.
