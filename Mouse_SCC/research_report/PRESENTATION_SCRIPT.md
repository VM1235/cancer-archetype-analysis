# Presentation script: Mouse SCC archetype analysis

Deck: "Mouse SCC archetype analysis" (claude.ai artifact). Main talk = title + 8 numbered slides + conclusion, about **12–14 minutes**. The appendix (A1–A10) is for questions only.

Speak in your own words; the lines below are a guide, not something to read out. Numbers in **bold** are the ones worth saying aloud.

---

## Title slide (≈1 min)

> This weekend I tested whether the archetype pipeline from my breast-cancer chapter works on a second, very different dataset: mouse skin squamous cell carcinoma from Pérez-González et al., Nature Communications 2026. It has about **24,000 tumour cells from 12 tumours**, with the paper's four EMT states: epithelial, early hybrid, hybrid and late EMT.
>
> The question was simple: do these cells fill a Pareto triangle, and are the EMT states the corners? The short answer is on the slide: there's a shared EMT path, but a robust three-corner triangle isn't supported yet. Let me show you how I got there, because the route matters.

**Transition:** "First, the pipeline exactly as in the breast chapter."

## 1 · First fit, all cells (≈1.5 min)

> PCA on the cells, with at most 800 per tumour so that no big tumour dominates, then PCHA for k = 3 to 7. Each k is tested against 200 shuffled datasets. The rule for choosing k was written down before fitting: the smallest k with p below 0.01, not beyond the ESV elbow.
>
> k = 3 wins, **p = 0.005**: the real data beat all 200 shuffles. Two corners make sense: late EMT, with collagen genes, and dividing epithelial cells, with basal keratins.
>
> But look at the right panel. The third corner is coloured almost entirely purple, the multiome samples. **79%** of its cells are from the nuclei assay. That corner is technology, not biology.

**Transition:** "So the obvious fix: fit on one technology only."

## 2 · scRNA cells only (≈1.5 min)

> Same pipeline, fitted only on the 8 scRNA tumours. The triangle is still significant, p = 0.005. The late and dividing-epithelial corners are unchanged.
>
> The odd corner is still there, though: this thin spike on the left. Coloured by tumour, it's almost all one tumour, **SC15: 94%** of its scRNA specialists. Leaving SC15 out doesn't remove it, because a second tumour, SC46, carries the same profile.

**Transition:** "Before chasing that, I checked what happens if we only use EMT genes."

## 3 · EMT genes only (≈1 min)

> These are the 176 KS EMT genes from Tan 2014 that exist in mouse, the same panel as in my breast work. The cells line up beautifully from epithelial to late EMT, but **no triangle beats the shuffled data: p = 0.83**.
>
> The lesson: neat label separation is not evidence for a simplex. On EMT genes alone, this looks like a continuum, a line, not a triangle.

**Transition:** "Back to the all-gene triangle: who counts as a specialist?"

## 4 · Groves-style specialists (≈1 min)

> Groves 2022 defines specialists as the 5% of cells closest to each corner. Using that rule on the scRNA triangle, the tips are each dominated by one or two tumours: **SC15** at the odd corner, **SC2b1** at the late corner (92%), SC32 and SC33 at the epithelial corner. Hybrid cells are almost never specialists: **1.2%**.
>
> This is the same problem as patient TN-0135 in the breast data: the most extreme cells at a corner tend to come from one sample.

## 5 · Specialists under three rules (≈1 min)

> The specialist fraction depends heavily on the rule. Hybrid specialists are 18% with my purity-gap rule, 1.2% with Groves' 5% rule, and 0% with Groves' weight-above-0.95 rule. So I don't trust any single percentage; I trust the ranking, which holds under all three: **hybrid cells are the least specialised state**.

**Transition:** "Now, what is that odd corner?"

## 6 · Odd axis removed (≈2 min)

> I looked at the odd corner's genes: long genes like Camk1d, Pcdh7 and Maml2. That score is high in every multiome sample and in SC15 and SC46, across all EMT states. It's a whole-sample shift, not a cell state, most likely technical: how the RNA was captured or how reads were counted. I can't prove the cause without the raw reads (appendix A6).
>
> So I regressed that score out of every gene and refitted. The odd corner disappears, and the triangle I predicted before seeing any data appears: **late EMT, proliferating hybrid, and epithelial / early hybrid**. The hybrid corner is not just dividing cells: dividing epithelial cells don't go there.
>
> But here's the catch: the triangle is now weak. **p = 0.04**, so 7 of 200 shuffles do as well. That passes Groves' 0.05 bar, but not my pre-registered 0.01. Most of the original strong p-value came from the technical axis. And this analysis was designed after seeing the data, so it's exploratory.

**Transition:** "So is this triangle real biology, or does it come from mixing different tumours?"

## 7 · Triangles inside single tumours (≈1.5 min)

> I fitted a triangle inside each of the three tumours that contain all the states: SC15, SC32 and SC46. Only **SC46** has a real triangle (p = 0.005), and it's the only one with a hybrid corner. SC15 and SC32 show no triangle (p ≈ 0.8).
>
> What every tumour does share is a late-EMT corner with the same genes: gene-signature correlation **0.66–0.71** with the pooled late corner (appendix A9).

## 8 · Every tumour in one tumour's triangle (≈1.5 min)

> Finally, I froze SC32's triangle and placed all the other tumours' cells in it, without refitting. In every one of the 12 tumours, including the multiome ones, you see the same order: late EMT at one end, hybrid next to it, then early hybrid, then epithelial. **67–98%** of each tumour's late-EMT cells land at SC32's late corner. Hybrid cells sit between the ends, never at a corner.

**Transition:** "Putting it together."

## Conclusion (≈1.5 min)

> What the data support: one EMT path shared across all 12 tumours and both technologies, with hybrid cells as intermediates rather than specialists. Also, that the first strong triangle was driven by a technical axis. That's a useful warning for anyone running archetype analysis on integrated single-cell and nuclei data.
>
> What isn't supported yet: a robust three-task triangle, and any statement about tasks or transitions, since these are snapshots.
>
> The proposition I'd like your view on: EMT here behaves like a continuum with tumour-specific side directions, rather than a Pareto trade-off. The cheapest decisive test is to compare a line (k = 2) against the triangle directly. After that: settle the technical axis with gene-length and raw-read checks, look more closely at SC46, and run the MAGIC and Harmony versions.

---

## Likely questions, with short answers

**"How much variance do the PCs you used explain?"**
The triangle lives in PC1–PC2. On all genes that's **14.2%** of total variance (both versions A and B). On EMT genes it's 33.7%; after removing the V1 axis, 10.9%; inside single tumours, 10–14%. That's low, but not unusual: Hausser's bulk METABRIC tumours had 19.2% in PC1–PC2 (11.1% + 8.1%), while Groves' MAGIC-imputed mouse TKO single-cell data needed 9 PCs to reach 85%. Imputation concentrates variance into fewer PCs, and we haven't imputed yet. Within the PC1–PC2 plane, the triangle explains about 99% of the variance (ESV 0.989–0.994).

**"Isn't the triangle just tumour structure?"**
Partly, yes; that's why I did slides 7 and 8. Leaving out any one tumour keeps 90–99% of cells at the same corner (A8), but corners can be held up by groups of similar tumours. Inside single tumours a triangle holds only in SC46.

**"Why p < 0.01 and not 0.05 like Groves?"**
It was fixed in advance to be conservative with thousands of cells. For the first three versions the choice made no difference; after removing the V1 axis it decides the answer (p = 0.04), which is why I call that result weak and exploratory.

**"Did you just tune k to get 3?"**
No. All k from 3 to 7 were fitted and tested, and the rule was fixed before fitting. The EMT-gene version ended with no k at all.

**"Is V1 definitely technical?"**
Not proven. It's a whole-sample shift present in every nuclei sample and two scRNA tumours, independent of EMT state and of depth. The leading explanation is intronic reads boosting long genes. The definitive test needs spliced/unspliced counts from the raw reads.

**"Has archetype analysis been applied to tumours before?"**
Yes, many times: Hausser 2019 on bulk TCGA and METABRIC tumours; Groves 2022 on SCLC cell lines, two patients' tumours (single-cell) and three mouse TKO tumours (single-cell); Blais & Jeukens 2024 on 481 DLBCL patients (bulk). What's specific here is the EMT-state question in this mouse SCC dataset and the within-tumour and cross-tumour checks. I haven't done a systematic literature search for EMT-specific archetype work, so I wouldn't claim it's the first.

**"Why not MAGIC, like Groves?"**
To see the raw shape and keep the shuffled null fair first. MAGIC smooths cells toward their neighbours, which can create apparent shape. It's the next planned run.
