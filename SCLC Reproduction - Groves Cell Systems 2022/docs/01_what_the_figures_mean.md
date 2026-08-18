# What Figure 1A–C is saying

These three panels are one argument, not three unrelated plots. Groves et al. (Cell Systems 2022) claim that human SCLC cell lines do not fill expression space as a blob. They sit inside a low-dimensional **polytope** whose corners are **archetypes**: extreme expression programs. Each real cell line is a mixture of those extremes. The same geometry then shows up in bulk tumors.

The method is Pareto Task Inference (ParTI / archetypal analysis), following Hausser et al. (Nature Communications 2019). The biological idea is older still: if a cell cannot maximize every function at once, the data cloud is shaped by **trade-offs**. The vertices are the pure tasks; the interior is compromise.

---

## The object being fitted

Start with a matrix: **120 SCLC cell lines × ~16,000 genes**, log-expression, already batch-corrected for Minna vs CCLE. That matrix is too wide to fit a polytope in gene space, so it is reduced to its leading principal components (here, 12 PCs, which capture about half the variance). Archetypal analysis (PCHA) then finds \(k\) points — the archetypes — such that:

1. every cell line is a **convex combination** of the archetypes (mixture weights \(\ge 0\) that sum to 1);
2. every archetype is a **convex combination** of real cell lines (the corners are not imaginary points outside the data).

Geometrically, the archetypes are the vertices of a simplex that tries to wrap the cloud. If \(k=5\), that simplex is a 4-dimensional tetrahedron-like object (a \(k\)-vertex simplex lives in \(k-1\) dimensions).

The names later attached to those vertices (“SCLC-A”, “proliferation”, and so on) are **not** produced by PCHA. PCHA only returns coordinates and mixture weights. Names come from matching vertices to known subtypes (Panel B) and from gene-set analysis (out of scope here).

---

## Panel A — how many corners?

Two questions, answered together.

### 1. Explained sample variance (ESV)

For each candidate \(k = 2, \ldots, 15\), fit PCHA and ask what fraction of the 12-PC cloud is captured by mixtures of \(k\) archetypes. Plot ESV vs \(k\), and the increment \(\mathrm{ESV}(k)-\mathrm{ESV}(k-1)\).

If the cloud is really a simplex with a small number of vertices, extra archetypes after the true \(k\) buy little new variance: an **elbow**. An elbow is a hint, not a proof. You can always raise ESV by adding vertices.

### 2. The t-ratio permutation test

A simplex can wrap noise. The t-ratio asks whether the fitted simplex fills an unusually large fraction of the data’s own convex hull:

\[
t = \frac{\mathrm{volume}(\text{archetype simplex})}{\mathrm{volume}(\text{convex hull of the samples})}
\]

in the natural dimension of a \(k\)-vertex simplex, i.e. the first **\(k-1\)** PCs.

- \(t\) near 1: the archetypes sit near the true hull — the cloud is simplex-shaped.
- \(t\) near 0: the fitted simplex is a small object inside a rounder cloud.

PCHA always returns some \(t < 1\). Significance comes from a **null**: independently shuffle each PC across samples (destroy correlations between axes, keep each PC’s marginal distribution), refit PCHA, recompute \(t\). Repeat many times.

\[
p = \Pr(t_{\text{shuffle}} \ge t_{\text{real}})
\]

A small \(p\) means: you do not get a simplex this “full” after breaking the geometry.

**Why not just take the smallest \(p\)?** Adding vertices almost always improves the fit, so \(p\) tends to fall as \(k\) grows. The paper’s rule is: take the **smallest \(k\) that is significant** (and that is also supported by the elbow and by biology). In Groves et al. that is **\(k=5\)**: \(k=4\) is not significant (\(p \approx 0.059\)), \(k=5\) is (\(p \approx 0.034\)).

Our official reproduction uses the same calling convention as their MATLAB ParTI code (see the companion implementation note). We recover their t-ratios to three decimals. With that protocol, \(k=4\) is not significant and \(k=5\) is. That is the result we report.

The 2-D polytope drawing is a **cartoon**. Fitting happens in 4-D for \(k=5\) (or in 12-D for the ESV curve). The scatter is PC1 vs PC2 so a human can see five vertices and the cloud inside them. Edges between all pairs are drawn because a simplex is the convex hull of its vertices; crossing edges in 2-D are expected when a 4-D object is flattened.

---

## Panel B — are the five corners the five known subtypes?

SCLC is already typed by master transcription factors into SCLC-A (ASCL1), SCLC-A2 (a related ASCL1-high group), SCLC-N (NEUROD1), SCLC-P (POU2F3), and SCLC-Y (YAP1). Those labels were assigned **independently of PCHA**, by clustering and TF expression, and published with the Groves data.

Panel B is a validation, not a second archetype fit.

For each of the five archetypes:

1. measure Euclidean distance of every cell line to that archetype in PC space;
2. sort the 120 lines into **10 equal bins** (12 lines each; bin 0 = closest);
3. test whether a subtype is over-represented in that bin (hypergeometric test, Benjamini–Hochberg FDR \(q < 0.1\));
4. call a match only if the enrichment **peaks in bin 0** (the lines nearest the vertex), not somewhere in the middle of the simplex.

If the geometry is biological, each subtype should spike at bin 0 for **exactly one** archetype, and the five spikes should land on five different archetypes. That is a 1-to-1 correspondence: the mathematical corners are the known transcriptional extremes.

A line in the interior of the polytope is then interpreted as a **mixture** of subtypes, which is the paper’s bridge to plasticity. Panel B does not prove plasticity; it only shows that the vertices align with the classical types.

---

## Panel C — do real tumors live in the same polytope?

Cell lines are the clean system (clonal, no stroma). The scientific worry is that the five-corner geometry is a culture artifact.

Panel C brings in **81 human SCLC tumors** (George / “Thomas” bulk RNA-seq), ComBat-corrected **together** with the 120 lines so that “cell line vs tumor” is the batch, not a biological axis you want to keep.

Two checks:

### Left: scatter in a **new** PCA of all 201 samples

The five archetypes were **not refit** on tumors. They are the cell-line vertices, mapped into this combined PCA and held fixed. If tumors fall inside that polytope, the tumor transcriptomes are mixtures of the same extremes. If they shot out past a vertex, tumors would need extra programs that cell lines do not have.

In the paper (and here) tumors sit inside, often toward the SCLC-A region rather than stretched all the way to every vertex. That is expected if bulk tumors are mixtures and/or biased toward one subtype.

### Right: variance explained

- Fit PCA on **tumors only**. Cumulative variance vs number of components is the **ceiling**: the best a linear basis can do on tumors.
- Take PCs from the **combined** 201-sample PCA and ask how much tumor variance those axes capture.
- Shuffle the combined matrix and repeat, as a **null** (unrelated axes should explain little tumor variance).

The claim is not “5 PCs explain 80% of all tumor variance.” It is: **5 combined PCs capture ~80% of what a tumor-only PCA captures with 5 components.** The tumor cloud is largely the same subspace the cell lines already span. Grey shuffle curves stay low; the match is not accidental.

---

## How the three panels lock together

| Panel | Question | Answer we are reproducing |
|---|---|---|
| A | How many extremes? | **Five**, smallest significant simplex |
| B | What are those extremes? | The five **known SCLC subtypes** |
| C | Is this only cell lines? | **No** — tumors occupy the same polytope / subspace |

If A failed, B and C would be decorations around the wrong \(k\). If B failed, the polytope would be a geometric curiosity. If C failed, the story would not leave the dish. All three are required for the claim that SCLC is organized by a small number of trade-off tasks that also describe human tumors.

---

## What these figures are not

- They are not a clustering of cell lines. Clustering partitions; archetypes are **vertices**. A sample can be 40% A and 60% N.
- They do not assign the later task labels (proliferation, secretion, etc.). That is gene-set enrichment on genes that peak near each vertex (Figure 1E in the paper, not reproduced here).
- They do not, by themselves, prove that cells move between archetypes. Plasticity is a later claim (single-cell and perturbation work in the same paper).
- The 2-D pictures are projections. Significance was tested in \(k-1\) dimensions, not in the cartoon.
