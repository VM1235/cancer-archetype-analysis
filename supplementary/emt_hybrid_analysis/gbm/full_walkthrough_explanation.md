# Full walkthrough: from the papers to the breast result to the GBM result

This document connects every piece: the two papers your work is built on,
what "sir" asked for, how the previous Claude session interpreted and built
it for breast cancer, and how I adapted it for glioblastoma (GBM) — with
the reasoning, math, code, and actual numbers at each step, starting from
first principles.

---

## Part 1 — The two papers underneath everything

Your repo (`cancer-archetype-analysis`) actually rests on **three** papers,
not one. It's important to keep them separate because they answer
different questions and Panel D is where they get combined.

### 1a. Groves et al., *Cell Systems* 2022 — "what shape is a population of cells?"

This is the paper the whole repo (SCLC, breast, GBM "Panel A/B/C") is
reproducing the *method* of. The idea, from basic to advanced:

- **Basic idea.** If you measure the expression of thousands of genes in
  many cell lines, each cell line is a point in a very high-dimensional
  space (one axis per gene). Most of that variation is redundant, so you
  first compress it with **PCA (Principal Component Analysis)** — a
  standard technique that finds the handful of directions (principal
  components, PCs) along which the cells actually differ the most, and
  re-expresses each cell line as coordinates along those PCs instead of
  along thousands of genes.

- **The core hypothesis (Pareto Task Inference).** Cells often have to be
  good at several conflicting "tasks" at once (e.g. "migrate fast" vs.
  "proliferate fast" vs. "resist stress"). If a task trade-off is sharp
  enough, evolution/selection tends to push the population toward a small
  number of extreme, highly-specialized strategies, with everyone else
  being some **mixture** of those strategies. Geometrically, this predicts
  the cloud of points (cell lines, in PC space) should look like a
  **polytope** — a simplex (triangle in 2D, tetrahedron in 3D, etc.) —
  whose corners ("archetypes") are the pure task-specialists, with
  everybody else sitting inside as a weighted blend.

- **PCHA (Principal Convex Hull Analysis)**, also called **archetypal
  analysis**, is the algorithm (`py_pcha` in your code, wrapped in
  `src/archetypes.py` / `pcha_bounded.py`) that, given a chosen number of
  corners *k*, finds the best-fit *k*-vertex polytope around the data and
  returns two things per cell line:
  - its **distance** to each vertex (used for "which archetype is this
    sample closest to"), and
  - its **weight vector S** — a set of *k* non-negative numbers that sum
    to 1, one per archetype, saying how much of a "blend" of each
    specialist strategy that cell line is. A cell line that's 100% one
    archetype has weight [1, 0, 0, ...]; a cell line that's an even blend
    of all *k* has weight [1/k, 1/k, ..., 1/k].

- **Is the polytope real, or just noise?** Any point cloud can be
  "fit" with a polytope — the question is whether it's a *significantly
  better* fit than you'd get from random data with the same PCA structure
  but no real corner structure. Groves' test for this is the **t-ratio**:
  the actual polytope's volume divided by the volume of the smallest
  ellipsoid that contains all the data. You compute this ratio on the
  real data, then **shuffle** the data many times (permutation test,
  usually 500–1000 shuffles) and compute the same ratio on each shuffle
  to build a null distribution. If the real t-ratio sits far in the tail
  of that null distribution, p is small and you can trust the polytope
  shape as real, not an artifact of fitting a shape to noise.
  - This is exactly the number reported as `t_ratio_parti_500.csv` in
    every Panel A folder in your repo, and it turns out to matter a lot
    for interpreting Panel D (see Part 4 and Part 5).

- **Enrichment (Panel B).** Once you trust (or at least have) a polytope,
  you can ask: does an *independent* label — one that had nothing to do
  with PCA or PCHA, e.g. a published subtype call — cluster near one
  particular vertex more than chance would predict? This is tested with
  a **hypergeometric test** (`src/enrichment.hypergeometric_enrichment`):
  bin the cell lines by distance-to-vertex into deciles, and ask "is this
  subtype over-represented in bin 0 (closest to the vertex) more than
  you'd expect if subtype assignment were random?" Multiple-testing
  correction (Benjamini-Hochberg FDR) is applied because you're running
  many such tests (every subtype × every vertex × every bin).

### 1b. George, Jolly, Xu, Somarelli & Levine, *Cancer Res.* 2017 — "how E, how M, how hybrid is this sample?"

This is the second paper — the PDF you attached — and it answers a
completely different kind of question: not "what shape is the population,"
but "for a single sample, where does it sit on the epithelial ↔
mesenchymal spectrum, on a continuous scale?"

- **Biological background (basic).** Epithelial-to-mesenchymal transition
  (EMT) is a process cancer cells use to become more mobile/invasive —
  losing epithelial cell-adhesion machinery (E-cadherin, tight junctions)
  and gaining mesenchymal machinery (vimentin, motility). It was long
  treated as binary (a cell either "is" epithelial or "is" mesenchymal),
  but there's good evidence cells can sit in a **stable hybrid E/M**
  state, not just transiently passing through the middle.

- **The model (intermediate).** The authors trained an **ordinal
  multinomial logistic regression** on the NCI-60 cancer cell line panel
  (60 well-characterized lines with known E / hybrid-E/M / M status from
  prior CDH1/VIM protein measurements). "Ordinal" matters here — unlike a
  generic 3-class classifier, it encodes the belief that the three
  classes are *ordered* (E < E/M < M), so hybrid is modeled as literally
  sitting *between* the two extremes rather than being an unrelated third
  category.
  - They searched combinatorially over ~480 candidate genes/gene-ratios
    and found the best 2-predictor model uses **X1 = log2(CLDN7)** and
    **X2 = log2(VIM) − log2(CDH1)**, with fitted coefficients
    β = (α1=−7.87, α2=0.0413, β1=1.36, β2=−1.96) (their Eq. 2).
  - Eq. 2–4 turn (X1, X2) into three probabilities: P(E), P(hybrid E/M),
    P(M).
  - **Eq. 5 (the punchline)** collapses those three probabilities into a
    single continuous number, **μ ∈ [0, 2]**: 0 = pure epithelial, 2 =
    pure mesenchymal, 1 = maximally hybrid. Thresholds: μ<0.5 → "E",
    0.5≤μ≤1.5 → "E/M", μ>1.5 → "M".
  - They validated this against real experiments (SNAIL-induced EMT,
    flow cytometry, immunofluorescence) across multiple cancer types
    (lung, prostate, colorectal — Table 2, Figures 3–4 of the paper), and
    showed μ correlates with patient survival, though the direction of
    that correlation is tissue- and subtype-specific (Figure 5).

- **Advanced/important nuance — the sign ambiguity.** The paper's Eq. 2
  literally *subtracts* the linear predictor: `log(π/(1-π)) = α − (β1·X1
  + β2·X2)`. Taken at face value with the reported β's, that implies
  *more* CLDN7 (an epithelial gene) pushes the score toward
  *mesenchymal* — backwards. This is a known hazard of re-implementing an
  ordinal/cumulative-logit model from a paper without the original fitted
  object (different software packages flip this sign silently). Because
  of this, any reproduction of this model needs an independent sanity
  check against a trusted label that has nothing to do with CLDN7/CDH1/VIM,
  before trusting μ's direction.

---

## Part 2 — What "sir" actually asked for, and how it was decoded

You gave the previous Claude session ("Chat1"/"Chat2", the two transcripts
you uploaded) an instruction referring to **"KS-epi" and "KS-mes,"** with
no other context, pointing at a private repo. Nothing in the codebase
matched those terms verbatim — Claude had to reverse-engineer the intent
by actually reading your repo and the George et al. paper you'd also
supplied. Here's the reasoning chain (visible across Chat1 and Chat2):

1. Your `src/enrichment.py` already had a **hypergeometric test** for
   *categorical* labels vs. distance bins (Panel B's machinery).
   Categorical tests don't apply to a **continuous** score like μ.
2. The George et al. paper's μ is exactly such a continuous score, and
   the standard way to compare a continuous distribution between two
   groups is the **two-sample Kolmogorov-Smirnov test**
   (`scipy.stats.ks_2samp`) — which measures the largest gap between two
   groups' cumulative distributions and gives a p-value for "these two
   samples were drawn from the same underlying distribution."
3. So "KS-epi" / "KS-mes" was decoded as: *run a KS test comparing some
   continuous score between an epithelial-leaning group and a
   mesenchymal-leaning group* — and the natural pairing with "are hybrid
   cells generalists" (a running theme "sir" seems interested in) is to
   test the **archetype generalist score** (not μ itself) across the
   E / hybrid / M groups defined by μ.
4. This produced the concrete pipeline: score EMT with George et al.'s
   model → get a generalist score from the existing PCHA weights → run
   `KS-epi` (E vs hybrid) and `KS-mes` (M vs hybrid) on the generalist
   score, to directly test "do hybrid cells behave like archetype
   generalists?"

This is genuinely a reasonable, well-justified decoding of an ambiguous
instruction — not a guess dressed up as certainty; Chat1 states this
reasoning explicitly and flags it as "my best reading."

---

## Part 3 — Concepts and tools used, defined plainly

| Concept / tool | What it is | Where it's used here |
|---|---|---|
| **log2(TPM+1)** | A standard way to represent gene expression: TPM (transcripts per million) is a normalized read-count; +1 avoids log(0); log2 compresses the huge dynamic range. | Every expression matrix in this repo. |
| **PCA** | Compresses thousands of gene axes into a handful of PCs capturing most of the variance. | `src/pca.py`, first step of every Panel A. |
| **PCHA / archetypal analysis** | Fits the smallest-volume polytope (k corners) that contains (or nearly contains) the data in PC space; returns each sample's distance-to-vertex and mixture weights. | `src/archetypes.py`, `pcha_bounded.py`. |
| **t-ratio permutation test** | Tests whether the fitted polytope's shape is statistically real vs. an artifact, by comparing to shuffled-data nulls. | `t_ratio_parti_500.csv` in every `panel_a/`. |
| **Hypergeometric enrichment test** | Tests whether a categorical label (e.g. PAM50 subtype, Verhaak subtype) is over-represented near a specific archetype vertex, vs. chance. | `src/enrichment.hypergeometric_enrichment`, Panel B. |
| **Benjamini-Hochberg FDR** | Multiple-testing correction, controls the expected proportion of false positives among all "significant" results when you run many tests at once. | `_bh_qvalues` in `enrichment.py`. |
| **Ordinal multinomial logistic regression** | A classifier for ordered categories (E < E/M < M) that outputs calibrated probabilities per category. | George et al.'s trained model, reproduced in `emt_score.py`. |
| **Shannon entropy (normalized)** | A measure of how "spread out" a probability-like vector is; 0 = all weight on one item, 1 (after normalizing by log(k)) = spread perfectly evenly over k items. | `archetype_weight_entropy()` — this *is* the "generalist score." |
| **Two-sample Kolmogorov-Smirnov (KS) test** | Non-parametric test of whether two samples come from the same continuous distribution, based on the max gap between their empirical CDFs. Makes no assumption of normality — appropriate for small, possibly skewed samples like n=10–20 cell lines. | `ks_two_group()` / `ks_epi_mes()`; the actual "KS-epi"/"KS-mes" tests. |
| **Spearman correlation (rho)** | Rank-based correlation, robust to non-linear-but-monotonic relationships and outliers; used here as an unthresholded sanity check. | Both Panel D scripts, μ (or PN-MES score) vs. generalist score. |
| **Z-score** | `(x − mean) / std`, puts a gene's expression on a common, comparable scale across samples, needed before averaging different genes together into one program score. | Verhaak marker scoring (`02_assign_verhaak.py`), used by GBM Panel D. |
| **Tools**: Python 3, `pandas`/`numpy` (data handling), `scipy.stats` (hypergeometric, KS, Spearman, `expit`/sigmoid), `scipy.cluster` (clustering, SCLC only), `matplotlib` (figures), `py_pcha`/custom PCHA (archetypal fitting), `git`/GitHub (version control — this is literally how I read your repo), bash/Python scripting environment (where I actually executed everything below). | | |

---

## Part 4 — What the previous Claude session built for breast (Chat1 + Chat2 → your uploaded files)

**Inputs it worked from:** your live GitHub repo (cloned), the George et
al. PDF, and your existing Breast Cancer Panel A/B outputs (already-fit
PCHA on 63 DepMap invasive-breast cell lines, k=4, plus PAM50 subtype
labels).

**Four new/edited files** (all in your upload, all verified to actually
run in that session):

1. **`src/emt_score.py`** — implements George et al.'s Eq. 2–5 exactly,
   using their published β. Crucially, it computes μ under **both** sign
   conventions ("paper" literal, and "flipped") rather than guessing, and
   provides `sanity_check_against_labels()` to pick the convention that
   agrees with an independent PAM50-based prior (Basal/TNBC = more
   mesenchymal, LumB/Her2 = more epithelial — chosen because that's
   standard breast cancer biology, and PAM50 calls have nothing to do
   with CLDN7/CDH1/VIM, so this check isn't circular).
2. **`src/enrichment.py` additions** — `archetype_weight_entropy()` (the
   generalist score), `ks_two_group()` (generic 2-sample KS wrapper), and
   `ks_epi_mes()` (runs KS-epi and KS-mes together, returns a tidy table).
3. **`Breast Cancer/codes/05_emt_generalist_panelD.py`** — the
   orchestration script: loads the existing Panel A matrix + PCHA
   weights, scores μ, picks the sign convention, bins into E/E-M/M, gets
   the generalist score, runs the two KS tests, does the Spearman
   robustness check, cross-references Panel B's PAM50-archetype
   enrichment, and saves a two-panel figure.
4. **`Breast Cancer/docs/panel_d_emt_generalist.md`** — the write-up.

**What it found (real numbers, this run, n=63 breast lines):**

- **Sign check:** "flipped" convention agreed with biology (mean μ[Basal]
  = 0.955 vs mean μ[LumB+Her2] = 0.172) — the paper's literal sign did
  not. This is exactly the kind of ambiguity flagged in Part 1b above,
  caught and resolved rather than silently assumed.
- **Category counts (flipped, breast):** roughly 42 E, 10 hybrid E/M, 11 M.
- **Generalist score** = entropy of each line's k=4 PCHA weights
  (`S_k4_parti.npy`), already computed by your Panel A — nothing was
  refit.
- **KS-epi (E vs hybrid):** D=0.371, **p=0.166 — not significant**
  (likely underpowered at n=10 hybrid lines).
- **KS-mes (M vs hybrid):** D=0.609, **p=0.030 — significant.** Hybrid
  E/M lines have a higher median generalist score (0.599) than pure
  mesenchymal lines (0.387).
- **Spearman rho(μ, generalist score) = −0.01, p=0.95** — flat, which
  the doc correctly notes is the *expected* result if the true
  relationship is an inverted-U (both pure extremes → specialist, hybrid
  middle → generalist) rather than a straight line; Spearman only
  detects monotonic relationships, so a null Spearman doesn't contradict
  the categorical KS finding.
- **Interpretation, stated carefully in the doc:** partial support for
  "hybrid E/M breast cells behave like archetype generalists" — clearly
  stronger against the mesenchymal pole than the epithelial pole, and
  explicitly caveated: breast's own Panel A t-ratio was **not
  significant** (best p=0.082 at k=4), so the polytope itself isn't
  confirmed as a real simplex — the generalist/specialist framing is
  mechanically well-defined regardless, but its biological force is
  weaker without a confirmed simplex.

---

## Part 5 — What I built for glioblastoma (GBM), and why it's not a copy-paste

You asked me to "do similar analysis for GBM data." The word **similar**
is doing real work here — I did not literally rerun the breast script on
GBM's matrix, because doing so would have been scientifically wrong, and
I want to walk through exactly why and what I did instead, since that's
the crux of the whole exercise.

### 5a. First, I checked whether the breast recipe even applies

The George et al. model's entire discriminating power comes from **CDH1**
(E-cadherin) being present in epithelial cells and absent in mesenchymal
ones. Glioblastoma cells are **not epithelial** — they arise from
glial/neural lineage and never had E-cadherin-based junctions in the
first place. I didn't just assert this; I pulled the actual numbers from
your GBM expression matrix (54 DepMap GBM cell lines):

| gene | median log2(TPM+1), 54 GBM lines | interpretation |
|---|---|---|
| CDH1 | 0.37 | essentially off, in *every* line, PN or MES |
| VIM | 11.21 (scale ~0–13) | essentially maxed out, in *almost every* line |

I then actually ran George et al.'s model on this data anyway (both sign
conventions), to demonstrate the failure with real numbers rather than
just reasoning about it:

| sign | mean μ(Mesenchymal-labeled) | mean μ(Proneural-labeled) | "agrees with biology"? | usable? |
|---|---|---|---|---|
| paper | 0.038 | 0.095 | No | No — 51/54 lines land in "E" |
| flipped | 1.9999 | 1.9048 | Yes, barely | **No — 53/54 lines pinned at μ=2.000, essentially zero variance** |

Even the "correct-direction" convention is useless here: it isn't
measuring Proneural-vs-Mesenchymal biology, it's just measuring "this
cell has no E-cadherin," which is true of the whole dataset regardless of
subtype. Applying the breast formula to GBM would have produced numbers
that *look* like a real analysis but carry no real signal — I treated
this as a hard stop, not a caveat to footnote.

### 5b. Finding the right substitute axis — reusing the repo's own logic, not inventing new biology

Your repo already names the correct axis for GBM, in `Glioblastoma/README.md`:
*"Verhaak's four subtypes are not four exclusive states; the antagonistic
pair is Proneural–Mesenchymal"* (citing Jolly lab, *iScience* 2024). This
is structurally the same shape as E ↔ hybrid ↔ M — two opposed
transcriptional programs with a hybrid middle — just biologically the
correct pair *for this tissue*.

Even better: your `Glioblastoma/codes/02_assign_verhaak.py` had **already
computed** exactly the ingredients needed, independent of any archetype
fitting: for each of the 54 GBM cell lines, a **z-score-based marker
score** for three Wang et al. 2017 programs (Classical: EGFR, NES,
NOTCH3, SMO, GAS1, GLI2, NFKBIA; Mesenchymal: CD44, CHI3L1, RELB, TRADD,
TNFRSF1A, NFKB1, STAT3; Proneural: OLIG2, PDGFRA, DLL3, NKX2-2, SOX10,
ASCL1, NKX2-1), sitting in `results/panel_b/verhaak_labels_panelA.csv`.
I reused this directly rather than recomputing it or inventing a new
scoring scheme.

**PN-MES score** (my continuous stand-in for μ) = `score_Mesenchymal −
score_Proneural`, one number per cell line.

**Category.** George et al. had a published threshold (μ 0.5/1.5) from
an NCI-60-fitted model — there's no equivalent published cutoff for this
z-score difference, so I used **data-driven tertiles** instead (bottom
third = PN, middle third = Hybrid, top third = MES) and said so
explicitly, rather than implying it was a validated clinical threshold.
This also happened to give much better-balanced groups than breast's:
**18 PN / 18 Hybrid / 18 MES**, vs. breast's 42/10/11.

I sanity-checked this tertile split against your existing argmax
`verhaak_subtype` label (unchanged, computed independently): 17/18
tertile-PN lines were argmax-Proneural, 15/18 tertile-MES lines were
argmax-Mesenchymal — consistent, as expected.

### 5c. Same generalist-score machinery, reused without modification

`archetype_weight_entropy()` — the exact same function from
`src/enrichment.py` used for breast — computed on GBM's own PCHA weights
(`S_k7_parti.npy`, k=7 is GBM Panel A's elbow-selected k, same file
Panel B already uses).

`ks_epi_mes()` — the exact same function, called with relabeled arguments
(`epi_label="PN"`, `mes_label="MES"`, `hybrid_label="Hybrid"`). No code
changes needed; this function was written generically enough in the
breast session to be reused as-is.

### 5d. The actual GBM result (real numbers, this run, n=54 lines)

| test | D | p | median(group) | median(Hybrid) | breast's equivalent |
|---|---|---|---|---|---|
| KS-PN (PN vs Hybrid) | 0.167 | **0.972** | 0.555 (PN) | 0.516 | KS-epi: p=0.166 |
| KS-MES (MES vs Hybrid) | 0.222 | **0.781** | 0.514 (MES) | 0.516 | KS-mes: p=**0.030** (significant) |

Spearman rho(PN-MES score, generalist score) = **−0.183, p=0.184** across
all 54 lines — flat/weak, and even trending the opposite direction from
"hybrids are more spread-out."

**Neither KS test is significant, and the three groups' generalist scores
overlap almost completely** (see the figure — all three boxplots sit
around 0.51–0.55 median, with heavily overlapping ranges). This is a
**genuine negative result**, not an artifact of small sample size —
GBM's groups (n=18 each) are actually larger than breast's hybrid group
(n=10), and the KS statistics themselves (D=0.17–0.22) are much smaller
than breast's KS-mes (D=0.61), so this isn't "underpowered to detect
something real," it's "the effect isn't there in this data."

### 5e. The caveat that's *stronger* for GBM than it was for breast

Breast's Panel A polytope had a borderline, non-significant t-ratio (best
p=0.082 at k=4) — already a caveat in the original write-up. GBM's is
weaker across the board: **every k tested (3 through 7) has p > 0.05**,
best is p=0.068 at k=7, and — unlike breast, where Panel B found a
significant Basal-subtype peak near one archetype — **GBM's Panel B
enrichment found zero significant subtype-archetype peaks at k=7.** There
is currently no archetype in this fit that can be confidently called
"the Proneural vertex" or "the Mesenchymal vertex." I stated this
plainly in the GBM doc rather than letting the negative KS result stand
without that context — a negative result on top of an unconfirmed
polytope is weaker evidence than a negative result on top of a confirmed
one, and the doc says so.

---

## Part 6 — Putting the two results side by side: what this actually tells you

| | Breast | GBM |
|---|---|---|
| Phenotype axis | George et al. 2017 CLDN7/CDH1/VIM EMT score (μ) — biologically valid, epithelial tissue | Custom PN-MES marker z-score (reusing Verhaak/Wang 2017 markers) — George et al. model is *invalid* here (no E-cadherin biology) |
| Hybrid group size | n=10 | n=18 |
| Polytope significance (best t-ratio p) | 0.082 (k=4) — borderline | 0.068 (k=7) — borderline, slightly worse |
| Panel B subtype-archetype peaks | Basal peaks significantly near one vertex | none significant |
| KS-epi/KS-PN (weaker pole) | p=0.166, not significant | p=0.972, not significant |
| KS-mes/KS-MES (stronger pole) | **p=0.030, significant** | p=0.781, not significant |
| Conclusion | Partial support: hybrid cells look more generalist than mesenchymal specialists, not clearly vs. epithelial specialists | No support: hybrid PN/MES cells don't look different from either pole |

**The honest synthesis:** "hybrid-phenotype cells are archetype
generalists" is **not a universal law that replicates across cancer
types** in this repo — it showed up (partially) in breast and did not
show up in GBM, on the same statistical machinery. That's a genuinely
useful, reportable finding in itself (and arguably more scientifically
interesting than if GBM had just mechanically reproduced breast's
result) — but it also means neither result should be over-claimed as
strong evidence of a real biological effect, given both polytopes sit at
the edge of statistical significance to begin with.

## Part 7 — Where everything now lives

- `Breast Cancer/codes/05_emt_generalist_panelD.py`,
  `src/emt_score.py`, `src/enrichment.py` (additions), 
  `Breast Cancer/docs/panel_d_emt_generalist.md` — from your previous
  Claude session, ready to commit to your repo as-is.
- `Glioblastoma/codes/03_pnmes_generalist_panelD.py`,
  `Glioblastoma/docs/panel_d_pnmes_generalist.md`,
  `results/panel_d_pnmes/*.csv`, `figures/Figure_1D_gbm_pnmes_generalist.png`
  — built and run by me this session, following the same conventions
  (reads only existing Panel A/B outputs, writes only to its own
  `panel_d_pnmes/` folder, never touches Panel A/B/C).
