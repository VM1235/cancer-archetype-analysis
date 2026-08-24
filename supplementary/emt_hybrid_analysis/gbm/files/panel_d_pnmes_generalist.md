# Panel D (GBM) — Proneural/Mesenchymal score vs archetype generalist score

Added 2026-08-19, mirrors `Breast Cancer/docs/panel_d_emt_generalist.md`.
Answers the GBM version of the same question breast Panel D asked: *are
"hybrid" glioblastoma cell lines — with mixed Proneural/Mesenchymal
transcriptional identity — archetype "generalists" rather than
specialists?* Does not touch Panel A/B/C.

## Why this is not just "rerun the breast script on GBM data"

The breast Panel D used George, Jolly, Xu, Somarelli & Levine's 2017
CLDN7 / CDH1 / VIM epithelial-EMT model to get a continuous phenotype
score. That model needs **CDH1 (E-cadherin)** as a live, varying
predictor — it is the epithelial half of the epithelial↔mesenchymal
axis the model was built to describe.

Glioblastoma is not an epithelial tumor. It arises from glial/neural
lineage; GBM cells never had E-cadherin junctions to lose. We checked
this directly rather than assuming it: running George et al.'s model
literally on this repo's 54 GBM cell lines gives

| gene | median log2(TPM+1) across 54 GBM lines |
|---|---|
| CDH1 | 0.37 (near the assay floor) |
| VIM  | 11.21 (near the ceiling; scale is ~0–13) |

CDH1 is barely expressed in *any* GBM line — Proneural or
Mesenchymal — and VIM is high in nearly all of them, because glial
cells are intrinsically vimentin-positive regardless of subtype, not
because they've undergone an epithelial-to-mesenchymal transition.
Feeding these into the model's Eq. 2–5 pins the score at the
mesenchymal extreme almost everywhere:

| sign convention | mean mu(Mesenchymal) | mean mu(Proneural) | agrees with biology? | discriminates at all? |
|---|---|---|---|---|
| paper | 0.038 | 0.095 | No | No — 51/54 lines land in "E" |
| flipped | 1.9999 | 1.9048 | Yes (barely) | **No — 53/54 lines are pinned at mu=2.000** |

Even the sign that "agrees with biology" on average does so with
essentially zero variance: it isn't measuring Proneural-vs-Mesenchymal
biology, it's measuring "this cell doesn't express E-cadherin," which
is true of the whole dataset. Using this score for GBM would be a
category error, not a minor caveat — so this Panel D does not use it.

## What we used instead

This repo's own `Glioblastoma/README.md` already names the correct
axis for GBM: *"Verhaak's four subtypes are not four exclusive
states; the antagonistic pair is Proneural–Mesenchymal"* (citing the
Jolly lab, iScience 2024) — which is exactly the same kind of
"two-poles-plus-a-hybrid-middle" structure that E/hybrid-E-M/M is for
epithelial cancers. So instead of forcing an epithelial model onto
glial cells, Panel D reuses the marker-gene scores this repo *already
computes independently of PCHA* in `02_assign_verhaak.py`
(Wang et al. 2017 marker sets: 7 Classical, 7 Mesenchymal, 7
Proneural genes; each cell line's score per program = mean z-score of
its marker genes across the 54 lines).

**PN-MES score** = `score_Mesenchymal − score_Proneural` (continuous,
one number per cell line, positive = more Mesenchymal-leaning,
negative = more Proneural-leaning). Classical-program genes are not
part of this axis (Classical is a separate, non-antagonistic program
in Wang's scheme) — Classical-labeled lines are simply wherever their
Mesenchymal/Proneural z-scores happen to put them, same as everyone
else.

**Category (PN / Hybrid / MES).** George et al. had a published
threshold (mu 0.5 / 1.5) from their NCI-60-trained ordinal model. This
marker z-score difference has no equivalent published cutoff, so we
used **data-driven tertiles** instead (bottom third = PN, middle third
= Hybrid, top third = MES) — stated explicitly so it isn't mistaken
for a validated clinical threshold. This choice also gives much better
balanced group sizes than breast's n=10 hybrid line: **18 PN / 18
Hybrid / 18 MES** out of 54 lines.

Cross-tab against this repo's existing argmax `verhaak_subtype` label
(Panel B, unchanged) shows the tertile split is consistent with it —
17/18 tertile-PN lines are argmax-Proneural, 15/18 tertile-MES lines
are argmax-Mesenchymal — with most of the disagreement in the Hybrid
tertile, as expected for cells scored near the argmax decision
boundary. Classical-labeled lines spread across all three tertiles (9
Hybrid, 3 MES, 1 PN), consistent with Classical being a separate
program rather than a third pole of this specific axis.

Output: `results/panel_d_pnmes/pnmes_scores.csv` (per-line Classical/
Mesenchymal/Proneural marker z-scores, pnmes_score, category).

## Step 2 — generalist score (`src/enrichment.archetype_weight_entropy`, unchanged)

Same function as breast, same idea: normalized Shannon entropy of each
cell line's PCHA archetype-mixture weights (`S_k7_parti.npy` — GBM
Panel A's k, `results/panel_a/suggested_k.txt`). 0 = concentrated on
one archetype (specialist), 1 = spread evenly over all k=7 archetypes
(generalist).

**Important caveat, stronger than breast's.** Breast's k=4 polytope
had a borderline, non-significant t-ratio (best p=0.082). GBM's is
weaker still: `results/panel_a/t_ratio_parti_500.csv` shows **p > 0.05
for every k tested (k=3–7; best p=0.068 at k=7)**, and Panel B's own
subtype-vs-archetype hypergeometric enrichment (`enrichment_verhaak.csv`)
found **zero** significant bin-0 peaks at this k — i.e. this repo
cannot currently point to "the Proneural archetype" or "the
Mesenchymal archetype" the way breast Panel B could point to a
Basal-leaning vertex. Read everything below as more exploratory than
the (already-caveated) breast result, not as a confirmed finding.

Output: `results/panel_d_pnmes/archetype_weights_and_generalist_score.csv`.

## Step 3 — KS-PN / KS-MES (`src/enrichment.ks_epi_mes`, unchanged — just relabeled)

Same function breast used, called with `epi_label="PN"`,
`mes_label="MES"`, `hybrid_label="Hybrid"` (it's generic; no code
changes needed). Two pre-specified two-sample KS tests:

- **KS-PN**: generalist score, PN lines (n=18) vs Hybrid lines (n=18)
- **KS-MES**: generalist score, MES lines (n=18) vs Hybrid lines (n=18)

This run:

| test | D | p | median(group) | median(Hybrid) |
|---|---|---|---|---|
| KS-PN  | 0.167 | 0.972 | 0.555 (PN)  | 0.516 |
| KS-MES | 0.222 | 0.781 | 0.514 (MES) | 0.516 |

Neither is significant, and the medians barely differ (all three
groups sit around generalist score ≈0.51–0.55). Robustness check with
no thresholding at all: Spearman rho(pnmes_score, generalist_score) =
**-0.183, p = 0.184** across all 54 lines — also not significant, and
if anything trends the opposite direction from an "extremes are
specialists, hybrids are generalists" hypothesis (weak negative
correlation, not the inverted-U breast's raw scatter suggested; see
`Figure_1D_gbm_pnmes_generalist.png`).

## Reading these numbers (GBM, k=7, n=54, this run)

- **No evidence that hybrid PN/MES glioblastoma lines are archetype
  generalists.** Both KS tests are far from significance (p=0.97,
  p=0.78) and the three groups' generalist-score distributions
  visibly overlap almost completely in the figure.
- This is a genuine negative result on this dataset/method, not
  "inconclusive because of the caveats" — the caveats (non-significant
  polytope, no confirmed PN/MES archetype identity) mean we should be
  cautious about *over-interpreting a positive* finding, but they
  don't manufacture this particular negative one; the effect sizes
  (KS D=0.17–0.22, rho=-0.18) are just genuinely small next to
  breast's KS-mes D=0.61.
- Contrast with breast: breast found a significant KS-mes (p=0.030,
  hybrid E/M > mesenchymal specialists) with a directionally
  consistent (if underpowered) KS-epi. GBM shows neither. Put
  together, "hybrid-phenotype cells are archetype generalists" does
  **not** replicate as a general cross-cancer pattern in this repo —
  it may be specific to how breast's epithelial-EMT axis relates to
  DepMap breast's PCHA polytope, not a universal rule about hybrid
  transcriptional states.
- Two candidate (untested here) explanations, offered as hypotheses
  only: (1) GBM's polytope itself is weaker/less real than breast's
  (worse t-ratio, no significant subtype peaks), so there may be
  nothing coherent for a PN/MES axis to align with; (2) the
  Proneural-Mesenchymal axis may genuinely relate to archetype
  structure differently than epithelial-EMT does biologically — e.g.
  if it aligns with one specific archetype direction rather than with
  "distance from every vertex," a directional test (per-archetype,
  see `distance_bins()` + `ks_two_group()`) might behave differently
  from this pooled entropy test. Not explored in this run.

## Caveats carried over from the rest of this repo

- GBM's Panel A t-ratio permutation test was **not significant at any
  k tested** (k=3–7, best p=0.068 at k=7) — weaker than breast's
  already-non-significant best (p=0.082 at k=4). The polytope itself
  is not confirmed to be a real simplex; "specialist/generalist" is a
  correspondingly weaker claim here than in breast, and much weaker
  than in the Groves SCLC reproduction (which was significant).
- Panel B's Verhaak-vs-archetype enrichment found no significant peaks
  at k=7, so there is no independently-validated "this archetype =
  Proneural" / "this archetype = Mesenchymal" identity to sanity-check
  against, unlike breast Panel B's Basal/archetype-4 peak.
- Tertile-based PN/Hybrid/MES thresholds are a data-driven convenience,
  not a validated clinical cutoff (contrast with George et al.'s
  paper-defined 0.5/1.5 mu thresholds for breast/NCI-60).
- Marker-based Mesenchymal/Proneural scores use only 7 genes per
  program (Wang et al. 2017's compact set), same limitation the breast
  analysis flagged for its 2-gene EMT score — cross-check against a
  fuller signature (e.g. the original Verhaak 2010 or Wang 2017
  full gene sets, or ssGSEA) before treating this as a precise
  per-line phenotype.
- This repo's own README explicitly separates the Groves-style
  multi-archetype k=7 fit used here from the hypothesis-driven k=2
  PN-MES-axis fit (`results/panel_a_k2_12/`), which was built
  specifically to test the PN-MES antagonism directly in PC space.
  This Panel D deliberately uses the *same* k as Panel B (k=7, for
  apples-to-apples comparison with breast's k-matches-Panel-B
  convention), not the k=2 axis fit; re-running the generalist score
  from `S_k2_12pc.npy` instead is a natural follow-up, flagged but not
  done here to keep this analysis parallel to breast's.

## How to rerun

From the repository root, after Panels A and B already exist for GBM:

```bash
.venv/bin/python -u "Glioblastoma/codes/03_pnmes_generalist_panelD.py"
```

Reads only existing `results/panel_a/` and `results/panel_b/` outputs;
writes only to `results/panel_d_pnmes/` and the two
`Figure_1D_gbm_pnmes_generalist...` image paths. Safe to rerun any
time; it does not call PCHA again.
