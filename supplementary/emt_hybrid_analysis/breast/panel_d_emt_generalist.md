# Panel D — EMT score vs archetype generalist score

Added 2026-08-19. Answers: *are hybrid E/M breast cell lines archetype
"generalists" rather than specialists?* Does not touch Panel A/B/C.

## Inputs (nothing recomputed)

- `results/panel_a/S_k{k}_parti.npy` — PCHA mixture weights, one row per
  cell line, one column per archetype, rows sum to 1 (already fit).
- `results/panel_a/pc_scores.npy`, `archetypes_k{k}_parti.npy` — same
  polytope Panel B already uses.
- `data/processed/input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv`
  — same 16,500 × 63 matrix as Panel A; we only read three rows
  (CLDN7, CDH1, VIM) out of it.
- `results/panel_b/pam50_labels_panelA.csv` — used once, only to pick a
  sign convention (see below), not part of the archetype test itself.

## Step 1 — EMT score (`src/emt_score.py`)

Implements George, Jolly, Xu, Somarelli & Levine, *Cancer Res*
77(22):6415–6428 (2017) (`papers/` does not yet have a copy — see
citation above; PDF supplied by the user, not re-hosted here). Their
model uses exactly two predictors — log2(CLDN7) and log2(VIM/CDH1) —
in a 2-threshold ordinal logistic regression fit on the NCI-60 panel,
with published coefficients (their Results section):

```
alpha1 = -7.87, alpha2 = 0.0413, beta1 = 1.36, beta2 = -1.96
```

fed into their Equations 2–5 to give each sample P(E), P(hybrid E/M),
P(M), and a single score mu in [0, 2] (0 = epithelial, 1 = maximal
hybrid, 2 = mesenchymal), thresholded at 0.5 / 1.5 into three
categories.

**Sign-convention caveat.** The paper's Eq. 2 subtracts the linear
predictor from each threshold. Taken completely literally, that makes
higher CLDN7 (a tight-junction / epithelial gene) push mu *toward*
mesenchymal — backwards. This is a known hazard reproducing a
cumulative-logit model from a paper without the original fitted
object; MATLAB's `mnrfit` and textbook GLM notation disagree on this
sign in general. `emt_score.py` computes mu under both the literal
("paper") and negated ("flipped") convention and does not pick one
blindly.

`05_emt_generalist_panelD.py` picks between them using an **independent**
check that has nothing to do with archetypes: PAM50 Basal/TNBC lines
are a standard, well-established proxy for high EMT activity in
breast cancer; LumB/Her2 are more epithelial. On this run:

| sign | mean mu(Basal) | mean mu(LumB+Her2) | agrees with biology? |
|---|---|---|---|
| paper | 1.304 | 1.926 | No |
| **flipped** | **0.955** | **0.172** | **Yes** |

So the script uses the flipped convention. **If you rerun this with
different cell lines or a different cancer type, re-check this table
before trusting mu** — the script prints it every run and raises an
error rather than proceeding if neither convention agrees.

Output: `results/panel_d_emt/emt_scores.csv` (per-line CLDN7, CDH1,
VIM, X1, X2, mu, category, PAM50 label — both sign conventions are
also kept in the intermediate table for audit).

## Step 2 — generalist score (`src/enrichment.archetype_weight_entropy`)

PCHA's `S` matrix already gives each cell line's fractional membership
across the k=4 archetypes (a convex combination — this is exactly the
specialist/generalist readout used in Hausser 2019 and Groves 2022's
own Pareto Task Inference framework, not a new distance metric we
invented). We take the normalized Shannon entropy of each row:

- entropy → 0: weight concentrated on one archetype (specialist)
- entropy → 1: weight spread evenly over all k=4 archetypes (generalist)

Output: `results/panel_d_emt/archetype_weights_and_generalist_score.csv`.

## Step 3 — KS-epi / KS-mes (`src/enrichment.ks_epi_mes`)

Two pre-specified two-sample Kolmogorov–Smirnov tests (no multiple
testing correction needed for two planned comparisons):

- **KS-epi**: generalist score, E lines (n=42) vs E/M lines (n=10)
- **KS-mes**: generalist score, M lines (n=11) vs E/M lines (n=10)

This run:

| test | D | p | median(group) | median(E/M) |
|---|---|---|---|---|
| KS-epi | 0.371 | 0.166 | 0.491 (E) | 0.599 |
| KS-mes | 0.609 | **0.030** | 0.387 (M) | 0.599 |

Plus a non-thresholded robustness check: Spearman rho(mu,
generalist_score) across all 63 lines was -0.01 (p=0.95) — expected,
since a monotone correlation is the wrong shape to test an "extremes
are specialists, middle is generalist" (inverted-U) hypothesis; see
`Figure_1D_breast_emt_generalist.png` panel 1 for the actual shape
instead of relying on this number alone.

## Reading these numbers (breast, k=4, n=63, this run)

- E/M lines have a **higher median generalist score than the
  mesenchymal lines, and the difference is statistically significant**
  (KS-mes p=0.030).
- E/M lines also have a higher median generalist score than epithelial
  lines, but **that difference is not significant** at n=10 vs n=42
  (KS-epi p=0.166) — likely underpowered (only 10 hybrid lines), not
  necessarily a true null; do not over-read this as "no relationship
  to epithelial specialists."
- Interpretation: on this dataset, **hybrid E/M lines behave more like
  generalists relative to mesenchymal specialists**, with a directionally
  consistent but not-yet-significant relationship to epithelial
  specialists. This is a real, testable partial answer to "are EMT
  hybrids generalists" — not a settled yes, and n=10 hybrid lines is
  small.
- Cross-reference with Panel B: LumB peaks at archetype 3, Basal peaks
  at archetype 4 (`enrichment_pam50.csv`) — i.e. archetype 4 is this
  polytope's best candidate for a "mesenchymal-leaning" vertex, useful
  if you want to look at KS-mes broken down per-archetype instead of
  pooled across all k=4 (not done here; `distance_bins()` +
  `ks_two_group()` from `src/enrichment.py` support that directly).

## Caveats carried over from the rest of this repo

- Panel A's t-ratio permutation test was **not significant** for
  breast at any k tested (`results/panel_a/t_ratio_parti_500.csv`,
  best p=0.082 at k=4) — i.e. the polytope itself is not a
  statistically confirmed simplex. Panel D's generalist score still
  works mechanically (PCHA always returns weights for any k), but
  "generalist/specialist" is a much weaker claim on a
  not-significantly-simplex-shaped point cloud than it would be on a
  confirmed one (e.g. SCLC in the Groves reproduction). State this
  caveat wherever Panel D results are reported.
- n=10 hybrid E/M lines is small; KS tests are exact but low-powered
  here. Do not report KS-epi's non-significance as a negative result
  without that caveat.
- The EMT model was trained on NCI-60 (pan-cancer), not on breast
  specifically, and CLDN7/CDH1/VIM alone is a coarse 2-gene score.
  Cross-check against a fuller signature (e.g. Tan et al. 2014's KS-
  based EMT score, or ssGSEA on a larger EMT gene set) before treating
  mu as a precise per-line phenotype.

## How to rerun

From the repository root, after Panels A and B already exist for breast:

```bash
.venv/bin/python -u "Breast Cancer/codes/05_emt_generalist_panelD.py"
```

Reads only existing `results/panel_a/` and `results/panel_b/` outputs
plus the Panel A input matrix; writes only to `results/panel_d_emt/`
and the two `Figure_1D...` image paths. Safe to rerun any time; it
does not call PCHA again.
