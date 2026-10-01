# Primer — Anderson Fig 3: what the paper did, what we tried, what’s next

**Audience:** you, starting from scratch.  
**Paper:** Anderson et al., PML archetypes in LUAD premalignant lesions (MCR / bioRxiv).  
**Folder:** `Anderson_PML_Archetypes/`

---

# Part I — What Figure 3 is (basic → advanced)

## 1. Basic idea (no jargon)

The authors took **137 premalignant lung lesions** (early abnormal spots, not full tumors yet) and asked:

> Do these lesions sit in a few extreme “personalities,” or is every lesion a unique snowflake?

They found **four extreme types** (archetypes):

| Name | Rough meaning |
|---|---|
| **A1 Normal-like** | Looks more like normal airway biology |
| **A2 Inflammation** | Immune / inflammatory programs high |
| **A3 Cell adhesion** | Adhesion / matrix-ish programs |
| **A4 Proliferation** | Cell-cycle / growth programs |

**Figure 3** is the picture of that: every lesion as a point on a 2D map, with the four extremes as corners of a shape.

- **Fig 3A:** same map four times; each panel colors “how close am I to this one extreme?”  
- **Fig 3B:** one map; points near one extreme get that color; everyone else is grey (“Secondary”); circles = AAH, triangles = AIS/MIA.

**Intuition:** think of four corners of a room. Most people stand inside the room (mixtures). A few stand near one corner (almost pure type). Grey = middle of the room.

---

## 2. Intermediate: what data goes into the map

### Samples
137 PMLs from 4 studies (cohorts):

| Cohort | Approx. n | Histology |
|---|---:|---|
| Lung PCA (GSE319666) | 80 | mostly AAH + some AIS/MIA |
| GSE102511 | 17 | AAH |
| GSE166720 | 17 | AIS/MIA |
| GSE193725 | 23 | AIS/MIA |

### From genes → 9 numbers per lesion
They don’t plot all ~20,000 genes. They:

1. Clean expression (normalize, remove study effects, patient correlation, etc.).  
2. Score each lesion on **9 gene modules** (bundles of genes that move together) using **GSVA**.  
3. So each lesion becomes a point in **9-dimensional** space (one score per module).

### From 9D → 4 corners (archetypes)
**Archetypal analysis (PCHA / ParetoTI)** finds the tightest “tetrahedron” (4 corners in high-D) that wraps the cloud of lesions.

Those 4 corners are A1–A4 (after they rename them by biology).

### From high-D → 2D plot
**PCA** compresses the 9 module scores to PC1 and PC2 for the eye.  
The black shape in Fig 3 is the 4 corners drawn in that 2D view.

### Who is “primary”?
For each corner, take the closest **12.5%** of lesions (an “octile”).  
If a lesion is uniquely in exactly one corner’s octile → that color.  
Otherwise → Secondary (grey).

---

## 3. Advanced: exact author pipeline (from their public code)

Repo: `github.com/kelley27/PML_Archetypes` (mirrored at `reference/author_repo/`).

| Step | Script | What it does |
|---|---|---|
| Residuals | `Step0_compute_residuals.R` | **Counts** → TMM → **voom** → design `~ Study + TopGenes` → `duplicateCorrelation(Patient)` (two-pass) → limma residuals |
| Modules | `Step1`–`Step2` | Build the 9 modules (already published as `pml_modules.rds`) |
| Choose k=4 | `Step3` | Elbow + t-ratio tests |
| Fig 3 core | `Step4_compute_archetypes.R` | GSVA on residuals → ParetoTI bootstrap (n=5000, 80% subsample, seed 42792) → average vertices → distances in **module space** → PCA for plot → octile bins → **hardcoded rename** Arch1→A1, Arch2→A4, Arch3→A2, Arch4→A3 |
| External data | `Step5_…` | Don’t re-fit corners; **project** new GSVA onto saved `PML_arcA4_pml.rds` |

**Critical details:**

1. **`voom` needs raw counts**, not TPM/FPKM.  
2. Distances for coloring/assignment are in **9-module space**, not “eyeball distance on the PNG.”  
3. PCA axis signs are arbitrary (left/right can flip with no biology change).  
4. Biological names are a **manual relabel** of one specific fit — not a universal algorithm.

**What is public vs gated:**

| Public | Request from authors |
|---|---|
| Code (Step0–5) | `PML_resid_duplicateCor.rds` (the residual matrix) |
| `pml_modules.rds` | Per-cohort count/residual objects |
| `PML_arcA4_pml.rds` (fitted corners!) | |
| `PML_Archetypes.rds` (published labels/distances) | |

So: **code is open; the sample-level residual matrix that Fig 3 was built on is not.**

---

# Part II — What we tried (in order) and what we got

Goal: rebuild Fig 3A/3B for all 137 lesions, paper colors, as close as possible.

## Attempt 1 — v1 public-table rebuild (`results/fig3/`)

**What we did**

- Lung PCA: counts-ish log.  
- Other 3 cohorts: public **TPM/FPKM** → log2(x+1).  
- Harmonize → limma residuals `~ Study + TopGenes` (**no** patient blocking yet).  
- GSVA with author modules.  
- Fit our own PCHA in Python → plot.

**Result**

- Figures existed; colors matched.  
- PC variance roughly in the ballpark (~20% / 16% vs paper ~22% / 15%).  
- Label overlap with authors OK-ish; **shape of the figure already wrong** (corners not where paper puts them).

---

## Attempt 2 — SRA → real counts (parked)

**Idea:** download FASTQs for the 57 external samples, Salmon-quant to counts, then run true Step0.

**Result:** downloads stalled repeatedly (0/57 finished). You asked to stop blocking on this.

---

## Attempt 3 — Path B / v2 (`results/fig3_v2/`)

**What we did**

- Same public expression matrix.  
- Add **`duplicateCorrelation(Patient)`** (paper-like residual design).  
- GSVA again.  
- ParetoTI wouldn’t install (Bioc deps) → **500× bootstrap PCHA** instead.  
- **Bug (later found):** fitted/distances in **PC1–3**, not full 9D GSVA.

**QC**

| Metric | Value |
|---|---|
| Primary agreement (both primary) | ~82% |
| Signed-distance correlation | ~0.54–0.83 |

**Visual:** still starkly different from paper (e.g. A1 on opposite side of the plot).

---

## Attempt 4 — Author repo read + project onto their corners (`results/fig3_v2_author_arc/`)

**Insight:** we already had `PML_arcA4_pml.rds` = their **real** fitted tetrahedron.  
Like Step5: don’t re-fit; measure distance from our GSVA points to **their** corners.

**QC (big jump)**

| Metric | Re-fit v2 | Author corners |
|---|---|---|
| Primary agreement | ~82% | **~94%** |
| Distance r | 0.54–0.83 | **0.91–0.95** |

**Visual:** labels/biology much better; **plot still looks wrong** — many points outside the black shape; axes still not paper-like until cosmetic flip.

---

## Attempt 5 — Fix re-fit code (`results/fig3_v2_refit9d/`)

**Bug:** `08` used PC subspace; paper uses 9D GSVA.  
**Fix:** `10_refit_pcha_gsva9d.py` — re-fit in 9D.

**Result:** method closer to paper; QC **not** better than locking author corners; figure **still** doesn’t look like Fig 3.

---

## Attempt 6 — Cosmetic PC1 flip (`results/fig3_v2_author_arc_pc1flip/`)

Flip PC1 so A1 is on the **right** like the paper.

**Result:** sides look closer; **points still spill outside the hull** → still “looks bad.”

---

# Part III — Why it still looks bad (the bottleneck)

## The one-sentence bottleneck

> We never had the same **input matrix** the paper’s Fig 3 was built from (count→voom residuals for all 137). Everything downstream is then a different cloud — so corners, PCA, and “points inside the shape” cannot match the screenshot.

## Unpack that

1. **~42% of samples** (3 cohorts) only have public TPM/FPKM.  
2. Paper’s Step0 is **voom on counts**. Feeding log(TPM+1) into a limma-like residual step is **not** a small error; it changes the geometry.  
3. Archetypes are **extremes**. Warped clouds move corners a lot.  
4. Fig 3 is a **discovery plot** of that cloud. Unlike Breast/Hausser-style figures where you **project into a fixed triangle**, here the map *is* the result. Wrong map → wrong picture.  
5. Code bugs (PC vs 9D) and PCA flips were real but **secondary**. Fixing them improved honesty of the method; they did not conjure the paper’s cloud.

## Two questions people mix up

| Question | Our best answer with public data |
|---|---|
| “Do we recover **who** is near which archetype?” | Pretty well if we use **their** corners (~94% / r≈0.9) |
| “Does the **PNG** look like Fig 3?” | **No** — points don’t sit in their shape |

You were judging (rightly) by the second. We kept reporting the first. That’s why it felt like gaslighting.

---

# Part IV — Insights (what we learned)

1. **Method instinct was mostly right** — Study + TopGenes + patient blocking matches Step0; GSVA settings match; octile logic matches.  
2. **Data gate ≠ code gate** — repo is public; residuals are not.  
3. **Re-fitting on bad GSVA makes the picture worse** than projecting onto published corners.  
4. **Author names are a hardcoded map** from one ParetoTI run — so “recover A1” means recover that biology, not a universal vertex index.  
5. **Breast Fig 4 ≠ this problem** — fixed reference space vs rebuild-the-discovery-map.  
6. **“Looks bad” = geometric mismatch**, not mainly wrong color names.

---

# Part V — Next steps (priority order)

1. **Email authors** for `PML_resid_duplicateCor.rds` (and/or per-cohort residuals). Draft lives in `docs/HANDOFF.md`. This is the highest-ceiling fix.  
2. **Re-check GEO** for raw count tables on GSE102511 / GSE166720 / GSE193725 (sometimes both TPM and counts exist).  
3. **Keep** `fig3_v2_author_arc` (+ optional PC1 flip) as the honest “best public-data” product — report QC numbers, don’t claim pixel match.  
4. **SRA/Salmon** only if GEO has no counts and authors don’t share.  
5. **ParetoTI install** only after real residuals/counts exist.

---

# Part VI — File cheat sheet

| What | Where |
|---|---|
| This primer | `docs/PRIMER_FIG3_FROM_SCRATCH.md` |
| Full session log | `docs/SESSION_2026-09-09_FIG3.md` |
| Author-repo addendum | `docs/ADDENDUM_2026-09-09_AUTHOR_REPO.md` |
| Agent handoff | `docs/HANDOFF.md` |
| Paper crop | `results/fig3/paper_fig3AB_crop.png` |
| Best public-data figs | `results/fig3_v2_author_arc/` (and `_pc1flip/`) |
| Corrected re-fit (still not pretty) | `results/fig3_v2_refit9d/` |
| Author code mirror | `reference/author_repo/` |

---

# One paragraph to remember

Figure 3 is a map of 137 lesions in a space built from **count-based residuals → 9 module scores → four extreme corners → PCA view**. We rebuilt that path with **public TPM/FPKM for many samples**, so we built a **different map**. Fixing code and locking their published corners improved **who matches whom**; it cannot force our points to sit inside their shape. Until we get their residual matrix (or true counts), a screenshot-faithful Fig 3 is out of reach — and that’s expected, not a failure of understanding.
