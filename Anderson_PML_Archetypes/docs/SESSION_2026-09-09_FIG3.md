# Session note — Anderson Fig 3A/3B (2026-09-09)

**Workspace:** `/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes`  
**Paper:** Anderson et al., *Mol Cancer Res* (MCR-25-1093) — `papers/mcr-25-1093.pdf`  
**Author GitHub:** https://github.com/kelley27/PML_Archetypes  
**Goal:** Reproduce **Figure 3A** (signed distance to each archetype in PC1–PC2 of 9-module GSVA) and **Figure 3B** (primary top-12.5% octile vs secondary; shapes = histology).

This note covers what was built today, where files live, how to re-run, QC numbers, and why the replica still looks starkly different from the paper.

---

## 1. Executive summary

| Item | Outcome |
|---|---|
| Target | Fig 3A + 3B for all **137** PMLs across **4** cohorts |
| Best current rebuild | **Path B / v2** → `results/fig3_v2/` |
| Earlier rebuild | **v1** → `results/fig3/` |
| Count-level rebuild via SRA | **Parked** (downloads stall; 0/57 finished) |
| ParetoTI (paper’s R package) | **Not installed** (Bioconductor Depends failed) → used **500× bootstrap PCHA** in Python instead |
| Looks like paper corners? | **No** — layout of A1–A4 in PC space differs strongly |
| Labels / biology roughly OK? | **Partially** — primary-label agreement **82.4%** (n=34 both-primary); distance Pearson r **0.54–0.83** |
| Honest description | “Improved public-data rebuild,” **not** a pixel-faithful Fig 3 |

**Next step that would actually fix the map:** obtain author `PML_resid_duplicateCor.rds` (Path A), or finish true counts for the three external cohorts.

---

## 2. Paper method (what we are aiming at)

1. Expression for 137 PMLs (Lung PCA + GSE102511 + GSE166720 + GSE193725).  
2. Counts → **TMM / voom** → limma residuals `~ Study + TopGenes` with **`duplicateCorrelation(Patient)`**.  
3. **GSVA** on 9 published PML modules (`pml_modules.rds`).  
4. **ParetoTI** bootstrap archetypes (k=4; paper uses large bootstrap, seed `42792`).  
5. PCA of GSVA scores → plot PC1–PC2.  
6. **3A:** color by signed (−) distance to each archetype vertex.  
7. **3B:** primary = unique membership in closest **12.5%** octile per archetype; else secondary; AAH = circle, AIS/MIA = triangle.

**Paper PC variance (approx.):** PC1 **21.7%**, PC2 **14.6%**.

**Paper colors (matched in our plots):**

| Archetype | Hex |
|---|---|
| A1 Normal-like | `#416c3e` |
| A2 Inflammation | `#67b3bb` |
| A3 Cell adhesion | `#843385` |
| A4 Proliferation | `#ddb137` |
| Secondary | `#d0d0d0` |

---

## 3. Data used today

### 3.1 Cohorts (n = 137)

| Study | Histology | n | Public expression on GEO |
|---|---|---:|---|
| Lung_PCA (GSE319666) | AAH | 68 | **Counts** ✓ |
| Lung_PCA (GSE319666) | AIS/MIA | 12 | **Counts** ✓ |
| GSE102511 | AAH | 17 | TPM only |
| GSE166720 | AIS/MIA | 17 | FPKM only |
| GSE193725 | AIS/MIA | 23 | FPKM only |

**Blocker:** Paper pipeline needs **counts** for all cohorts. Only Lung PCA has public counts. The other three cohorts (~57 samples) are TPM/FPKM tables only. Author residual RDS is request-only.

### 3.2 Author assets (on disk)

Under `data/author/`:

| File | Role |
|---|---|
| `pml_modules.rds` | 9 gene modules for GSVA |
| `PML_Archetypes.csv` / `.rds` | Published distances, bins, `Combined` primary/secondary labels |
| `PML_arcA4_pml.rds` | Author archetype object (reference) |
| `ENS_SYM_ID.rds`, `ENS_ENT_ID.rds` | ID maps |

### 3.3 Processed matrices we built

Under `data/processed/` (mostly from **v1** build):

| File | Description |
|---|---|
| `PML_expr_log_harmonized.rds` | Harmonized log expression: Lung PCA log-CPM-like + log2(TPM/FPKM+1) for others |
| `PML_sample_info.csv` | sample_id, study, histology, patient, TopGenes, author_combined |
| `PML_resid.rds` / `.csv` | **v1** limma residuals (`~ Study + TopGenes`, **no** patient blocking) |
| `public_PML_gsm_map.csv` | GSM ↔ cohort map |
| `modules/` | Module gene lists exported for convenience |

**Path B residuals** live under results (not re-copied to `data/processed/`):

| File | Description |
|---|---|
| `results/fig3_v2/PML_resid_v2.rds` | limma residuals with **`duplicateCorrelation(Patient)`** on the harmonized log matrix |
| `results/fig3_v2/gsva_scores.csv` | 9-module GSVA (Gaussian) on v2 residuals |

### 3.4 SRA / Salmon attempt (parked)

- Mapped 57 external GSMs → SRA/ENA FASTQs.  
- Built Salmon index (GENCODE v44).  
- Scripts: `codes/04_quant_sra_cohort.sh`, `codes/04b_quant_ena_resumable.py`, `data/sra/run_all_public_pml.sh`.  
- **Status:** downloads repeatedly stalled; **0/57** quantified. Deprioritized by user preference (unblock Path B instead of waiting on SRA).

---

## 4. What we did (pipeline chronology)

### Phase A — v1 approximate rebuild (`results/fig3/`)

1. Build harmonized expression matrix from GEO tables.  
2. limma residuals `~ Study + TopGenes` (**no** `duplicateCorrelation`).  
3. GSVA on author modules.  
4. Python multi-start **PCHA** (project `src.archetypes`), not ParetoTI bootstrap.  
5. Plot Fig 3A/3B with paper colors; QC vs author CSV.

**Scripts:** `codes/01_build_pml_matrix.R` → `02_compute_residuals.R` → `02b_gsva.R` → `03_fit_archetypes_plot_fig3.py`

**v1 QC (approx.):** primary agreement ~88% (earlier dual-primary definition); distance r ~0.53–0.83; PC ~20.5% / 15.9%. Shape still off vs paper crop.

### Phase B — SRA counts (attempted, stopped)

Tried to replace TPM/FPKM with Salmon counts for the three external cohorts. Blocked on network/download stalls. Left in place but not required for current figures.

### Phase C — Path B / v2 (`results/fig3_v2/`) — **current best**

1. Recompute residuals on `PML_expr_log_harmonized.rds` with:
   - design `~ Study + TopGenes`
   - **`duplicateCorrelation` + `block = Patient`**
2. GSVA (Gaussian) on author `pml_modules.rds`.  
3. Tried ParetoTI → **install failed** (hard Depends: BioQC, AnnotationHub, GO.db, AUCell, etc.).  
4. Fallback: **500 bootstrap PCHA** fits (80% subsample, seed `42792`), average aligned vertices (`codes/08_…`).  
5. Map unordered PCHA vertices → paper names A1–A4 by overlap with author primary samples.  
6. Primary = unique top **12.5%** octile; plot 3A/3B; QC.

**Scripts:** `codes/07_pathB_resid_gsva_archetypes.R` → `codes/08_pathB_bootstrap_pcha_plot.py`

**Also written but not the live Path B path:** `05_residuals_paperlike.R`, `06_gsva_paretoti.R` (assume count/voom-style inputs / ParetoTI).

---

## 5. Codes (inventory)

| Script | Role | Used for current figs? |
|---|---|---|
| `01_build_pml_matrix.R` | GEO → harmonized log matrix + sample info | Yes (upstream of both v1/v2) |
| `02_compute_residuals.R` | v1 limma residuals | v1 only |
| `02b_gsva.R` | v1 GSVA | v1 only |
| `03_fit_archetypes_plot_fig3.py` | v1 PCHA + Fig 3A/B | v1 only |
| `04_quant_sra_cohort.sh` | SRA/Salmon wrapper | Parked |
| `04b_quant_ena_resumable.py` | Resumable ENA download + Salmon | Parked |
| `05_residuals_paperlike.R` | Count/voom-style residuals draft | Not run end-to-end today |
| `06_gsva_paretoti.R` | GSVA + ParetoTI draft | Blocked (ParetoTI) |
| `07_pathB_resid_gsva_archetypes.R` | **v2** patient-blocked residuals + GSVA | **Yes** |
| `08_pathB_bootstrap_pcha_plot.py` | **v2** bootstrap PCHA + Fig 3A/B + QC | **Yes** |

Shared Python archetype engine: repo `src/archetypes.py` (`fit_pcha`, `simplex_volume`).

### Re-run Path B

```bash
export KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 MPLBACKEND=Agg
cd /Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes

Rscript codes/07_pathB_resid_gsva_archetypes.R
../.venv/bin/python -u codes/08_pathB_bootstrap_pcha_plot.py
```

R packages: local lib `Anderson_PML_Archetypes/rlib` (limma, GSVA, etc.).  
Python: repo root `.venv`.

---

## 6. Figures and outputs

### 6.1 Paper reference crops

| File | Notes |
|---|---|
| `results/fig3/paper_fig3AB_crop.png` | Cropped paper Fig 3A/B for side-by-side comparison |
| `results/fig3/paper_fig3_page.png` | Broader page capture |

**Paper layout (approximate vertex positions in PC1–PC2):**

| Archetype | Region on paper Fig 3B |
|---|---|
| A1 Normal-like (green) | **Bottom-right** |
| A2 Inflammation (teal) | **Top** |
| A3 Cell adhesion (purple) | **Far left** |
| A4 Proliferation (gold) | **Bottom-left** |

### 6.2 v1 outputs — `results/fig3/`

- `fig3A_distance_panels.png` / `.pdf`
- `fig3B_primary_assignments.png` / `.pdf`
- `gsva_scores.csv`, `sample_assignments.csv`, `archetypes_pc.csv`, `pc_scores.csv`

### 6.3 v2 outputs (current) — `results/fig3_v2/`

| File | Description |
|---|---|
| `fig3A_distance_panels.png` / `.pdf` | Four panels, signed distance coloring |
| `fig3B_primary_assignments.png` / `.pdf` | Primary vs secondary + histology shapes |
| `PML_resid_v2.rds` / `.csv` | Patient-blocked residuals |
| `gsva_scores.csv` / `.rds` | Module scores |
| `sample_assignments.csv` | Our primary, author Combined, distances, PC1/PC2 |
| `archetypes_pc.csv` | Vertex coordinates in PC space |
| `sample_info_v2.csv` | Sample annotation copy used in v2 |
| `fit_method.txt` | `pcha_python` |
| `plot_log.txt` | Bootstrap + QC log |

**v2 archetype PC coordinates (after label mapping):**

| Archetype | PC1 | PC2 |
|---|---:|---:|
| A1 Normal-like | −0.90 | −0.59 |
| A2 Inflammation | −1.11 | +0.39 |
| A3 Cell adhesion | +1.02 | −0.69 |
| A4 Proliferation | +0.49 | +0.95 |

**v2 layout on the page:** A1 bottom-left, A2 top-left, A3 bottom-right, A4 top — **not** the paper layout.

**v2 PC variance:** PC1 **20.6%**, PC2 **16.1%** (paper ~21.7% / 14.6%).

---

## 7. QC vs authors (Path B / v2)

From `results/fig3_v2/plot_log.txt`:

| Metric | Value |
|---|---|
| Bootstrap fits kept | 500 / 500 |
| Primary agreement (both call primary) | **82.4%** (n=34) |
| Signed-distance Pearson r vs author | A1 **0.67**, A2 **0.54**, A3 **0.83**, A4 **0.56** |
| Our primary counts | 18 each archetype + 65 Secondary |
| Author Combined counts | A1 12, A2 13, A3 16, A4 15, Secondary 81 |

### Histology composition (sanity check)

Trend **directionally** matches paper / author labels (A1/A2 more AAH; A4 more AIS/MIA), but is weaker / noisier than the published figure:

| Set | A1 %AAH | A2 %AAH | A4 %AIS/MIA |
|---|---:|---:|---:|
| Author Combined primaries | 83% | 77% | 53% |
| Our primaries | 61% | 78% | 67% |

Author-A1 samples in **our** PC space sit near mean PC ≈ **(−0.9, −0.3)** (left), whereas the paper draws A1 on the **right**. Names follow samples; page corners do not match.

---

## 8. Discrepancies (detailed)

### 8.1 Visual layout (main user concern)

| | Paper Fig 3B | Our Fig 3B (v2) |
|---|---|---|
| A1 | Bottom-**right** | Bottom-**left** |
| A2 | **Top** | Top-**left** |
| A3 | **Far left** | Bottom-**right** |
| A4 | Bottom-**left** | **Top** |

This is stark. It is **not** explained by a single PC1 sign flip alone: after flipping PC1, A1/A3 move toward paper sides, but **A2 vs A4 vertical roles still differ** (paper top vertex = A2; ours = A4). So the fitted tetrahedron geometry itself differs, not only axis polarity.

### 8.2 Why “labels OK” can coexist with “figure looks wrong”

We label vertices by **sample overlap** with author primaries, not by forcing corners to paper PC coordinates.

- **Who:** sample X is near our green vertex and author-A1 → call it A1. (~82% when both primary.)  
- **Where:** that green vertex can still sit on the opposite side of the plot from the paper.

Comparing **page regions** (paper right vs our right) mixes different biologies. Comparing **named colors** is the fairer check — and even that is only approximate.

### 8.3 Root causes (data / method gaps)

1. **Expression units wrong for ~42% of samples**  
   GSE102511 / 166720 / 193725 entered as TPM/FPKM, not counts→TMM/voom. Harmonized log matrices are a substitute, not equivalent.

2. **Residuals computed on that substitute matrix**  
   Path B adds patient blocking (good), but still cannot recover the paper’s residual space without count-level (or author residual) input.

3. **Archetypes are extremes**  
   PCHA vertices sit on the edge of the GSVA cloud. Warping the cloud moves corners a lot even if many mid-cloud sample–sample relationships are vaguely preserved → large visual change in 3A/3B.

4. **2D PCA projection**  
   Different 9-D clouds → different PC1/PC2 orientations (including sign flips) → tetrahedron can look rotated/reflected on the page.

5. **ParetoTI vs bootstrap PCHA**  
   We approximated paper bootstrap with project PCHA (500 fits, not ParetoTI n=5000). Secondary difference relative to (1)–(3).

6. **Primary calling counts differ**  
   We assign 18 primaries per archetype (tight octile uniqueness); authors have fewer primaries / more Secondary (81). Thresholding / uniqueness rules or distance scales differ.

### 8.4 What is *not* the main bug

- Not “we randomly swapped color names with no sample basis.” Vertex→name mapping is data-driven from author Combined labels.  
- Not comparable failure mode to Breast Cancer / Hausser-style panels that **project new cells into a fixed reference triangle** (DepMap / METABRIC). Those lock corner positions by construction. Anderson Fig 3 **is** the discovery map; wrong upstream matrix moves the whole map.

### 8.5 Contrast with `Breast Cancer/figures/Figure_4_ks_panelA_gse173634_sc_magic.png`

| | Breast / Hausser-style Fig 4 check | Anderson Fig 3 |
|---|---|---|
| Space | Fit once on **available** reference (DepMap / METABRIC) | Must rebuild paper’s GSVA/PC space |
| New data | **Projected** into fixed PCs / triangle | Used to **re-fit** cloud + archetypes |
| Success criterion | Cells fall in triangle / subtype near vertices | Corners + distances match published figure |
| Why Breast looks “fine” | Corners locked | Corners free → sensitive to missing counts |

---

## 9. Decisions made with the user today

1. Keep **all 4 cohorts / 137 PMLs** (no tiny pilot).  
2. Use **paper colors** for comparison.  
3. After SRA pain: **prefer Path B unblock** over waiting on counts.  
4. Do **not** claim pixel-identical Fig 3 without author residuals or much higher distance correlation / matching layout.  
5. Clarified: page-corner mismatch is expected given current inputs; named-archetype biology is only partially recovered.

---

## 10. Recommended next steps

**Superseded by** [`ADDENDUM_2026-09-09_AUTHOR_REPO.md`](ADDENDUM_2026-09-09_AUTHOR_REPO.md) §7 after the author repo went public. Short version:

1. **Try first:** project v2 GSVA onto fixed vertices in `data/author/PML_arcA4_pml.rds` (no re-fit).  
2. Recheck GEO for raw counts on the three external cohorts.  
3. **Still highest ceiling:** email for `PML_resid_duplicateCor.rds`.  
4. ParetoTI install / SRA: lower priority until counts or author residuals land.

---

## 11. Related docs

| Doc | Purpose |
|---|---|
| `docs/HANDOFF.md` | Agent handoff / next actions / email draft |
| `docs/HIGH_FIDELITY_STATUS.md` | Short note that SRA is parked |
| `docs/README.md` | Folder overview / run commands |
| `docs/SESSION_2026-09-09_FIG3.md` | **This file** — today’s full record |
| `docs/ADDENDUM_2026-09-09_AUTHOR_REPO.md` | Author repo public: confirmed methods, data still gated, projection shortcut |

---

## 12. Quick file map

```
Anderson_PML_Archetypes/
  codes/
    01–03          # v1 pipeline
    04 / 04b       # SRA (parked)
    05–06          # paperlike drafts (ParetoTI path)
    07–08          # Path B / v2 (CURRENT)
  data/
    author/        # modules + published assignments
    processed/     # harmonized expr + sample info (+ v1 resid)
    sra/           # partial download/quant machinery
  results/
    fig3/          # v1 figs + paper crops
    fig3_v2/       # CURRENT Path B figs + resid/GSVA/QC
  docs/
    HANDOFF.md
    SESSION_2026-09-09_FIG3.md
  rlib/            # local R packages
```
