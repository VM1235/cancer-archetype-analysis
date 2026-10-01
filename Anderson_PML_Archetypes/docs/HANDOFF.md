# Handoff: Anderson et al. Fig 3A/3B

**Updated:** 2026-09-09 (~03:30 IST)  
**Workspace:** `/Users/apple/Desktop/thesis_iisc/Project_1`  
**Folder:** `Anderson_PML_Archetypes/`  
**Paper:** `papers/mcr-25-1093.pdf`  
**Author code:** https://github.com/kelley27/PML_Archetypes  

**Full session write-up (data, codes, figures, discrepancies):**  
[`docs/SESSION_2026-09-09_FIG3.md`](SESSION_2026-09-09_FIG3.md)

**From-scratch primer (paper Fig 3 → our attempts → bottleneck):**  
[`docs/PRIMER_FIG3_FROM_SCRATCH.md`](PRIMER_FIG3_FROM_SCRATCH.md)

**Addendum (author repo now public — methods confirmed; residuals still request-only; try projecting onto `PML_arcA4_pml.rds`):**  
[`docs/ADDENDUM_2026-09-09_AUTHOR_REPO.md`](ADDENDUM_2026-09-09_AUTHOR_REPO.md)

**Projection experiment (done):** `results/fig3_v2_author_arc/` — v2 GSVA → author XC; primary agreement **93.5%**, distance r **0.91–0.95**. Page layout still PC1-flipped vs paper.

---

## 1. Status (read this first)

**Path B (v2) is DONE** — improved public-data rebuild without SRA:

- `results/fig3_v2/fig3A_distance_panels.png`
- `results/fig3_v2/fig3B_primary_assignments.png`
- Patient-blocked limma residuals + GSVA + **500× bootstrap PCHA** (ParetoTI install failed on Bioc deps)
- QC vs authors: primary agreement **82.4%** (n=34); distance r ≈ **0.54–0.83**
- PC var **20.6% / 16.1%** (paper 21.7% / 14.6%)

**v1** still in `results/fig3/` (earlier PCHA, no bootstrap).

**SRA path remains DEPRIORITIZED** (downloads stall). Next meaningful upgrade = **author residual RDS** (Path A).

---

## 2. Recommended ways forward (no SRA)

### Path A — Best for “looks like the paper” (recommended next)
**Email authors for residual RDS**, then run their exact Step4.

Request from corresponding author (Jennifer Beane, `jbeane@bu.edu`) / Kelley Anderson GitHub (`kelley27`):

- `PML_resid_duplicateCor.rds`
- optionally `dge_PML.rds`, `PML_arcA4_pml.rds` (arc already on GitHub)

Then:

```bash
# put RDS in Anderson_PML_Archetypes/data/author_input/ or similar
# reuse their Step4 logic / our 06_gsva_paretoti.R adapted to that resid
```

Author GitHub already has: `pml_modules.rds`, `PML_Archetypes.rds`, `PML_arcA4_pml.rds`.

### Path B — Best progress *now* without waiting (do this while emailing)
Improve the **existing public-table rebuild** (no new downloads):

1. Install **ParetoTI** and fit with paper settings (`n=5000`, `sample_prop=0.8`, `seed=42792`) on current GSVA scores  
2. Add `duplicateCorrelation(Patient)` to limma residuals on current harmonized matrix  
3. Replot with paper colors  
4. QC vs `data/author/PML_Archetypes.csv` (distance Pearson r + label agreement)

Inputs already on disk:
- `data/processed/PML_resid.rds` / `PML_expr_log_harmonized.rds`
- `data/processed/PML_sample_info.csv`
- `results/fig3/gsva_scores.csv`
- author modules/labels in `data/author/`

```bash
export KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 MPLBACKEND=Agg
cd Anderson_PML_Archetypes

# Install ParetoTI once
Rscript -e '.libPaths("rlib"); remotes::install_github("vitkl/ParetoTI", lib="rlib", upgrade="never")'

# Recompute residuals with patient blocking (extend 02 or 05 for non-count expr)
# Then GSVA + ParetoTI + plot (adapt 06 + 03)
```

Honest label for user: **“improved public-data rebuild”**, still not count/voom-identical.

### Path C — Optional later (only if user insists)
Resume SRA quant (`data/sra/run_all_public_pml.sh` / `codes/04b_quant_ena_resumable.py`).  
Expect overnight+ and frequent stalls. Partial FASTQ resume works with curl `-C -` if work dir is kept. **Do not treat as critical path.**

---

## 3. What the user wants

1. Fig **3A** (distance-colored PC panels) + **3B** (primary octile + histology shapes)  
2. Paper **colors** (green / teal / purple / gold)  
3. Methods as close as possible  
4. All **4 cohorts / 137 PMLs**  
5. After SRA pain: **prefer unblocking over perfect counts**

---

## 4. v1 deliverables (done)

| Item | Path |
|---|---|
| Fig 3A | `results/fig3/fig3A_distance_panels.png` |
| Fig 3B | `results/fig3/fig3B_primary_assignments.png` |
| Assignments | `results/fig3/sample_assignments.csv` |
| Paper color reference crop | `results/fig3/paper_fig3AB_crop.png` |

**v1 method:** Lung PCA counts + 3× TPM/FPKM → log → limma `~ Study + TopGenes` → GSVA → Python PCHA → plot.

**QC vs authors:** primary agreement ~88% (dual-primary); distance r ~0.53–0.83 → shape still off.

**Paper colors:**
- A1 `#416c3e` · A2 `#67b3bb` · A3 `#843385` · A4 `#ddb137` · Secondary `#d0d0d0`
- 3A: sequential per-archetype (not RdYlBu); 3B: AAH `o`, AIS/MIA `^`

---

## 5. Why perfect match needs more than GEO tables

Paper: counts → TMM/voom → duplicateCorrelation → GSVA → ParetoTI bootstrap.

| Cohort | On GEO | Blocker |
|---|---|---|
| GSE319666 Lung PCA | **counts** ✓ | none |
| GSE102511 / 166720 / 193725 | TPM/FPKM only | no public count matrix |

Author residuals = request-only. SRA FASTQs exist but downloads stall on this machine.

---

## 6. Suggested next-agent plan (in order)

1. **Stop** any stuck `04b`/`curl`/`run_all_public` jobs.  
2. Draft short **author data-request email** (Path A) for the user to send.  
3. **Path B now:** install ParetoTI; tighten residuals + refit + replot into `results/fig3_v2/` or `results/fig3_paperlike/`.  
4. Report distance-correlation / label agreement vs author CSV.  
5. Only revisit SRA if user explicitly wants it again.

---

## 7. Draft email (Path A)

Subject: Request for residual expression matrix — PML_Archetypes (MCR-25-1093)

> Dear Dr. Beane / Dr. Anderson,  
> I’m reproducing Figure 3 from your MCR 2026 PML archetype paper using your public GitHub (`kelley27/PML_Archetypes`). GEO provides TPM/FPKM for GSE102511/GSE166720/GSE193725 rather than count matrices. Could you share `PML_resid_duplicateCor.rds` (and sample annotations if needed) so we can run the published archetype/plotting steps?  
> Thank you — [Name], IISc

---

## 8. Environment notes

- Always: `export KMP_DUPLICATE_LIB_OK=TRUE`  
- R packages: `Anderson_PML_Archetypes/rlib` (edgeR, limma, GSVA present; **ParetoTI not installed yet**)  
- Python: repo `.venv`  
- Scripts: `codes/01`–`03` = v1; `04`/`04b` = SRA (optional); `05`/`06` = paperlike (written for count matrices — adapt `05` for Path B continuous expr or skip to GSVA on existing resid)

---

## 9. Do / Don’t

**Do:** Path A + Path B; keep author modules; keep paper colors; QC vs author assignments.  
**Don’t:** Block the chat on SRA; claim pixel-identical Fig 3 without residuals or high distance r; commit FASTQs.

---

## 10. File map

```
Anderson_PML_Archetypes/
  codes/01–03     # v1 working pipeline + paper-colored plots
  codes/04b       # SRA quant (OPTIONAL / flaky)
  codes/05–06     # paperlike residuals/ParetoTI (adapt for Path B)
  data/author/    # modules + published assignments
  data/processed/ # v1 matrices + sample info
  data/sra/       # maps/index/partial downloads (ignorable for now)
  results/fig3/   # CURRENT figures
  docs/HANDOFF.md # this file
```
