# SCLC Panel B/C fix (18 Aug 2026)

## FLAG — containment does not support “mostly contained”

Under the ParTI-matched k=5 simplex, **4/81 tumors (4.9%)** have all barycentric weights ≥ 0 in the combined-PCA 4-D fitting space.

That is **not** “mostly contained.” The old claim was never a computed fraction; it was a PC1–PC2 scatter impression. The 2-D picture can look filled while points sit outside the 4-simplex.

**Panel B 5/5 subtype↔vertex match still holds** (same folds; vertex *indices* permute). The scientific headline that changes is Panel C containment.

---

## What changed

Scoped to SCLC `codes/run_panel_b.py`, `codes/run_panel_c.py`, then `export_sclc_reproduction_figures.py` for **Figure 1B and 1C only**. Figure 1A bytes were restored after the exporter rewrote all three PNGs. Panel A scripts and breast/GBM were not modified.

| Item | Before | After |
|---|---|---|
| Archetype file (B and C) | `results/panel_a/archetypes_k5.npy` **(5, 12)** first-pass, δ=0.1 | `results/panel_a/archetypes_k5_parti.npy` **(5, 4)** ParTI, δ=0 |
| Panel B distances | all 12 PCs | first **k−1 = 4** PCs (`scores[:, :4]`), same pattern as breast/GBM |
| Panel C gene-space inverse | 12-D scores into 12-PC PCA | pad (5,4) → (5,12) with zeros on PC5–12, then inverse |
| Containment | not computed | barycentric lstsq in combined PCA first 4 PCs; `w ≥ -1e-6` |

Snapshots: `Aug 18/before/` (pre-fix outputs + old 1B/1C) and `Aug 18/after/` (new outputs + new 1B/1C). Live results remain in `results/panel_b/`, `results/panel_c/`, `figures/`.

---

## Panel B — before vs after (author `NEW_10_2020` labels)

Equal-count **10** bins, hypergeometric + BH FDR 0.1, hit = significant **and** peak at bin 0. Unchanged machinery; only distances changed.

| Subtype | Before: arc, fold, q | After: arc, fold, q | Bin-0 peak? | 1-to-1? |
|---|---|---|---|---|
| SCLC-N | Arc **1**, 4.800, 4.11e-08 | Arc **5**, 4.800, 4.11e-08 | yes / yes | yes / yes |
| SCLC-A2 | Arc **2**, 4.286, 1.80e-07 | Arc **3**, 4.286, 1.80e-07 | yes / yes | yes / yes |
| SCLC-A | Arc **3**, 2.500, 1.84e-04 | Arc **2**, 2.500, 2.07e-04 | yes / yes | yes / yes |
| SCLC-Y | Arc **4**, 10.000, 2.63e-09 | Arc **4**, 10.000, 2.63e-09 | yes / yes | yes / yes |
| SCLC-P | Arc **5**, 10.000, 1.42e-10 | Arc **1**, 10.000, 1.42e-10 | yes / yes | yes / yes |

Folds and raw p-values for the five hits are the same; BH q for A moved slightly (table-wide correction). Vertex numbering is a permutation except Y (still Arc 4). **No subtype lost its bin-0 match; no two subtypes share a vertex.**

**Verdict: PASS** on the claim “clean 5/5 subtype-to-archetype match.” Arc labels in any write-up that said N=1, A2=2, A=3, Y=4, P=5 must be updated to **P=1, A=2, A2=3, Y=4, N=5**.

---

## Panel C — containment and variance

| Quantity | Before | After |
|---|---|---|
| Containment | visual “mostly contained”; **no number** | **4/81 (4.9%)** inside the 4-simplex |
| Tumors inside | — | `t.S00825`, `t.S00837`, `t.S01861`, `t.S02249` |
| 5 combined PCs / tumor-only ceiling | 80.1% | **80.1%** (unchanged; this statistic does not use archetype coordinates) |

CSV: `results/panel_c/containment_summary.csv`, per-tumor flags in `results/panel_c/tcga_sample_subtype_archetype.csv`.

**Verdict: CHANGED.** “Tumors mostly contained” **does not hold** as a 4-D simplex statement. PC1–PC2 can still look overlapping (Figure 1C). Groves’ qualitative Panel C story is not reproduced by the barycentric test used for breast (0/1097) and GBM (0/154). SCLC is 4/81 — better than those, not “mostly.”

---

## PASS / CHANGED

| Original claim | Verdict |
|---|---|
| 5/5 subtype match at bin 0 | **PASS** (vertex indices permuted) |
| Tumors mostly contained in the cell-line polytope | **CHANGED — FAIL as stated.** Computed containment **4/81**. |
