## 9. Cleaning audit and reruns (17 October 2026)

### 9.1 Why this section exists
The results in Sections 3 to 7 depend on a cleaning step (Section 8.4): it adds within-sample triangles in PDAC and GSE114687 and removes them in 4T1. A reader can fairly ask whether the five datasets were cleaned in a comparable way, and whether the step did what the report says. This section answers that. It was done in three parts: **read the cleaning code of every dataset** (what is actually applied), **measure what the rules keep**, and **rerun the fits the audit says could change**. The decision taken was not to force one cleaning rule on all five datasets (they are very different systems) but to keep each dataset's own rule, fix what is wrong, and test how much the findings depend on the two free choices (the UMI cutoff and the 4T1 contamination cutoff).

### 9.2 What each dataset's cleaning actually does (read from the code)

| | Basic cell QC (before any fit) | Cleaning arm used for within-sample fits | Contamination handling | Technical genes dropped |
|---|---|---|---|---|
| **SCC** | authors' own (more than 700 genes, under 20% mitochondrial) | none (its problem is the assay, not depth) | none | no |
| **Cook** | none added here | at least 10,000 UMI **and** mitochondrial fraction at most 0.07 | none (cell lines) | **no** (the docstring says yes; the code never applies it) |
| **PDAC** | 1,000 to 7,000 genes, mitochondrial at most 15% | at least 10,000 UMI | tumour-cell cluster selection (Section 5.1) | yes, mouse list |
| **4T1** | at least 500 genes, mitochondrial at most 15% | at least 10,000 UMI | immune and stromal panel score at most 0.3 (cut set by eye) | yes, mouse list |
| **GSE114687** | at least 500 genes, mitochondrial at most 15% | at least 10,000 UMI | none (cultured cells) | **no** (mouse list, matches nothing) |

### 9.3 What the audit found

**Finding 1: the technical-gene filter did nothing in the two human datasets.** The list of technical genes (Malat1, Neat1, haemoglobin, liver ambient-RNA genes, Gm and AY transcripts) was written with mouse gene symbols and matched by exact spelling. Human data uses upper case (`MALAT1`). In the old GSE114687 cleaned fits the filter printed `technical genes dropped: []`, and `MALAT1` was still in the data; Cook's cleaning script never applied the filter. The statement "technical genes dropped" in Section 2.2 was therefore wrong for Cook and GSE114687.

**Finding 2: the same 10,000-UMI cutoff keeps very different shares of cells** (measured on the five processed expression files).

![Figure A1. Share of cells kept by 10,000 UMI](../cleaning_audit/figures/figA1_umi_keep_fraction.png)

**Figure A1 explained.** Blue bars: the share of all cells in each dataset with at least 10,000 UMI (Cook 89%, PDAC 80%, SCC 71%, 4T1 69%, GSE114687 45%). Orange dots: the lowest single sample (SCC tumour MU8 30%, PDAC library iv-100cell 52%, 4T1 clone E 21%, GSE114687 HuMEC Mock 35%, Cook OVCA420 83%). "Cleaned" therefore means "most cells" in Cook and "under half" in GSE114687, and a whole clone (4T1 E, 142 cells left) can drop out of the analysis. In SCC the multiome assay keeps 50% and the scRNA assay 86%.

**Finding 3: SCC has no cleaning arm, and pooled cleaned fits exist only for 4T1 and PDAC** (as far as the scripts and result folders show). The Cook, SCC and GSE114687 pooled fits appear to use the lighter QC. This limits how directly pooled results compare across datasets (Section 9.8).

**Finding 4: two free choices were never tested.** The 4T1 contamination cutoff (0.3) was set by eye, and the 10,000-UMI cutoff was borrowed from Cook and PDAC (no scan). PDAC cleaned fits on the KS genes were never run.

### 9.4 Code changes (code only; no old result touched)
All new options default to the old behaviour. The old settings were checked to reproduce the old observed t-ratio exactly on one small unit of each dataset (Cook OVCA420 0.4446, PDAC M1-Met 0.4586, 4T1 EM2 0.4721, GSE114687 MCF10A Mock 0.6310) before anything was queued. Non-default settings write to **new** output folders (suffixes below), so old and new can be compared side by side.

| Script | New option | Output suffix |
|---|---|---|
| `Breast Cancer/GSE114687_pseudospace/codes/08_within_unit.py` | `--drop-human-tech` (MALAT1, NEAT1, HBB, HBA1, HBA2) | `_ht` (UMI sweep adds `_t5000`, `_t15000` through the existing `--min-umi` and `--out-suffix`) |
| `Mouse_SCC/Cook/codes/07_within_line_clean.py` | `--drop-human-tech` (same list) | `_ht` |
| `Breast Cancer/Mouse_4T1_GSE303505/codes/08_within_unit.py` | `--min-umi`, `--dc-cut` | `_u<value>`, `_c<value>` when not 10,000 and 0.3 |
| `Mouse_SCC/Mouse_PDAC_GSE173958/codes/08_within_unit.py` | `--min-umi` | `_u<value>` when not 10,000 |

The human list is deliberately short. The mouse list also contains liver ambient-RNA genes (Alb, Apoa1, Saa1, Serpina1...) that are contamination in PDAC but real expression in cultured breast cells; a first version dropped SAA1, SAA2 and SERPINA1 in the human data and was corrected before any fit was run.

### 9.5 Flow of the reruns
1. Edit and check the four scripts (Section 9.4).
2. Write `docs/cleaning_audit/gen_jobs.py`, which lists every fit to run, in priority order, as `jobs.tsv` (id, expected result file, command). 104 fits in five groups:
   - Cook, cleaned per-line fits with the human fix: 18 fits (4 lines, all genes and KS, k = 3 and 4, plus A549 time-balanced).
   - GSE114687, cleaned within-sample fits with the human fix: 24 fits (6 units, all genes and KS, k = 3 and 4).
   - 4T1 contamination cutoff 0.2 and 0.5: 20 fits (5 clones, all genes and KS, k = 3).
   - PDAC cleaned KS fits: 18 fits (9 units, k = 3 and 4).
   - GSE114687 UMI sweep at 5,000 and 15,000: 24 fits (6 units, all genes and KS, k = 3).
3. Run them with `docs/cleaning_audit/run_queue.py`, four fits at a time (the machine has 16 GB; one fit uses about 350 MB and one core). A fit whose result file already exists is skipped, so the queue can be stopped and restarted without losing work. Each fit writes its own log; `logs/queue.log` records start, done and fail lines.
4. Compare old and new with `compare.py` (writes `COMPARE.md`), track with `progress.py` (writes `PROGRESS.md`), draw the figures with `make_figures.py`.
5. Score the Hallmark sets at corners with `hallmark_scc_4t1.py` (SCC and 4T1, the two datasets that never had this check).

All 104 fits finished with no failures (started 01:58, finished 07:00 on 17 October; one fit takes 1.5 to 25 minutes, median 12). The fits themselves are unchanged: 200 shuffles, 15 PCHA starts for the real fit, 5 per shuffle, seed fixed.

### 9.6 Results

#### 9.6.1 The human technical-gene fix changes nothing

![Figure A2. p before versus after the human technical-gene fix](../cleaning_audit/figures/figA2_human_gene_fix.png)

**Figure A2 explained.** Each dot is one fit (orange Cook, 18 fits; blue GSE114687, 24 fits). Horizontal axis: the permutation p with the old filter (which matched no gene); vertical axis: the p after dropping MALAT1, NEAT1 and haemoglobin genes. The diagonal is "no change"; the dotted lines mark p = 0.05. Every dot sits on the diagonal (a few are jittered apart for visibility); none crosses a dotted line.
- **Cook:** all 18 fits keep their pass or fail status. A549 on all genes at k = 3 still passes (p = 0.010 before and after) and still fails when time-balanced (0.224 to 0.269). MCF7 on KS genes at k = 3 still passes (0.030). KS fits are identical to three decimals, as expected, since MALAT1 and NEAT1 are not KS genes.
- **GSE114687:** 18 of 24 cleaned fits pass before and after (all genes 11 of 12, KS genes 7 of 12). The largest change is one fit (HuMEC both arms, all genes, k = 3: 0.925 to 0.886), a fail either way.
- Conclusion: the figure "18 of 24" in Section 7.5 stands. Only the wording "technical genes dropped" needed correcting.

#### 9.6.2 The 4T1 contamination cutoff (0.2, 0.3, 0.5)

![Figure A4. 4T1 within-clone p by contamination cutoff](../cleaning_audit/figures/figA4_4t1_cutoff.png)

**Figure A4 explained.** Rows are within-clone fits at k = 3 (all genes, then KS genes; the five clones with enough cells after cleaning), columns are the cutoff on the immune and stromal panel score (a stricter 0.2, the 0.3 used in Section 6, a looser 0.5). Numbers are permutation p; darker is smaller; a star marks p below 0.05. Looser cutoffs keep more cells (for example EM3: 1,208, 1,507 and 1,582 cells at 0.2, 0.3 and 0.5).
- **All genes:** no clone passes at any cutoff (0 of 5 at each). "Within a single clone there is almost nothing" does not depend on the 0.3 chosen by eye.
- **KS genes:** clone EM2 passes at all three cutoffs (0.035, 0.015 and 0.050; the last is borderline, p = 0.0498). The parental 4T1 clone passes only at 0.5 (p = 0.005, against 0.100 at 0.3 and 0.2). A looser cutoff leaves more immune and stromal cells, which can form a corner of their own; the corner genes of this fit were not checked, so this one pass is unexplained.

#### 9.6.3 PDAC cleaned KS-gene fits (never run before)

| Unit | KS k = 3 | KS k = 4 | All genes k = 3 (Section 5.4) |
|---|---|---|---|
| M1-SS | 0.333 | **0.005** | **0.005** |
| M1-PT | 0.990 | 0.154 | **0.020** |
| M1-Liver | 0.995 | 0.970 | **0.020** |
| M1-Met | 0.627 | 0.597 | 0.517 |
| M2-PT | 0.622 | 0.910 | **0.005** |
| M2-Liver | 0.856 | 0.746 | **0.005** |
| M2-Met | 0.124 | **0.040** | **0.005** |
| Mouse 1 (all libraries) | 1.000 | 0.985 | 0.124 |
| Mouse 2 (all libraries) | 0.866 | **0.025** | **0.005** |

On the EMT genes alone **no unit passes at k = 3** (0 of 9, against 7 of 9 on all genes). Three of nine pass at k = 4 (M1-SS, M2-Met, mouse 2); k = 4 is beyond the smallest k that passes in the all-genes fits, so these are weaker evidence. The within-library triangles in PDAC therefore come from the other programmes (proliferation, hypoxia, epithelial identity), not from the EMT genes, which agrees with Section 5.5.

#### 9.6.4 GSE114687 UMI cutoff (5,000, 10,000, 15,000)

![Figure A3. GSE114687 within-sample p by UMI cutoff](../cleaning_audit/figures/figA3_gse_umi_sweep.png)

**Figure A3 explained.** Rows are the six within-sample units for all genes and KS genes at k = 3 (human technical genes dropped); columns the UMI cutoff. Numbers are permutation p, a star marks p below 0.05. Cells kept differ a lot by column: at 5,000 every unit is capped at the 4,000-cell fit limit (a random subsample); at 10,000 the units keep 2,536 to 4,000 cells; at 15,000 only 284 to 1,838.

| Passing fits out of 6 | 5,000 | 10,000 | 15,000 |
|---|---|---|---|
| All genes | 5 | 5 | 6 |
| KS genes | 1 | 3 | 4 |

- **All genes:** the three MCF10A units (Mock, TGF-β, both arms) pass at every cutoff with p = 0.005. HuMEC is the fragile cell line: HuMEC TGF-β fails at 5,000 (p = 1.000), only just passes at 10,000 (0.025) and passes at 15,000; HuMEC (both arms) fails at 10,000 (0.886) but passes at 5,000 and 15,000.
- **KS genes:** a KS-only triangle exists in MCF10A Mock, HuMEC Mock and HuMEC TGF-β at 10,000 and 15,000 but **not at 5,000**. So the claim that the EMT genes alone give a triangle needs cells with at least about 10,000 UMI. This agrees with the depth control in Section 7.8, which tested the 10,000-cell version directly, and means the claim is tied to a depth condition.
- At 15,000 so few cells remain that a pass is a smaller sample, not automatically stronger evidence.

#### 9.6.5 Hallmark check of corner names for SCC and 4T1 (the two datasets without one)
Method as in Sections 5 and 7 (code `docs/cleaning_audit/hallmark_scc_4t1.py`): the 10% of cells nearest each corner are scored on MSigDB Hallmark sets (mean z of the set's genes; score = nearest 10% minus the rest; mouse symbols are the capitalised human symbols); the highest set is reported.

![Figure A5. Hallmark scores at SCC and 4T1 corners](../cleaning_audit/figures/figA5_hallmark_scc_4t1.png)

**Figure A5 explained.** Rows are corners (with the share of cycling cells, S or G2M, among the cells nearest the corner); columns the seven Hallmark sets; red means the set is higher near the corner, blue lower.
- **4T1:** in the pooled cleaned fit and every within-clone fit the cycling corner (94 to 97% cycling cells; 64 to 80% in the weaker EM1 and EM3 fits) is top on E2F targets (+0.43 to +0.81) and the G2-M checkpoint. Clone M's corner in the pooled fit scores highest on TGF-β (+0.25) and hypoxia (+0.15). The epithelial-leaning corner has no strong Hallmark set (EMT +0.05). **No 4T1 corner is named by an EMT set.** This confirms the Section 6 reading: the corners are cell division and clone identity.
- **SCC:** the proliferating corner is top on E2F and G2-M (E2F +0.42, +0.47 and +0.58 in the three fits), the late-EMT corner top on EMT (+0.36 to +0.43), as named. The corner called technical has no set above +0.25 (TGF-β +0.23 to +0.25, a weak match), which is consistent with an axis that is not a biological programme.
- The Hallmark check now covers all five datasets (PDAC, Cook and GSE114687 were done earlier).

### 9.7 What holds and what depends on cleaning

| Finding | Holds under the cleaning variants tried? |
|---|---|
| Pooled over samples, EMT genes give no simplex beyond sample identity | Yes (unchanged by every rerun) |
| Corners that appear are programmes, not EMT states | Yes (Hallmark now confirms for all five datasets) |
| Hybrid cells are mostly generalists | Not re-tested here; checked under light QC and cleaned in 4T1 (72 to 85% and 83 to 87%) |
| 4T1: almost nothing within a single clone | Yes (0 of 5 on all genes at every cutoff; one unexplained KS pass at 0.5) |
| PDAC: within-library triangles are not made by EMT genes | Yes (KS k = 3: 0 of 9) |
| GSE114687: a triangle inside a sample after cleaning | Yes for MCF10A on all genes at 5,000 to 15,000; HuMEC depends on the cutoff |
| GSE114687: the EMT genes alone give a triangle | Only at 10,000 UMI or more |
| Cook: A549 triangle | Depends on cleaning and on time balance (unchanged) |

### 9.8 Limits of the audit and what was not done
- The datasets still differ in how they are cleaned (Section 9.2): SCC has no cleaning arm, pooled cleaned fits appear to exist only for 4T1 and PDAC, and the contamination step differs by tissue. These differences are real properties of the datasets and were left as they are; results are compared as patterns that hold in each dataset, not as numbers across datasets.
- The UMI sweep was run for GSE114687 only. The 4T1 and PDAC UMI cutoffs, and the SCC cleaned arm, were not scanned.
- The path-versus-triangle test was not run for SCC, PDAC and the new fits; the planted-triangle power check was not run.
- Shuffle p-values floor at 0.005 (200 shuffles); a "star" on a 15,000-UMI fit uses far fewer cells.
- The one unexplained pass (4T1 parental, KS genes, cutoff 0.5) was not followed up.

### 9.9 Files and reproduction
Everything is in `docs/cleaning_audit/`: `SUMMARY.md` (short results), `COMPARE.md` (old versus new p for every fit), `HANDOFF.md` (how to continue), `PROGRESS.md`, `jobs.tsv`, `gen_jobs.py`, `run_queue.py`, `run_job.sh`, `progress.py`, `compare.py`, `make_figures.py`, `hallmark_scc_4t1.py`, `hallmark_scc_4t1.csv`, `figures/` and `logs/`. To rerun (from `Project_1/`, with `.venv`; do not use `python -I`, which breaks the local `scc_common` import):
```
nohup python3 docs/cleaning_audit/run_queue.py > docs/cleaning_audit/logs/queue_stdout.log 2>&1 &   # resumable, 4 at a time
python3 docs/cleaning_audit/compare.py                 # COMPARE.md
.venv/bin/python docs/cleaning_audit/hallmark_scc_4t1.py
.venv/bin/python docs/cleaning_audit/make_figures.py
```
New result folders sit next to the old ones in each dataset's `results/within_unit/` (or `within_line/`) with the suffixes in Section 9.4.

---

