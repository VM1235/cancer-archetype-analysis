# Cleaning-rule audit: handoff note

Started 2026-10-17. Read this first if you are a new session continuing the work.

## Why this exists
The five-dataset EMT archetype report (`docs/emt_single_cell/EMT_archetype_analysis_FINAL_REPORT.md`) says "cleaning" (>= 10,000 UMI, technical genes dropped) changes results: it adds within-sample triangles in PDAC and GSE114687 and removes them in 4T1. The user (Vedashree) emailed her supervisor (Prof. Mohit Jolly) on 2026-10-17 listing next steps: Hallmark check of corner names across datasets, cleaning-threshold sensitivity, then cross-dataset interpretation. She chose **option 2**: keep dataset-appropriate cleaning rules, audit the gaps, fix only what is needed (no full rerun with one common rule).

## What the audit found (all checked in code and data)
1. **Technical-gene filter did nothing in the two human datasets (GSE114687, Cook).** The list is mouse symbols (`Malat1`...); human data uses `MALAT1`. Confirmed: old GSE114687 `--clean` printed `technical genes dropped: []`. Cook's cleaning script never applied it at all. The report's "technical genes dropped" is therefore wrong for these two.
2. **10,000 UMI keeps very different fractions:** Cook 89%, PDAC 80%, SCC 71% (multiome 50%), 4T1 69% (clone E 21%), GSE114687 45% (HuMEC 35-40%).
3. SCC has no cleaning arm (its issue is the assay, not UMI). Pooled cleaned fits exist only for 4T1 and PDAC.
4. 4T1 contamination cutoff 0.3 was set by eye (no sensitivity check). PDAC cleaned KS fits were never run. Cook alone uses mito <= 0.07.

## What was changed (code only; old results untouched)
All new options default to the old behaviour; old settings reproduce the old observed t-ratio exactly (checked on Cook OVCA420, PDAC M1-Met, 4T1 EM2, GSE114687 MCF10A_Mock). Non-default settings write to NEW folders.
- `Breast Cancer/GSE114687_pseudospace/codes/08_within_unit.py`: `--drop-human-tech` (MALAT1 NEAT1 HBB HBA1 HBA2), output suffix `_ht`. Already had `--min-umi` and `--out-suffix`.
- `Mouse_SCC/Cook/codes/07_within_line_clean.py`: `--drop-human-tech`, suffix `_ht`.
- `Breast Cancer/Mouse_4T1_GSE303505/codes/08_within_unit.py`: `--min-umi`, `--dc-cut` (suffix `_u<value>`, `_c<value>` when non-default).
- `Mouse_SCC/Mouse_PDAC_GSE173958/codes/08_within_unit.py`: `--min-umi` (suffix `_u<value>`).
I deliberately did NOT drop SAA1/SAA2/SERPINA1 in the human list: they are liver ambient genes in PDAC but real expression in cultured breast cells.

## The queue (104 fits, 4 at a time, resumable)
- `gen_jobs.py` writes `jobs.tsv` (priority order). `run_queue.py` runs it; a fit whose result json exists is skipped, so rerunning is safe. `run_job.sh` runs one fit and logs to `logs/<id>.log`; `logs/queue.log` has START/DONE/FAIL lines.
- `progress.py` writes `PROGRESS.md`. `compare.py` writes `COMPARE.md` (old vs new p-values, with a pass/fail column).
- Groups, in priority order: (1) Cook + GSE114687 reruns with the human fix (could change the emailed "18 of 24"); (2) 4T1 cutoff 0.2/0.5; (3) PDAC cleaned KS fits; (4) GSE114687 UMI sweep 5,000 and 15,000.
- Start/resume: `cd /Users/apple/Desktop/thesis_iisc/Project_1 && nohup python3 docs/cleaning_audit/run_queue.py > docs/cleaning_audit/logs/queue_stdout.log 2>&1 &`

## STATUS (2026-10-17 07:00): all 104 fits finished, no failures. Results are in SUMMARY.md and COMPARE.md.
Hallmark for SCC and 4T1 done (hallmark_scc_4t1.csv). Remaining: cross-dataset interpretation (robust vs cleaning-dependent table), correcting the report wording, optional path test on SCC/PDAC/4T1.

## (Original plan, now mostly done)
- Run `compare.py`; tell the user plainly whether GSE114687's emailed "18 of 24 pass" (11 of 12 all-genes, 7 of 12 KS) changes, and whether Cook A549 still passes.
- Hallmark corner check for 4T1 (script exists only in `Breast Cancer/GSE114687_pseudospace/codes/09_hallmark_corners.py`; SCC has `Mouse_SCC/hallmark_check/`, check whether it already covers SCC).
- Then the comparison across datasets: sort each finding as robust (holds with light QC, cleaning, and across thresholds) or cleaning-dependent; interpret from robust ones only.
- Correct the report wording ("technical genes dropped") if the fix matters.

## Notes
- Run from `Project_1/` with `.venv/bin/python -u` (do not use `-I`: it breaks the local `scc_common` import).
- 16 GB machine: no more than 4 fits at once.
- Nothing is committed to git; the user decides when.
- The user prefers plain-language explanations and dislikes predictions-first on new datasets.
