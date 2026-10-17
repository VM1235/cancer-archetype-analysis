# Cleaning-audit progress
_Updated 2026-10-17 07:00_

**104 of 104 fits finished.**

| Group | Done | Total |
|---|---|---|
| Cook (human tech genes) | 18 | 18 |
| GSE114687 (human tech genes) | 24 | 24 |
| 4T1 contamination cutoff 0.2/0.5 | 20 | 20 |
| PDAC cleaned KS fits | 18 | 18 |
| GSE114687 UMI sweep | 24 | 24 |

Failed jobs: 0

## To resume (safe to rerun; finished fits are skipped)
```
cd /Users/apple/Desktop/thesis_iisc/Project_1
nohup python3 docs/cleaning_audit/run_queue.py > docs/cleaning_audit/logs/queue_stdout.log 2>&1 &
```
Check it is running: `ps aux | grep 08_within_unit | grep -v grep | wc -l` (and `tail docs/cleaning_audit/logs/queue.log`).

Then: `python3 docs/cleaning_audit/compare.py` writes COMPARE.md (old vs new p-values).
