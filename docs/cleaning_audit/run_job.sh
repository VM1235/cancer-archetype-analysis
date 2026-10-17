#!/bin/bash
# run_job.sh <id> <expected_json> <command...>  -- skips if the result already exists (so the queue can be restarted safely)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONWARNINGS=ignore
ID="$1"; JSON="$2"; shift 2
L=docs/cleaning_audit/logs
if [ -f "$JSON" ]; then echo "$(date +%T) SKIP $ID" >> $L/queue.log; exit 0; fi
echo "$(date +%T) START $ID" >> $L/queue.log
S=$(date +%s)
if eval "$@" > "$L/$ID.log" 2>&1 && [ -f "$JSON" ]; then
  echo "$(date +%T) DONE  $ID ($(( $(date +%s) - S ))s)" >> $L/queue.log
else
  echo "$(date +%T) FAIL  $ID (see logs/$ID.log)" >> $L/queue.log
fi
