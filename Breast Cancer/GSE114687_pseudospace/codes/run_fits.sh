#!/bin/bash
# Step 3 for one gene set: ESV curve + each k as its own process, then a summary.
# Usage (from Project_1/):  bash "Breast Cancer/GSE114687_pseudospace"/codes/run_fits.sh all in_vivo 200
GENESET=${1:-all}; TAG=${2:-in_vivo}; NPERM=${3:-200}
PY=.venv/bin/python
D="Breast Cancer/GSE114687_pseudospace"
S="$D/codes/03_fit_and_test.py"
OUT="$D/results/${GENESET}_${TAG}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
$PY -u "$S" --geneset $GENESET --tag $TAG --esv > "$OUT"/log_03_esv.txt 2>&1 &
for K in 3 4 5 6 7; do
  $PY -u "$S" --geneset $GENESET --tag $TAG --k $K --n-perm $NPERM > "$OUT"/log_03_k$K.txt 2>&1 &
done
wait
$PY -u "$S" --geneset $GENESET --tag $TAG --summarize | tee "$OUT"/log_03_summary.txt
