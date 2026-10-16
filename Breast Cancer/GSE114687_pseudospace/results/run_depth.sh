#!/bin/bash
# usage: run_depth.sh variant unit geneset   (variant rd = regress depth on cleaned cells; dm = depth-matched matrix)
cd /Users/apple/Desktop/thesis_iisc/Project_1
S="Breast Cancer/GSE114687_pseudospace/codes"; R="Breast Cancer/GSE114687_pseudospace/results"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
if [ "$1" = "rd" ]; then
  .venv/bin/python -u "$S/08_within_unit.py" --unit $2 --geneset $3 --k 3 --n-perm 200 --clean --regress-depth > "$R/log_08rd_$2_$3_k3.txt" 2>&1
else
  EXPR_NAME=expr_tumour_depthmatched.h5ad .venv/bin/python -u "$S/08_within_unit.py" --unit $2 --geneset $3 --k 3 --n-perm 200 --out-suffix _dm > "$R/log_08dm_$2_$3_k3.txt" 2>&1
fi
