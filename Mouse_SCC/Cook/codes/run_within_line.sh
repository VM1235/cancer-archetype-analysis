#!/bin/bash
# Cook per-line fits after cleaning (and time balancing). From Project_1/.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONWARNINGS=ignore
D=Mouse_SCC/Cook
mkdir -p $D/results/within_line
for V in clean clean_tbal; do for K in 3 4; do for G in all ks; do for L in A549 DU145 MCF7 OVCA420; do echo "$L $G $K $V"; done; done; done; done | xargs -P 6 -L 1 sh -c \
 '.venv/bin/python -u '$D'/codes/07_within_line_clean.py --line $0 --geneset $1 --k $2 --variant $3 > '$D'/results/within_line/log_$1_$0_$3_k$2.txt 2>&1'
echo done > $D/results/within_line/DONE
