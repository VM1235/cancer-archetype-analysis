#!/bin/bash
# All within-unit fits, 6 at a time. From Project_1/:  bash Mouse_SCC/Mouse_PDAC_GSE173958/codes/run_within.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONWARNINGS=ignore
D=Mouse_SCC/Mouse_PDAC_GSE173958
for K in 3 4; do for G in ks all; do for U in M2-Met M2-PT M2-Liver M1-PT M1-Met M1-Liver M1-SS M1 M2 invitro; do
  echo "$U $G $K"; done; done; done | xargs -P 6 -L 1 sh -c \
 '.venv/bin/python -u '$D'/codes/08_within_unit.py --unit $0 --geneset $1 --k $2 > '$D'/results/within_unit/log_$1_$0_k$2.txt 2>&1'
echo done > $D/results/within_unit/DONE
