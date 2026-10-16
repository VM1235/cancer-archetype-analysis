#!/bin/bash
# Within-unit fits again, after removing low-UMI cells and technical genes. All genes only. From Project_1/.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONWARNINGS=ignore
D=Mouse_SCC/Mouse_PDAC_GSE173958
for K in 3 4; do for U in M2-Met M2-PT M2-Liver M1-PT M1-Met M1-Liver M1-SS M1 M2; do echo "$U all $K"; done; done | xargs -P 6 -L 1 sh -c \
 '.venv/bin/python -u '$D'/codes/08_within_unit.py --clean --unit $0 --geneset $1 --k $2 > '$D'/results/within_unit/log_clean_$1_$0_k$2.txt 2>&1'
echo done > $D/results/within_unit/DONE_clean
