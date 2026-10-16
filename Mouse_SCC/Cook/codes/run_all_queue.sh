#!/bin/bash
# Runs the remaining Cook fits one after another (each uses ~6 processes). From Project_1/.
export PYTHONWARNINGS=ignore
for CFG in "ks pooled" "all A549" "ks A549" "all DU145" "ks DU145" "all MCF7" "ks MCF7" "all OVCA420" "ks OVCA420"; do
  set -- $CFG
  nice bash Mouse_SCC/Cook/codes/run_fits.sh $1 $2 200 > Mouse_SCC/Cook/results/queue_$1_$2.log 2>&1
done
echo done > Mouse_SCC/Cook/results/queue_done.flag
