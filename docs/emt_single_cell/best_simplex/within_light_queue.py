#!/usr/bin/env python3
"""Hallmark EMT within-unit fits on the light-QC arm (the cleaned arm is run by the first Hallmark queue). 4 at a time, 200 shuffles, k = 3 and 4.
Run from Project_1/:  python3 docs/emt_single_cell/best_simplex/within_light_queue.py [file_to_wait_for]"""
import os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]; os.chdir(ROOT)
if len(sys.argv) > 1:
    while not Path(sys.argv[1]).exists(): time.sleep(60)
PY = ".venv/bin/python"
env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONWARNINGS="ignore")
LOG = "/private/tmp/claude-501/-Users-apple-Desktop-thesis-iisc-Project-1/4353b0d5-02a1-4798-8b62-cbbe4d13b418/scratchpad/hm_logs"
jobs = []
for u in ("M1-PT", "M1-SS", "M1-Liver", "M1-Met", "M2-PT", "M2-Liver", "M2-Met", "M1", "M2"):
    jobs += [f'{PY} -u Mouse_SCC/Mouse_PDAC_GSE173958/codes/08_within_unit.py --unit {u} --geneset hallmark --k {k} > {LOG}/L_pdac_{u}_{k}.log 2>&1' for k in (3, 4)]
for u in ("EM1", "EM2", "EM3", "M", "4T1"):   # clone E light-QC Hallmark is already in the first Hallmark queue
    jobs += [f'{PY} -u "Breast Cancer/Mouse_4T1_GSE303505/codes/08_within_unit.py" --unit {u} --geneset hallmark --k {k} > {LOG}/L_t4_{u}_{k}.log 2>&1' for k in (3, 4)]
for u in ("MCF10A", "MCF10A_Mock", "MCF10A_TGFB", "HuMEC", "HuMEC_Mock", "HuMEC_TGFB"):
    jobs += [f'{PY} -u "Breast Cancer/GSE114687_pseudospace/codes/08_within_unit.py" --unit {u} --geneset hallmark --k {k} > {LOG}/L_gse_{u}_{k}.log 2>&1' for k in (3, 4)]
for l in ("A549", "DU145", "MCF7", "OVCA420"):   # Cook light QC = per-line fit on the forward arm (02_pca --lines, then 03)
    c = "Mouse_SCC/Cook/codes"
    jobs.append(f'{PY} -u {c}/02_pca.py --geneset hallmark --lines {l} --tag {l} > {LOG}/L_cook_{l}_pca.log 2>&1 && {PY} -u {c}/03_fit_and_test.py --geneset hallmark --tag {l} --k 3 4 --n-perm 200 > {LOG}/L_cook_{l}_fit.log 2>&1')
with ThreadPoolExecutor(4) as ex: list(ex.map(lambda c: subprocess.run(c, shell=True, env=env).returncode, jobs))
print("LIGHT WITHIN DONE", flush=True)
