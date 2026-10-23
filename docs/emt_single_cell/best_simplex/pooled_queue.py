#!/usr/bin/env python3
"""Queue of the new pooled fits (Hallmark EMT genes; cleaned pooled arms; SCC KS scRNA-only / six-tumour arms).
Each run: PCA -> ESV curve -> PCHA + permutation test for k = 3, 4, 5 in parallel (100 shuffles, to fit the time budget; the
real fits elsewhere used 200). Resumable: a run whose fit_k5.json exists is skipped. Two runs at a time.
Run from Project_1/:  python3 docs/emt_single_cell/best_simplex/pooled_queue.py [file_to_wait_for]
"""
import os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]; os.chdir(ROOT)
if len(sys.argv) > 1:
    while not Path(sys.argv[1]).exists(): time.sleep(60)
PY = ".venv/bin/python"; NPERM = 100; KS = (3, 4, 5)
BASE = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONWARNINGS="ignore")
D = {"SCC": "Mouse_SCC", "COOK": "Mouse_SCC/Cook", "PDAC": "Mouse_SCC/Mouse_PDAC_GSE173958", "T4": "Breast Cancer/Mouse_4T1_GSE303505", "GSE": "Breast Cancer/GSE114687_pseudospace"}
# (dataset, geneset, tag, extra args for 02_pca.py, extra env)
RUNS = []
for gs in ("hallmark",):
    RUNS += [("SCC", gs, "main", [], {}), ("COOK", gs, "pooled", ["--lines", "all", "--tag", "pooled"], {}), ("PDAC", gs, "in_vivo", ["--assay", "in_vivo"], {}),
             ("T4", gs, "main", [], {}), ("T4", gs, "clean", ["--tag", "clean"], {"EXPR_NAME": "expr_tumour_clean.h5ad"}), ("GSE", gs, "main", [], {})]
for gs in ("all", "ks", "hallmark"):
    RUNS += [("COOK", gs, "pooled_clean", ["--lines", "all", "--tag", "pooled_clean"], {"EXPR_NAME": "expr_log_filtered_clean.h5ad"}),
             ("GSE", gs, "clean", ["--tag", "clean"], {"EXPR_NAME": "expr_tumour_clean.h5ad"}),
             ("PDAC", gs, "in_vivo_clean", ["--assay", "in_vivo", "--tag", "in_vivo_clean"], {"EXPR_NAME": "expr_tumour_clean.h5ad"})]
for gs in ("ks", "hallmark"):
    RUNS += [("SCC", gs, "scrna", ["--assay", "scRNA"], {}), ("SCC", gs, "scrna6", ["--assay", "scRNA", "--exclude-tumours", "SC15", "SC46", "--tag", "scrna6"], {})]
RES = {"SCC": "Mouse_SCC/results", "COOK": "Mouse_SCC/Cook/results", "PDAC": "Mouse_SCC/Mouse_PDAC_GSE173958/results", "T4": "Breast Cancer/Mouse_4T1_GSE303505/results", "GSE": "Breast Cancer/GSE114687_pseudospace/results"}

def run(r):
    ds, gs, tag, extra, env = r
    out = Path(RES[ds]) / f"{gs}_{tag}"
    if (out / "fit_k5.json").exists(): return
    e = dict(BASE, **env); code = f"{D[ds]}/codes"
    log = lambda n: open(f"/private/tmp/claude-501/-Users-apple-Desktop-thesis-iisc-Project-1/4353b0d5-02a1-4798-8b62-cbbe4d13b418/scratchpad/pooled_logs/pooled_{ds}_{gs}_{tag}_{n}.log", "w")
    subprocess.run([PY, "-u", f"{code}/02_pca.py", "--geneset", gs] + extra, env=e, stdout=log("pca"), stderr=subprocess.STDOUT)
    tg = [PY, "-u", f"{code}/03_fit_and_test.py", "--geneset", gs, "--tag", tag]
    ps = [subprocess.Popen(tg + ["--esv"], env=e, stdout=log("esv"), stderr=subprocess.STDOUT)]
    ps += [subprocess.Popen(tg + ["--k", str(k), "--n-perm", str(NPERM)], env=e, stdout=log(f"k{k}"), stderr=subprocess.STDOUT) for k in KS]
    for p in ps: p.wait()
    print(time.strftime("%H:%M"), "done", ds, gs, tag, flush=True)

with ThreadPoolExecutor(2) as ex: list(ex.map(run, RUNS))
print("ALL POOLED DONE", flush=True)
