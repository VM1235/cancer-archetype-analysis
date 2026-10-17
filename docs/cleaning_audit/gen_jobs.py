#!/usr/bin/env python3
"""Write jobs.tsv: one line per fit, in priority order. Columns: id <TAB> expected_result_json <TAB> command.
Run from Project_1/:  python3 docs/cleaning_audit/gen_jobs.py
Priority: (1) Cook + GSE114687 reruns with the human technical-gene fix (could change emailed numbers),
          (2) 4T1 contamination-cutoff sweep, PDAC cleaned KS fits, (3) GSE114687 UMI sweep.
"""
PY = ".venv/bin/python -u"
COOK = "Mouse_SCC/Cook"; GSE = "Breast Cancer/GSE114687_pseudospace"; T4 = "Breast Cancer/Mouse_4T1_GSE303505"; PD = "Mouse_SCC/Mouse_PDAC_GSE173958"
jobs = []
def add(i, js, cmd): jobs.append((i, js, cmd))

# 1a. Cook, cleaned per-line fits with human technical genes dropped
for gs in ("all", "ks"):
    for line in ("A549", "DU145", "MCF7", "OVCA420"):
        for k in (3, 4):
            add(f"cook_{gs}_{line}_k{k}_ht", f"{COOK}/results/within_line/{gs}_{line}_clean_ht/fit_k{k}.json",
                f'{PY} {COOK}/codes/07_within_line_clean.py --line {line} --geneset {gs} --k {k} --variant clean --drop-human-tech')
for gs in ("all", "ks"):                       # A549's result depended on time balance
    add(f"cook_{gs}_A549_k3_tbal_ht", f"{COOK}/results/within_line/{gs}_A549_clean_tbal_ht/fit_k3.json",
        f'{PY} {COOK}/codes/07_within_line_clean.py --line A549 --geneset {gs} --k 3 --variant clean_tbal --drop-human-tech')
# 1b. GSE114687, cleaned within-sample fits with human technical genes dropped
UNITS = ("MCF10A_Mock", "HuMEC_Mock", "MCF10A_TGFB", "HuMEC_TGFB", "MCF10A", "HuMEC")
for k in (3, 4):
    for gs in ("all", "ks"):
        for u in UNITS:
            add(f"gse_{gs}_{u}_k{k}_ht", f"{GSE}/results/within_unit/{gs}_{u}_clean_ht/fit_k{k}.json",
                f'{PY} "{GSE}/codes/08_within_unit.py" --unit {u} --geneset {gs} --k {k} --clean --drop-human-tech')
# 2a. 4T1 contamination cutoff 0.2 and 0.5 (0.3 already done), k=3, within clone
for cut in (0.2, 0.5):
    for gs in ("all", "ks"):
        for u in ("EM1", "EM2", "EM3", "M", "4T1"):
            add(f"t4_{gs}_{u}_k3_c{cut}", f"{T4}/results/within_unit/{gs}_{u}_clean_dc_c{cut}/fit_k3.json",
                f'{PY} "{T4}/codes/08_within_unit.py" --unit {u} --geneset {gs} --k 3 --clean --decontam --dc-cut {cut}')
# 2b. PDAC cleaned KS fits (never run), k=3,4
for k in (3, 4):
    for u in ("M2-Met", "M2-PT", "M2-Liver", "M1-PT", "M1-Met", "M1-Liver", "M1-SS", "M1", "M2"):
        add(f"pdac_ks_{u}_k{k}", f"{PD}/results/within_unit/ks_{u}_clean/fit_k{k}.json",
            f'{PY} {PD}/codes/08_within_unit.py --clean --unit {u} --geneset ks --k {k}')
# 3. GSE114687 UMI sweep at 5,000 and 15,000 (human fix on), k=3
for t in (5000, 15000):
    for gs in ("all", "ks"):
        for u in UNITS:
            add(f"gse_{gs}_{u}_k3_t{t}", f"{GSE}/results/within_unit/{gs}_{u}_clean_ht_t{t}/fit_k3.json",
                f'{PY} "{GSE}/codes/08_within_unit.py" --unit {u} --geneset {gs} --k 3 --clean --drop-human-tech --min-umi {t} --out-suffix _t{t}')
with open("docs/cleaning_audit/jobs.tsv", "w") as f:
    for i, js, cmd in jobs: f.write(f"{i}\t{js}\t{cmd}\n")
print(len(jobs), "jobs")
