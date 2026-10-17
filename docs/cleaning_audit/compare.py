#!/usr/bin/env python3
"""Write COMPARE.md: old versus new permutation p for every rerun fit. From Project_1/: python3 docs/cleaning_audit/compare.py"""
import os, json
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
COOK="Mouse_SCC/Cook/results/within_line"; GSE="Breast Cancer/GSE114687_pseudospace/results/within_unit"
T4="Breast Cancer/Mouse_4T1_GSE303505/results/within_unit"; PD="Mouse_SCC/Mouse_PDAC_GSE173958/results/within_unit"
def rd(p):
    try: d = json.load(open(p)); return d["p_value"], d["n_cells"], d["t_ratio"]
    except Exception: return None
def fmt(x): return "n/a" if x is None else f"{x[0]:.3f} ({x[1]} cells)"
sections = []
def sec(title, note, rows):
    out = [f"## {title}", note, "", "| Fit | Before | After | Passes (p<0.05) before -> after |", "|---|---|---|---|"]
    for name, a, b in rows:
        ra, rb = rd(a), rd(b)
        pa = "-" if ra is None else ("yes" if ra[0] < .05 else "no"); pb = "-" if rb is None else ("yes" if rb[0] < .05 else "no")
        out.append(f"| {name} | {fmt(ra)} | {fmt(rb)} | {pa} -> {pb} |")
    sections.append("\n".join(out))
rows = [(f"{gs} {l} k{k}", f"{COOK}/{gs}_{l}_clean/fit_k{k}.json", f"{COOK}/{gs}_{l}_clean_ht/fit_k{k}.json") for gs in ("all","ks") for l in ("A549","DU145","MCF7","OVCA420") for k in (3,4)]
rows += [(f"{gs} A549 time-balanced k3", f"{COOK}/{gs}_A549_clean_tbal/fit_k3.json", f"{COOK}/{gs}_A549_clean_tbal_ht/fit_k3.json") for gs in ("all","ks")]
sec("Cook: cleaned, before/after dropping human technical genes (MALAT1, NEAT1, haemoglobin)", "Before = old cleaned fit (filter matched nothing).", rows)
U=("MCF10A_Mock","HuMEC_Mock","MCF10A_TGFB","HuMEC_TGFB","MCF10A","HuMEC")
rows=[(f"{gs} {u} k{k}", f"{GSE}/{gs}_{u}_clean/fit_k{k}.json", f"{GSE}/{gs}_{u}_clean_ht/fit_k{k}.json") for k in (3,4) for gs in ("all","ks") for u in U]
sec("GSE114687: cleaned, before/after dropping human technical genes", "Before = old cleaned fit (filter matched nothing).", rows)
for t in (5000, 15000):
    rows=[(f"{gs} {u} k3", f"{GSE}/{gs}_{u}_clean_ht/fit_k3.json", f"{GSE}/{gs}_{u}_clean_ht_t{t}/fit_k3.json") for gs in ("all","ks") for u in U]
    sec(f"GSE114687: UMI cutoff 10,000 -> {t:,}", "Before = 10,000 with human technical genes dropped.", rows)
for c in (0.2, 0.5):
    rows=[(f"{gs} {u} k3", f"{T4}/{gs}_{u}_clean_dc/fit_k3.json", f"{T4}/{gs}_{u}_clean_dc_c{c}/fit_k3.json") for gs in ("all","ks") for u in ("EM1","EM2","EM3","M","4T1")]
    sec(f"4T1: contamination cutoff 0.3 -> {c}", "Within-clone, cleaned (>=10,000 UMI).", rows)
rows=[(f"ks {u} k{k}", "none", f"{PD}/ks_{u}_clean/fit_k{k}.json") for k in (3,4) for u in ("M2-Met","M2-PT","M2-Liver","M1-PT","M1-Met","M1-Liver","M1-SS","M1","M2")]
sec("PDAC: cleaned KS-gene fits (never run before)", "No earlier result to compare.", rows)
open("docs/cleaning_audit/COMPARE.md","w").write("# Cleaning audit: old vs new\n\n" + "\n\n".join(sections) + "\n")
print("wrote docs/cleaning_audit/COMPARE.md")
