# Processed cell-line matrix

Written by `codes/01_build_input_panelA.py`.

| File | Role |
|---|---|
| `input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv` | Genes × 63 invasive breast carcinoma cell lines |
| `input_panelA_models_used.csv` | Metadata rows that went into that matrix |
| `input_panelA_entrez_mapped.csv` | Same expression, Entrez IDs (PAM50 / genefu) |
| `input_panelA_build_report.txt` | Filter counts from the last build |

Codes load the first CSV as Panel A/B input. Do not log-transform again.
