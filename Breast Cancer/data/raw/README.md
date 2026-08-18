# Raw DepMap downloads

Expected files:

1. `Model.csv` — model metadata (in git).  
2. `OmicsExpressionProteinCodingGenesTPMLogp1.csv` — protein-coding log2(TPM+1). **Gitignored** (~500 MB).

Download from DepMap / [download.depmap.org](https://depmap.org). Do not transform; the table is already log2(TPM+1).

Then from the repo root:

```bash
.venv/bin/python -u "Breast Cancer/codes/01_build_input_panelA.py"
```
