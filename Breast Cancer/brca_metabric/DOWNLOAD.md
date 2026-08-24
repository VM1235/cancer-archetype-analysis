# METABRIC-BRCA (cBioPortal)

Raw study files for `prepare_metabric_brca.py` are **not in git** (~1.3 GB).

## Download

1. Open [METABRIC in cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric).
2. Use **Download** → **Study** (or sync via the cBioPortal Datahub).
3. Extract so this folder contains at least:

| File | Role |
|------|------|
| `data_mrna_illumina_microarray.txt` | log2 HT-12 expression (~657 MB) |
| `data_clinical_sample.txt` | Sample-level clinical |
| `data_clinical_patient.txt` | Patient-level (`CLAUDIN_SUBTYPE`) |
| `meta_*.txt` | cBioPortal metadata (optional) |

Expected layout (same as cBioPortal export):

```
Breast Cancer/brca_metabric/
  data_mrna_illumina_microarray.txt
  data_clinical_sample.txt
  data_clinical_patient.txt
  meta_mrna_illumina_microarray.txt
  ...
```

## After download

From the **repository root**:

```bash
.venv/bin/python -u "Breast Cancer/codes/prepare_metabric_brca.py"
```

Writes prepared tables to `Breast Cancer/results/panel_c_metabric_ks_genelist/` (gene-matched METABRIC matrix + clinical metadata). See `Breast Cancer/docs/KS_METABRIC_PIPELINE.md` for the full KS → METABRIC workflow.
