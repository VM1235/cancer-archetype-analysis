# GSE161529 (Pal et al. 2021) — local data

Was `Breast Cancer/untitled folder`. Raw GEO + Figshare objects stay **local** (gitignored). Small tables (`SampleStats_GSE161529.txt`, `metadata_gse161529.csv`, `features.tsv`) can go in git.

## Place here

| File / dir | Source | Git |
|---|---|---|
| `GSE161529_RAW.tar` + extracted `GSE161529_RAW/` | GEO | ignored |
| `SeuratObject_{TNBC,HER2,ERTotal}Tum.rds` | Figshare 10.6084/m9.figshare.17058077 | ignored |
| `features.tsv` | GEO | ok |
| `SampleStats_GSE161529.txt` | Chen et al. companion | ok |
| `metadata_gse161529.csv` | from `00_extract_metadata_gse161529.R` | ok |
| `processed/*.npz` | from `00_build_matrix_gse161529.py` | ignored |

## Build

From the **repo root**:

```bash
Rscript "Breast Cancer/codes/00_extract_metadata_gse161529.R"
.venv/bin/python -u "Breast Cancer/codes/00_build_matrix_gse161529.py"
.venv/bin/python -u "Breast Cancer/codes/1run_panelA_ks_project_gse161529_sc.py"
```
