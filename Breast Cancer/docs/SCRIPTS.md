# Breast `codes/` — what to run, in order

Always from the **repository root**. Shared math is repo-root `src/`.

`1run_*.py` = MAGIC version of the older `run_panelA_ks_project_*` scripts. Prefer `1run_*` for numbers in the Pal reverse story.

## 1. DepMap KS Panel A / B / C

| Script | Role |
|---|---|
| `01_build_input_panelA.py` | 63-line matrix from DepMap |
| `00_build_ks_genelist_input.py` | Restrict to Tan 2014 KS (206 genes) |
| `run_panelA_ks_genelist.py` | PCHA k=3 + 500-shuffle t-ratio |
| `run_panelA_ks_genelist_extendedk.py` | k=4–7 |
| `03_map_and_pam50.R` + `04_match_pam50_to_panelA.py` | PAM50 labels |
| `run_panelB_ks_genelist.py` | PAM50 enrichment on DepMap vertices |
| `prepare_tcga_brca.py` / `prepare_metabric_brca.py` | Bulk tumors |
| `run_panelC_*.py` | Project tumors into **DepMap** space (ComBat) |

## 2. scRNA → DepMap KS triangle

Raw Wu / GSE173634 counts live in repo-root `Sahoo/` (gitignored dumps).

| Script | Role |
|---|---|
| `1run_panelA_ks_project_gse173634_sc.py` | Gambardella → DepMap k=3 |
| `1run_panelA_ks_project_gse176078_sc.py` | Wu → DepMap k=3 |
| `00_extract_metadata_gse161529.R` + `00_build_matrix_gse161529.py` | Build Pal matrix (`data/gse161529/`) |
| `1run_panelA_ks_project_gse161529_sc.py` | Pal → DepMap k=3 |

## 3. Reverse: fit Pal, project others in

| Script | Role |
|---|---|
| `17.9_run_panelA_ks_genelist_gse161529.py` | Pal PCHA + ESV (no permutation) |
| `20_t_ratio_gse161529.py` | 50-shuffle t-ratio at k=3 and k=4 |
| `18_project_gse176078_onto_gse161529.py` | Wu → Pal |
| `19_project_gse173634_onto_gse161529.py` | GSE173634 → Pal |
| `21_panelB_gse161529.py` | Name Pal vertices |
| `22_project_depmap_onto_gse161529.py` | 63 lines → Pal |

`magic_utils.py` is imported by the `1run_*` scripts; not run alone.
