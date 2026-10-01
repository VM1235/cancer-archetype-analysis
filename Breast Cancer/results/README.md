# Breast `results/` map

Do not mix these folders. Panel A **fit** vs sc **projection** vs reverse Pal fit are different geometries.

## DepMap KS cell lines (fit here)

| Folder | What |
|---|---|
| `panel_a/` | Full ~16k-gene Panel A (not KS) |
| `panel_a_genelist/` | Other genelist track |
| `panel_a_ks_genelist/` | **Official KS k=3** on 63 lines (`p=0.022`) |
| `panel_a_ks_genelist_extendedk/` | Same lines, k=4–7 |
| `panel_b/` `panel_b_ks_genelist/` | PAM50 enrichment on DepMap vertices |
| `panel_c_*` | TCGA / METABRIC bulk projected into **DepMap** space |

## scRNA projected **into DepMap** KS k=3 (Hausser Fig. 4)

| Folder | Dataset |
|---|---|
| `panel_a_ks_genelist_sc_gse173634[_magic/_k4]` | Gambardella cell-line sc |
| `panel_a_ks_genelist_sc_gse176078[_magic/_k4]` | Wu tumor sc |
| `panel_a_ks_genelist_sc_gse161529_magic` | Pal tumor sc → DepMap triangle (13.3% inside) |

## Pal GSE161529 **self-fit** (reverse)

| Folder | What |
|---|---|
| `panel_a_ks_genelist_gse161529/` | Pal PCHA + 50-shuffle t-ratio |
| `panel_b_gse161529/` | Name Pal vertices (IHC + KS Epi/Mes) |
| `gse176078_on_gse161529/` | Wu → Pal triangle |
| `gse173634_on_gse161529/` | GSE173634 → Pal triangle |
| `depmap_on_gse161529/` | 63 lines → Pal triangle |

MAGIC `*.npz` caches are gitignored (>100 MB).
