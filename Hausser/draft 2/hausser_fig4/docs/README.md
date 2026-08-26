# Hausser Fig. 4 — single cells in the METABRIC archetype shape

Reproduction of **Figure 4a / 4c / 4d** from Hausser et al. 2019
(*Nat. Commun.* 10:5423), following `Karaayvaz2018/scBC.R`.

## What the paper does

1. Fit PCA + 4-archetype tetrahedron on **METABRIC bulk** breast tumors (ParTI).
2. Filter Karaayvaz 2018 scRNA-seq → ~**650 cancer cells** from 6 tumors.
3. Keep ~**1964 genes** expressed in ≥50% of those cells and present in METABRIC.
4. **Fig 4a:** project cells onto METABRIC **PC1–3**; draw the tetrahedron
   (archetypes re-projected on the **same common genes**).
5. **Fig 4c:** project onto METABRIC **PC1, PC2, PC50**.
6. **Fig 4d:** METABRIC PCs explain ~25% of sc-PC variance; shuffled ~1.5%.

## Why the first draft looked wrong

It plotted bulk tumors + sc cells + archetypes fitted in *full-gene* PC space,
while projecting cells with only a gene subset. That collapsed the sc cloud to
a speck at the origin. The original code re-projects archetypes with the
**same common-gene loadings** as the cells, and Fig 4a shows **cells + tetrahedron
only** (no bulk cloud).

## Run

```bash
.venv/bin/python -u "Hausser/draft 2/hausser_fig4/codes/01_project_sc_onto_metabric_shape.py"
.venv/bin/python -u "Hausser/draft 2/hausser_fig4/codes/02_export_figure_4.py"
.venv/bin/python -u "Hausser/draft 2/hausser_fig4/codes/03_fev_panel.py"
```

Inputs come from the original release under
`Hausser_Original/.../Karaayvaz2018/` (`expMatrix.csv`, `arcsOrig_genes.csv`,
`GSE118389_counts_rsem.txt`).
