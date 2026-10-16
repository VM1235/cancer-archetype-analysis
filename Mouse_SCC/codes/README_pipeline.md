# Mouse SCC archetype pipeline (non-MAGIC)

Fits a Pareto simplex (archetypes) to the 23,959 mouse SCC cells from GSE288137
and asks where the paper's EMT states sit in it. Run everything from `Project_1/`.

Pre-registered predictions: `Mouse_SCC/results/PREREGISTERED_PREDICTIONS.md`
(written before any fit).

## Steps

| Step | Script | What it does | Time |
|---|---|---|---|
| 1 | `01_build_anndata.py` | counts -> log1p CP10k; keep genes in >=1% of cells; drop mito/ribo genes (14,502 genes) | ~5 min |
| 2 | `02_pca.py --geneset all` | tumour-balanced fit cells (<=800 per tumour = 9,600), exact PCA, project all cells, PC diagnostics | ~2 min |
| 3 | `run_fits.sh all main 200` | ESV curve; PCHA for k = 3..7; permutation test (200 shuffles) | hours (parallel per k) |
| 4 | `04_describe_simplex.py --geneset all --k K` | weights for every cell, label tables, tumour check, enrichment, marker genes, figures | ~2 min |
| 5 | `05_robustness.py --geneset all --k K --variant loto` | full refit leaving out each tumour; also `scrna`, `nophase`, `seed` | ~3 min per refit |

`--geneset ks` runs the same thing on the 176 Tan 2014 KS genes found in mouse
(the gene panel used for the Pal breast fits).

Outputs go to `Mouse_SCC/results/<geneset>_<tag>/`.

## Choices and the reason for each

- **Tumour-balanced fit cells (<=800 per tumour).** Four tumours supply 75% of late-EMT cells.
  Without a cap the biggest tumours decide where the vertices go. All cells are still placed
  in the fitted simplex afterwards.
- **No gene scaling.** Genes are centred only, as in Hausser. Scaling would give rare, noisy
  genes the same weight as strong ones.
- **No batch correction in the main fit.** We want to see the raw structure. Batch effects
  are measured instead (step 2 diagnostics) and tested (step 5: tumour left out, scRNA only).
- **Exact PCA via the cell-by-cell matrix ("gram").** Gives the same PCs as a full SVD
  (checked on this data: score correlation 1.0000000) but 7x faster, which makes the
  repeated refits in step 5 practical.
- **PCHA engine from `src/archetypes.py`.** Same code as the breast work, so numbers are
  comparable: fit on the first k-1 PCs, delta = 0, 15 starts, keep the largest volume.
- **Permutation null.** Each of the first k-1 PCs is shuffled independently, then refit with
  5 starts. p = (1 + shuffles at least as good) / (1 + shuffles). The +1 stops p from being
  reported as exactly 0.
- **Specialist rule.** Largest weight minus second largest >= 0.35, as in the Pal work.
- **Vertex names are not automatic.** Step 4 prints the evidence (states, clusters, tumours,
  genes). Names are given only after the tumour breakdown is checked.

## Caveats to keep in mind

- With ~10,000 cells, a small p only says the data is not a shuffled blob. Batch structure
  alone can do that. It is not proof of a Pareto front.
- On all genes, PC1 is explained more by tumour and assay (multiome nuclei vs whole cells)
  than by EMT state. Check step 5 (`scrna`, `loto`) before trusting any vertex.
- The ESV "elbow" is a heuristic, not a test.
