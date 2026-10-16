# Cook cell lines, TGF-β time course (GSE147405)

Four human cell lines (A549, DU145, MCF7, OVCA420) sampled by scRNA-seq at 0 d, 8 h, 1 d, 3 d and 7 d after TGF-β1, plus three removal arms (8 h, 1 d and 3 d after withdrawal). Time is an EMT axis that does not depend on any expression score. The removal-arm cells are held out of every fit and projected afterwards.

## Pipeline

Run from the repository root.

| Step | Script | Purpose |
|---|---|---|
| 1 | `codes/01_build_anndata.py` | merge the four matrices (11,250 shared genes), CP10k + log1p |
| 2 | `codes/02_pca.py --lines all\|A549\|...` | pooled (at most 550 cells per line) or single-line PCA |
| 3 | `codes/run_fits.sh` | PCHA k = 3 to 7, 200 shuffles |
| 4 | `codes/04_figures.py` | overview figures |
| 6 | `codes/06_groves_specialists.py` | closest-5% specialist rule |
| 7 | `codes/07_within_line_clean.py` | per-line fits after cleaning (at least 10,000 UMI, mitochondrial fraction at most 0.07), optional time balancing |
| 8 | `codes/08_path_vs_triangle.py` | triangle versus a pseudotime path |
| 9 to 12 | `09_describe_labels`, `10_pooled_tetrahedron`, `11_per_line_triangles`, `12_corner_programmes` | label tables, figures, gene-set scores per corner |

Raw files go in `data/raw/` (GEO supplementary files, not tracked).

## Results in brief

- Pooled: a simplex at every k from 3 to 7 (p = 0.005), but each vertex is one cell line.
- Per line, cleaned: A549 on all genes (p = 0.010) and MCF7 on KS genes (p = 0.03) pass at k = 3; no line passes at k = 4. A549 loses significance when time points are balanced (p = 0.224).
- Corners are proliferation and TGF-β response; KS epithelial and mesenchymal gene sets stay near zero at every corner.
- A pseudotime path explains 7 of 8 per-line fits as well as a triangle.
- In A549 the middle time point (1 d) is 89% generalists; removal-arm cells stay on the late side.

Details: [`../../docs/emt_single_cell/EMT_archetype_analysis_FINAL_REPORT.md`](../../docs/emt_single_cell/EMT_archetype_analysis_FINAL_REPORT.md), section 4.
