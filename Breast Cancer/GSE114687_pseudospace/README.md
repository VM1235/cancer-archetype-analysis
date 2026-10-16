# MCF10A and HuMEC pseudospace EMT (GSE114687)

Human mammary epithelial cells (MCF10A, immortalised; HuMEC, primary) grown as colonies, so EMT progresses with distance from the colony edge. Each line has an untreated (Mock) and a TGF-β arm, each with inner-colony and outer-colony samples (McFaline-Figueroa et al.). Colony position is an EMT axis that does not depend on an expression score. The CROP-seq screen in the same series is not used. Downloads are listed in `data/raw/DOWNLOAD_LOG.md` (not tracked).

## Pipeline

Run from the repository root.

| Step | Script | Purpose |
|---|---|---|
| 1 | `codes/01_build_anndata.py` | read the two Matrix Market files, keep the authors' cells, QC, CP10k + log1p |
| 1c | `codes/01c_emt_labels.py` | KS EMT score and state quartiles |
| 2 | `codes/02_pca.py` | PCA, at most 1,500 fit cells per sample |
| 3 | `codes/run_fits.sh ks\|all main 200` | pooled PCHA k = 3 to 7, 200 shuffles |
| 4, 6 | `04_describe_simplex`, `06_groves_specialists` | description, specialist rules |
| 8 | `codes/08_within_unit.py [--clean --regress-depth --out-suffix]` | fits inside one sample or one cell line |
| 8b, 8c, 9, 10 | `08b_describe_corners`, `08c_path_vs_triangle`, `09_hallmark_corners`, `10_describe_labels` | corner genes, path test, Hallmark scoring, state tables |
| 1e, 14 | `01e_depth_matched.py`, `14_depth_summary.py` | depth control |
| 13 | `codes/13_mail_figures.py` | slide figures in `results/mail_figures/` |

Memory: no more than about four `08_within_unit.py` jobs at once on a 16 GB machine.

## Results

See [`RESULTS.md`](RESULTS.md). In brief: pooled fits only recover the samples; inside samples, after cleaning, 18 of 24 fits pass (KS genes too in three samples) with three recurring corners (cycling, differentiated epithelial, EMT / matrix) aligned with colony position; the triangle beats a one-dimensional EMT path in 6 of 12 tests; the result survives depth regression (8 of 8) and depth matching (7 of 8); a hybrid state is the most generalist state in 7 of 8 fits (gap rule) and 8 of 8 (closest-5% rule).
