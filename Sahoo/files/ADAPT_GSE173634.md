# Wiring MAGIC into run_panelA_ks_project_gse173634_sc.py

**Superseded:** use the fixed scripts
`Breast Cancer/codes/1run_panelA_ks_project_gse173634_sc.py` and
`Breast Cancer/codes/magic_utils.py`.

Sahoo protocol (from their GitHub Rmagic code, not just the paper Methods
sentence):

```r
expr <- t(counts)
expr <- library.size.normalize(expr)
expr <- sqrt(expr)
magic(expr, genes=<genes of interest>, solver="approximate")
```

Python equivalent: lib-size normalize (`rescale="median"`) + sqrt on the
FULL gene × QC-cell matrix, then `magic.MAGIC(solver="approximate").fit_transform(..., genes=KS)`.

Do **not** use the earlier HVG + log1p CP10k draft — that was wrong relative
to Sahoo's deposited code.
