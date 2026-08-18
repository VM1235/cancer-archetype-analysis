# Glioblastoma

Placeholder for the same Groves / ParTI pipeline on glioblastoma. No analysis has been run yet.

When work starts, keep the same layout as SCLC and breast:

| Path | Put here |
|---|---|
| `data/` | Expression matrices the codes will load (document exact filenames in this README) |
| `codes/` | Dataset-specific drivers; reuse repo-root `src/` for PCHA, PCA, enrichment |
| `figures/` | Final panel figures |
| `results/` | Intermediate fits, tables, null distributions |

Expected data format (to match the other two projects): **genes × samples** CSV, log-expression, one column per sample. If two studies are merged, batch-correct before PCHA (ComBat), as in Groves Panel A (cell-line studies) or Panel C (lines vs tumors).

Do not copy SCLC or breast \(k\), PC count, or subtype names. Choose dimension from this matrix (e.g. PCs to ~50% variance) and \(k\) from ESV + t-ratio on these samples.

## How to run

Nothing to run yet. After scripts exist, invoke them from the **repository root**, for example:

```bash
.venv/bin/python -u "Glioblastoma/codes/run_panelA.py"
```
