#!/usr/bin/env Rscript
# Residualize harmonized PML expression for Study + TopGenes (limma), matching
# author Step0 design (without edgeR filterByExpr / duplicateCorrelation for
# cross-platform TPM+counts mix).

suppressPackageStartupMessages({
  .libPaths(c(
    "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes/rlib",
    .libPaths()
  ))
  library(limma)
})

ROOT <- "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes"
PROC <- file.path(ROOT, "data/processed")

expr <- readRDS(file.path(PROC, "PML_expr_log_harmonized.rds"))
info <- readRDS(file.path(PROC, "PML_sample_info.rds"))
stopifnot(identical(colnames(expr), info$sample_id))

# Drop very low-variance genes
rv <- apply(expr, 1, var, na.rm = TRUE)
keep <- is.finite(rv) & rv > 1e-8
expr <- expr[keep, , drop = FALSE]
message("Genes after variance filter: ", nrow(expr))

info$Study <- factor(info$Study)
design <- model.matrix(~ Study + TopGenes, data = info)
fit <- lmFit(expr, design)
fit <- eBayes(fit)
resid <- residuals(fit, expr)

saveRDS(resid, file.path(PROC, "PML_resid.rds"))
# CSV for Python (genes as first column)
out <- cbind(gene = rownames(resid), as.data.frame(resid, check.names = FALSE))
data.table::fwrite(out, file.path(PROC, "PML_resid.csv"))
message("Saved residuals: ", nrow(resid), " x ", ncol(resid))
