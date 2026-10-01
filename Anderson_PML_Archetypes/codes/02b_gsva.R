#!/usr/bin/env Rscript
# Compute GSVA scores for author PML modules on residual expression.

suppressPackageStartupMessages({
  .libPaths(c(
    "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes/rlib",
    .libPaths()
  ))
  library(GSVA)
})

ROOT <- "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes"
PROC <- file.path(ROOT, "data/processed")
AUTH <- file.path(ROOT, "data/author")
OUT <- file.path(ROOT, "results/fig3")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

resid <- readRDS(file.path(PROC, "PML_resid.rds"))
modules <- readRDS(file.path(AUTH, "pml_modules.rds"))

# Keep genes present in residual matrix
modules <- lapply(modules, function(g) intersect(g, rownames(resid)))
message("Module sizes present: ", paste(names(modules), lengths(modules), sep = "=", collapse = ", "))

# GSVA Gaussian kernel (continuous residuals), matching author Step4 settings
param <- gsvaParam(
  exprData = as.matrix(resid),
  geneSets = modules,
  kcdf = "Gaussian",
  maxDiff = TRUE,
  absRanking = FALSE
)
scores <- gsva(param, verbose = TRUE)
# scores: gene sets x samples
scores <- as.matrix(scores)

write.csv(t(scores), file.path(OUT, "gsva_scores.csv"), quote = TRUE)
saveRDS(scores, file.path(OUT, "gsva_scores.rds"))
message("GSVA done: ", nrow(scores), " modules x ", ncol(scores), " samples")
