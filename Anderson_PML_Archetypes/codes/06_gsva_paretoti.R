#!/usr/bin/env Rscript
# GSVA + ParetoTI bootstrap archetypes (Anderson Step4) + Fig 3A/B plots.
# Falls back to a message if ParetoTI is not installed.

suppressPackageStartupMessages({
  .libPaths(c(
    "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes/rlib",
    .libPaths()
  ))
  library(GSVA)
  library(ggplot2)
})

ROOT <- "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes"
AUTH <- file.path(ROOT, "data/author")
OUT <- file.path(ROOT, "results/fig3_paperlike")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42792)

resid <- as.matrix(data.table::fread(file.path(OUT, "PML_resid_paperlike.csv")), rownames = 1)
# fread matrix helper
resid_dt <- data.table::fread(file.path(OUT, "PML_resid_paperlike.csv"))
genes <- resid_dt[[1]]
resid <- as.matrix(resid_dt[, -1, with = FALSE])
rownames(resid) <- genes

modules <- readRDS(file.path(AUTH, "pml_modules.rds"))
modules <- lapply(modules, function(g) intersect(g, rownames(resid)))
message("Module sizes: ", paste(names(modules), lengths(modules), sep = "=", collapse = ", "))

param <- gsvaParam(
  exprData = resid,
  geneSets = modules,
  kcdf = "Gaussian",
  maxDiff = TRUE,
  absRanking = FALSE
)
scores <- as.matrix(gsva(param, verbose = TRUE))
write.csv(t(scores), file.path(OUT, "gsva_scores.csv"))

has_pareto <- requireNamespace("ParetoTI", quietly = TRUE)
if (!has_pareto) {
  message("ParetoTI not installed. Install with:")
  message('  remotes::install_github("vitkl/ParetoTI")')
  message("Writing GSVA scores only; run Python PCHA fallback or install ParetoTI.")
  quit(save = "no", status = 0)
}

library(ParetoTI)
GSVAScores <- as.data.frame(scores)
arc_data <- fit_pch_bootstrap(
  GSVAScores,
  n = 5000,
  sample_prop = 0.8,
  noc = as.integer(4),
  delta = 0,
  seed = 42792,
  conv_crit = 1e-6
)
arc_pml <- average_pch_fits(arc_data)
saveRDS(arc_pml, file.path(OUT, "PML_arcA4_pml.rds"))

# Distances + PC projection
data_attr <- merge_arch_dist(
  arc_data = arc_pml,
  data = GSVAScores,
  feature_data = GSVAScores,
  colData = data.frame(sample_id = colnames(GSVAScores), row.names = colnames(GSVAScores)),
  dist_metric = "euclidean",
  colData_id = "sample_id",
  rank = FALSE
)
pcs <- project_to_pcs(
  arc_pml,
  as.matrix(GSVAScores),
  n_dim = 6,
  pc_method = "svd",
  log2 = FALSE,
  zscore = FALSE
)
saveRDS(list(data_attr = data_attr, pcs = pcs, arc = arc_pml), file.path(OUT, "archetype_fit.rds"))
message("ParetoTI fit complete")
