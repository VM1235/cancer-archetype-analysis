#!/usr/bin/env Rscript
# Path B: improved public-data Fig 3 rebuild (no SRA).
# - limma residuals ~ Study + TopGenes with duplicateCorrelation(Patient)
# - GSVA (Gaussian) on author modules
# - ParetoTI bootstrap k=4 if available, else note for Python PCHA
# - export tables for plotting / QC

suppressPackageStartupMessages({
  .libPaths(c(
    "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes/rlib",
    .libPaths()
  ))
  library(limma)
  library(GSVA)
  library(data.table)
})

ROOT <- "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes"
PROC <- file.path(ROOT, "data/processed")
AUTH <- file.path(ROOT, "data/author")
OUT <- file.path(ROOT, "results/fig3_v2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42792)

expr <- readRDS(file.path(PROC, "PML_expr_log_harmonized.rds"))
info <- fread(file.path(PROC, "PML_sample_info.csv"))
auth <- fread(file.path(AUTH, "PML_Archetypes.csv"))

# Align to author sample order
ord <- match(auth$sample_id, info$sample_id)
stopifnot(!anyNA(ord))
info <- info[ord]
expr <- expr[, info$sample_id, drop = FALSE]
stopifnot(identical(colnames(expr), info$sample_id))

info$Study <- factor(info$study)
info$Patient <- factor(info$patient)
info$TopGenes <- as.numeric(info$TopGenes)

# Drop near-zero variance genes
rv <- apply(expr, 1, var, na.rm = TRUE)
expr <- expr[is.finite(rv) & rv > 1e-8, , drop = FALSE]
message("Genes: ", nrow(expr), " Samples: ", ncol(expr))

design <- model.matrix(~ Study + TopGenes, data = as.data.frame(info))
v <- expr  # already log-scale continuous (not counts) — use limma directly
fit0 <- lmFit(v, design)
# Patient blocking (paper uses duplicateCorrelation on voom; we apply same idea)
corfit <- duplicateCorrelation(v, design, block = info$Patient)
message("consensus correlation: ", signif(corfit$consensus, 4))
fit <- lmFit(v, design, block = info$Patient, correlation = corfit$consensus)
fit <- eBayes(fit)
resid <- residuals(fit, v)

fwrite(data.table(gene = rownames(resid), resid), file.path(OUT, "PML_resid_v2.csv"))
saveRDS(resid, file.path(OUT, "PML_resid_v2.rds"))

# GSVA
modules <- readRDS(file.path(AUTH, "pml_modules.rds"))
modules <- lapply(modules, function(g) intersect(g, rownames(resid)))
message("Module sizes: ", paste(names(modules), lengths(modules), sep = "=", collapse = ", "))

param <- gsvaParam(
  exprData = as.matrix(resid),
  geneSets = modules,
  kcdf = "Gaussian",
  maxDiff = TRUE,
  absRanking = FALSE
)
scores <- as.matrix(gsva(param, verbose = TRUE))
# scores: modules x samples
write.csv(t(scores), file.path(OUT, "gsva_scores.csv"))
saveRDS(scores, file.path(OUT, "gsva_scores.rds"))
message("GSVA done")

# ParetoTI if available
if (requireNamespace("ParetoTI", quietly = TRUE)) {
  library(ParetoTI)
  GSVAScores <- as.data.frame(scores)
  message("Fitting ParetoTI bootstrap n=5000 ...")
  arc_data <- fit_pch_bootstrap(
    GSVAScores,
    n = 5000L,
    sample_prop = 0.8,
    noc = 4L,
    delta = 0,
    seed = 42792,
    conv_crit = 1e-6
  )
  arc_pml <- average_pch_fits(arc_data)
  saveRDS(arc_pml, file.path(OUT, "PML_arcA4_v2.rds"))

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
    arc_pml, as.matrix(GSVAScores), n_dim = 6, pc_method = "svd",
    log2 = FALSE, zscore = FALSE
  )
  saveRDS(list(data_attr = data_attr, pcs = pcs, arc = arc_pml),
          file.path(OUT, "archetype_fit_v2.rds"))

  # Export PC + distances for Python plotting fallback / QC
  pc_mat <- t(pcs$data)
  # distance columns from data_attr
  dd <- as.data.frame(data_attr$data)
  fwrite(data.table(sample_id = rownames(pc_mat), pc_mat), file.path(OUT, "pc_scores_pareto.csv"))
  fwrite(dd, file.path(OUT, "arch_distances_pareto.csv"))
  writeLines("paretoti", file.path(OUT, "fit_method.txt"))
  message("ParetoTI fit saved")
} else {
  writeLines("pcha_python", file.path(OUT, "fit_method.txt"))
  message("ParetoTI not installed — use Python PCHA on gsva_scores.csv")
}

# Sample info for plotting
info$author_combined <- auth$Combined
fwrite(info, file.path(OUT, "sample_info_v2.csv"))
message("Wrote ", OUT)
