#!/usr/bin/env Rscript
# Groves Panel C batch correction: sva::ComBat, intercept-only model, reference batch.
# Matches Cell-line-tumor-batch-correction-and-clustering.Rmd:
#   bc2 <- ComBat(as.matrix(dfc), batch=..., mod=mod0, ref.batch='m')
# Here batches are cell_line vs tumor; ref.batch is cell_line.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: combat_cellline_tumor.R merged.csv batch.csv out.csv [ref_batch] [rlib]")
}
merged_path <- args[[1]]
batch_path <- args[[2]]
out_path <- args[[3]]
ref_batch <- if (length(args) >= 4) args[[4]] else "cell_line"
lib <- if (length(args) >= 5) args[[5]] else file.path(getwd(), "rlib")
.libPaths(c(lib, .libPaths()))

if (!requireNamespace("sva", quietly = TRUE)) {
  stop("sva is not installed in ", lib)
}
library(sva)

merged <- as.matrix(read.csv(merged_path, row.names = 1, check.names = FALSE))
batch_df <- read.csv(batch_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("sample", "batch") %in% names(batch_df))) {
  stop("batch.csv must have columns sample, batch")
}
batch_df$sample <- as.character(batch_df$sample)
colnames(merged) <- as.character(colnames(merged))
if (!identical(batch_df$sample, colnames(merged))) {
  m <- match(colnames(merged), batch_df$sample)
  if (any(is.na(m))) stop("batch.csv samples do not match merged columns")
  batch_df <- batch_df[m, ]
}
batch <- factor(batch_df$batch)
if (!(ref_batch %in% levels(batch))) {
  stop("ref_batch not in batch levels: ", paste(levels(batch), collapse = ", "))
}
phen <- data.frame(row.names = colnames(merged))
mod0 <- model.matrix(~ 1, phen)
message("ComBat genes=", nrow(merged), " samples=", ncol(merged),
        " ref.batch=", ref_batch)
bc <- ComBat(dat = merged, batch = batch, mod = mod0, par.prior = TRUE,
             ref.batch = ref_batch)
write.csv(bc, file = out_path, quote = TRUE)
message("Wrote ", out_path)
