#!/usr/bin/env Rscript
# Paper-faithful residuals for Fig 3 (Anderson et al. Step0-style):
# per-study counts -> combined PML matrix -> TMM -> voom ->
# limma ~ Study + TopGenes with duplicateCorrelation(Patient).
#
# Inputs (gene x sample CSVs, Ensembl IDs without version):
#   data/processed/counts_GSE319666_LungPCA.csv   (from public GEO counts)
#   data/processed/counts_GSE102511.csv           (Salmon)
#   data/processed/counts_GSE166720.csv
#   data/processed/counts_GSE193725.csv
#   data/processed/PML_sample_info.csv            (137 PMLs)

suppressPackageStartupMessages({
  .libPaths(c(
    "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes/rlib",
    .libPaths()
  ))
  library(edgeR)
  library(limma)
  library(data.table)
})

ROOT <- "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes"
PROC <- file.path(ROOT, "data/processed")
OUT <- file.path(ROOT, "results/fig3_paperlike")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

info <- fread(file.path(PROC, "PML_sample_info.csv"))
auth <- fread(file.path(ROOT, "data/author/PML_Archetypes.csv"))
# ensure author order
info <- info[match(auth$sample_id, info$sample_id)]
stopifnot(!anyNA(info$sample_id))

read_counts <- function(path) {
  dt <- fread(path)
  gene <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- sub("\\..*$", "", gene)
  storage.mode(mat) <- "double"
  mat
}

# Lung PCA from GEO counts, subset to PML samples in info
lung_raw <- file.path(ROOT, "data/raw/GSE319666_counts.tsv.gz")
lung <- fread(lung_raw)
setnames(lung, 1, "gene")
lung[, gene := sub("\\..*$", "", gene)]
lung <- lung[!duplicated(gene)]
lung_ids <- info$sample_id[info$Study == "Lung_PCA" | info$study == "Lung_PCA"]
# column name may be Study or study
if (!"study" %in% names(info) && "Study" %in% names(info)) info$study <- info$Study
lung_ids <- info$sample_id[info$study == "Lung_PCA"]
stopifnot(all(lung_ids %in% names(lung)))
lung_mat <- as.matrix(lung[, ..lung_ids])
rownames(lung_mat) <- lung$gene
fwrite(data.table(gene = rownames(lung_mat), lung_mat),
       file.path(PROC, "counts_GSE319666_LungPCA.csv"))

need <- c(
  GSE102511 = file.path(PROC, "counts_GSE102511.csv"),
  GSE166720 = file.path(PROC, "counts_GSE166720.csv"),
  GSE193725 = file.path(PROC, "counts_GSE193725.csv"),
  Lung_PCA = file.path(PROC, "counts_GSE319666_LungPCA.csv")
)
missing <- need[!file.exists(need)]
if (length(missing)) {
  stop("Missing count matrices:\n", paste(missing, collapse = "\n"))
}

mats <- lapply(need, read_counts)
# Align samples to GSM / library IDs in info
align <- function(mat, study) {
  ids <- info$sample_id[info$study == study]
  # Salmon outputs use sample_id (GSM); Lung uses library id
  miss <- setdiff(ids, colnames(mat))
  if (length(miss)) stop(study, " missing samples: ", paste(miss, collapse = ","))
  mat[, ids, drop = FALSE]
}
mats2 <- list(
  Lung_PCA = align(mats$Lung_PCA, "Lung_PCA"),
  GSE102511 = align(mats$GSE102511, "GSE102511"),
  GSE166720 = align(mats$GSE166720, "GSE166720"),
  GSE193725 = align(mats$GSE193725, "GSE193725")
)

common <- Reduce(intersect, lapply(mats2, rownames))
message("Common genes: ", length(common))
counts <- do.call(cbind, lapply(mats2, function(m) m[common, , drop = FALSE]))
counts <- counts[, info$sample_id, drop = FALSE]

# Round Salmon estimated counts for edgeR
counts <- round(counts)

sample_info <- data.frame(
  sample_id = info$sample_id,
  Histology = info$histology,
  Study = factor(info$study),
  Patient = factor(info$patient),
  row.names = info$sample_id,
  stringsAsFactors = FALSE
)

dge <- DGEList(counts = counts, samples = sample_info)
keep <- filterByExpr(dge, group = dge$samples$Histology)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
message("Genes after filterByExpr: ", nrow(dge))

# TopGenes (% of counts in top 1000 genes)
cs <- rowSums(dge$counts)
top <- names(sort(cs, decreasing = TRUE))[seq_len(min(1000L, length(cs)))]
top_sum <- colSums(dge$counts[top, , drop = FALSE])
dge$samples$TopGenes <- as.numeric(top_sum / colSums(dge$counts) * 100)

design <- model.matrix(~ Study + TopGenes, dge$samples)
v <- voom(dge, design = design)
corfit <- duplicateCorrelation(v, design, block = dge$samples$Patient)
v2 <- voom(dge, design, block = dge$samples$Patient, correlation = corfit$consensus)
corfit2 <- duplicateCorrelation(v2, design, block = dge$samples$Patient)
fit <- lmFit(v2, design, block = dge$samples$Patient, correlation = corfit2$consensus)
fit <- eBayes(fit)
resid <- residuals(fit, v2)

saveRDS(resid, file.path(OUT, "PML_resid_paperlike.rds"))
fwrite(data.table(gene = rownames(resid), resid), file.path(OUT, "PML_resid_paperlike.csv"))
fwrite(data.table(sample_id = rownames(dge$samples), dge$samples),
       file.path(OUT, "PML_dge_samples.csv"))
message("Saved residuals ", nrow(resid), " x ", ncol(resid), " to ", OUT)
