#!/usr/bin/env Rscript
# Count-level Fig 3 rebuild (no author residual RDS):
#   Lung PCA GEO counts + NCBI HISAT2/featureCounts for GSE102511/166720/193725
#   -> Ensembl intersect -> Step0 TMM/voom/duplicateCorrelation -> GSVA
#
# Outputs: results/fig3_ncbi_voom/

suppressPackageStartupMessages({
  .libPaths(c(
    "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes/rlib",
    .libPaths()
  ))
  library(data.table)
  library(edgeR)
  library(limma)
  library(GSVA)
})

ROOT <- "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes"
RAW <- file.path(ROOT, "data/raw")
PROC <- file.path(ROOT, "data/processed")
AUTH <- file.path(ROOT, "data/author")
OUT <- file.path(ROOT, "results/fig3_ncbi_voom")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42792)

strip_ens <- function(x) sub("\\..*$", "", x)

info <- fread(file.path(PROC, "PML_sample_info.csv"))
auth <- fread(file.path(AUTH, "PML_Archetypes.csv"))
ord <- match(auth$sample_id, info$sample_id)
stopifnot(!anyNA(ord))
info <- info[ord]

ens_ent <- as.data.frame(readRDS(file.path(AUTH, "ENS_ENT_ID.rds")))
ens_ent$entrezgene_id <- as.character(ens_ent$entrezgene_id)
ens_ent$ensembl_gene_id <- strip_ens(ens_ent$ensembl_gene_id)
ens_ent <- ens_ent[!is.na(ens_ent$entrezgene_id) & ens_ent$entrezgene_id != "" &
                     !is.na(ens_ent$ensembl_gene_id) & ens_ent$ensembl_gene_id != "", ]
# prefer first Ensembl per Entrez
ens_ent <- ens_ent[!duplicated(ens_ent$entrezgene_id), ]
ent2ens <- setNames(ens_ent$ensembl_gene_id, ens_ent$entrezgene_id)

read_ncbi_ensembl <- function(path, gsm_ids) {
  dt <- fread(path)
  setnames(dt, 1, "GeneID")
  dt[, GeneID := as.character(GeneID)]
  miss <- setdiff(gsm_ids, names(dt))
  if (length(miss)) stop("Missing GSMs in ", basename(path), ": ", paste(miss, collapse = ", "))
  ens <- ent2ens[dt$GeneID]
  keep <- !is.na(ens)
  dt <- dt[keep]
  ens <- ens[keep]
  mat <- as.matrix(dt[, ..gsm_ids])
  storage.mode(mat) <- "double"
  # collapse duplicate Ensembl (sum counts)
  split_idx <- split(seq_len(nrow(mat)), ens)
  out <- matrix(0, nrow = length(split_idx), ncol = ncol(mat))
  rownames(out) <- names(split_idx)
  colnames(out) <- gsm_ids
  for (i in seq_along(split_idx)) {
    rows <- split_idx[[i]]
    if (length(rows) == 1L) out[i, ] <- mat[rows, ]
    else out[i, ] <- colSums(mat[rows, , drop = FALSE])
  }
  round(out)
}

# ---- Lung PCA (Ensembl, GEO counts) ----
lung <- fread(file.path(RAW, "GSE319666_counts.tsv.gz"))
setnames(lung, 1, "gene")
lung[, gene := strip_ens(gene)]
lung <- lung[!duplicated(gene)]
lung_ids <- info$sample_id[info$study == "Lung_PCA"]
stopifnot(all(lung_ids %in% names(lung)))
lung_mat <- as.matrix(lung[, ..lung_ids])
rownames(lung_mat) <- lung$gene
storage.mode(lung_mat) <- "double"
lung_mat <- round(lung_mat)
message("Lung PCA: ", nrow(lung_mat), " genes x ", ncol(lung_mat), " PMLs")

# ---- NCBI counts, columns already GSM = sample_id ----
ncbi_files <- c(
  GSE102511 = "GSE102511_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  GSE166720 = "GSE166720_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  GSE193725 = "GSE193725_raw_counts_GRCh38.p13_NCBI.tsv.gz"
)
ncbi_mats <- lapply(names(ncbi_files), function(st) {
  ids <- info$sample_id[info$study == st]
  m <- read_ncbi_ensembl(file.path(RAW, ncbi_files[[st]]), ids)
  message(st, ": ", nrow(m), " Ensembl x ", ncol(m), " PMLs")
  m
})
names(ncbi_mats) <- names(ncbi_files)

common <- Reduce(intersect, c(list(rownames(lung_mat)), lapply(ncbi_mats, rownames)))
message("Common Ensembl genes: ", length(common))

mats <- c(list(Lung_PCA = lung_mat[common, , drop = FALSE]),
          lapply(ncbi_mats, function(m) m[common, , drop = FALSE]))
counts <- do.call(cbind, mats)
counts <- counts[, info$sample_id, drop = FALSE]
stopifnot(identical(colnames(counts), info$sample_id))
message("Combined counts: ", nrow(counts), " x ", ncol(counts))

fwrite(data.table(gene = rownames(counts), counts),
       file.path(PROC, "PML_counts_ncbi_combined.csv"))
saveRDS(counts, file.path(PROC, "PML_counts_ncbi_combined.rds"))

# ---- Step0 ----
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

cs <- rowSums(dge$counts)
top <- names(sort(cs, decreasing = TRUE))[seq_len(min(1000L, length(cs)))]
dge$samples$TopGenes <- as.numeric(colSums(dge$counts[top, , drop = FALSE]) /
                                     colSums(dge$counts) * 100)

design <- model.matrix(~ Study + TopGenes, dge$samples)
message("voom pass 1 ...")
v <- voom(dge, design = design)
message("duplicateCorrelation pass 1 ...")
corfit <- duplicateCorrelation(v, design, block = dge$samples$Patient)
message("consensus cor 1: ", signif(corfit$consensus, 4))
message("voom pass 2 ...")
v2 <- voom(dge, design, block = dge$samples$Patient, correlation = corfit$consensus)
message("duplicateCorrelation pass 2 ...")
corfit2 <- duplicateCorrelation(v2, design, block = dge$samples$Patient)
message("consensus cor 2: ", signif(corfit2$consensus, 4))
fit <- lmFit(v2, design, block = dge$samples$Patient, correlation = corfit2$consensus)
fit <- eBayes(fit)
resid <- residuals(fit, v2)
message("Residuals: ", nrow(resid), " x ", ncol(resid))

saveRDS(resid, file.path(OUT, "PML_resid_ncbi_voom.rds"))
fwrite(data.table(gene = rownames(resid), resid), file.path(OUT, "PML_resid_ncbi_voom.csv"))

info_out <- copy(info)
info_out$TopGenes <- dge$samples$TopGenes[match(info_out$sample_id, rownames(dge$samples))]
info_out$Histology <- info_out$histology
info_out$Study <- info_out$study
info_out$Patient <- info_out$patient
info_out$author_combined <- auth$Combined[match(info_out$sample_id, auth$sample_id)]
fwrite(info_out, file.path(OUT, "sample_info.csv"))

# ---- GSVA (author modules, Gaussian) ----
modules <- readRDS(file.path(AUTH, "pml_modules.rds"))
modules <- lapply(modules, function(g) intersect(strip_ens(g), rownames(resid)))
message("Module sizes after residual gene filter: ",
        paste(names(modules), lengths(modules), sep = "=", collapse = ", "))
if (any(lengths(modules) < 2)) {
  warning("Some modules have <2 genes after filtering")
}

param <- gsvaParam(
  exprData = as.matrix(resid),
  geneSets = modules,
  kcdf = "Gaussian",
  maxDiff = TRUE,
  absRanking = FALSE
)
scores <- as.matrix(gsva(param, verbose = TRUE))
write.csv(t(scores), file.path(OUT, "gsva_scores.csv"))
saveRDS(scores, file.path(OUT, "gsva_scores.rds"))

writeLines(c(
  "counts: GSE319666 GEO + NCBI GRCh38.p13 raw counts (HISAT2/featureCounts)",
  "IDs: Entrez->Ensembl via ENS_ENT_ID.rds; gene intersect then filterByExpr(Histology)",
  "residuals: Step0 TMM/voom two-pass duplicateCorrelation(Patient) ~ Study + TopGenes",
  "GSVA: kcdf=Gaussian on author pml_modules.rds"
), file.path(OUT, "method.txt"))

message("Wrote ", OUT)
