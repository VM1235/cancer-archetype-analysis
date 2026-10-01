#!/usr/bin/env Rscript
# Build a genes x 137-PML expression matrix + sample metadata for Fig 3A/B.
# Lung PCA: raw counts (GSE319666). Public cohorts: published TPM/FPKM matrices.
# Gene IDs harmonized to Ensembl (no version) using author ENS_SYM_ID.rds.

suppressPackageStartupMessages({
  .libPaths(c(
    file.path(dirname(normalizePath(".")), "Anderson_PML_Archetypes", "rlib"),
    "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes/rlib",
    .libPaths()
  ))
  library(data.table)
})

ROOT <- "/Users/apple/Desktop/thesis_iisc/Project_1/Anderson_PML_Archetypes"
RAW <- file.path(ROOT, "data/raw")
AUTH <- file.path(ROOT, "data/author")
PROC <- file.path(ROOT, "data/processed")
dir.create(PROC, showWarnings = FALSE, recursive = TRUE)

auth <- fread(file.path(AUTH, "PML_Archetypes.csv"))
ens_sym <- readRDS(file.path(AUTH, "ENS_SYM_ID.rds"))
ens_sym <- as.data.frame(ens_sym)
ens_sym <- ens_sym[!duplicated(ens_sym$hgnc_symbol) & ens_sym$hgnc_symbol != "", ]
sym2ens <- setNames(ens_sym$ensembl_gene_id, ens_sym$hgnc_symbol)

strip_ens <- function(x) sub("\\..*$", "", x)

# ---- Lung PCA counts ----
lung <- fread(file.path(RAW, "GSE319666_counts.tsv.gz"))
gene_col <- names(lung)[1]
setnames(lung, gene_col, "gene")
lung[, gene := strip_ens(gene)]
lung <- lung[!duplicated(gene)]
lung_meta <- fread(file.path(PROC, "GSE319666_sample_meta.csv"))
lung_pml_ids <- auth$sample_id[!grepl("^GSM", auth$sample_id)]
stopifnot(all(lung_pml_ids %in% names(lung)))
lung_mat <- as.matrix(lung[, ..lung_pml_ids])
rownames(lung_mat) <- lung$gene
# integer-ish counts
storage.mode(lung_mat) <- "double"

lung_info <- data.frame(
  sample_id = lung_pml_ids,
  matrix_col = lung_pml_ids,
  study = "Lung_PCA",
  histology = lung_meta$histology[match(lung_pml_ids, lung_meta$sample_id)],
  patient = sub("_[^_]+$", "", lung_pml_ids),
  stringsAsFactors = FALSE
)
lung_info$histology[lung_info$histology %in% c("AIS", "MIA")] <- "AIS/MIA"
lung_info$histology[lung_info$histology == "ADC"] <- "LUAD"

# ---- Public GEO matrices (symbol -> Ensembl) ----
pub_map <- fread(file.path(PROC, "public_PML_gsm_map.csv"))

read_symbol_mat <- function(path, gene_col_candidates = c("Unique ID", "Gene", "gene", "V1")) {
  dt <- fread(path)
  if (names(dt)[1] %in% c("V1", "") || names(dt)[1] == "") {
    setnames(dt, 1, "gene")
  } else if (!"gene" %in% names(dt)) {
    hit <- intersect(gene_col_candidates, names(dt))
    if (!length(hit)) hit <- names(dt)[1]
    setnames(dt, hit[1], "gene")
  }
  dt <- dt[!is.na(gene) & gene != ""]
  dt <- dt[!duplicated(gene)]
  dt
}

collapse_to_ensembl <- function(dt, sample_cols) {
  ens <- sym2ens[dt$gene]
  keep <- !is.na(ens)
  dt <- dt[keep]
  ens <- ens[keep]
  mat <- as.matrix(dt[, ..sample_cols])
  storage.mode(mat) <- "double"
  # average duplicate Ensembl mappings
  split_idx <- split(seq_len(nrow(mat)), ens)
  out <- matrix(NA_real_, nrow = length(split_idx), ncol = ncol(mat))
  rownames(out) <- names(split_idx)
  colnames(out) <- sample_cols
  for (i in seq_along(split_idx)) {
    rows <- split_idx[[i]]
    if (length(rows) == 1L) out[i, ] <- mat[rows, ]
    else out[i, ] <- colMeans(mat[rows, , drop = FALSE], na.rm = TRUE)
  }
  out
}

# GSE102511 TPM (columns like 001-AAH)
g102 <- read_symbol_mat(file.path(RAW, "GSE102511_Smruthy_etal_allsamples_TPM.txt.gz"))
m102 <- pub_map[study == "GSE102511"]
stopifnot(all(m102$matrix_col %in% names(g102)))
mat102 <- collapse_to_ensembl(g102, m102$matrix_col)
colnames(mat102) <- m102$sample_id
info102 <- data.frame(
  sample_id = m102$sample_id,
  matrix_col = m102$matrix_col,
  study = "GSE102511",
  histology = m102$histology,
  patient = sub("-.*$", "", m102$matrix_col),
  stringsAsFactors = FALSE
)

# GSE166720 FPKM
g166 <- read_symbol_mat(file.path(RAW, "GSE166720_NIHAD_FPKM_normalized.txt.gz"))
# first col may be unnamed gene
if (!"gene" %in% names(g166)) setnames(g166, 1, "gene")
m166 <- pub_map[study == "GSE166720"]
stopifnot(all(m166$matrix_col %in% names(g166)))
mat166 <- collapse_to_ensembl(g166, m166$matrix_col)
colnames(mat166) <- m166$sample_id
info166 <- data.frame(
  sample_id = m166$sample_id,
  matrix_col = m166$matrix_col,
  study = "GSE166720",
  histology = m166$histology,
  patient = m166$matrix_col,
  stringsAsFactors = FALSE
)

# GSE193725 FPKM
g193 <- read_symbol_mat(file.path(RAW, "GSE193725_GGO_pNS_Solid_RNAseq_FPKM.txt.gz"))
m193 <- pub_map[study == "GSE193725"]
# column names in file use dots: GGO_10 vs GGO.N_12; titles are GGO_10
# fread may keep GGO_10 as-is
miss <- setdiff(m193$matrix_col, names(g193))
if (length(miss)) {
  # try replacing _ with .
  alt <- gsub("_", ".", miss, fixed = TRUE)
  map_fix <- setNames(names(g193), names(g193))
  for (i in seq_along(miss)) {
    if (alt[i] %in% names(g193)) {
      # rename for join convenience
    }
  }
}
# Build column rename: matrix_col -> actual file col
resolve_col <- function(want, have) {
  if (want %in% have) return(want)
  w2 <- gsub("_", ".", want, fixed = TRUE)
  if (w2 %in% have) return(w2)
  stop("Missing column: ", want)
}
cols193 <- vapply(m193$matrix_col, resolve_col, character(1), have = names(g193))
mat193 <- collapse_to_ensembl(g193, unname(cols193))
colnames(mat193) <- m193$sample_id
info193 <- data.frame(
  sample_id = m193$sample_id,
  matrix_col = m193$matrix_col,
  study = "GSE193725",
  histology = m193$histology,
  patient = m193$matrix_col,
  stringsAsFactors = FALSE
)

# ---- Harmonize to common genes; expression scale ----
# Lung PCA: keep as counts for now (handled in residuals script).
# Public: log2(x+1) of TPM/FPKM then treat as continuous expression.
log_mat <- function(m) log2(m + 1)

mats_pub <- list(
  GSE102511 = log_mat(mat102),
  GSE166720 = log_mat(mat166),
  GSE193725 = log_mat(mat193)
)

common_genes <- Reduce(intersect, c(list(rownames(lung_mat)), lapply(mats_pub, rownames)))
message("Common Ensembl genes: ", length(common_genes))

lung_mat <- lung_mat[common_genes, , drop = FALSE]
mats_pub <- lapply(mats_pub, function(m) m[common_genes, , drop = FALSE])

# Combined continuous matrix for residualization:
# convert Lung PCA counts -> logCPM (library-size normalized) so all studies share a log scale.
lib <- colSums(lung_mat)
lung_logcpm <- sweep(log2(lung_mat + 0.5), 2, log2(lib / 1e6), "-")

expr <- cbind(lung_logcpm, mats_pub$GSE102511, mats_pub$GSE166720, mats_pub$GSE193725)
sample_info <- rbind(lung_info, info102, info166, info193)
# order to match author sample order
ord <- match(auth$sample_id, sample_info$sample_id)
stopifnot(!anyNA(ord))
sample_info <- sample_info[ord, , drop = FALSE]
expr <- expr[, sample_info$sample_id, drop = FALSE]

# TopGenes metric (% signal in top 1000 genes) on this expression matrix
rank_sum <- function(v) {
  o <- order(v, decreasing = TRUE)
  sum(v[o[seq_len(min(1000L, length(v)))]]) / sum(v)
}
# For log expression, use un-logged proxy via 2^x for TopGenes to mimic count concentration
top_genes <- vapply(seq_len(ncol(expr)), function(j) {
  x <- 2^expr[, j]
  rank_sum(x) * 100
}, numeric(1))
sample_info$TopGenes <- top_genes
sample_info$Histology <- sample_info$histology
sample_info$Study <- sample_info$study
sample_info$Patient <- sample_info$patient

# attach author labels for later comparison
sample_info$author_combined <- auth$Combined[match(sample_info$sample_id, auth$sample_id)]

saveRDS(expr, file.path(PROC, "PML_expr_log_harmonized.rds"))
saveRDS(sample_info, file.path(PROC, "PML_sample_info.rds"))
fwrite(sample_info, file.path(PROC, "PML_sample_info.csv"))
message("Saved ", nrow(expr), " genes x ", ncol(expr), " samples")
message("Histology:\n")
print(table(sample_info$Histology, sample_info$Study))
