#!/usr/bin/env Rscript
# Step 0b–0c: map gene symbols -> Entrez, then PAM50 via genefu::molecular.subtyping.
# Sample IDs are kept as the Panel A matrix column names (ACH-...).

trailing <- commandArgs(trailingOnly = TRUE)
if (length(trailing) >= 1) {
  breast <- normalizePath(trailing[1])
} else {
  # Prefer cwd if it already contains the Panel A matrix (run from Breast Cancer/).
  processed <- file.path(getwd(), "data", "processed", "input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv")
  if (file.exists(processed) || file.exists("input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv")) {
    breast <- normalizePath(getwd())
  } else {
    stop("Pass the Breast Cancer folder as the first argument.")
  }
}

in_matrix <- file.path(breast, "data", "processed", "input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv")
if (!file.exists(in_matrix)) {
  in_matrix <- file.path(breast, "input_panelA_invasive_breast_carcinoma_ccle_logtpm_filtered.csv")
}
out_entrez <- file.path(breast, "data", "processed", "input_panelA_entrez_mapped.csv")
out_map_log <- file.path(breast, "results", "panel_b", "gene_symbol_to_entrez_log.csv")
out_unmapped <- file.path(breast, "results", "panel_b", "gene_symbols_unmapped.csv")
out_pam50_raw <- file.path(breast, "results", "panel_b", "pam50_genefu_raw.csv")
out_report <- file.path(breast, "results", "panel_b", "step0bc_report.txt")

dir.create(file.path(breast, "results", "panel_b"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(breast, "data", "processed"), recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("genefu", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop(
    "R/genefu is not available. Need Bioconductor packages genefu and org.Hs.eg.db.\n",
    "Not substituting another classifier."
  )
}

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(genefu)
})
data(pam50.robust, package = "genefu")

cat("genefu version:", as.character(packageVersion("genefu")), "\n")
cat("org.Hs.eg.db version:", as.character(packageVersion("org.Hs.eg.db")), "\n")
print(args(molecular.subtyping))

expr <- read.csv(in_matrix, row.names = 1, check.names = FALSE)
symbols <- rownames(expr)
cat("Input matrix:", in_matrix, "\n")
cat("Genes x samples:", nrow(expr), "x", ncol(expr), "\n")
cat("Sample IDs kept as-is; n =", ncol(expr), "\n")

mapped <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = symbols,
  columns = "ENTREZID",
  keytype = "SYMBOL"
)
# One row per successful mapping; symbols with no Entrez are NA
map_first <- mapped[!duplicated(mapped$SYMBOL), ]
n_input <- length(symbols)
n_unmapped <- sum(is.na(map_first$ENTREZID))
unmapped_symbols <- map_first$SYMBOL[is.na(map_first$ENTREZID)]

# Multi-map: same symbol to multiple Entrez (select() can return several rows)
n_multi_symbol <- sum(duplicated(mapped$SYMBOL[!is.na(mapped$ENTREZID)]))

ok <- merge(
  data.frame(SYMBOL = symbols, stringsAsFactors = FALSE),
  mapped[!is.na(mapped$ENTREZID), ],
  by = "SYMBOL",
  all.x = FALSE
)
# If a symbol maps to several Entrez IDs, keep the first listed mapping and log it.
ok <- ok[!duplicated(ok$SYMBOL), ]

expr_m <- expr[ok$SYMBOL, , drop = FALSE]
rownames(expr_m) <- ok$ENTREZID

# Duplicate Entrez after mapping: keep the most variable gene (across samples)
dup <- duplicated(rownames(expr_m)) | duplicated(rownames(expr_m), fromLast = TRUE)
n_dup_entrez_genes <- sum(dup)
if (any(dup)) {
  vars <- apply(expr_m, 1, var)
  keep <- tapply(seq_len(nrow(expr_m)), rownames(expr_m), function(i) i[which.max(vars[i])])
  expr_m <- expr_m[unlist(keep), , drop = FALSE]
}

write.csv(expr_m, out_entrez)
write.csv(
  data.frame(
    symbol = map_first$SYMBOL,
    entrez = map_first$ENTREZID,
    mapped = !is.na(map_first$ENTREZID),
    stringsAsFactors = FALSE
  ),
  out_map_log,
  row.names = FALSE
)
write.csv(data.frame(symbol = unmapped_symbols, stringsAsFactors = FALSE), out_unmapped, row.names = FALSE)

cat("Step 0b mapping:\n")
cat("  input genes:", n_input, "\n")
cat("  failed to map (no Entrez):", n_unmapped, "\n")
cat("  mapped unique symbols:", nrow(ok), "\n")
cat("  genes collapsed as duplicate Entrez:", n_dup_entrez_genes, "\n")
cat("  Entrez genes in saved matrix:", nrow(expr_m), "x", ncol(expr_m), "samples\n")
cat("  wrote", out_entrez, "\n")

# PAM50: samples in rows, genes (Entrez) in columns
data_mat <- t(as.matrix(expr_m))
storage.mode(data_mat) <- "double"
annot <- data.frame(
  probe = colnames(data_mat),
  EntrezGene.ID = colnames(data_mat),
  stringsAsFactors = FALSE
)
rownames(annot) <- annot$probe

pam <- molecular.subtyping(
  sbt.model = "pam50",
  data = data_mat,
  annot = annot,
  do.mapping = TRUE,
  verbose = TRUE
)

subtype <- as.character(pam$subtype)
names(subtype) <- names(pam$subtype)
if (is.null(names(subtype)) || any(is.na(names(subtype))) || any(names(subtype) == "")) {
  names(subtype) <- rownames(data_mat)
}

proba <- as.data.frame(pam$subtype.proba)
if (is.null(rownames(proba))) {
  rownames(proba) <- rownames(data_mat)
}
conf <- apply(as.matrix(proba), 1, function(x) max(x, na.rm = TRUE))

out <- data.frame(
  cell_line_genefu = names(subtype),
  pam50_subtype = subtype,
  confidence_score = conf[names(subtype)],
  stringsAsFactors = FALSE
)
proba_aligned <- proba[out$cell_line_genefu, , drop = FALSE]
out <- cbind(out, proba_aligned)
write.csv(out, out_pam50_raw, row.names = FALSE)

failed <- is.na(out$pam50_subtype) | out$pam50_subtype == "" | out$pam50_subtype == "NA"
cat("\nStep 0c PAM50 (genefu::molecular.subtyping, sbt.model='pam50'):\n")
cat("  n samples classified:", nrow(out), "\n")
cat("  subtype counts:\n")
print(table(out$pam50_subtype, useNA = "always"))
cat("  failed/empty calls:", sum(failed), "\n")
if (any(failed)) {
  cat("  failed IDs:\n")
  print(out$cell_line_genefu[failed])
} else {
  cat("  no failed/empty subtype calls\n")
}
cat("  min confidence (max subtype.proba):", min(out$confidence_score, na.rm = TRUE), "\n")
cat("  wrote", out_pam50_raw, "\n")

report <- c(
  paste("genefu", packageVersion("genefu")),
  paste("n_input_genes", n_input),
  paste("n_unmapped", n_unmapped),
  paste("n_entrez_matrix", nrow(expr_m)),
  paste("n_pam50", nrow(out)),
  paste("n_failed_pam50", sum(failed))
)
writeLines(report, out_report)
cat("Wrote", out_report, "\n")
