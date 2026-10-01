#!/usr/bin/env Rscript
# Extract per-cell barcode / sample / cluster metadata from the Pal et al.
# 2021 (GSE161529) companion Seurat objects, published alongside:
#   Chen, Pal, Lindeman, Visvader, Smyth (2022) Scientific Data 9:96
#   "R code and downstream analysis objects for the scRNA-seq atlas of
#    normal and tumorigenic human breast tissue"
#   Figshare DOI: 10.6084/m9.figshare.17058077
#
# Output is shaped to play the same role GSE176078's metadata.csv played
# in your existing 1run_panelA_ks_project_gse176078_sc.py script: one row
# per cell, with a barcode you can match against the raw count matrix and
# a subtype/cluster label to filter and color by.
#
# BEFORE RUNNING: download these three files from the Figshare deposit and
# put them in RDS_DIR (see chat for the exact link):
#   SeuratObject_TNBCTum.rds
#   SeuratObject_HER2Tum.rds
#   SeuratObject_ERTotalTum.rds
# (these are already restricted to the epithelial/tumour cluster within
#  each subtype -- see Table 1 of the Scientific Data paper -- so no
#  separate "Cancer Epithelial" filter is needed here, unlike GSE176078.)
#
# Usage (repo root):
#   Rscript "Breast Cancer/codes/00_extract_metadata_gse161529.R"

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
codes_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
RDS_DIR <- normalizePath(file.path(codes_dir, "..", "data", "gse161529"), mustWork = TRUE)
OUT_CSV <- file.path(RDS_DIR, "metadata_gse161529.csv")

objects <- list(
  list(file = "SeuratObject_TNBCTum.rds",    subtype = "TNBC"),
  list(file = "SeuratObject_HER2Tum.rds",    subtype = "HER2+"),
  list(file = "SeuratObject_ERTotalTum.rds", subtype = "ER+")
)

extract_one <- function(path, subtype) {
  if (!file.exists(path)) {
    stop("Missing ", path, " -- download it from the Figshare deposit first.")
  }
  message("Reading ", path)
  so <- readRDS(path)
  # Objects were built under Seurat v3.1.1; UpdateSeuratObject() makes them
  # readable with whatever Seurat version you have installed now.
  so <- tryCatch(UpdateSeuratObject(so), error = function(e) so)

  md <- so@meta.data
  full_barcode <- rownames(md)

  # Cell names are "<SampleComb>_<10x barcode>", e.g.
  # "TN_B1_0554_AAACCTGAGCTAGTGG-1". Sample names themselves contain
  # underscores, so recover the barcode by matching the trailing 16bp
  # ACGT token (+ optional "-N" suffix) rather than splitting on "_".
  bc_pattern <- "[ACGT]{16}(-[0-9]+)?$"
  raw_barcode <- regmatches(full_barcode, regexpr(bc_pattern, full_barcode))
  sample_comb <- substr(full_barcode, 1, nchar(full_barcode) - nchar(raw_barcode) - 1)

  data.frame(
    full_barcode = full_barcode,
    raw_barcode  = raw_barcode,
    # back to the GEO / SampleStats.txt form, e.g. "TN-B1-0554"
    sample_name  = gsub("_", "-", sample_comb),
    subtype      = subtype,
    # cluster id *within this Tum object only* (resolution 0.1, 0-indexed).
    # There is no published LumA/LumB/Basal/Cycling-style label for
    # GSE161529 the way GSE176078 has celltype_minor -- treat this as a
    # coarse within-subtype grouping, not a named cell state.
    cluster      = as.integer(as.character(md$seurat_clusters)),
    group        = md$group,
    stringsAsFactors = FALSE
  )
}

all_md <- do.call(rbind, lapply(objects, function(o) {
  extract_one(file.path(RDS_DIR, o$file), o$subtype)
}))

message(sprintf("Total cells: %d", nrow(all_md)))
print(table(all_md$subtype))
print(table(all_md$sample_name, all_md$subtype))

dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(all_md, OUT_CSV, row.names = FALSE)
message("Wrote ", OUT_CSV)
