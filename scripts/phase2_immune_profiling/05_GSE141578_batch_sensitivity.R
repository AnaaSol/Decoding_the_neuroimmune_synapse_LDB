#!/usr/bin/env Rscript
################################################################################
# GSE141578 batch-covariate sensitivity analysis.
#
# The main pseudobulk DESeq2 result (04_GSE141578_disease_vs_HC_pseudobulk.R)
# uses design ~batch + disease_group. This script re-fits the SAME pseudobulk
# matrix (no re-QC/re-clustering -- loads the already-saved CD4+ Seurat object)
# under two alternative designs, to show the CXCR4/CXCL12 non-replication
# result is not an artifact of the batch-adjustment choice:
#   (a) ~disease_group only (no batch covariate at all)
#   (b) PD-only vs HC, restricted to samples where batch is not perfectly
#       confounded with the disease contrast (drops the 2 DLB samples, which
#       both fall in batch2 -- see 04's own documentation of this confound)
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(dplyr)
})

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
output_dir  <- file.path(project_dir, "data/GSE141578_processed")

cat("================================================================================\n")
cat("Phase 2c: GSE141578 batch-covariate sensitivity analysis\n")
cat("================================================================================\n\n")

cd4_cells <- readRDS(file.path(output_dir, "GSE141578_CD4_Tcells_seurat.rds"))
if (packageVersion("Seurat") >= "5.0.0" && length(Layers(cd4_cells[["RNA"]])) > 1) {
  cd4_cells <- JoinLayers(cd4_cells)
}
meta_cd4 <- cd4_cells@meta.data
raw_counts <- GetAssayData(cd4_cells, layer = "counts")

build_pseudobulk <- function(cells_subset_idx) {
  meta_sub <- meta_cd4[cells_subset_idx, ]
  counts_sub <- raw_counts[, cells_subset_idx, drop = FALSE]
  group_ids <- paste(meta_sub$sample_id, meta_sub$disease_group, sep = "__")
  unique_groups <- unique(group_ids)
  pb_matrix <- sapply(unique_groups, function(g) {
    idx <- which(group_ids == g)
    if (length(idx) == 1) counts_sub[, idx] else Matrix::rowSums(counts_sub[, idx])
  })
  pb_matrix <- as.matrix(pb_matrix)
  col_data <- data.frame(
    group_id = unique_groups,
    sample_id = sub("__.*", "", unique_groups),
    disease_group = sub(".*__", "", unique_groups),
    row.names = unique_groups, stringsAsFactors = FALSE
  )
  col_data$batch <- setNames(meta_sub$batch, meta_sub$sample_id)[col_data$sample_id]
  col_data$disease_group <- factor(col_data$disease_group, levels = c("HC", "disease"))
  col_data$batch <- factor(col_data$batch, levels = c("batch1", "batch2"))
  cell_counts <- table(group_ids)
  keep <- names(cell_counts)[cell_counts >= 5]
  list(counts = pb_matrix[, keep, drop = FALSE], coldata = col_data[keep, , drop = FALSE])
}

run_model <- function(pb, design_formula, contrast, label) {
  cat(sprintf("-- %s --\n", label))
  cat(sprintf("   Design: %s\n", deparse(design_formula)))
  cat(sprintf("   N samples: %d (%s)\n", nrow(pb$coldata),
              paste(table(pb$coldata$disease_group), collapse = " vs. ")))
  print(table(pb$coldata$batch, pb$coldata$disease_group))

  dds <- DESeqDataSetFromMatrix(countData = pb$counts, colData = pb$coldata, design = design_formula)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = contrast, alpha = 0.05, pAdjustMethod = "BH")
  res_df <- as.data.frame(res) %>% tibble::rownames_to_column("gene") %>% filter(!is.na(padj))
  n_sig <- sum(res_df$padj < 0.05)
  cat(sprintf("   %d genes tested, %d significant (padj<0.05)\n", nrow(res_df), n_sig))

  for (g in c("CXCR4", "CXCL12")) {
    row <- res_df[res_df$gene == g, ]
    if (nrow(row) == 0) cat(sprintf("   %-8s not detected/tested\n", g))
    else cat(sprintf("   %-8s log2FC=%+.3f padj=%.3g\n", g, row$log2FoldChange, row$padj))
  }
  cat("\n")
  res_df$model <- label
  res_df
}

results_all <- list()

# (a) No batch covariate, full 22-sample design
pb_full <- build_pseudobulk(seq_len(nrow(meta_cd4)))
results_all[["no_batch"]] <- run_model(pb_full, ~disease_group,
                                        c("disease_group", "disease", "HC"),
                                        "Sin covariable de batch (22 muestras)")

# (b) PD-only vs HC (drop the 2 DLB samples, both in batch2)
pd_hc_idx <- which(meta_cd4$condition %in% c("PD", "HC"))
pb_pd <- build_pseudobulk(pd_hc_idx)
batch_x_group <- table(pb_pd$coldata$batch, pb_pd$coldata$disease_group)
use_batch_pd <- all(batch_x_group >= 2)
results_all[["pd_only"]] <- run_model(
  pb_pd, if (use_batch_pd) ~batch + disease_group else ~disease_group,
  c("disease_group", "disease", "HC"),
  sprintf("Solo PD vs. HC (sin DLB, %s)", ifelse(use_batch_pd, "con batch", "sin batch -- insuficiente"))
)

combined <- do.call(rbind, results_all)
out_file <- file.path(output_dir, "pseudobulk_DEG_batch_sensitivity.csv")
write.csv(combined, out_file, row.names = FALSE)
cat(sprintf("Guardado: %s\n\n", out_file))

cat("--- R Session Info ---\n")
print(sessionInfo())
