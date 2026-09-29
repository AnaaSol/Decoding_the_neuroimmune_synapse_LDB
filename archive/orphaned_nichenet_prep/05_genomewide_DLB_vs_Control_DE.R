#!/usr/bin/env Rscript
################################################################################
# Genome-wide DLB vs Control differential expression in ACC neurons
# Author: Ana (session 2026-08-19)
# Phase: 3 - Neuronal Vulnerability Profiling (extension)
#
# Purpose:
#   The original Phase 3 pipeline only computed expression for a curated list
#   of 16 receptors. NicheNet ligand-activity prediction requires a genuine,
#   genome-wide "geneset of interest" (the transcriptional response to be
#   explained) and a "background" of all expressed genes in the receiver
#   population. This script computes that DE result directly from
#   ACC_neurons_only.rds so NicheNet is not run against a toy/curated gene
#   list.
#
# Input:  data/GSE178146_processed/ACC_neurons_only.rds
# Output: data/GSE178146_processed/neuron_DLB_vs_Control_genomewide_DE.csv
#         data/GSE178146_processed/neuron_background_expressed_genes.txt
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

set.seed(42)

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
phase3_dir  <- file.path(project_dir, "data/GSE178146_processed")

cat("Loading ACC_neurons_only.rds...\n")
neurons <- readRDS(file.path(phase3_dir, "ACC_neurons_only.rds"))
cat(sprintf("  %d neurons, %d genes\n", ncol(neurons), nrow(neurons)))
cat("  condition table:\n"); print(table(neurons$condition))

Idents(neurons) <- "condition"

cat("\nRunning genome-wide Wilcoxon DE: DLB vs Control...\n")
de <- FindMarkers(
  neurons,
  ident.1 = "DLB",
  ident.2 = "Control",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.01,
  only.pos = FALSE
)
de$gene <- rownames(de)
de <- de %>% arrange(p_val_adj)

n_sig <- sum(de$p_val_adj < 0.05, na.rm = TRUE)
cat(sprintf("  %d genes tested, %d significant (padj<0.05)\n", nrow(de), n_sig))

write.csv(de, file.path(phase3_dir, "neuron_DLB_vs_Control_genomewide_DE.csv"), row.names = FALSE)

# Background = expressed genes (detected in >=1% of neurons), for NicheNet
pct_expr <- rowSums(GetAssayData(neurons, layer = "counts") > 0) / ncol(neurons)
background_genes <- names(pct_expr[pct_expr >= 0.01])
writeLines(background_genes, file.path(phase3_dir, "neuron_background_expressed_genes.txt"))
cat(sprintf("  %d background expressed genes (>=1%% of neurons) written\n", length(background_genes)))

cat("\nTop 15 DE genes:\n")
print(head(de[, c("gene","avg_log2FC","p_val_adj")], 15))

cat("\nDone.\n")
