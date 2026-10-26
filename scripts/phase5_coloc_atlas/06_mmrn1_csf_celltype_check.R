#!/usr/bin/env Rscript
################################################################################
# Phase 5 / step 6: is MMRN1 actually expressed in CSF T cells?
# Checks the Phase 1c "MMRN1 colocalizes in (bulk BLUEPRINT) CD4+ T cells"
# result at single-cell resolution: % of cells with >=1 MMRN1 UMI per cluster
# of the GSE141578 CSF object, alongside T-cell/monocyte markers.
################################################################################
suppressPackageStartupMessages(library(Seurat))
project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
s <- readRDS(file.path(project_dir, "data/GSE141578_processed/GSE141578_CSF_merged_seurat.rds"))
s <- JoinLayers(s)
genes <- intersect(c("MMRN1", "SNCA", "CD3E", "CD4", "CD8A", "CD14"), rownames(s))
cnt <- GetAssayData(s, layer = "counts")[genes, ]
cl <- as.character(s$seurat_clusters)
pct <- sapply(split(seq_len(ncol(s)), cl), function(ix) 100 * Matrix::rowMeans(cnt[, ix, drop = FALSE] > 0))
out <- data.frame(cluster = colnames(pct), n_cells = as.vector(table(cl)[colnames(pct)]), t(round(pct, 2)))
write.csv(out, file.path(project_dir, "data/coloc_atlas/effector_expression/MMRN1_CSF_GSE141578_by_cluster.csv"), row.names = FALSE)
print(out)
cat(sprintf("MMRN1+ cells: %d of %d\n", sum(cnt["MMRN1", ] > 0), ncol(s)))
