#!/usr/bin/env Rscript
################################################################################
# scRNA-seq Processing with Seurat
# Author: Ana
# Date: 2025-11-25
# Phase: 2 - Immune Cell Profiling
#
# Purpose: Process CSF T cells from GSE161192 using Seurat
#   - Load 10X data
#   - QC and filtering
#   - Normalization and scaling
#   - Dimensionality reduction (PCA, UMAP)
#   - Clustering
#   - Cell type annotation
#
# Input: GSE161192 filtered feature-barcode matrices
# Output: Processed Seurat object with clusters and annotations
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(42)  # Reproducibility

cat("================================================================================\n")
cat("scRNA-seq Processing with Seurat - GSE161192 (CSF T cells)\n")
cat("================================================================================\n\n")

# Set paths
project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
data_dir <- file.path(project_dir, "data/GSE161192_raw")
output_dir <- file.path(project_dir, "data/GSE161192_processed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
samples <- data.frame(
  sample_id = c("P1U", "P1S", "P2U", "P2S"),
  patient = c("P1", "P1", "P2", "P2"),
  condition = c("Unstimulated", "Stimulated", "Unstimulated", "Stimulated"),
  gsm = c("GSM4905012", "GSM4905013", "GSM4905015", "GSM4905014")
)

cat("Samples to process:\n")
print(samples)
cat("\n")

################################################################################
# Step 1: Load 10X data
################################################################################

cat("Step 1: Loading 10X Genomics data...\n")

seurat_objects <- list()

for (i in 1:nrow(samples)) {
  sample_name <- samples$sample_id[i]
  gsm_id <- samples$gsm[i]

  cat(sprintf("  Loading %s (%s)...\n", sample_name, gsm_id))

  # Path to filtered matrix
  data_path <- file.path(data_dir,
                         sprintf("%s_%s_filtered_feature_bc_matrix",
                                gsm_id, sample_name))

  # Load 10X data
  counts <- Read10X(data.dir = data_path)

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = "CSF_Tcells",
    min.cells = 3,
    min.features = 200
  )

  # Add metadata
  seurat_obj$sample <- sample_name
  seurat_obj$patient <- samples$patient[i]
  seurat_obj$condition <- samples$condition[i]

  seurat_objects[[sample_name]] <- seurat_obj

  cat(sprintf("    %d cells, %d features\n",
              ncol(seurat_obj), nrow(seurat_obj)))
}

cat("\n")

################################################################################
# Step 2: QC Metrics and Filtering
################################################################################

cat("Step 2: Calculating QC metrics...\n")

for (sample_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_name]]

  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  # Calculate ribosomal percentage
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

  seurat_objects[[sample_name]] <- seurat_obj
}

# QC thresholds for CSF T cells
# nFeature_max raised from 2500 to 6000: activated/expanded T cells routinely
# express 3000-6000 genes. A ceiling of 2500 removes the most activated clones —
# exactly the population of interest — and produced cell counts ~3x lower than
# reported in Gate et al. 2021 (the source study for this dataset).
nFeature_min <- 200
nFeature_max <- 6000
nCount_max <- 10000
percent_mt_max <- 15

cat(sprintf("\nQC thresholds:\n"))
cat(sprintf("  nFeature: %d - %d\n", nFeature_min, nFeature_max))
cat(sprintf("  nCount: < %d\n", nCount_max))
cat(sprintf("  percent.mt: < %d%%\n", percent_mt_max))
cat("\n")

# Filter cells
for (sample_name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[sample_name]]

  cells_before <- ncol(seurat_obj)

  seurat_obj <- subset(seurat_obj, subset =
                        nFeature_RNA > nFeature_min &
                        nFeature_RNA < nFeature_max &
                        nCount_RNA < nCount_max &
                        percent.mt < percent_mt_max)

  cells_after <- ncol(seurat_obj)

  cat(sprintf("  %s: %d → %d cells (%.1f%% retained)\n",
              sample_name, cells_before, cells_after,
              (cells_after/cells_before)*100))

  seurat_objects[[sample_name]] <- seurat_obj
}

cat("\n")

################################################################################
# Step 3: Merge samples
################################################################################

cat("Step 3: Merging samples...\n")

# Merge all samples
merged_seurat <- merge(
  seurat_objects[[1]],
  y = seurat_objects[2:length(seurat_objects)],
  add.cell.ids = names(seurat_objects),
  project = "CSF_Tcells"
)

cat(sprintf("  Total cells after merge: %d\n", ncol(merged_seurat)))
cat(sprintf("  Total features: %d\n\n", nrow(merged_seurat)))

################################################################################
# Step 4: Normalization and Scaling
################################################################################

cat("Step 4: Normalizing and scaling data...\n")

# Normalize
merged_seurat <- NormalizeData(merged_seurat,
                               normalization.method = "LogNormalize",
                               scale.factor = 10000,
                               verbose = FALSE)

# Find variable features
merged_seurat <- FindVariableFeatures(merged_seurat,
                                      selection.method = "vst",
                                      nfeatures = 2000,
                                      verbose = FALSE)

# Scale data
all_genes <- rownames(merged_seurat)
merged_seurat <- ScaleData(merged_seurat,
                           features = all_genes,
                           vars.to.regress = c("percent.mt", "nCount_RNA"),
                           verbose = FALSE)

cat("  ✓ Normalization complete\n\n")

################################################################################
# Step 5: Dimensionality Reduction
################################################################################

cat("Step 5: Dimensionality reduction...\n")

# PCA
merged_seurat <- RunPCA(merged_seurat,
                        features = VariableFeatures(object = merged_seurat),
                        npcs = 50,
                        verbose = FALSE)

cat("  ✓ PCA complete (50 PCs)\n")

# Determine optimal number of PCs (use elbow plot heuristic)
# For T cells, typically 15-30 PCs is sufficient
n_pcs <- 20
cat(sprintf("  Using %d PCs for downstream analysis\n", n_pcs))

# UMAP
merged_seurat <- RunUMAP(merged_seurat,
                         dims = 1:n_pcs,
                         verbose = FALSE)

cat("  ✓ UMAP complete\n\n")

################################################################################
# Step 6: Clustering
################################################################################

cat("Step 6: Clustering cells...\n")

# Find neighbors
merged_seurat <- FindNeighbors(merged_seurat,
                               dims = 1:n_pcs,
                               verbose = FALSE)

# Find clusters (resolution 0.5 is good for T cell subsets)
merged_seurat <- FindClusters(merged_seurat,
                              resolution = 0.5,
                              verbose = FALSE)

n_clusters <- length(unique(Idents(merged_seurat)))
cat(sprintf("  ✓ Identified %d clusters\n\n", n_clusters))

################################################################################
# Step 7: Cell Type Annotation (T cell markers)
################################################################################

cat("Step 7: Annotating T cell subtypes...\n")

# Key T cell markers
tcell_markers <- list(
  "CD4+ T cells" = c("CD4", "IL7R"),
  "CD8+ T cells" = c("CD8A", "CD8B"),
  "Regulatory T cells" = c("FOXP3", "IL2RA"),
  "Memory T cells" = c("IL7R", "CCR7"),
  "Effector T cells" = c("GZMA", "GZMB", "PRF1"),
  "Activated T cells" = c("CD69", "CD25"),
  "Proliferating" = c("MKI67", "TOP2A")
)

cat("  Calculating marker expression scores...\n")

for (cell_type in names(tcell_markers)) {
  markers <- tcell_markers[[cell_type]]
  # Only use markers that exist in the dataset
  markers_present <- markers[markers %in% rownames(merged_seurat)]

  if (length(markers_present) > 0) {
    score_name <- gsub(" ", "_", gsub("\\+", "pos", cell_type))
    merged_seurat <- AddModuleScore(
      merged_seurat,
      features = list(markers_present),
      name = score_name,
      ctrl = 50
    )
  }
}

cat("  ✓ Cell type scores calculated\n\n")

################################################################################
# Step 8: Save processed object
################################################################################

cat("Step 8: Saving processed Seurat object...\n")

output_file <- file.path(output_dir, "CSF_Tcells_seurat_processed.rds")
saveRDS(merged_seurat, file = output_file)

cat(sprintf("  ✓ Saved to: %s\n\n", output_file))

################################################################################
# Step 9: Generate summary plots
################################################################################

cat("Step 9: Generating summary plots...\n")

# Create plots directory
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE)

# UMAP by cluster
p1 <- DimPlot(merged_seurat, reduction = "umap", label = TRUE) +
  ggtitle("UMAP - Clusters") +
  theme_minimal()

# UMAP by sample
p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "sample") +
  ggtitle("UMAP - Sample") +
  theme_minimal()

# UMAP by condition
p3 <- DimPlot(merged_seurat, reduction = "umap", group.by = "condition") +
  ggtitle("UMAP - Condition") +
  theme_minimal()

# Save plots
ggsave(file.path(plots_dir, "umap_clusters.pdf"), p1, width = 8, height = 6)
ggsave(file.path(plots_dir, "umap_sample.pdf"), p2, width = 8, height = 6)
ggsave(file.path(plots_dir, "umap_condition.pdf"), p3, width = 8, height = 6)

# QC violin plots
p_qc <- VlnPlot(merged_seurat,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0.1)
ggsave(file.path(plots_dir, "qc_metrics.pdf"), p_qc, width = 12, height = 4)

cat("  ✓ Plots saved to:", plots_dir, "\n\n")

################################################################################
# Summary
################################################################################

cat("================================================================================\n")
cat("Processing Complete!\n")
cat("================================================================================\n\n")

cat("Summary:\n")
cat(sprintf("  Total cells: %d\n", ncol(merged_seurat)))
cat(sprintf("  Total features: %d\n", nrow(merged_seurat)))
cat(sprintf("  Number of clusters: %d\n", n_clusters))
cat(sprintf("  Samples: %s\n", paste(unique(merged_seurat$sample), collapse = ", ")))
cat("\nOutput files:\n")
cat(sprintf("  - Seurat object: %s\n", output_file))
cat(sprintf("  - Plots: %s\n", plots_dir))
cat("\nNext steps:\n")
cat("  → Run 02_integrate_tcr_data.R to add TCR information\n")
cat("  → Identify expanded clones in CD4+ T cells\n")
cat("================================================================================\n")
