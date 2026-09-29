#!/usr/bin/env Rscript
################################################################################
# GSE141578 — CSF T Cell Disease vs. HC Pseudobulk DESeq2 Analysis
# Author: Ana
# Date: 2026-07-17
# Phase: 2 - Immune Cell Profiling (extended dataset)
#
# Purpose:
#   Gate et al. 2021 (Science 374:868) released a larger CSF cohort (GSE141578)
#   alongside the smaller stimulation dataset (GSE161192). GSE141578 provides
#   22 CSF samples: 9 PD + 2 DLB + 11 HC — sufficient for pseudobulk DE that
#   failed on GSE161192 (N=2 donors, defect D7/D8).
#
#   This script:
#     1. Loads all CSF samples from data/GSE141578_raw/csf_10x/
#     2. QC, normalization, integration, UMAP, clustering
#     3. CD4+ T cell identification via canonical markers
#     4. Pseudobulk DESeq2: disease (PD + DLB combined) vs. HC
#     5. Exports top DEGs and compares to manuscript claims
#
# Input:  data/GSE141578_raw/csf_10x/{CSF*}/ — standard 10x directories
# Output: data/GSE141578_processed/
#           pseudobulk_DEG_disease_vs_HC.csv
#           cd4_disease_upregulated_genes.txt
#           plots/
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(42)

cat("================================================================================\n")
cat("GSE141578 CSF T cells — Disease vs. HC Pseudobulk Analysis\n")
cat("================================================================================\n\n")

project_dir   <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
data_dir      <- file.path(project_dir, "data/GSE141578_raw/csf_10x")
output_dir    <- file.path(project_dir, "data/GSE141578_processed")
plots_dir     <- file.path(output_dir, "plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)

################################################################################
# Sample metadata — verified per-sample against individual GEO GSM records
# (GSE141578, Gate et al. 2021 Science; diagnosis: field on each GSM page)
# HC: CSF1 CSF3 CSF4 CSF8 CSF12 CSF13 CSF16 CSF20 CSF23 CSF28 CSF32 (n=11)
# PD: CSF7 CSF10 CSF14 CSF21 CSF24 CSF26 CSF27 CSF29 CSF30 (n=9)
# DLB: CSF25 CSF31 (n=2)
# Series design confirms n=11 HC / n=11 disease (PD+DLB) across 22 CSF samples.
# NOTE: CSF32 was previously excluded here as "not reported" — that was wrong;
# GSM4404059 explicitly lists "diagnosis: healthy control". Corrected 2026-08-18.
################################################################################

sample_meta <- data.frame(
  sample_id = c(
    "CSF1",  "CSF3",  "CSF4",  "CSF8",  "CSF12",
    "CSF13", "CSF16", "CSF20", "CSF23", "CSF28", "CSF32",
    "CSF7",  "CSF10", "CSF14", "CSF21", "CSF24",
    "CSF26", "CSF27", "CSF29", "CSF30",
    "CSF25", "CSF31"
  ),
  condition = c(
    rep("HC", 11),
    rep("PD", 9),
    rep("DLB", 2)
  ),
  stringsAsFactors = FALSE
)

################################################################################
# Technical batch covariate — GEO has no age/sex fields for this series, but
# Sample_submission_date splits cleanly into two library-prep batches:
#   batch1 = GSM4208xxx series, submitted Dec 06 2019 (14 samples)
#   batch2 = GSM4404xxx series, submitted Mar 09 2020 (8 samples)
# Verified per-GSM via individual GEO records (2026-08-19). Not perfectly
# confounded with disease_group (each batch has HC + disease samples), but
# both DLB samples (CSF25, CSF31) fall in batch2 — DLB-specific effects
# cannot be fully separated from batch2 with this design.
################################################################################

batch1_samples <- c("CSF14","CSF10","CSF21","CSF24","CSF7",
                     "CSF1","CSF3","CSF4","CSF8","CSF12","CSF13","CSF23","CSF20","CSF16")
sample_meta$batch <- ifelse(sample_meta$sample_id %in% batch1_samples, "batch1", "batch2")

# For pseudobulk DE: combine PD + DLB into "disease" group
sample_meta$disease_group <- ifelse(sample_meta$condition == "HC", "HC", "disease")

cat(sprintf("Samples: %d HC | %d PD | %d DLB\n",
            sum(sample_meta$condition == "HC"),
            sum(sample_meta$condition == "PD"),
            sum(sample_meta$condition == "DLB")))
cat(sprintf("Pseudobulk design: %d disease vs. %d HC\n\n",
            sum(sample_meta$disease_group == "disease"),
            sum(sample_meta$disease_group == "HC")))

################################################################################
# Step 1: Load samples
################################################################################

cat("Step 1: Loading samples...\n")

seurat_list <- list()

for (i in seq_len(nrow(sample_meta))) {
  sid   <- sample_meta$sample_id[i]
  spath <- file.path(data_dir, sid)

  if (!dir.exists(spath)) {
    cat(sprintf("  SKIP %s — directory not found\n", sid))
    next
  }

  counts <- tryCatch(
    Read10X(data.dir = spath),
    error = function(e) { cat(sprintf("  ERROR reading %s: %s\n", sid, e$message)); NULL }
  )
  if (is.null(counts)) next

  obj <- CreateSeuratObject(counts = counts, project = sid,
                            min.cells = 3, min.features = 200)
  obj$sample_id    <- sid
  obj$condition    <- sample_meta$condition[i]
  obj$disease_group <- sample_meta$disease_group[i]
  obj$batch        <- sample_meta$batch[i]

  seurat_list[[sid]] <- obj
  cat(sprintf("  %s (%s): %d cells\n", sid, sample_meta$condition[i], ncol(obj)))
}

cat(sprintf("\n  Loaded %d/%d samples\n\n", length(seurat_list), nrow(sample_meta)))

################################################################################
# Step 2: QC and filtering
################################################################################

cat("Step 2: QC filtering...\n")

qc_thresholds <- list(
  nFeature_min  = 200,
  nFeature_max  = 6000,
  nCount_max    = 50000,
  percent_mt_max = 15   # CSF cells; slightly more permissive than brain nuclei
)

for (sid in names(seurat_list)) {
  obj <- seurat_list[[sid]]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  n_before <- ncol(obj)
  obj <- subset(obj,
                subset = nFeature_RNA >= qc_thresholds$nFeature_min &
                         nFeature_RNA <= qc_thresholds$nFeature_max &
                         nCount_RNA   <= qc_thresholds$nCount_max   &
                         percent.mt   <  qc_thresholds$percent_mt_max)
  n_after <- ncol(obj)
  cat(sprintf("  %s: %d → %d cells (%.1f%% kept)\n",
              sid, n_before, n_after, 100 * n_after / n_before))
  seurat_list[[sid]] <- obj
}

# Drop samples with fewer than 10 cells post-QC
n_before_drop <- length(seurat_list)
seurat_list <- Filter(function(obj) ncol(obj) >= 10, seurat_list)
if (length(seurat_list) < n_before_drop) {
  cat(sprintf("\n  Dropped %d sample(s) with <10 cells post-QC\n",
              n_before_drop - length(seurat_list)))
}
cat(sprintf("\n  Retained %d samples for downstream analysis\n\n", length(seurat_list)))

################################################################################
# Step 3: Normalize, find variable features, then merge
################################################################################

cat("Step 3: Normalizing and merging...\n")

seurat_list <- lapply(seurat_list, function(obj) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 3000, verbose = FALSE)
  obj
})

merged <- merge(
  seurat_list[[1]],
  y       = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list),
  project = "GSE141578_CSF"
)

cat(sprintf("  Merged: %d cells, %d features\n\n", ncol(merged), nrow(merged)))

################################################################################
# Step 4: Scale, PCA, UMAP, clustering
################################################################################

cat("Step 4: Dimensionality reduction and clustering...\n")

# Use variable features from merged object
merged <- FindVariableFeatures(merged, selection.method = "vst",
                               nfeatures = 3000, verbose = FALSE)
merged <- ScaleData(merged,
                    vars.to.regress = c("percent.mt", "nCount_RNA"),
                    verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)

# Harmony batch correction by sample (prevents sample-driven clustering)
harmony_available <- requireNamespace("harmony", quietly = TRUE)
if (harmony_available) {
  library(harmony)
  merged <- RunHarmony(merged, group.by.vars = "sample_id", verbose = FALSE)
  reduction_use <- "harmony"
  cat("  ✓ Harmony batch correction applied\n")
} else {
  cat("  NOTE: harmony not installed — using PCA (sample effects not corrected)\n")
  reduction_use <- "pca"
}

n_dims <- 20
merged <- RunUMAP(merged, reduction = reduction_use, dims = 1:n_dims, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = reduction_use, dims = 1:n_dims, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.5, verbose = FALSE)

n_clusters <- length(unique(Idents(merged)))
cat(sprintf("  ✓ %d clusters identified\n\n", n_clusters))

################################################################################
# Step 5: CD4+ T cell identification
################################################################################

cat("Step 5: Identifying CD4+ T cells...\n")

# Score each cluster for CD4+ T cell markers
cd4_markers  <- c("CD4", "IL7R", "CCR7", "CD3D", "CD3E")
cd8_markers  <- c("CD8A", "CD8B")
nk_markers   <- c("GNLY", "NKG7", "KLRD1")
mono_markers <- c("LYZ", "CD14", "FCGR3A")

marker_lists <- list(
  CD4pos   = cd4_markers[cd4_markers   %in% rownames(merged)],
  CD8pos   = cd8_markers[cd8_markers   %in% rownames(merged)],
  NK       = nk_markers[nk_markers     %in% rownames(merged)],
  Monocyte = mono_markers[mono_markers %in% rownames(merged)]
)

for (ct in names(marker_lists)) {
  if (length(marker_lists[[ct]]) > 0) {
    merged <- AddModuleScore(merged, features = list(marker_lists[[ct]]),
                             name = ct, ctrl = 50)
  }
}

# Identify CD4+ clusters: higher CD4pos score than CD8pos, and CD4pos > 0
cluster_scores <- merged@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    cd4_score  = mean(get("CD4pos1"), na.rm = TRUE),
    cd8_score  = mean(get("CD8pos1"), na.rm = TRUE),
    .groups = "drop"
  )

cd4_clusters <- cluster_scores %>%
  filter(cd4_score > 0 & cd4_score > cd8_score) %>%
  pull(seurat_clusters)

cat(sprintf("  CD4+ clusters: %s\n", paste(cd4_clusters, collapse = ", ")))

Idents(merged) <- "seurat_clusters"
if (length(cd4_clusters) == 0) {
  cat("  WARNING: No CD4+ clusters identified by score threshold.\n")
  cat("  Falling back to all T cells (CD4pos score > 0 in any cluster).\n")
  cd4_clusters <- cluster_scores %>%
    filter(cd4_score > cd8_score) %>%
    pull(seurat_clusters)
}

cd4_cells <- subset(merged, idents = cd4_clusters)
cat(sprintf("  CD4+ T cells: %d / %d total cells (%.1f%%)\n\n",
            ncol(cd4_cells), ncol(merged), 100 * ncol(cd4_cells) / ncol(merged)))

################################################################################
# Step 6: Pseudobulk DESeq2 — disease vs. HC in CD4+ T cells
################################################################################

cat("Step 6: Pseudobulk DESeq2 (disease vs. HC in CD4+ T cells)...\n")

pseudobulk_results <- NULL

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  cat("  WARNING: DESeq2 not installed.\n")
  cat("  Install: BiocManager::install('DESeq2')\n\n")
} else {
  suppressPackageStartupMessages(library(DESeq2))

  meta_cd4 <- cd4_cells@meta.data
  required_cols <- c("sample_id", "disease_group")
  missing_cols <- setdiff(required_cols, colnames(meta_cd4))

  if (length(missing_cols) > 0) {
    cat(sprintf("  ERROR: Missing columns: %s\n\n", paste(missing_cols, collapse=", ")))
  } else {

    # ------------------------------------------------------------------
    # Build pseudobulk matrix: sum raw counts per (sample × disease_group)
    # ------------------------------------------------------------------
    # merge() keeps one "counts" layer per sample in SeuratObject v5;
    # GetAssayData() cannot pull from multiple layers at once, so join first.
    cd4_cells   <- JoinLayers(cd4_cells)
    raw_counts  <- GetAssayData(cd4_cells, layer = "counts")
    group_ids   <- paste(meta_cd4$sample_id, meta_cd4$disease_group, sep = "__")
    unique_groups <- unique(group_ids)

    pb_matrix <- sapply(unique_groups, function(g) {
      cell_idx <- which(group_ids == g)
      if (length(cell_idx) == 1) raw_counts[, cell_idx]
      else Matrix::rowSums(raw_counts[, cell_idx])
    })
    pb_matrix <- as.matrix(pb_matrix)

    # Build colData
    sample_to_batch <- setNames(meta_cd4$batch, meta_cd4$sample_id)
    col_data <- data.frame(
      group_id     = unique_groups,
      sample_id    = sub("__.*", "", unique_groups),
      disease_group = sub(".*__", "", unique_groups),
      row.names    = unique_groups,
      stringsAsFactors = FALSE
    )
    col_data$batch <- sample_to_batch[col_data$sample_id]
    col_data$disease_group <- factor(col_data$disease_group,
                                     levels = c("HC", "disease"))
    col_data$batch <- factor(col_data$batch, levels = c("batch1", "batch2"))

    # Drop groups with < 5 cells
    cell_counts_per_group <- table(group_ids)
    keep_groups <- names(cell_counts_per_group)[cell_counts_per_group >= 5]
    pb_matrix <- pb_matrix[, keep_groups, drop = FALSE]
    col_data  <- col_data[keep_groups, , drop = FALSE]

    n_disease <- sum(col_data$disease_group == "disease")
    n_hc      <- sum(col_data$disease_group == "HC")

    cat(sprintf("  Pseudobulk groups: %d disease, %d HC (min 5 cells/group)\n",
                n_disease, n_hc))

    if (n_disease < 2 || n_hc < 2) {
      cat(sprintf("  Insufficient replication (need ≥2 per group) — skipping DESeq2.\n\n"))
    } else {

      # Use batch as a covariate only if both batches have >=2 samples in
      # each disease_group cell after the min-cell filter, otherwise the
      # model is unidentifiable (perfect confounding).
      batch_x_group <- table(col_data$batch, col_data$disease_group)
      use_batch <- all(batch_x_group >= 2)
      design_formula <- if (use_batch) ~ batch + disease_group else ~ disease_group
      cat(sprintf("  Batch x disease_group table:\n"))
      print(batch_x_group)
      cat(sprintf("  Design: %s (%s)\n",
                  deparse(design_formula),
                  ifelse(use_batch, "batch covariate included",
                         "batch dropped — insufficient cells per batch x group cell")))

      tryCatch({
        dds <- DESeqDataSetFromMatrix(
          countData = pb_matrix,
          colData   = col_data,
          design    = design_formula
        )

        # Minimal pre-filtering: keep genes with at least 10 total reads
        dds <- dds[rowSums(counts(dds)) >= 10, ]

        dds <- DESeq(dds, quiet = TRUE)
        res <- results(dds,
                       contrast   = c("disease_group", "disease", "HC"),
                       alpha      = 0.05,
                       pAdjustMethod = "BH")
        res_df <- as.data.frame(res) %>%
          tibble::rownames_to_column("gene") %>%
          filter(!is.na(padj)) %>%
          arrange(padj, desc(abs(log2FoldChange)))

        n_sig   <- sum(res_df$padj < 0.05, na.rm = TRUE)
        n_up    <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
        n_down  <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)

        cat(sprintf("  ✓ DESeq2 complete: %d significant genes (padj<0.05)\n", n_sig))
        cat(sprintf("     Upregulated in disease: %d\n", n_up))
        cat(sprintf("     Downregulated in disease: %d\n\n", n_down))

        if (n_sig > 0) {
          cat("  Top 15 upregulated in disease (CD4+ T cells):\n")
          top_up <- res_df %>% filter(padj < 0.05, log2FoldChange > 0) %>% head(15)
          for (j in seq_len(nrow(top_up))) {
            cat(sprintf("    %-15s log2FC=%.2f padj=%.2e\n",
                        top_up$gene[j], top_up$log2FoldChange[j], top_up$padj[j]))
          }
          cat("\n")
        }

        # Check manuscript claims
        cat("  Manuscript gene claims (CXCR4, key cytokines):\n")
        claim_genes <- c("CXCR4", "CXCL12", "IL17A", "IFNG", "TNF",
                         "IL2", "GZMB", "PRF1", "LAG3", "HAVCR2")
        for (g in claim_genes) {
          row <- res_df[res_df$gene == g, ]
          if (nrow(row) == 0) {
            cat(sprintf("    %-10s not detected\n", g))
          } else {
            sig_flag <- ifelse(!is.na(row$padj) & row$padj < 0.05, "***", "ns")
            cat(sprintf("    %-10s log2FC=%+.2f padj=%.2e %s\n",
                        g, row$log2FoldChange, row$padj, sig_flag))
          }
        }
        cat("\n")

        # Save results
        write.csv(res_df,
                  file.path(output_dir, "pseudobulk_DEG_disease_vs_HC.csv"),
                  row.names = FALSE)

        # Upregulated gene list for downstream use (Phase 4 ligand list)
        up_genes <- res_df %>%
          filter(padj < 0.05, log2FoldChange > 0) %>%
          pull(gene)
        writeLines(up_genes,
                   file.path(output_dir, "cd4_disease_upregulated_genes.txt"))
        cat(sprintf("  ✓ Saved: pseudobulk_DEG_disease_vs_HC.csv (%d genes)\n", nrow(res_df)))
        cat(sprintf("  ✓ Saved: cd4_disease_upregulated_genes.txt (%d upregulated)\n\n", length(up_genes)))

        pseudobulk_results <- res_df

      }, error = function(e) {
        cat(sprintf("  ERROR in DESeq2: %s\n\n", e$message))
      })
    }
  }
}

################################################################################
# Step 7: Save Seurat objects and plots
################################################################################

cat("Step 7: Saving objects and plots...\n")

saveRDS(merged,    file.path(output_dir, "GSE141578_CSF_merged_seurat.rds"))
saveRDS(cd4_cells, file.path(output_dir, "GSE141578_CD4_Tcells_seurat.rds"))
cat("  ✓ Seurat objects saved\n")

# UMAP plots
p1 <- DimPlot(merged, reduction = "umap", label = TRUE) +
  ggtitle("UMAP — Clusters (all CSF cells)") + theme_minimal()
p2 <- DimPlot(merged, reduction = "umap", group.by = "disease_group") +
  ggtitle("UMAP — Disease group") + theme_minimal()
p3 <- DimPlot(merged, reduction = "umap", group.by = "condition") +
  ggtitle("UMAP — Condition (HC/PD/DLB)") + theme_minimal()
ggsave(file.path(plots_dir, "umap_overview.pdf"),
       p1 | p2 | p3, width = 18, height = 5)

p4 <- DimPlot(cd4_cells, reduction = "umap", group.by = "disease_group") +
  ggtitle("UMAP — CD4+ T cells") + theme_minimal()
ggsave(file.path(plots_dir, "umap_cd4_cells.pdf"), p4, width = 7, height = 6)

# Volcano plot if DESeq2 ran
if (!is.null(pseudobulk_results)) {
  vol_df <- pseudobulk_results %>%
    mutate(
      sig = padj < 0.05 & abs(log2FoldChange) > 0.5,
      label = ifelse(sig & rank(padj) <= 20, gene, "")
    )
  p_vol <- ggplot(vol_df, aes(x = log2FoldChange, y = -log10(padj),
                              color = sig, label = label)) +
    geom_point(alpha = 0.6, size = 1) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("grey60", "#e74c3c")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.4) +
    labs(title = "Pseudobulk DESeq2: disease vs. HC (CD4+ T cells)",
         subtitle = "GSE141578 — n=11 disease, n=11 HC",
         x = "log2 fold change", y = "-log10(padj)") +
    theme_minimal() + theme(legend.position = "none")
  tryCatch(
    ggsave(file.path(plots_dir, "volcano_disease_vs_HC.pdf"), p_vol, width = 8, height = 6),
    error = function(e) cat("  NOTE: ggrepel not available, skipping volcano labels\n")
  )
}

cat("  ✓ Plots saved\n\n")

################################################################################
# Summary
################################################################################

cat("================================================================================\n")
cat("Analysis Complete\n")
cat("================================================================================\n\n")
cat(sprintf("  Samples loaded:    %d\n", length(seurat_list)))
cat(sprintf("  Total cells:       %d\n", ncol(merged)))
cat(sprintf("  CD4+ T cells:      %d\n", ncol(cd4_cells)))
cat(sprintf("  Pseudobulk DESeq2: %s\n",
            ifelse(is.null(pseudobulk_results), "NOT RUN", "COMPLETE")))
cat("\nOutput files:\n")
cat(sprintf("  %s/\n", output_dir))
cat("    GSE141578_CSF_merged_seurat.rds\n")
cat("    GSE141578_CD4_Tcells_seurat.rds\n")
cat("    pseudobulk_DEG_disease_vs_HC.csv\n")
cat("    cd4_disease_upregulated_genes.txt\n")
cat("    plots/\n")
cat("\nNext steps:\n")
cat("  → Phase 4: use cd4_disease_upregulated_genes.txt as ligand list\n")
cat("  → Compare DEGs with GSE161192 Wilcoxon results\n")
cat("================================================================================\n")
