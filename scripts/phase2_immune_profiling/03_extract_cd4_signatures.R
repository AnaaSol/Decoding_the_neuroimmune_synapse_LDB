#!/usr/bin/env Rscript
################################################################################
# Extract Gene Signatures from Expanded CD4+ T Cells (CLEAN REBUILD FIX)
# Author: Ana
# Date: 2025-11-25
# Phase: 2 - Immune Cell Profiling
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Matrix) # Necesario para manejo de matrices dispersas
})

set.seed(42)  # Reproducibility

cat("================================================================================\n")
cat("Gene Signature Extraction - Expanded CD4+ T Cells\n")
cat("================================================================================\n\n")

# Set paths
project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
processed_dir <- file.path(project_dir, "data/GSE161192_processed")
phase1_dir <- file.path(project_dir, "data")
output_dir <- processed_dir
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Step 1: Load data
################################################################################

cat("Step 1: Loading integrated Seurat object...\n")

seurat_file <- file.path(processed_dir, "CSF_Tcells_with_TCR.rds")
if (!file.exists(seurat_file)) {
  stop("Error: Integrated Seurat object not found. Run 02_integrate_tcr_data.R first.")
}

seurat_obj <- readRDS(seurat_file)
cat(sprintf("  ✓ Loaded %d cells\n\n", ncol(seurat_obj)))

# Load LBD risk genes from Phase 1
cat("Loading LBD risk genes from Phase 1...\n")
lbd_genes_file <- file.path(phase1_dir, "LBD_risk_genes_ranked.tsv")

if (file.exists(lbd_genes_file)) {
  lbd_genes <- read.table(lbd_genes_file, header = TRUE, sep = "\t")
  cat(sprintf("  ✓ Loaded %d LBD risk genes\n", nrow(lbd_genes)))
} else {
  cat("  Warning: LBD risk genes file not found\n")
  lbd_genes <- NULL
}

cat("\n")

################################################################################
# Step 2: Subsetting to CD4+ T cells (MANUAL STITCHING WITH GENE ALIGNMENT)
################################################################################

cat("Step 2: Subsetting to CD4+ T cells...\n")

# 1. Identificar células y hacer subset
cd4_cells_ids <- colnames(seurat_obj)[seurat_obj$is_CD4]
seurat_subset <- subset(seurat_obj, cells = cd4_cells_ids)

cat("  ✓ Initial subset complete. Aligning and stitching layers manually...\n")

# A. Asegurar que miramos el ensayo RNA
DefaultAssay(seurat_subset) <- "RNA"
rna_assay <- seurat_subset[["RNA"]]

# B. Buscar capas de conteo
all_layers <- Layers(rna_assay)
count_layers <- all_layers[grep("counts", all_layers)]

cat(sprintf("    Found %d count layers: %s\n", length(count_layers), paste(count_layers, collapse=", ")))

if (length(count_layers) == 0) {
    stop("CRITICAL ERROR: No 'counts' layers found in RNA assay.")
}

# C. Extraer matrices y encontrar genes comunes
list_of_matrices <- list()
all_gene_lists <- list()

for (layer_name in count_layers) {
    # Extraer matriz
    mat <- LayerData(seurat_subset, assay = "RNA", layer = layer_name)
    list_of_matrices[[layer_name]] <- mat
    all_gene_lists[[layer_name]] <- rownames(mat)
}

# Encontrar la intersección de genes (genes presentes en TODAS las capas)
common_genes <- Reduce(intersect, all_gene_lists)

cat(sprintf("    Genes in layer 1: %d\n", length(all_gene_lists[[1]])))
cat(sprintf("    Common genes across all layers: %d\n", length(common_genes)))

if (length(common_genes) < 100) {
    stop("Error: Too few common genes found (<100). Check if samples are from same organism/annotation.")
}

# D. Recortar matrices a los genes comunes y Unir
# Esto soluciona el error "number of rows must match"
cat("    Aligning matrices to common genes...\n")
aligned_matrices <- lapply(list_of_matrices, function(x) {
    return(x[common_genes, ])
})

# Unir columnas (Ahora sí funcionará porque las filas son idénticas)
clean_counts <- do.call(cbind, aligned_matrices)

# E. Alinear Metadatos
clean_meta <- seurat_subset@meta.data
cells_in_meta <- rownames(clean_meta)
cells_in_counts <- colnames(clean_counts)
common_cells <- intersect(cells_in_meta, cells_in_counts)

# Reordenar para asegurar coincidencia perfecta
clean_counts <- clean_counts[, common_cells]
clean_meta <- clean_meta[common_cells, ]

# 3. Crear objeto Seurat NUEVO y LIMPIO
seurat_cd4 <- CreateSeuratObject(counts = clean_counts, meta.data = clean_meta)
DefaultAssay(seurat_cd4) <- "RNA"

cat("  ✓ Object fully rebuilt with aligned genes.\n")

# 4. Normalizar
cat("  ✓ Running Normalization and Scaling on clean object...\n")
seurat_cd4 <- NormalizeData(seurat_cd4, verbose = FALSE)
seurat_cd4 <- FindVariableFeatures(seurat_cd4, verbose = FALSE)
seurat_cd4 <- ScaleData(seurat_cd4, verbose = FALSE)

cat(sprintf("  ✓ %d CD4+ T cells ready for analysis\n", ncol(seurat_cd4)))

# Count expanded vs non-expanded
n_expanded <- sum(seurat_cd4$is_expanded, na.rm = TRUE)
n_not_expanded <- sum(!seurat_cd4$is_expanded, na.rm = TRUE)

cat(sprintf("  Expanded: %d cells\n", n_expanded))
cat(sprintf("  Non-expanded: %d cells\n\n", n_not_expanded))

################################################################################
# Step 3: Differential expression analysis
################################################################################

cat("Step 3: Differential expression: Expanded vs Non-expanded CD4+ T cells...\n")

# Set identity to expansion status
Idents(seurat_cd4) <- seurat_cd4$is_expanded

deg_markers <- NULL

if (n_expanded >= 3 && n_not_expanded >= 3) {
  
  tryCatch({
      cat("  Running FindMarkers (Wilcoxon)...\n")
      
      deg_markers <- FindMarkers(
        seurat_cd4,
        ident.1 = "TRUE",   # Expanded
        ident.2 = "FALSE",  # Non-expanded
        test.use = "wilcox",
        logfc.threshold = 0.1, 
        min.pct = 0.05,        
        verbose = FALSE
      )
      
      if (nrow(deg_markers) == 0) {
          cat("  ⚠ No genes passed thresholds.\n")
          deg_markers <- NULL
      } else {
          deg_markers$gene <- rownames(deg_markers)
          deg_markers <- deg_markers[order(deg_markers$p_val_adj), ]
          
          n_sig <- sum(deg_markers$p_val_adj < 0.05, na.rm=TRUE)
          cat(sprintf("  ✓ Found %d differentially expressed genes (padj < 0.05)\n", n_sig))
    
          # Top upregulated genes in expanded cells
          top_up <- head(deg_markers[deg_markers$avg_log2FC > 0, ], 20)
          cat("\n  Top 10 upregulated genes in expanded CD4+ cells:\n")
          if(nrow(top_up) > 0) {
             print(top_up[1:min(10, nrow(top_up)), c("gene", "avg_log2FC", "p_val_adj")])
          }
    
          # Save DEG results
          write.csv(deg_markers,
                    file.path(output_dir, "DEG_expanded_vs_nonexpanded_CD4.csv"),
                    row.names = FALSE)
      }

  }, error = function(e) {
      cat("  Error during FindMarkers: ", e$message, "\n")
      deg_markers <- NULL
  })

} else {
  cat("  Warning: Not enough cells for differential expression analysis\n")
}

cat("\n")

################################################################################
# Step 3b: Pseudobulk DE (DESeq2) — donor-aware, stimulation-deconfounded
#
# Rationale: Step 3 above uses Wilcoxon on individual cells. With N=2 donors
# (P1, P2) this is pseudo-replication: p-values treat ~4000 cells as independent
# observations, grossly inflating significance (padj ~ 1e-45 for common genes).
#
# Pseudobulk approach: aggregate raw counts per (donor × stimulation × expansion),
# run DESeq2 with design ~stimulation + expansion_status to control for the
# α-syn stimulation arm in GSE161192.
#
# CAVEAT: N=2 donors → 8 pseudobulk groups max, ~5 df for residuals.
# Statistical power is very low; treat these results as directional only and
# validate in an independent cohort.
################################################################################

cat("Step 3b: Pseudobulk DE (DESeq2) — donor-aware, stimulation-deconfounded...\n")

pseudobulk_markers <- NULL

tryCatch({

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    cat("  WARNING: DESeq2 not installed — skipping pseudobulk analysis.\n")
    cat("  Install with: BiocManager::install('DESeq2')\n\n")

  } else {
    suppressPackageStartupMessages(library(DESeq2))

    required_cols <- c("patient", "condition", "is_expanded")
    missing_cols <- setdiff(required_cols, colnames(seurat_cd4@meta.data))
    if (length(missing_cols) > 0) {
      cat(sprintf("  WARNING: Missing metadata columns for pseudobulk: %s\n",
                  paste(missing_cols, collapse = ", ")))
      cat("  Ensure 01_process_scrna_seurat.R was run with patient/condition metadata.\n\n")

    } else {

      # --- 1. Build pseudobulk groups ---
      valid_mask <- !is.na(seurat_cd4$is_expanded)
      meta_pb    <- seurat_cd4@meta.data[valid_mask, ]

      meta_pb$expansion_status <- ifelse(meta_pb$is_expanded, "Expanded", "NonExpanded")
      # 'condition' here = stimulation status (Stimulated / Unstimulated)
      meta_pb$pb_group <- paste(meta_pb$patient,
                                sub("stimulated", "Stim",
                                    sub("Unstimulated", "Unstim", meta_pb$condition)),
                                meta_pb$expansion_status, sep = "_")

      # --- 2. Aggregate raw counts per group ---
      raw_counts <- GetAssayData(seurat_cd4, assay = "RNA", layer = "counts")
      raw_counts_valid <- raw_counts[, valid_mask]

      pb_groups <- sort(unique(meta_pb$pb_group))
      cat(sprintf("  Pseudobulk groups (%d):\n    %s\n",
                  length(pb_groups), paste(pb_groups, collapse = "\n    ")))

      pb_mat <- sapply(pb_groups, function(g) {
        idx <- rownames(meta_pb)[meta_pb$pb_group == g]
        if (length(idx) == 1) raw_counts_valid[, idx]
        else rowSums(raw_counts_valid[, idx, drop = FALSE])
      })
      colnames(pb_mat) <- pb_groups

      # --- 3. Build colData for DESeq2 ---
      # Parse group names back into factor columns
      pb_coldata <- data.frame(
        row.names    = pb_groups,
        patient      = factor(sub("_(Unstim|Stim)_(Expanded|NonExpanded)$", "", pb_groups)),
        stimulation  = factor(sub("^P[12]_", "", sub("_(Expanded|NonExpanded)$", "", pb_groups))),
        expansion    = factor(sub("^P[12]_(Unstim|Stim)_", "", pb_groups),
                              levels = c("NonExpanded", "Expanded"))
      )

      # Drop groups with fewer than 5 cells (unreliable pseudobulk)
      n_cells_per_group <- sapply(pb_groups, function(g)
        sum(meta_pb$pb_group == g))
      keep_groups <- pb_groups[n_cells_per_group >= 5]
      dropped     <- pb_groups[n_cells_per_group < 5]

      if (length(dropped) > 0) {
        cat(sprintf("  Groups dropped (<5 cells): %s\n", paste(dropped, collapse = ", ")))
      }

      if (length(keep_groups) < 4) {
        cat(sprintf("  Only %d groups with ≥5 cells — insufficient for DESeq2. Skipping.\n\n",
                    length(keep_groups)))

      } else {

        pb_mat_f    <- pb_mat[, keep_groups, drop = FALSE]
        pb_coldata_f <- pb_coldata[keep_groups, , drop = FALSE]

        # Check both levels of expansion are represented
        if (length(unique(pb_coldata_f$expansion)) < 2) {
          cat("  WARNING: Only one expansion level in filtered groups. Cannot run DE. Skipping.\n\n")

        } else {

          # --- 4. Run DESeq2 ---
          dds <- DESeqDataSetFromMatrix(
            countData = pb_mat_f,
            colData   = pb_coldata_f,
            design    = ~stimulation + expansion
          )

          # Gene-level count filter: ≥5 counts in ≥2 samples
          dds <- dds[rowSums(counts(dds) >= 5) >= 2, ]
          cat(sprintf("  Genes retained after count filter: %d\n", nrow(dds)))
          cat(sprintf("  Pseudobulk samples used: %d  (donor df = %d)\n",
                      ncol(dds), ncol(dds) - 3))

          dds <- DESeq(dds, quiet = TRUE)

          res <- results(dds,
                         contrast = c("expansion", "Expanded", "NonExpanded"),
                         alpha    = 0.05)
          res_df       <- as.data.frame(res)
          res_df$gene  <- rownames(res_df)
          res_df       <- res_df[order(res_df$padj, na.last = TRUE), ]

          n_sig_pb  <- sum(res_df$padj < 0.05, na.rm = TRUE)
          n_up_pb   <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
          n_down_pb <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)

          cat(sprintf("  ✓ DESeq2 complete: %d sig genes (padj<0.05): %d up, %d down\n",
                      n_sig_pb, n_up_pb, n_down_pb))
          cat("  CAVEAT: N=2 donors → very low power. Use as directional guide only.\n")

          if (n_sig_pb > 0) {
            cat("\n  Top pseudobulk DEGs:\n")
            top_pb <- head(res_df[!is.na(res_df$padj), c("gene","log2FoldChange","pvalue","padj")], 10)
            print(top_pb)
          }

          write.csv(res_df,
                    file.path(output_dir, "DEG_pseudobulk_expanded_vs_nonexpanded_CD4.csv"),
                    row.names = FALSE)
          cat(sprintf("\n  ✓ Pseudobulk results → DEG_pseudobulk_expanded_vs_nonexpanded_CD4.csv\n\n"))

          pseudobulk_markers <- res_df

          # Export pseudobulk upregulated genes (for comparison with Wilcoxon list)
          if (n_up_pb > 0) {
            pb_up_genes <- res_df$gene[!is.na(res_df$padj) &
                                        res_df$padj < 0.05 &
                                        res_df$log2FoldChange > 0]
            write.table(pb_up_genes,
                        file.path(output_dir, "expanded_CD4_upregulated_genes_pseudobulk.txt"),
                        row.names = FALSE, col.names = FALSE, quote = FALSE)
          }
        }
      }
    }
  }

}, error = function(e) {
  cat(sprintf("  ERROR in pseudobulk analysis: %s\n\n", e$message))
})

cat("\n")

################################################################################
# Step 4: Overlap with LBD risk genes
################################################################################

overlap_genes <- c()

if (!is.null(lbd_genes) && !is.null(deg_markers) && nrow(deg_markers) > 0) {
  cat("Step 4: Checking overlap with LBD risk genes...\n")

  # Get significant DEGs
  sig_degs <- deg_markers$gene[deg_markers$p_val_adj < 0.05]

  # Check overlap with LBD risk genes
  lbd_symbols <- lbd_genes$SYMBOL
  overlap_genes <- intersect(sig_degs, lbd_symbols)

  cat(sprintf("  Significant DEGs: %d\n", length(sig_degs)))
  cat(sprintf("  LBD risk genes: %d\n", nrow(lbd_genes)))
  cat(sprintf("  Overlap: %d genes\n", length(overlap_genes)))

  if (length(overlap_genes) > 0) {
    cat("\n  ⚠ Overlapping genes (LBD risk + expanded CD4+):\n")
    overlap_df <- merge(
      deg_markers[deg_markers$gene %in% overlap_genes, ],
      lbd_genes[, c("SYMBOL", "ZSTAT", "P")],
      by.x = "gene", by.y = "SYMBOL"
    )
    colnames(overlap_df)[colnames(overlap_df) == "P"] <- "LBD_P"
    colnames(overlap_df)[colnames(overlap_df) == "ZSTAT"] <- "LBD_ZSTAT"

    print(overlap_df[, c("gene", "avg_log2FC", "p_val_adj", "LBD_ZSTAT", "LBD_P")])

    write.csv(overlap_df,
              file.path(output_dir, "LBD_risk_genes_in_expanded_CD4.csv"),
              row.names = FALSE)
  }

  cat("\n")
}

################################################################################
# Step 5: Pathway analysis prep
################################################################################

cat("Step 5: Preparing gene lists for pathway analysis...\n")

if (!is.null(deg_markers) && nrow(deg_markers) > 0) {
  # Export gene lists
  up_genes <- deg_markers$gene[deg_markers$avg_log2FC > 0 & deg_markers$p_val_adj < 0.05]
  down_genes <- deg_markers$gene[deg_markers$avg_log2FC < 0 & deg_markers$p_val_adj < 0.05]

  write.table(up_genes,
              file.path(output_dir, "expanded_CD4_upregulated_genes.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  write.table(down_genes,
              file.path(output_dir, "expanded_CD4_downregulated_genes.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  cat(sprintf("  ✓ Upregulated genes: %d\n", length(up_genes)))
  cat(sprintf("  ✓ Downregulated genes: %d\n\n", length(down_genes)))
}

################################################################################
# Step 6: Visualization (FIXED: REMOVE NAs)
################################################################################

cat("Step 6: Creating visualizations...\n")

if (!is.null(deg_markers) && nrow(deg_markers) > 0) {
  # --- 1. Volcano Plot (Sin cambios) ---
  deg_markers$significance <- ifelse(
    deg_markers$p_val_adj < 0.05 & abs(deg_markers$avg_log2FC) > 0.5,
    "Significant", "Not significant"
  )

  p_volcano <- ggplot(deg_markers, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                       color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(title = "Differential Expression: Expanded vs Non-expanded CD4+ T cells",
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme_minimal()

  ggsave(file.path(plots_dir, "volcano_plot_expanded_cd4.pdf"),
         p_volcano, width = 8, height = 6)

  # --- 2. Heatmap (FIX: Limpiar NAs) ---
  
  # Definir genes top
  top_genes <- c(
    head(deg_markers$gene[deg_markers$avg_log2FC > 0], 20),
    head(deg_markers$gene[deg_markers$avg_log2FC < 0], 20)
  )
  top_genes <- unique(top_genes)
  
  # A. Crear un objeto limpio SOLO para el plot (Sin NAs)
  # Identificar células que NO son NA en is_expanded
  valid_cells <- colnames(seurat_cd4)[!is.na(seurat_cd4$is_expanded)]
  seurat_plot <- subset(seurat_cd4, cells = valid_cells)
  
  cat(sprintf("  Filtered %d cells with 'NA' status for heatmap plotting.\n", 
              ncol(seurat_cd4) - ncol(seurat_plot)))

  # B. Escalar datos en el objeto limpio
  # Es necesario re-escalar porque es un subset nuevo
  cat("  Scaling data for heatmap...\n")
  seurat_plot <- ScaleData(seurat_plot, features = top_genes, verbose = FALSE)

  # C. Generar Heatmap
  # group.bar = TRUE dibuja la barra de colores arriba indicando los grupos
  p_heatmap <- DoHeatmap(seurat_plot, 
                         features = top_genes, 
                         group.by = "is_expanded", 
                         size = 4,
                         draw.lines = TRUE) + # Dibuja líneas entre grupos
    scale_fill_gradientn(colors = c("blue", "white", "red")) +
    ggtitle("Top DEGs in Expanded CD4+ T Cells")

  ggsave(file.path(plots_dir, "heatmap_top_degs_cd4.pdf"),
         p_heatmap, width = 10, height = 8)

  cat("  ✓ Plots saved\n\n")
}
################################################################################
# Summary
################################################################################

cat("================================================================================\n")
cat("Gene Signature Extraction Complete!\n")
cat("================================================================================\n\n")

if (!is.null(deg_markers)) {
  cat("Summary:\n")
  cat(sprintf("  CD4+ T cells analyzed: %d\n", ncol(seurat_cd4)))
  cat(sprintf("  Expanded cells: %d\n", n_expanded))
  cat(sprintf("  DEGs identified: %d\n", nrow(deg_markers)))
  cat(sprintf("  Significant DEGs (padj<0.05): %d\n",
              sum(deg_markers$p_val_adj < 0.05, na.rm=TRUE)))

  if (length(overlap_genes) > 0) {
    cat(sprintf("  Overlap with LBD risk genes: %d\n", length(overlap_genes)))
  }
}

cat("\nOutput files:\n")
cat(sprintf("  - %s  [Wilcoxon — cells as units, inflated significance]\n",
            file.path(output_dir, "DEG_expanded_vs_nonexpanded_CD4.csv")))
cat(sprintf("  - %s  [DESeq2 pseudobulk — donor-aware, preferred]\n",
            file.path(output_dir, "DEG_pseudobulk_expanded_vs_nonexpanded_CD4.csv")))
cat(sprintf("  - %s\n", file.path(output_dir, "expanded_CD4_upregulated_genes.txt")))
cat(sprintf("  - %s  [pseudobulk version]\n",
            file.path(output_dir, "expanded_CD4_upregulated_genes_pseudobulk.txt")))
cat(sprintf("  - Plots: %s\n", plots_dir))

cat("\n✓ Phase 2 Complete!\n")
cat("\nNext: Phase 3 - Neuronal vulnerability analysis (GSE178146)\n")
cat("================================================================================\n")

# Emit session info for reproducibility record
cat("\n--- R Session Info ---\n")
print(sessionInfo())
