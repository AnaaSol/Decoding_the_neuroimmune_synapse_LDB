#!/usr/bin/env Rscript
################################################################################
# TCR Repertoire Analysis with scRepertoire
# Author: Ana
# Date: 2025-11-25
# Phase: 2 - Immune Cell Profiling
#
# Purpose: Integrate TCR sequencing data with scRNA-seq
#   - Load TCR contig annotations from 10X VDJ
#   - Calculate clonality metrics
#   - Identify expanded clones
#   - Integrate with Seurat object
#   - Focus on CD4+ T cells with clonal expansion
#
# Input:
#   - TCR contig annotations (filtered_contig_annotations.csv)
#   - Processed Seurat object from step 01
# Output:
#   - Seurat object with TCR information
#   - Clonality analysis results
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(scRepertoire)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)  # Reproducibility

cat("================================================================================\n")
cat("TCR Repertoire Analysis - GSE161192\n")
cat("================================================================================\n\n")

# Set paths
project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
data_dir <- file.path(project_dir, "data/GSE161192_raw")
processed_dir <- file.path(project_dir, "data/GSE161192_processed")
output_dir <- processed_dir
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Step 1: Load processed Seurat object
################################################################################

cat("Step 1: Loading processed Seurat object...\n")

seurat_file <- file.path(processed_dir, "CSF_Tcells_seurat_processed.rds")
if (!file.exists(seurat_file)) {
  stop("Error: Processed Seurat object not found. Run 01_process_scrna_seurat.R first.")
}

seurat_obj <- readRDS(seurat_file)
cat(sprintf("  ✓ Loaded Seurat object with %d cells\n\n", ncol(seurat_obj)))

################################################################################
# Step 2: Load TCR data
################################################################################

cat("Step 2: Loading TCR contig annotations...\n")

# Sample metadata
samples <- data.frame(
  sample_id = c("P1U", "P1S", "P2U", "P2S"),
  gsm = c("GSM4891327", "GSM4891328", "GSM4891329", "GSM4891330")
)

tcr_list <- list()

for (i in 1:nrow(samples)) {
  sample_name <- samples$sample_id[i]
  gsm_id <- samples$gsm[i]

  cat(sprintf("  Loading %s...\n", sample_name))

  # Path to TCR contig file
  tcr_file <- file.path(data_dir,
                        sprintf("%s_%s_filtered_contig_annotations.csv.gz",
                               gsm_id, sample_name))

  if (file.exists(tcr_file)) {
    # Read TCR data
    contig_data <- read.csv(gzfile(tcr_file), stringsAsFactors = FALSE)

    tcr_list[[sample_name]] <- contig_data

    cat(sprintf("    %d TCR contigs loaded\n", nrow(contig_data)))
  } else {
    cat(sprintf("    Warning: TCR file not found: %s\n", tcr_file))
  }
}

cat("\n")

################################################################################
# Step 3: Process TCR data with scRepertoire (FIXED VERSION)
################################################################################

cat("Step 3: Processing TCR data with scRepertoire...\n")

# 1. Combinar (Esto generará prefijos dobles como P1U_P1U_..., pero es predecible)
combined_tcr <- combineTCR(
  tcr_list,
  samples = names(tcr_list),
  ID = names(tcr_list), 
  cells = "T-AB"
)

# 2. CIRUGÍA: Corregir los Barcodes manualmente para que coincidan con Seurat
# Queremos cambiar "P1U_P1U_Barcode-1" a "P1U_Barcode-1"
cat("  ✓ Fixing barcodes to match Seurat format exactly...\n")

for (i in seq_along(combined_tcr)) {
    # Obtener el prefijo actual (ej. P1U)
    sample_prefix <- names(tcr_list)[i]
    
    # Definir el patrón erróneo (P1U_P1U_)
    wrong_pattern <- paste0(sample_prefix, "_", sample_prefix, "_")
    
    # Definir el patrón correcto (P1U_)
    correct_pattern <- paste0(sample_prefix, "_")
    
    # Reemplazar en la columna barcode
    combined_tcr[[i]]$barcode <- gsub(wrong_pattern, correct_pattern, combined_tcr[[i]]$barcode)
}

# 3. VERIFICACIÓN FINAL (Ojos que ven, script que no falla)
cat("\n  --- QUALITY CHECK ---\n")
tcr_example <- combined_tcr[[1]]$barcode[1]
seurat_example <- colnames(seurat_obj)[1]

cat(sprintf("  TCR Barcode format:    %s\n", tcr_example))
cat(sprintf("  Seurat Barcode format: %s\n", seurat_example))

# Verificar si el formato del TCR está contenido en Seurat (ignorando orden)
if (tcr_example %in% colnames(seurat_obj)) {
    cat("  ✓ MATCH CONFIRMED! Los barcodes coinciden.\n")
} else {
    # Si sigue fallando, probamos si es culpa del sufijo -1
    if (!grepl("-1", tcr_example) && grepl("-1", seurat_example)) {
        cat("  ! ALERTA: Falta el sufijo '-1' en el TCR. Agregándolo...\n")
        for (i in seq_along(combined_tcr)) {
            combined_tcr[[i]]$barcode <- paste0(combined_tcr[[i]]$barcode, "-1")
        }
        cat("  ✓ Sufijo -1 agregado.\n")
    } else {
        cat("  ⚠ ADVERTENCIA: Aún parecen diferentes. Revisa los prints de arriba.\n")
    }
}
cat("\n")

cat("  ✓ TCR data combined across samples\n")

# Calculate clonotype diversity (Código original)
cat("\n  Clonotype diversity per sample:\n")
for (sample_name in names(combined_tcr)) {
  tcr_data <- combined_tcr[[sample_name]]
  n_clonotypes <- length(unique(tcr_data$CTaa))
  n_cells <- nrow(tcr_data)
  cat(sprintf("    %s: %d cells, %d unique clonotypes\n",
              sample_name, n_cells, n_clonotypes))
}

cat("\n")

################################################################################
# Step 4: Integrate TCR with Seurat
################################################################################

cat("Step 4: Integrating TCR data with Seurat object...\n")

# Combine into single object for Seurat integration
seurat_obj <- combineExpression(
  combined_tcr,
  seurat_obj,
  cloneCall = "aa",  # Use amino acid sequence for clonotype definition
  proportion = FALSE
)

cat("  ✓ TCR data integrated into Seurat object\n")

# Count cells with TCR information
cells_with_tcr <- sum(!is.na(seurat_obj$CTaa))
cat(sprintf("  Cells with TCR: %d / %d (%.1f%%)\n\n",
            cells_with_tcr, ncol(seurat_obj),
            (cells_with_tcr / ncol(seurat_obj)) * 100))

################################################################################
# Step 5: Identify clonal expansion
################################################################################

cat("Step 5: Analyzing clonal expansion...\n")

# Calculate clonal frequency
clonotype_freq <- table(seurat_obj$CTaa[!is.na(seurat_obj$CTaa)])
clonotype_freq <- sort(clonotype_freq, decreasing = TRUE)

# Expansion threshold: ≥2 cells with identical CDR3 amino acid sequence (cloneCall="aa")
# Note: amino-acid identity is less stringent than nucleotide identity; different
# DNA rearrangements can encode the same CDR3 aa. Report both aa and nt thresholds.
EXPANSION_THRESHOLD_AA <- 2  # primary (aa-level)
EXPANSION_THRESHOLD_STRICT <- 3  # secondary / sensitivity check

expanded_clones <- names(clonotype_freq[clonotype_freq >= EXPANSION_THRESHOLD_AA])
n_expanded <- length(expanded_clones)

# Stricter threshold
expanded_clones_strict <- names(clonotype_freq[clonotype_freq >= EXPANSION_THRESHOLD_STRICT])
n_expanded_strict <- length(expanded_clones_strict)
cells_expanded_strict <- sum(clonotype_freq[expanded_clones_strict])

cat(sprintf("  Total unique clonotypes: %d\n", length(clonotype_freq)))
cat(sprintf("  Expanded clones (≥%d cells, aa-level): %d (%d cells)\n",
            EXPANSION_THRESHOLD_AA, n_expanded, sum(clonotype_freq[expanded_clones])))
cat(sprintf("  Expanded clones (≥%d cells, stricter): %d (%d cells)\n",
            EXPANSION_THRESHOLD_STRICT, n_expanded_strict, cells_expanded_strict))
cat(sprintf("  Largest clone size: %d cells\n", max(clonotype_freq)))

# Clone-size distribution statistics
cat("\n  Clone size distribution:\n")
size_breaks <- c(1, 2, 3, 5, 10, 50, 100, Inf)
size_labels <- c("singleton (1)", "2 cells", "3-4 cells", "5-9 cells",
                 "10-49 cells", "50-99 cells", "100+ cells")
size_counts <- cut(as.numeric(clonotype_freq), breaks = size_breaks,
                   labels = size_labels, right = FALSE)
size_table <- table(size_counts)
for (i in seq_along(size_table)) {
  cat(sprintf("    %s: %d clonotypes\n", names(size_table)[i], size_table[i]))
}
cat(sprintf("  Gini index (clonal skewness): %.3f\n", {
  f <- sort(as.numeric(clonotype_freq))
  n <- length(f)
  2 * sum(seq_len(n) * f) / (n * sum(f)) - (n + 1) / n
}))

cat("\n  Top 10 most expanded clones:\n")
top_clones <- head(clonotype_freq, 10)
for (i in 1:length(top_clones)) {
  cat(sprintf("    %d. %s: %d cells\n", i, names(top_clones)[i], top_clones[i]))
}
cat("\n")

# Save clone-size summary
clone_size_df <- data.frame(
  clonotype = names(clonotype_freq),
  n_cells = as.numeric(clonotype_freq),
  is_expanded_aa2 = as.numeric(clonotype_freq) >= EXPANSION_THRESHOLD_AA,
  is_expanded_strict = as.numeric(clonotype_freq) >= EXPANSION_THRESHOLD_STRICT
)
write.csv(clone_size_df,
          file.path(output_dir, "clonotype_size_table.csv"),
          row.names = FALSE)

# Add expansion status to metadata (primary: ≥2 cells, aa-level)
seurat_obj$is_expanded <- seurat_obj$CTaa %in% expanded_clones
seurat_obj$is_expanded[is.na(seurat_obj$CTaa)] <- NA

# Also store strict expansion status for sensitivity analysis
seurat_obj$is_expanded_strict <- seurat_obj$CTaa %in% expanded_clones_strict
seurat_obj$is_expanded_strict[is.na(seurat_obj$CTaa)] <- NA

################################################################################
# Step 6: Focus on CD4+ T cells (UPDATED FOR SEURAT V5)
################################################################################

cat("Step 6: Identifying CD4+ T cells...\n")

# 1. FIX PARA SEURAT V5: Unir capas antes de acceder a los datos
# Esto soluciona el error "GetAssayData doesn't work for multiple layers"
if (packageVersion("Seurat") >= "5.0.0") {
  cat("  ✓ Detected Seurat v5: Joining layers for data access...\n")
  seurat_obj <- JoinLayers(seurat_obj)
}

# 2. Usar FetchData en lugar de GetAssayData (Es más robusto)
# Extraemos CD4 y CD8A directamente a un data.frame
gene_expression <- FetchData(seurat_obj, vars = c("CD4", "CD8A"), layer = "data")

# 3. Definir células CD4+ (CD4 > 0.5 y CD8A < 0.5)
# Usamos las columnas del data.frame que acabamos de crear
is_cd4_pos <- gene_expression$CD4 > 0.5 & gene_expression$CD8A < 0.5

n_cd4 <- sum(is_cd4_pos)
cat(sprintf("  CD4+ T cells identified: %d (%.1f%%)\n",
            n_cd4, (n_cd4 / ncol(seurat_obj)) * 100))

# Add to metadata
seurat_obj$is_CD4 <- is_cd4_pos

# CD4+ cells with expansion
cd4_expanded <- is_cd4_pos & !is.na(seurat_obj$is_expanded) & seurat_obj$is_expanded
n_cd4_expanded <- sum(cd4_expanded)

cat(sprintf("  CD4+ cells with clonal expansion: %d\n", n_cd4_expanded))
cat(sprintf("  These cells are likely α-synuclein reactive!\n\n"))

# Save list of expanded CD4+ clonotypes
if (n_cd4_expanded > 0) {
    expanded_cd4_clones <- unique(seurat_obj$CTaa[cd4_expanded])
    write.csv(
      data.frame(clonotype = expanded_cd4_clones),
      file.path(output_dir, "expanded_CD4_clonotypes.csv"),
      row.names = FALSE
    )
} else {
    cat("  ⚠ No expanded CD4+ clones found.\n")
}

################################################################################
# Step 7: Save integrated object
################################################################################

cat("Step 7: Saving integrated Seurat object...\n")

output_file <- file.path(output_dir, "CSF_Tcells_with_TCR.rds")
saveRDS(seurat_obj, file = output_file)

cat(sprintf("  ✓ Saved to: %s\n\n", output_file))

################################################################################
# Step 8: Generate clonality plots
################################################################################

cat("Step 8: Generating clonality visualizations...\n")

# UMAP with clonal expansion
p1 <- DimPlot(seurat_obj,
              group.by = "is_expanded",
              cols = c("TRUE" = "red", "FALSE" = "gray80"),
              order = TRUE) +
  ggtitle("Clonally Expanded T Cells") +
  theme_minimal()

# UMAP highlighting CD4+ cells
p2 <- DimPlot(seurat_obj,
              group.by = "is_CD4",
              cols = c("TRUE" = "blue", "FALSE" = "gray80"),
              order = TRUE) +
  ggtitle("CD4+ T Cells") +
  theme_minimal()

# CD4 expression
p3 <- FeaturePlot(seurat_obj, features = "CD4") +
  ggtitle("CD4 Expression") +
  theme_minimal()

# Clonal frequency histogram
clone_sizes <- as.data.frame(table(clonotype_freq))
colnames(clone_sizes) <- c("Clone_Size", "Frequency")
clone_sizes$Clone_Size <- as.numeric(as.character(clone_sizes$Clone_Size))

p4 <- ggplot(clone_sizes, aes(x = Clone_Size, y = Frequency)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  scale_x_continuous(breaks = seq(1, max(clone_sizes$Clone_Size), by = 1)) +
  labs(title = "Clonotype Size Distribution",
       x = "Number of cells per clonotype",
       y = "Number of clonotypes") +
  theme_minimal()

# Save plots
ggsave(file.path(plots_dir, "umap_clonal_expansion.pdf"), p1, width = 8, height = 6)
ggsave(file.path(plots_dir, "umap_cd4_cells.pdf"), p2, width = 8, height = 6)
ggsave(file.path(plots_dir, "cd4_expression.pdf"), p3, width = 8, height = 6)
ggsave(file.path(plots_dir, "clonotype_size_distribution.pdf"), p4, width = 8, height = 6)

cat("  ✓ Plots saved\n\n")

################################################################################
# Step 9: Export results table (FIXED)
################################################################################

cat("Step 9: Exporting results table...\n")

# Asegurarnos de tener los datos de expresión disponibles
# (Recuperamos los datos usando FetchData, igual que en el paso 6)
data_for_export <- FetchData(seurat_obj, vars = c("CD4", "CD8A"), layer = "data")

# Create results table
results_df <- data.frame(
  cell_barcode = colnames(seurat_obj),
  cluster = Idents(seurat_obj),
  # Nota: Si 'sample' o 'condition' no existen en tu metadata, usa 'orig.ident'
  sample = if("sample" %in% colnames(seurat_obj@meta.data)) seurat_obj$sample else seurat_obj$orig.ident,
  condition = if("condition" %in% colnames(seurat_obj@meta.data)) seurat_obj$condition else NA,
  clonotype = seurat_obj$CTaa,
  is_CD4 = seurat_obj$is_CD4,
  is_expanded = seurat_obj$is_expanded,
  
  # --- FIX: Usar la nueva tabla de datos ---
  CD4_expression = data_for_export$CD4,
  CD8A_expression = data_for_export$CD8A,
  # -----------------------------------------
  
  UMAP_1 = Embeddings(seurat_obj, "umap")[, 1],
  UMAP_2 = Embeddings(seurat_obj, "umap")[, 2]
)

# Filter for CD4+ expanded cells
# (Aseguramos que no haya NAs en is_expanded para el filtro)
keep_expanded <- !is.na(results_df$is_CD4) & results_df$is_CD4 & 
                 !is.na(results_df$is_expanded) & results_df$is_expanded

cd4_expanded_table <- results_df[keep_expanded, ]

write.csv(results_df,
          file.path(output_dir, "all_cells_with_tcr_annotations.csv"),
          row.names = FALSE)

write.csv(cd4_expanded_table,
          file.path(output_dir, "CD4_expanded_clones.csv"),
          row.names = FALSE)

cat(sprintf("  ✓ Results tables saved\n\n"))

################################################################################
# Summary
################################################################################

cat("================================================================================\n")
cat("TCR Integration Complete!\n")
cat("================================================================================\n\n")

cat("Summary:\n")
cat(sprintf("  Total cells: %d\n", ncol(seurat_obj)))
cat(sprintf("  Cells with TCR: %d (%.1f%%)\n",
            cells_with_tcr, (cells_with_tcr/ncol(seurat_obj))*100))
cat(sprintf("  Unique clonotypes: %d\n", length(clonotype_freq)))
cat(sprintf("  Expanded clones (≥2 cells): %d\n", n_expanded))
cat(sprintf("  CD4+ T cells: %d (%.1f%%)\n",
            n_cd4, (n_cd4/ncol(seurat_obj))*100))
cat(sprintf("  CD4+ expanded clones: %d\n", n_cd4_expanded))
cat("\nKey findings:\n")
cat(sprintf("  ⚠ %d CD4+ T cells show clonal expansion\n", n_cd4_expanded))
cat("  → These are likely α-synuclein reactive T cells!\n")
cat("\nOutput files:\n")
cat(sprintf("  - Integrated object: %s\n", output_file))
cat(sprintf("  - CD4+ expanded clones: %s\n",
            file.path(output_dir, "CD4_expanded_clones.csv")))
cat(sprintf("  - Plots: %s\n", plots_dir))
cat("\nNext steps:\n")
cat("  → Extract gene signatures from CD4+ expanded clones\n")
cat("  → Compare with LBD risk genes from Phase 1\n")
cat("  → Proceed to Phase 3 (neuronal vulnerability)\n")
cat("================================================================================\n")
