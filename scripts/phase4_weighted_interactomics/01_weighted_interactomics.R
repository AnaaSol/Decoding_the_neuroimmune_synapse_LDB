#!/usr/bin/env Rscript
################################################################################
# Phase 4: Weighted Interactomics (FIX: Graph Layout)
# Integration of Transcriptomics (scRNA) + Genetics (GWAS)
# Author: Ana
# Date: 2025-11-29
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(scales)
})

cat("================================================================================\n")
cat("Phase 4: Weighted Interactomics - Integrative Prioritization of Neuro-Immune Pathways\n")
cat("================================================================================\n\n")

# Set paths
project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
phase1_dir <- file.path(project_dir, "data") 
phase3_dir <- file.path(project_dir, "data/GSE178146_processed")
output_dir <- file.path(project_dir, "results_final")
plots_dir <- file.path(output_dir, "plots")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Step 1: Load Data Layers
################################################################################

cat("Step 1: Loading Data Layers...\n")

# A. Active Interactions
interactions_file <- file.path(phase3_dir, "active_ligand_receptor_pairs.csv")
if (!file.exists(interactions_file)) stop("Phase 3 results (interactions) not found.")
interactions <- read.csv(interactions_file)
cat(sprintf("  ✓ Loaded %d active ligand-receptor pairs from Phase 3\n", nrow(interactions)))

# B. Receptor Expression Levels
exp_file <- file.path(phase3_dir, "receptor_expression_by_celltype.csv")
if (!file.exists(exp_file)) stop("Phase 3 results (expression) not found.")
receptor_exp <- read.csv(exp_file, row.names = 1)
cat("  ✓ Loaded receptor expression data\n")

# C. Genetic Risk Scores
gwas_file <- file.path(phase1_dir, "LBD_risk_genes_ranked.tsv")
if (!file.exists(gwas_file)) stop("Phase 1 results (GWAS) not found.")
gwas_data <- read.table(gwas_file, header = TRUE, sep = "\t")
cat(sprintf("  ✓ Loaded %d GWAS risk genes from Phase 1\n", nrow(gwas_data)))

################################################################################
# Step 2: The "Weighting" Engine
################################################################################

cat("\nStep 2: Calculating Genetic Vulnerability Scores...\n")

# 1. GWAS lookup
gwas_scores <- setNames(gwas_data$ZSTAT, gwas_data$SYMBOL)
gwas_pvals <- setNames(gwas_data$P, gwas_data$SYMBOL)

# 2. Enrich Interactions Table
weighted_network <- interactions

# Add Genetic Risk Score for the RECEPTOR
weighted_network$Receptor_Risk_Z <- sapply(weighted_network$receptor, function(r) {
  if (r %in% names(gwas_scores)) return(gwas_scores[[r]]) else return(0)
})

# Add Expression Level
weighted_network$Max_Neuronal_Expression <- sapply(weighted_network$receptor, function(r) {
  if (r %in% rownames(receptor_exp)) return(max(receptor_exp[r, ])) else return(0)
})

# 3. Calculate Integrative Prioritization Score
# Formula: neuronal expression × (1 + max(0, GWAS Z)) — a heuristic weighting,
# NOT causal inference. Higher score = stronger expression + genetic support.
weighted_network$Prioritization_Score <- weighted_network$Max_Neuronal_Expression * (1 + pmax(0, weighted_network$Receptor_Risk_Z))

# Sort
weighted_network <- weighted_network %>% arrange(desc(Prioritization_Score))

# Identify Top Hits
top_hits <- weighted_network %>% filter(Receptor_Risk_Z > 2.0)

cat(sprintf("  ✓ Interactions processed: %d\n", nrow(weighted_network)))
cat(sprintf("  ★ PRIORITY HITS (Genetically Supported): %d\n", nrow(top_hits)))

# Save
write.csv(weighted_network, file.path(output_dir, "weighted_neuroimmune_network.csv"), row.names = FALSE)

################################################################################
# Step 3: Visualization - Scatter Plot
################################################################################

cat("\nStep 3: Generating visualizations...\n")

p_volcano_like <- ggplot(weighted_network, aes(x = Max_Neuronal_Expression, y = Receptor_Risk_Z)) +
  geom_point(aes(color = interaction_type, size = Prioritization_Score), alpha = 0.7) +
  geom_text(aes(label = ifelse(Prioritization_Score > quantile(Prioritization_Score, 0.8),
                               paste0(ligand, "->", receptor), "")),
            vjust = -0.5, size = 3, check_overlap = TRUE) +
  geom_hline(yintercept = 1.96, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(
    title = "Prioritization of Neuro-Immune Interactions",
    x = "Neuronal Receptor Expression (scRNA-seq)",
    y = "Receptor Genetic Risk Z-Score (GWAS)",
    size = "Prioritization Score"
  ) +
  theme_minimal()

ggsave(file.path(plots_dir, "prioritization_scatter_plot.pdf"), p_volcano_like, width = 10, height = 8)

################################################################################
# Step 4: Network Visualization (FIXED LAYOUT)
################################################################################

cat("Step 4: Generating Network Graph...\n")

# Filter for plotting (Top 20 interactions to avoid clutter)
plot_data <- weighted_network %>% head(20)

edges <- plot_data[, c("ligand", "receptor")]

# Create Nodes list explicitly distinguishing Source vs Target type for coloring
nodes_ligand <- data.frame(name = unique(edges$ligand), node_class = "T_Cell_Ligand", risk = 0)
nodes_receptor <- data.frame(
  name = unique(edges$receptor), 
  node_class = "Neuron_Receptor",
  risk = sapply(unique(edges$receptor), function(r) {
    val <- weighted_network$Receptor_Risk_Z[weighted_network$receptor == r][1]
    if(is.na(val)) 0 else val
  })
)

# Combine and remove duplicates
# If a node is both (e.g. CD38), we keep it as Receptor to show the risk
nodes <- rbind(nodes_ligand, nodes_receptor)
nodes <- nodes %>%
  arrange(desc(node_class)) %>% # Sort so Receptor comes first (if duplicates)
  distinct(name, .keep_all = TRUE)

# Create Graph
graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

# Plot using 'nicely' layout (handles self-loops like CD38->CD38 perfectly)
p_net <- ggraph(graph, layout = 'nicely') + 
  geom_edge_fan(aes(alpha = ..index..), show.legend = FALSE, 
                 arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) + 
  geom_node_point(aes(color = node_class, size = ifelse(node_class == "Neuron_Receptor", risk + 2, 2))) +
  geom_node_text(aes(label = name), vjust = 1.8, size = 3.5, fontface = "bold") +
  scale_color_manual(values = c("T_Cell_Ligand" = "#3498db", "Neuron_Receptor" = "#e74c3c")) +
  scale_size(range = c(3, 8)) +
  labs(
    title = "Integrative Neuro-Immune Interactome (Prioritization Network)",
    subtitle = "Nodes colored by source. Red node size = GWAS Z-score. Score = Expression × (1 + GWAS Z).",
    color = "Cell Origin",
    size = "GWAS Risk"
  ) +
  theme_void()

ggsave(file.path(plots_dir, "integrative_prioritization_network.pdf"), p_net, width = 12, height = 8)

cat("  ✓ Network plot saved.\n")

################################################################################
# Summary
################################################################################

cat("\n================================================================================\n")
cat("Phase 4 Complete: Analysis Finished.\n")
cat("================================================================================\n\n")

cat("Top Prioritized Interactions (By Expression + Genetics):\n")
print(head(weighted_network[, c("ligand", "receptor", "Prioritization_Score")], 10))

cat(sprintf("\nResults directory: %s\n", output_dir))
cat("================================================================================\n")

# Emit session info for reproducibility record
cat("\n--- R Session Info ---\n")
print(sessionInfo())
