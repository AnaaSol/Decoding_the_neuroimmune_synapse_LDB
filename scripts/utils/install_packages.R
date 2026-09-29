#!/usr/bin/env Rscript
################################################################################
# Install and pin all R packages required by the LBD Neuroimmune Synapse pipeline
# Run once before executing any analysis scripts:
#   Rscript scripts/utils/install_packages.R
#
# For full reproducibility, initialise renv after running this script:
#   renv::init()   # snapshots installed versions to renv.lock
################################################################################

cat("================================================================================\n")
cat("Installing R package dependencies\n")
cat("================================================================================\n\n")

# CRAN packages
cran_packages <- c(
  "Seurat",        # v5 — Hao et al., Nat Biotechnol 2024;42:293-304
  "ggplot2",       # plotting
  "dplyr",         # data manipulation
  "patchwork",     # plot composition
  "tidyr",         # data tidying
  "Matrix",        # sparse matrix support
  "scales",        # ggplot axis scales
  "pheatmap",      # heatmaps
  "igraph",        # network analysis
  "ggraph",        # network visualisation
  "harmony",       # batch correction — Korsunsky et al., Nat Methods 2019
  "ggrepel",       # volcano plot labels (GSE141578 pseudobulk script)
  "yaml"           # read pipeline_config.yaml
)

# Bioconductor packages
bioc_packages <- c(
  "scRepertoire",  # TCR repertoire analysis
  "DESeq2",        # pseudobulk differential expression
  "scDblFinder",   # doublet detection
  "BiocParallel"   # parallelisation
)

# Install CRAN packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  } else {
    cat(sprintf("  ✓ %s already installed (v%s)\n",
                pkg, as.character(packageVersion(pkg))))
  }
}

# Install Bioconductor packages
cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
  } else {
    cat(sprintf("  ✓ %s already installed (v%s)\n",
                pkg, as.character(packageVersion(pkg))))
  }
}

# Print session info for reproducibility record
cat("\n================================================================================\n")
cat("Installed package versions:\n")
cat("================================================================================\n")
all_pkgs <- c(cran_packages, bioc_packages)
for (pkg in all_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %-20s %s\n", pkg, as.character(packageVersion(pkg))))
  } else {
    cat(sprintf("  %-20s FAILED TO INSTALL\n", pkg))
  }
}

cat("\nR version:\n")
cat(sprintf("  %s\n", R.version.string))

cat("\n✓ Package installation complete.\n")
cat("Run renv::init() to snapshot versions to renv.lock for full reproducibility.\n")
cat("================================================================================\n")
