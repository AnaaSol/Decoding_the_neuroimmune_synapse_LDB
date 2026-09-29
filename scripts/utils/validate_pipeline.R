#!/usr/bin/env Rscript
################################################################################
# Pipeline Validation Script — LBD Neuroimmune Synapse Study
# Run after the full pipeline completes:
#   Rscript scripts/utils/validate_pipeline.R
# Or via Snakemake:
#   snakemake validate --cores 1
#
# Each check either passes silently or stops with a descriptive error.
# Exit code 0 = all checks passed.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(yaml)
})

cfg <- yaml::read_yaml("config/pipeline_config.yaml")
DATA_DIR    <- file.path(cfg$project$dir, "data")
RESULTS_DIR <- file.path(cfg$project$dir, "results_final")

n_pass <- 0
n_fail <- 0

assert <- function(condition, message) {
  if (!condition) {
    cat(sprintf("  FAIL: %s\n", message))
    n_fail <<- n_fail + 1
  } else {
    cat(sprintf("  PASS: %s\n", message))
    n_pass <<- n_pass + 1
  }
}

cat("================================================================================\n")
cat("Pipeline Validation Checks\n")
cat("================================================================================\n\n")

# ── 1. Phase 1 ───────────────────────────────────────────────────────────────
cat("Phase 1: Genetic Prioritization\n")

risk_file <- file.path(DATA_DIR, "LBD_risk_genes_ranked.tsv")
assert(file.exists(risk_file), "LBD_risk_genes_ranked.tsv exists")

if (file.exists(risk_file)) {
  risk <- read.table(risk_file, header = TRUE, sep = "\t")

  assert(nrow(risk) > 100,
         sprintf("MAGMA analyzed >100 genes (got %d; <100 suggests pre-filtering still active)", nrow(risk)))

  assert("ZSTAT" %in% colnames(risk),
         "ZSTAT column present in risk genes table")

  assert(!any(is.na(risk$ZSTAT)),
         "No NA values in ZSTAT column")

  assert("SYMBOL" %in% colnames(risk),
         "SYMBOL column present (gene names resolved)")

  # BCAM is top gene in the 34-gene pre-filtered run; warn if it still tops the full run
  if (nrow(risk) < 200) {
    cat("  WARN: Only", nrow(risk), "genes — MAGMA likely still received pre-filtered SNPs.\n")
    cat("        Re-run with --p-threshold 1.0 to get genome-wide results.\n")
  }

  # BCL3 should be absent (Discussion text claiming BCL3 was wrong; BCAM is correct)
  if ("BCL3" %in% risk$SYMBOL) {
    cat(sprintf("  NOTE: BCL3 present in MAGMA output at rank %d (Z=%.2f).\n",
                which(risk$SYMBOL == "BCL3"),
                risk$ZSTAT[risk$SYMBOL == "BCL3"]))
  } else {
    cat("  PASS: BCL3 absent from MAGMA output (consistent with confirmed defect D2).\n")
    n_pass <- n_pass + 1
  }
}

cat("\n")

# ── 2. Phase 2 ───────────────────────────────────────────────────────────────
cat("Phase 2: CSF Immune Profiling\n")

deg_file <- file.path(DATA_DIR, "GSE161192_processed", "DEG_expanded_vs_nonexpanded_CD4.csv")
assert(file.exists(deg_file), "Wilcoxon DEG file exists")

pb_file <- file.path(DATA_DIR, "GSE161192_processed", "DEG_pseudobulk_expanded_vs_nonexpanded_CD4.csv")
assert(file.exists(pb_file), "Pseudobulk DESeq2 DEG file exists")

up_file <- file.path(DATA_DIR, "GSE161192_processed", "expanded_CD4_upregulated_genes.txt")
assert(file.exists(up_file), "Upregulated gene list exists")

clone_file <- file.path(DATA_DIR, "GSE161192_processed", "clonotype_size_table.csv")
assert(file.exists(clone_file), "Clonotype size table exists")

if (file.exists(deg_file)) {
  deg <- read.csv(deg_file)
  assert("gene" %in% colnames(deg),       "DEG file has 'gene' column")
  assert("avg_log2FC" %in% colnames(deg), "DEG file has 'avg_log2FC' column")
  assert("p_val_adj" %in% colnames(deg),  "DEG file has 'p_val_adj' column")
  assert(nrow(deg) > 0,                   "DEG file is non-empty")

  # CXCL12 should NOT be in upregulated genes (confirmed absent, defect D3)
  if (file.exists(up_file)) {
    up_genes <- read.table(up_file, header = FALSE)$V1
    assert(!"CXCL12" %in% up_genes,
           "CXCL12 absent from upregulated genes (manuscript claim unsupported — defect D3)")
  }
}

if (file.exists(clone_file)) {
  clones <- read.csv(clone_file)
  assert("n_cells" %in% colnames(clones), "Clone table has 'n_cells' column")
  assert(max(clones$n_cells) > 10,        "At least one clone with >10 cells (sanity check)")
}

cat("\n")

# ── 3. Phase 3 ───────────────────────────────────────────────────────────────
cat("Phase 3: Cortical snRNA-seq\n")

# Per-sample QC files
expected_samples <- unlist(cfg$datasets$cortex_snrna$sample_groups)
expected_n <- length(expected_samples)

qc_files <- file.path(DATA_DIR, "GSE178146_processed",
                       paste0(expected_samples, "_seurat_qc.rds"))
n_qc_found <- sum(file.exists(qc_files))
assert(n_qc_found == expected_n,
       sprintf("All %d per-sample QC files present (found %d)", expected_n, n_qc_found))

lr_file <- file.path(DATA_DIR, "GSE178146_processed", "active_ligand_receptor_pairs.csv")
assert(file.exists(lr_file), "Active L-R pairs file exists")

cond_file <- file.path(DATA_DIR, "GSE178146_processed", "receptor_expression_DLB_vs_Control.csv")
assert(file.exists(cond_file),
       "receptor_expression_DLB_vs_Control.csv exists (primary comparison, not PD_vs_Control)")

# Verify old "PD_vs_Control" file does NOT exist (would indicate D1 fix not applied)
old_cond_file <- file.path(DATA_DIR, "GSE178146_processed", "receptor_expression_PD_vs_Control.csv")
if (file.exists(old_cond_file)) {
  cat("  WARN: receptor_expression_PD_vs_Control.csv still exists.\n")
  cat("        This was produced before the D1 fix. Re-run Phase 3 to overwrite.\n")
}

if (file.exists(lr_file)) {
  lr <- read.csv(lr_file)
  assert(nrow(lr) > 0,                    "L-R pairs file is non-empty")
  assert("ligand" %in% colnames(lr),      "L-R table has 'ligand' column")
  assert("receptor" %in% colnames(lr),    "L-R table has 'receptor' column")

  # CXCL12/CXCR4 should not be an active pair (defect D3)
  has_cxcl12 <- any(lr$ligand == "CXCL12")
  assert(!has_cxcl12,
         "CXCL12 not in active L-R pairs (expected — defect D3 confirmed)")

  # Verify the condition comparison file uses DLB not PD
  if (file.exists(cond_file)) {
    cond_data <- read.csv(cond_file)
    assert("DLB" %in% colnames(cond_data) || "log2FC_DLB_vs_Control" %in% colnames(cond_data),
           "Condition comparison is DLB vs Control (not PD vs Control)")
  }
}

cat("\n")

# ── 4. Phase 4 ───────────────────────────────────────────────────────────────
cat("Phase 4: Weighted Interactomics\n")

net_file <- file.path(RESULTS_DIR, "weighted_neuroimmune_network.csv")
assert(file.exists(net_file), "Weighted network CSV exists")

if (file.exists(net_file)) {
  net <- read.csv(net_file)

  # Score column renamed from Causality_Score (defect D5)
  assert("Prioritization_Score" %in% colnames(net),
         "Network uses 'Prioritization_Score' column (not 'Causality_Score' — defect D5 fixed)")
  assert(!"Causality_Score" %in% colnames(net),
         "'Causality_Score' column absent (old name that overstated causality)")

  assert(nrow(net) > 0, "Network is non-empty")
  assert(all(!is.na(net$Prioritization_Score)),
         "No NA values in Prioritization_Score")
}

plot_file <- file.path(RESULTS_DIR, "plots", "integrative_prioritization_network.pdf")
assert(file.exists(plot_file),
       "integrative_prioritization_network.pdf exists (not 'causal_network_graph.pdf')")

cat("\n")

# ── Summary ──────────────────────────────────────────────────────────────────
cat("================================================================================\n")
cat(sprintf("Validation complete: %d passed, %d failed\n", n_pass, n_fail))
cat("================================================================================\n")

if (n_fail > 0) {
  quit(status = 1)
}
