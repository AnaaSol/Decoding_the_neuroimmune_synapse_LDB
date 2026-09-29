#!/usr/bin/env Rscript
################################################################################
# eQTL colocalization (coloc::coloc.abf): LBD GWAS vs an eQTL dataset.
#
# Usage: Rscript 09_run_coloc.R [SOURCE]
#   SOURCE = "brain_acc" (default) -- GTEx anterior cingulate cortex, tissue-
#            matched to the GSE178146 snRNA-seq dataset used in Phase 3.
#   SOURCE = "blueprint_tcell" -- BLUEPRINT CD4+ naive T cell, matched to the
#            study's T-cell-mediated mechanistic hypothesis (a more direct
#            test of "does LBD genetic risk act through T-cell expression"
#            than the brain-tissue coloc alone).
#
# Addresses REPORT.md Sec 5.8/12.2 gap #2. Standard single-tissue coloc.abf
# (Wakefield ABF, default priors p1=p2=1e-4, p12=1e-5) -- NOT SuSiE fine-
# mapping, NOT Mendelian randomization. High PP4 = consistent with a shared
# causal variant; it is not proof of causality. Low/ambiguous PP4 (e.g. most
# mass on PP3, distinct causal variants) is a real, reportable result.
################################################################################

suppressPackageStartupMessages({
  library(coloc)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
source_name <- if (length(args) >= 1) args[1] else "brain_acc"

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
cfg <- yaml::read_yaml(file.path(project_dir, "config/pipeline_config.yaml"))
coloc_cfg <- cfg$phase1c_coloc

if (source_name == "brain_acc") {
  src_cfg <- coloc_cfg$eqtl_source
} else {
  matches <- Filter(function(s) s$name == source_name, coloc_cfg$eqtl_sources_extra)
  if (length(matches) == 0) stop(sprintf("Unknown eQTL source: %s", source_name))
  src_cfg <- matches[[1]]
}

cat("================================================================================\n")
cat(sprintf("Phase 1c: eQTL Colocalization (coloc.abf) vs %s (%s)\n", source_name, src_cfg$tissue))
cat("================================================================================\n\n")

gwas_dir <- file.path(project_dir, "data/coloc/gwas")
eqtl_dir <- file.path(project_dir, "data/coloc/eqtl", source_name)
results_dir <- file.path(project_dir, "data/coloc/results", source_name)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

n_total <- coloc_cfg$gwas$n_total
s_cases <- coloc_cfg$gwas$n_cases / n_total
n_eqtl <- src_cfg$sample_size

run_coloc_one_gene <- function(gene) {
  symbol <- gene$symbol
  ensembl <- gene$ensembl_id
  cat(sprintf("-- %s (%s) --\n", symbol, ensembl))

  gwas_file <- file.path(gwas_dir, paste0(symbol, "_gwas.tsv"))
  eqtl_file <- file.path(eqtl_dir, paste0(symbol, "_eqtl_raw.tsv"))

  if (!file.exists(gwas_file) || !file.exists(eqtl_file)) {
    cat("   SKIPPED: missing input file(s)\n\n")
    return(NULL)
  }

  gwas <- read.delim(gwas_file, stringsAsFactors = FALSE)
  eqtl <- read.delim(eqtl_file, stringsAsFactors = FALSE)

  if (nrow(gwas) == 0 || nrow(eqtl) == 0) {
    cat(sprintf("   SKIPPED: empty input (GWAS rows=%d, eQTL rows=%d)\n\n", nrow(gwas), nrow(eqtl)))
    return(data.frame(symbol = symbol, ensembl_id = ensembl, n_snps = 0,
                       PP0 = NA, PP1 = NA, PP2 = NA, PP3 = NA, PP4 = NA,
                       note = "empty eQTL input for this gene in this tissue"))
  }

  # eQTL Catalogue stores one row per rsid per variant; a single physical
  # variant can map to >1 rsid in dbSNP. Per their own FAQ/tutorial guidance,
  # drop duplicate rsids before colocalization.
  eqtl <- eqtl[!duplicated(eqtl$rsid), ]
  gwas <- gwas[!duplicated(gwas$rsid), ]

  merged <- merge(gwas, eqtl, by = "rsid", suffixes = c(".gwas", ".eqtl"))

  # Sanity QC: keep only SNPs where the GWAS and eQTL allele pairs match
  # (as an unordered set) -- excludes rare cross-file annotation errors.
  allele_pair_matches <- mapply(function(ea, oa, ref, alt) {
    setequal(c(toupper(ea), toupper(oa)), c(toupper(ref), toupper(alt)))
  }, merged$effect_allele, merged$other_allele, merged$ref, merged$alt)
  merged <- merged[allele_pair_matches, ]

  n_snps <- nrow(merged)
  cat(sprintf("   %d GWAS SNPs, %d eQTL SNPs, %d overlapping after rsid+allele match\n",
              nrow(gwas), nrow(eqtl), n_snps))

  if (n_snps < 20) {
    cat("   SKIPPED: fewer than 20 overlapping SNPs -- coloc result would be unreliable\n\n")
    return(data.frame(symbol = symbol, ensembl_id = ensembl, n_snps = n_snps,
                       PP0 = NA, PP1 = NA, PP2 = NA, PP3 = NA, PP4 = NA,
                       note = "insufficient overlapping SNPs"))
  }

  gwas_maf <- pmin(merged$eaf, 1 - merged$eaf)

  dataset_gwas <- list(
    snp = merged$rsid,
    beta = merged$beta.gwas,
    varbeta = merged$se.gwas^2,
    type = "cc",
    s = s_cases,
    N = n_total,
    MAF = gwas_maf
  )

  dataset_eqtl <- list(
    snp = merged$rsid,
    beta = merged$beta.eqtl,
    varbeta = merged$se.eqtl^2,
    type = "quant",
    N = n_eqtl,
    MAF = merged$maf
  )

  res <- tryCatch(
    coloc.abf(dataset1 = dataset_gwas, dataset2 = dataset_eqtl),
    error = function(e) {
      cat(sprintf("   ERROR running coloc.abf: %s\n\n", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(res)) {
    return(data.frame(symbol = symbol, ensembl_id = ensembl, n_snps = n_snps,
                       PP0 = NA, PP1 = NA, PP2 = NA, PP3 = NA, PP4 = NA,
                       note = "coloc.abf error"))
  }

  pp <- res$summary
  cat(sprintf("   PP4 (shared causal variant) = %.4f | PP3 (distinct variants) = %.4f\n\n",
              pp["PP.H4.abf"], pp["PP.H3.abf"]))

  write.csv(res$results, file.path(results_dir, paste0(symbol, "_coloc_per_snp.csv")), row.names = FALSE)

  data.frame(
    symbol = symbol, ensembl_id = ensembl, n_snps = n_snps,
    PP0 = pp["PP.H0.abf"], PP1 = pp["PP.H1.abf"], PP2 = pp["PP.H2.abf"],
    PP3 = pp["PP.H3.abf"], PP4 = pp["PP.H4.abf"], note = NA
  )
}

results <- do.call(rbind, lapply(coloc_cfg$target_genes, run_coloc_one_gene))
rownames(results) <- NULL
results$eqtl_source <- source_name

summary_file <- file.path(results_dir, "coloc_summary.csv")
write.csv(results, summary_file, row.names = FALSE)

cat("================================================================================\n")
cat(sprintf("Colocalization summary (GWAS vs %s eQTL, coloc.abf default priors)\n", source_name))
cat("================================================================================\n")
print(results[, c("symbol", "n_snps", "PP3", "PP4")])
cat(sprintf("\nFull results: %s\n", summary_file))
cat(sprintf("Per-SNP posterior tables: %s/{symbol}_coloc_per_snp.csv\n\n", results_dir))

cat("Interpretation guide (standard coloc.abf convention):\n")
cat("  PP4 > 0.8: strong evidence for a single shared causal variant\n")
cat("  PP3 > 0.8: strong evidence for two DISTINCT causal variants (genetic\n")
cat("             signal and eQTL signal are both real, but not the same variant)\n")
cat("  PP4 and PP3 both low/moderate: usually underpowered (too few SNPs, or\n")
cat("             weak eQTL/GWAS signal in this window) -- not evidence against\n")
cat("             colocalization, just an inconclusive test\n\n")

cat("--- R Session Info ---\n")
print(sessionInfo())
