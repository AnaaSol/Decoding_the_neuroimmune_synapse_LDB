#!/usr/bin/env Rscript
################################################################################
# Two-sample Mendelian randomization: genetically-predicted gene expression
# (eQTL, brain ACC or CD4+ T cell) -> LBD risk (GWAS).
#
# Completes the causal-inference chain that coloc.abf (Phase 1c) deliberately
# does not provide (coloc tests "consistent with a shared variant", not
# "X causes Y"). Two designs, depending on how many independent instruments
# are available for a gene:
#
#   (a) Lead-SNP Wald ratio: beta_MR = beta_GWAS / beta_eQTL at the single
#       most-significant eQTL SNP. Always computable if an instrument exists,
#       but cannot test for horizontal pleiotropy and is sensitive to weak-
#       instrument bias from that one SNP.
#
#   (b) Multi-SNP IVW, chr19 locus only: uses the 3 SNPs GCTA-COJO already
#       established as independent genome-wide-significant signals at this
#       locus (Phase 1b) as instruments for PVRL2/APOE/TOMM40/APOC1 -- a
#       principled instrument set (independence already proven), not an
#       ad hoc clump. Reported only where >=2 of the 3 SNPs are present with
#       a usable eQTL effect for that specific gene.
#
# Interpretation: beta_MR is in log-odds of LBD per unit of genetically-
# predicted normalized gene expression change. This is NOT proof of
# causality on its own (still subject to pleiotropy at the single-instrument
# level); it is one additional, standard piece of evidence layered on top of
# coloc.abf (Phase 1c), which already establishes which gene/tissue pairs are
# worth testing this way (i.e. this script targets MMRN1 specifically because
# Phase 1c found PP4=0.97 for it in T cells, not because it was fished for).
################################################################################

suppressPackageStartupMessages({
  library(MendelianRandomization)
  library(yaml)
})

cat("================================================================================\n")
cat("Phase 1e: Mendelian Randomization (eQTL -> LBD GWAS)\n")
cat("================================================================================\n\n")

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
cfg <- yaml::read_yaml(file.path(project_dir, "config/pipeline_config.yaml"))
coloc_cfg <- cfg$phase1c_coloc

gwas_dir <- file.path(project_dir, "data/coloc/gwas")
results_dir <- file.path(project_dir, "data/mr/results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

n_total <- coloc_cfg$gwas$n_total

# The 3 GCTA-COJO-verified independent signals at the chr19 APOE locus
# (Phase 1b, data/gcta_cojo/chr19_stepwise.jma.cojo) -- used as a principled
# multi-instrument set for chr19-locus genes only.
chr19_cojo_snps <- c("rs10405693", "rs769449", "rs1081105")
chr19_genes <- c("PVRL2", "APOE", "TOMM40", "APOC1")

load_merged <- function(symbol, source_name) {
  gwas_file <- file.path(gwas_dir, paste0(symbol, "_gwas.tsv"))
  eqtl_file <- file.path(project_dir, "data/coloc/eqtl", source_name, paste0(symbol, "_eqtl_raw.tsv"))
  if (!file.exists(gwas_file) || !file.exists(eqtl_file)) return(NULL)

  gwas <- read.delim(gwas_file, stringsAsFactors = FALSE)
  eqtl <- read.delim(eqtl_file, stringsAsFactors = FALSE)
  if (nrow(gwas) == 0 || nrow(eqtl) == 0) return(NULL)

  eqtl <- eqtl[!duplicated(eqtl$rsid), ]
  gwas <- gwas[!duplicated(gwas$rsid), ]
  merged <- merge(gwas, eqtl, by = "rsid", suffixes = c(".gwas", ".eqtl"))
  if (nrow(merged) == 0) return(NULL)

  allele_ok <- mapply(function(ea, oa, ref, alt) {
    setequal(c(toupper(ea), toupper(oa)), c(toupper(ref), toupper(alt)))
  }, merged$effect_allele, merged$other_allele, merged$ref, merged$alt)
  merged <- merged[allele_ok, ]

  # Align GWAS effect allele to the eQTL ALT (effect) allele -- flip GWAS
  # beta sign if the eQTL's alt allele is the GWAS's "other" allele.
  flip <- toupper(merged$effect_allele) != toupper(merged$alt)
  merged$beta.gwas[flip] <- -merged$beta.gwas[flip]
  merged
}

wald_ratio_mr <- function(beta_gwas, se_gwas, beta_eqtl, se_eqtl) {
  b <- beta_gwas / beta_eqtl
  se <- abs(se_gwas / beta_eqtl)
  p <- 2 * pnorm(-abs(b / se))
  list(beta = b, se = se, p = p, method = "Wald ratio (1 SNP)", n_snps = 1)
}

run_gene_tissue <- function(symbol, source_name) {
  cat(sprintf("-- %s / %s --\n", symbol, source_name))
  merged <- load_merged(symbol, source_name)
  if (is.null(merged)) {
    cat("   SKIPPED: no merged GWAS+eQTL data\n\n")
    return(NULL)
  }

  # Primary instrument: lead eQTL SNP (smallest eQTL p-value), require at
  # least suggestive significance (p<1e-5) to avoid a near-null instrument.
  # ("pvalue" appears in both input files -> merge() suffixed it to
  # pvalue.gwas/pvalue.eqtl; same for beta/se, handled elsewhere as .gwas/.eqtl)
  merged_sorted <- merged[order(merged$pvalue.eqtl), ]
  lead <- merged_sorted[1, ]
  if (lead$pvalue.eqtl >= 1e-5) {
    cat(sprintf("   SKIPPED: no eQTL SNP reaches p<1e-5 (best p=%.2e)\n\n", lead$pvalue.eqtl))
    return(data.frame(symbol = symbol, eqtl_source = source_name, method = "none",
                       n_snps = 0, beta_MR = NA, se_MR = NA, p_MR = NA,
                       lead_snp = NA, note = "no instrument at p<1e-5"))
  }

  wald <- wald_ratio_mr(lead$beta.gwas, lead$se.gwas, lead$beta.eqtl, lead$se.eqtl)
  cat(sprintf("   Wald ratio (lead SNP %s, eQTL p=%.2e): beta=%.3f, se=%.3f, p=%.3g\n",
              lead$rsid, lead$pvalue.eqtl, wald$beta, wald$se, wald$p))

  result_row <- data.frame(
    symbol = symbol, eqtl_source = source_name, method = wald$method,
    n_snps = 1, beta_MR = wald$beta, se_MR = wald$se, p_MR = wald$p,
    lead_snp = lead$rsid, note = NA
  )

  # Additional multi-SNP IVW for chr19-locus genes using the GCTA-COJO
  # independent signals, where available with a usable eQTL effect.
  if (symbol %in% chr19_genes) {
    cojo_rows <- merged[merged$rsid %in% chr19_cojo_snps, ]
    if (nrow(cojo_rows) >= 2) {
      mr_obj <- mr_input(
        bx = cojo_rows$beta.eqtl, bxse = cojo_rows$se.eqtl,
        by = cojo_rows$beta.gwas, byse = cojo_rows$se.gwas,
        snps = cojo_rows$rsid
      )
      ivw <- tryCatch(mr_ivw(mr_obj), error = function(e) NULL)
      if (!is.null(ivw)) {
        cat(sprintf("   IVW (%d GCTA-COJO instruments: %s): beta=%.3f, se=%.3f, p=%.3g\n",
                    nrow(cojo_rows), paste(cojo_rows$rsid, collapse = ","),
                    ivw@Estimate, ivw@StdError, ivw@Pvalue))
        result_row <- rbind(result_row, data.frame(
          symbol = symbol, eqtl_source = source_name,
          method = sprintf("IVW (%d GCTA-COJO instruments)", nrow(cojo_rows)),
          n_snps = nrow(cojo_rows), beta_MR = ivw@Estimate, se_MR = ivw@StdError,
          p_MR = ivw@Pvalue, lead_snp = NA, note = NA
        ))
      }
    }
  }
  cat("\n")
  result_row
}

all_results <- list()
for (gene in coloc_cfg$target_genes) {
  for (source_name in c("brain_acc", "blueprint_tcell")) {
    res <- run_gene_tissue(gene$symbol, source_name)
    if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
  }
}

results <- do.call(rbind, all_results)
results$padj_BH <- p.adjust(results$p_MR, method = "BH")
results <- results[order(results$p_MR), ]

out_file <- file.path(results_dir, "mendelian_randomization_summary.csv")
write.csv(results, out_file, row.names = FALSE)

cat("================================================================================\n")
cat("Resumen MR (ordenado por p sin corregir)\n")
cat("================================================================================\n")
print(results[, c("symbol", "eqtl_source", "method", "n_snps", "beta_MR", "p_MR", "padj_BH")])
cat(sprintf("\nCompleto: %s\n\n", out_file))

cat("--- R Session Info ---\n")
print(sessionInfo())
