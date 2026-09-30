#!/usr/bin/env Rscript
################################################################################
# SuSiE-coloc for MMRN1: does the PP4=0.97 coloc.abf result (Phase 1c, CD4+ T
# cell eQTL) hold up under a method that allows >1 causal variant per locus?
#
# coloc.abf (Phase 1c) assumes exactly one causal variant per dataset per
# locus. If MMRN1's region actually has 2 independent GWAS or eQTL signals,
# coloc.abf's PP4 could be a blend across them rather than reflecting one
# clean shared signal. coloc::runsusie() + coloc::coloc.susie() fit SuSiE
# (multi-variant) fine-mapping separately on each dataset using an LD matrix
# from the 1000G EUR reference panel, then test colocalization credible-set
# by credible-set.
################################################################################

suppressPackageStartupMessages({
  library(coloc)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
source_name <- if (length(args) >= 1) args[1] else "blueprint_tcell"

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
cfg <- yaml::read_yaml(file.path(project_dir, "config/pipeline_config.yaml"))
coloc_cfg <- cfg$phase1c_coloc
src_cfg <- if (source_name == "brain_acc") coloc_cfg$eqtl_source else
  Filter(function(s) s$name == source_name, coloc_cfg$eqtl_sources_extra)[[1]]

out_dir <- file.path(project_dir, "data/susie_coloc")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("================================================================================\n")
cat(sprintf("Phase 1f: SuSiE-coloc, MMRN1 vs %s\n", source_name))
cat("================================================================================\n\n")

gwas <- read.delim(file.path(project_dir, "data/coloc/gwas/MMRN1_gwas.tsv"), stringsAsFactors = FALSE)
eqtl <- read.delim(file.path(project_dir, "data/coloc/eqtl", source_name, "MMRN1_eqtl_raw.tsv"), stringsAsFactors = FALSE)
eqtl <- eqtl[!duplicated(eqtl$rsid), ]
gwas <- gwas[!duplicated(gwas$rsid), ]

merged <- merge(gwas, eqtl, by = "rsid", suffixes = c(".gwas", ".eqtl"))
allele_ok <- mapply(function(ea, oa, ref, alt) {
  setequal(c(toupper(ea), toupper(oa)), c(toupper(ref), toupper(alt)))
}, merged$effect_allele, merged$other_allele, merged$ref, merged$alt)
merged <- merged[allele_ok, ]
flip <- toupper(merged$effect_allele) != toupper(merged$alt)
merged$beta.gwas[flip] <- -merged$beta.gwas[flip]
cat(sprintf("Merged SNP set (rsid+allele matched, +/-1Mb scan window): %d SNPs\n", nrow(merged)))

# SuSiE fine-mapping is computationally heavy in the number of variants (LD
# matrix ops); restrict to a tighter +/-200kb window around the MMRN1 gene
# body (chr4:89,894,876-89,954,629) for tractability. This is standard
# practice for fine-mapping (vs. the wider scan window used for coloc.abf)
# and still comfortably covers the credible set region.
susie_window <- 200000
gene_start <- 89894876; gene_end <- 89954629
merged <- merged[merged$position >= (gene_start - susie_window) &
                  merged$position <= (gene_end + susie_window), ]
cat(sprintf("After restricting to +/-%dkb SuSiE window: %d SNPs\n", susie_window / 1000, nrow(merged)))

# --- Step 1: extract 1000G EUR genotypes for this exact SNP set (by rsID) ---
rsid_file <- file.path(out_dir, "mmrn1_rsids.txt")
writeLines(merged$rsid, rsid_file)

ref_prefix <- file.path(out_dir, "mmrn1_ref")
system2("plink", c("--bfile", "data/reference/g1000_eur", "--chr", "4",
                    "--extract", rsid_file, "--make-bed", "--out", ref_prefix),
        stdout = file.path(out_dir, "plink_extract.log"), stderr = file.path(out_dir, "plink_extract.log"))

bim <- read.table(paste0(ref_prefix, ".bim"), stringsAsFactors = FALSE)
ref_snps <- bim$V2
merged <- merged[merged$rsid %in% ref_snps, ]
merged <- merged[match(ref_snps, merged$rsid), ]  # order == bim/LD order
merged <- merged[!is.na(merged$rsid), ]
cat(sprintf("SNPs with reference-panel genotypes (LD-computable): %d\n", nrow(merged)))

# --- Step 2: signed LD correlation matrix (order matches merged/bim) ---
ld_prefix <- file.path(out_dir, "mmrn1_ld")
system2("plink", c("--bfile", ref_prefix, "--r", "square", "--out", ld_prefix),
        stdout = file.path(out_dir, "plink_ld.log"), stderr = file.path(out_dir, "plink_ld.log"))
LD <- as.matrix(read.table(paste0(ld_prefix, ".ld")))
rownames(LD) <- colnames(LD) <- merged$rsid
cat(sprintf("LD matrix: %d x %d\n\n", nrow(LD), ncol(LD)))

n_total <- coloc_cfg$gwas$n_total
s_cases <- coloc_cfg$gwas$n_cases / n_total
n_eqtl <- src_cfg$sample_size

d_gwas <- list(beta = merged$beta.gwas, varbeta = merged$se.gwas^2, type = "cc",
               s = s_cases, N = n_total, MAF = pmin(merged$eaf, 1 - merged$eaf),
               LD = LD, snp = merged$rsid, position = merged$position)
d_eqtl <- list(beta = merged$beta.eqtl, varbeta = merged$se.eqtl^2, type = "quant",
               N = n_eqtl, MAF = merged$maf, LD = LD, snp = merged$rsid,
               position = merged$position)

cat("-- Fitting SuSiE on GWAS dataset --\n")
susie_gwas <- runsusie(d_gwas, maxit = 500, repeat_until_convergence = FALSE)
cat(sprintf("   converged=%s, %d credible set(s) found\n",
            susie_gwas$converged, length(susie_gwas$sets$cs)))

cat("-- Fitting SuSiE on eQTL dataset --\n")
susie_eqtl <- runsusie(d_eqtl, maxit = 500, repeat_until_convergence = FALSE)
cat(sprintf("   converged=%s, %d credible set(s) found\n\n",
            susie_eqtl$converged, length(susie_eqtl$sets$cs)))

cat("-- coloc.susie: per-credible-set-pair colocalization --\n")
res <- coloc.susie(susie_gwas, susie_eqtl)
print(res$summary)

out_file <- file.path(out_dir, sprintf("MMRN1_susie_coloc_%s.csv", source_name))
write.csv(res$summary, out_file, row.names = FALSE)
cat(sprintf("\nGuardado: %s\n\n", out_file))

cat("--- R Session Info ---\n")
print(sessionInfo())
