#!/usr/bin/env Rscript
################################################################################
# Phase 5 / step 7: SuSiE fine-mapping + coloc.susie for every strong (PP4>0.8)
# atlas colocalization.
#
# coloc.abf assumes one causal variant per trait per region. SuSiE allows
# several (up to L), so coloc.susie can test colocalization credible set (CS)
# by credible set: does the specific GWAS signal that defines the locus share
# its causal variant with the eQTL, even if other signals exist nearby?
#
# Key implementation points:
#  - window: lead SNP +/- 250 kb (fine-mapping is quadratic in #SNPs)
#  - LD: 1000G EUR (n=503), signed r from PLINK. Both betas are aligned to the
#    reference panel's counted allele (bim A1) -- the signed LD must refer to
#    the same allele as the z-scores. (The earlier MMRN1-only script,
#    phase1_genetic_prior/14, did NOT align to A1; its inflated CS count is
#    consistent with that mismatch.)
#  - L = 5 signals max per trait
#  - diagnostic: susieR::estimate_s_rss(z, R, n); s near 0 = z-scores and LD
#    consistent, large s (> ~0.3) = out-of-sample LD mismatch, CSs unreliable
#  - if SuSiE finds no CS for a trait, coloc.susie cannot run -> reported as
#    such (the coloc.abf result then stands alone).
# Output: data/coloc_atlas/susie/susie_coloc_summary.tsv
################################################################################
suppressPackageStartupMessages({ library(coloc); library(susieR); library(data.table); library(yaml) })

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
cfg <- read_yaml(file.path(project_dir, "config/pipeline_config.yaml"))
acfg <- cfg$phase5_coloc_atlas
gw <- cfg$phase1c_coloc$gwas
out_dir <- file.path(project_dir, "data/coloc_atlas/susie")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ref <- file.path(project_dir, "data/reference/g1000_eur")
window <- 250000
L_max <- 5

strong <- fread(file.path(project_dir, "data/coloc_atlas/results/coloc_atlas_strong_PP4gt0.8.tsv"))
loci <- setNames(acfg$loci, sapply(acfg$loci, `[[`, "locus"))
ds_n <- setNames(sapply(acfg$datasets, `[[`, "sample_size"), sapply(acfg$datasets, `[[`, "id"))

run_pair <- function(row) {
  l <- loci[[row$locus]]
  g <- fread(file.path(project_dir, "data/coloc_atlas/gwas", paste0(row$locus, ".tsv")))
  e <- fread(cmd = paste("zcat", shQuote(file.path(project_dir, "data/coloc_atlas/eqtl", row$dataset,
                                                     paste0(row$locus, ".tsv.gz")))))
  e <- e[gene_id == row$gene_id & se > 0][!duplicated(pos)]
  g <- g[se > 0 & eaf > 0 & eaf < 1][!duplicated(pos)]
  m <- merge(e, g, by = "pos", suffixes = c(".e", ".g"))
  m <- m[abs(pos - l$lead_pos) <= window]
  m <- m[(alt == effect_allele & ref == other_allele) | (alt == other_allele & ref == effect_allele)]
  m[, beta.g := ifelse(alt == effect_allele, beta.g, -beta.g)]           # GWAS beta now per eQTL alt allele
  m <- m[grepl("^rs", rsid.e) & !duplicated(rsid.e)]

  # reference panel: extract by rsID, align alt allele to bim A1
  tag <- paste(row$locus, row$dataset, row$symbol, sep = "_")
  snp_file <- file.path(out_dir, paste0(tag, ".snps"))
  writeLines(m$rsid.e, snp_file)
  pre <- file.path(out_dir, tag)
  system2("plink", c("--bfile", ref, "--extract", snp_file, "--make-bed", "--out", pre), stdout = FALSE, stderr = FALSE)
  bim <- fread(paste0(pre, ".bim"), header = FALSE, col.names = c("chr", "rsid", "cm", "bp", "A1", "A2"))
  m <- merge(m, bim[, .(rsid.e = rsid, A1, A2)], by = "rsid.e")
  m <- m[(alt == A1 & ref == A2) | (alt == A2 & ref == A1)]
  m[, sgn := ifelse(alt == A1, 1, -1)]                                     # both betas -> per A1 allele
  m[, `:=`(bg = sgn * beta.g, be = sgn * beta.e)]
  m <- m[match(bim$rsid, rsid.e)][!is.na(rsid.e)]                           # order = bim order
  writeLines(m$rsid.e, snp_file)
  system2("plink", c("--bfile", pre, "--extract", snp_file, "--r", "square", "--out", pre), stdout = FALSE, stderr = FALSE)
  R <- as.matrix(fread(paste0(pre, ".ld"), header = FALSE))
  bim2 <- fread(paste0(pre, ".bim"), header = FALSE)$V2
  keep_bim <- bim2 %in% m$rsid.e
  R <- R[keep_bim, keep_bim]
  stopifnot(nrow(R) == nrow(m), all(bim2[keep_bim] == m$rsid.e))
  R[is.na(R)] <- 0
  dimnames(R) <- list(m$rsid.e, m$rsid.e)
  unlink(list.files(out_dir, pattern = paste0("^", gsub("([+.])", "\\\\\\1", tag), "\\.(bed|bim|fam|nosex|log|ld)$"), full.names = TRUE))

  n_e <- ds_n[[row$dataset]]
  z_g <- m$bg / m$se.g; z_e <- m$be / m$se.e
  s_g <- estimate_s_rss(z_g, R, n = gw$n_total)
  s_e <- estimate_s_rss(z_e, R, n = n_e)

  d1 <- list(beta = m$bg, varbeta = m$se.g^2, type = "cc", s = gw$n_cases / gw$n_total, N = gw$n_total,
             MAF = pmin(m$eaf, 1 - m$eaf), LD = R, snp = m$rsid.e, position = m$pos)
  d2 <- list(beta = m$be, varbeta = m$se.e^2, type = "quant", N = n_e, MAF = m$maf,
             LD = R, snp = m$rsid.e, position = m$pos)
  fit <- function(d) tryCatch(suppressWarnings(suppressMessages(
    runsusie(d, maxit = 1000, repeat_until_convergence = FALSE, L = L_max))), error = function(err) NULL)
  s1 <- fit(d1); s2 <- fit(d2)
  ncs <- function(s) if (is.null(s) || is.null(s$sets$cs)) 0L else length(s$sets$cs)
  res <- data.table(locus = row$locus, gene = row$symbol, dataset = row$dataset, study = row$study,
                    group = row$group, tier = row$tier, abf_PP4 = row$PP4, abf_top_snp = row$top_snp,
                    n_snps = nrow(m), lambda_s_gwas = s_g, lambda_s_eqtl = s_e,
                    n_cs_gwas = ncs(s1), n_cs_eqtl = ncs(s2),
                    lead_snp_in_panel = l$lead_snp %in% m$rsid.e)
  # which GWAS CS contains the locus lead SNP (or the coloc.abf top SNP)?
  cs_has <- function(s, snp) {
    if (ncs(s) == 0) return(NA_integer_)
    hit <- which(sapply(s$sets$cs, function(ix) snp %in% m$rsid.e[ix]))
    if (length(hit)) as.integer(sub("L", "", names(s$sets$cs)[hit[1]])) else NA_integer_
  }
  res[, gwas_cs_with_lead := cs_has(s1, l$lead_snp)]
  if (ncs(s1) > 0 && ncs(s2) > 0) {
    cs <- as.data.table(coloc.susie(s1, s2)$summary)
    fwrite(cs, file.path(out_dir, paste0(tag, "_coloc_susie_pairs.tsv")), sep = "\t")
    best <- cs[which.max(PP.H4.abf)]
    lead_pairs <- if (!is.na(res$gwas_cs_with_lead)) cs[idx1 == res$gwas_cs_with_lead] else cs[0]
    res[, `:=`(susie_max_PP4 = best$PP.H4.abf, susie_best_hit_gwas = best$hit1, susie_best_hit_eqtl = best$hit2,
               susie_PP4_leadCS = if (nrow(lead_pairs)) max(lead_pairs$PP.H4.abf) else NA_real_)]
  } else {
    res[, `:=`(susie_max_PP4 = NA_real_, susie_best_hit_gwas = NA_character_, susie_best_hit_eqtl = NA_character_,
               susie_PP4_leadCS = NA_real_)]
  }
  res
}

out <- rbindlist(lapply(seq_len(nrow(strong)), function(i) {
  r <- strong[i]
  cat(sprintf("[%d/%d] %s %s %s:%s ... ", i, nrow(strong), r$locus, r$symbol, r$study, r$group))
  x <- run_pair(r)
  cat(sprintf("nSNP=%d, CS gwas/eqtl=%d/%d, s=%.2f/%.2f, SuSiE max PP4=%s\n", x$n_snps, x$n_cs_gwas, x$n_cs_eqtl,
              x$lambda_s_gwas, x$lambda_s_eqtl, format(round(x$susie_max_PP4, 3))))
  x
}), fill = TRUE)
out[, verdict := fifelse(is.na(susie_max_PP4),
                         fifelse(n_cs_gwas == 0, "no GWAS credible set (SuSiE underpowered) - coloc.abf only",
                                 "no eQTL credible set - coloc.abf only"),
                         fifelse(susie_max_PP4 > 0.8, "confirmed by SuSiE-coloc",
                                 fifelse(susie_max_PP4 > 0.5, "suggestive under SuSiE-coloc", "not supported by SuSiE-coloc")))]
fwrite(out, file.path(out_dir, "susie_coloc_summary.tsv"), sep = "\t")
print(out[, .(locus, gene, group, tier, abf_PP4 = round(abf_PP4, 3), n_snps, cs = paste0(n_cs_gwas, "/", n_cs_eqtl),
              s_gwas = round(lambda_s_gwas, 2), s_eqtl = round(lambda_s_eqtl, 2),
              susie_PP4 = round(susie_max_PP4, 3), best = paste(susie_best_hit_gwas, susie_best_hit_eqtl), verdict)])
cat("Done.\n")
