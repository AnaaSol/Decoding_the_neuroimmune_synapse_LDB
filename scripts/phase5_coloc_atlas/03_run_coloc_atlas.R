#!/usr/bin/env Rscript
################################################################################
# Phase 5 / step 3: colocalization atlas.
#
# For every (dataset x locus x gene) with >=1 eQTL at p < min_eqtl_p inside the
# locus window, run coloc.abf (default priors p1=p2=1e-4, p12=1e-5, identical
# to Phase 1c so results are directly comparable) between the LBD GWAS and the
# eQTL. Alleles are matched on position + ref/alt; strand-ambiguous mismatches
# are dropped, swapped alleles are sign-flipped.
#
# Besides PP0..PP4, records for each test the "direction": sign of the eQTL
# effect of the GWAS RISK-increasing allele at the top H4 SNP (+1 = risk allele
# raises expression, -1 = lowers it). This is what makes the atlas interpretable
# as "risk acts by up/down-regulating gene X in cell type Y".
#
# Output: data/coloc_atlas/results/coloc_atlas_all.tsv (one row per test)
################################################################################
suppressPackageStartupMessages({ library(coloc); library(data.table); library(yaml) })

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
cfg <- yaml::read_yaml(file.path(project_dir, "config/pipeline_config.yaml"))
acfg <- cfg$phase5_coloc_atlas
gw <- cfg$phase1c_coloc$gwas
s_cc <- gw$n_cases / gw$n_total

out_dir <- file.path(project_dir, "data/coloc_atlas/results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gwas <- lapply(setNames(acfg$loci, sapply(acfg$loci, `[[`, "locus")), function(l) {
  g <- fread(file.path(project_dir, "data/coloc_atlas/gwas", paste0(l$locus, ".tsv")))
  g <- g[se > 0 & eaf > 0 & eaf < 1]
  g[!duplicated(pos)]
})

run_one <- function(e, g) {
  m <- merge(e, g, by = "pos", suffixes = c(".e", ".g"))
  same <- m$alt == m$effect_allele & m$ref == m$other_allele
  flip <- m$alt == m$other_allele & m$ref == m$effect_allele
  m <- m[same | flip]
  if (nrow(m) < 50) return(NULL)
  m[, beta.g := ifelse(alt == effect_allele, beta.g, -beta.g)]   # GWAS beta now for eQTL alt allele
  m <- m[!duplicated(pos) & se.e > 0 & maf > 0]
  m[, snp := paste0("chr", pos, "_", ref, "_", alt)]
  d1 <- list(beta = m$beta.g, varbeta = m$se.g^2, type = "cc", s = s_cc, snp = m$snp, position = m$pos)
  d2 <- list(beta = m$beta.e, varbeta = m$se.e^2, type = "quant", N = n_eqtl, MAF = m$maf,
             snp = m$snp, position = m$pos)
  r <- tryCatch(suppressMessages(coloc.abf(d1, d2)), error = function(err) { message("coloc error: ", conditionMessage(err)); NULL })
  if (is.null(r)) return(NULL)
  top <- r$results[which.max(r$results$SNP.PP.H4), ]
  # prior sensitivity: more conservative p12 = 5e-6 (Wallace 2020 recommendation)
  pp_cons <- coloc:::combine.abf(r$results$lABF.df1, r$results$lABF.df2, 1e-4, 1e-4, 5e-6, quiet = TRUE)
  k <- match(top$snp, m$snp)
  # direction for the risk-increasing allele: sign(beta_gwas * beta_eqtl), both on alt allele
  data.table(nsnps = nrow(m),
             PP0 = r$summary[["PP.H0.abf"]], PP1 = r$summary[["PP.H1.abf"]],
             PP2 = r$summary[["PP.H2.abf"]], PP3 = r$summary[["PP.H3.abf"]],
             PP4 = r$summary[["PP.H4.abf"]], PP4_p12_5e6 = pp_cons[["PP.H4.abf"]],
             top_snp = m$rsid.e[k], top_pos = m$pos[k], top_snp_PPH4 = top$SNP.PP.H4,
             top_gwas_p = m$pvalue.g[k], top_eqtl_p = m$pvalue.e[k],
             eqtl_beta_risk_allele = sign(m$beta.g[k]) * m$beta.e[k],
             min_gwas_p = min(m$pvalue.g), min_eqtl_p = min(m$pvalue.e))
}

res <- list()
for (ds in acfg$datasets) {
  n_eqtl <- ds$sample_size
  for (l in acfg$loci) {
    f <- file.path(project_dir, "data/coloc_atlas/eqtl", ds$id, paste0(l$locus, ".tsv.gz"))
    if (!file.exists(f)) next
    e <- fread(cmd = paste("zcat", shQuote(f)))
    if (nrow(e) == 0) next
    e <- e[se > 0]
    genes <- e[, .(minp = min(pvalue)), by = gene_id][minp < as.numeric(acfg$min_eqtl_p), gene_id]
    for (gid in genes) {
      r <- run_one(e[gene_id == gid][!duplicated(pos)], gwas[[l$locus]])
      if (is.null(r)) next
      res[[length(res) + 1]] <- cbind(data.table(dataset = ds$id, study = ds$study, group = ds$group,
                                                 class = ds$class, n = ds$sample_size, locus = l$locus,
                                                 gene_id = gid), r)
    }
  }
  cat(sprintf("%s %s:%s done (%d tests so far)\n", ds$id, ds$study, ds$group, length(res)))
}
all <- rbindlist(res)
fwrite(all, file.path(out_dir, "coloc_atlas_all.tsv"), sep = "\t")
cat(sprintf("Wrote %d coloc tests; %d with PP4 > 0.8\n", nrow(all), sum(all$PP4 > 0.8)))
