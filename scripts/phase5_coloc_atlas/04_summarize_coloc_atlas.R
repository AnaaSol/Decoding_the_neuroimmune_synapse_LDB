#!/usr/bin/env Rscript
################################################################################
# Phase 5 / step 4: annotate + summarize the colocalization atlas.
#
#  - maps Ensembl gene IDs -> symbols (GENCODE v39)
#  - calls a colocalization "strong" at PP4 > 0.8 (primary, default priors) and
#    reports whether it survives the conservative prior (p12 = 5e-6)
#  - heatmap: every locus-gene reaching PP4 > 0.5 in any dataset x all datasets
#  - per-locus best hit and immune-vs-brain breakdown
################################################################################
suppressPackageStartupMessages({ library(data.table); library(ggplot2); })

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
res_dir <- file.path(project_dir, "data/coloc_atlas/results")
cfg <- yaml::read_yaml(file.path(project_dir, "config/pipeline_config.yaml"))$phase5_coloc_atlas
d <- fread(file.path(res_dir, "coloc_atlas_all.tsv"))

# ---- gene symbols: GENCODE v39 (the annotation eQTL Catalogue uses) --------
# built by streaming the GTF's gene records only:
#   curl .../gencode.v39.basic.annotation.gtf.gz | zcat | awk '$3=="gene"' -> data/reference/gencode_v39_genes.tsv
sym <- fread(file.path(project_dir, "data/reference/gencode_v39_genes.tsv"))
sym <- sym[!duplicated(gene_id)]
d <- merge(d, sym, by = "gene_id", all.x = TRUE)
d[is.na(symbol), symbol := gene_id]
d[, label := paste0(locus, ": ", symbol)]
fwrite(d[order(-PP4)], file.path(res_dir, "coloc_atlas_annotated.tsv"), sep = "\t")

# ---- strong colocalizations -------------------------------------------------
strong <- d[PP4 > 0.8][order(locus, -PP4)]
strong[, robust_conservative_prior := PP4_p12_5e6 > 0.8]
# Does the colocalizing variant tag the GWAS locus lead signal? coloc.abf can
# return high PP4 for a weaker secondary association in the window (e.g. at
# GBA, where the rare lead rs2230288/E326K is absent from most eQTL panels).
# r2 to the lead from 1000G EUR (same reference panel as Phase 1b/1f).
ref <- file.path(project_dir, "data/reference/g1000_eur")
ld_tmp <- tempfile()
lead_ld <- rbindlist(lapply(cfg$loci, function(l) {
  system2("plink", c("--bfile", ref, "--r2", "--ld-snp", l$lead_snp, "--ld-window-kb", "1000",
                     "--ld-window", "99999", "--ld-window-r2", "0", "--out", ld_tmp), stdout = FALSE, stderr = FALSE)
  x <- fread(paste0(ld_tmp, ".ld"))
  data.table(locus = l$locus, top_snp = x$SNP_B, r2_to_lead = x$R2)
}))
strong <- merge(strong, lead_ld, by = c("locus", "top_snp"), all.x = TRUE)
strong[, tags_lead_signal := !is.na(r2_to_lead) & r2_to_lead >= 0.5]
strong[, tier := fifelse(tags_lead_signal & robust_conservative_prior, "Tier 1: lead signal, prior-robust",
                  fifelse(tags_lead_signal, "Tier 2: lead signal", "Secondary / non-lead signal"))]
setorder(strong, locus, -PP4)
strong[, direction := ifelse(eqtl_beta_risk_allele > 0, "risk allele UP", "risk allele DOWN")]
fwrite(strong, file.path(res_dir, "coloc_atlas_strong_PP4gt0.8.tsv"), sep = "\t")
cat(sprintf("Tests: %d (datasets=%d, loci=%d, genes=%d)\n", nrow(d), uniqueN(d$dataset), uniqueN(d$locus), uniqueN(d$gene_id)))
cat(sprintf("Strong (PP4>0.8): %d tests, %d distinct locus-genes\n", nrow(strong), uniqueN(strong$label)))
print(strong[, .(label, study, group, class, PP4 = round(PP4, 3), r2 = round(r2_to_lead, 2), tier,
                 PP4_cons = round(PP4_p12_5e6, 3), top_snp, direction)])

# per locus-gene: which tissue classes colocalize
by_gene <- d[, .(n_tested = .N, n_immune_tested = sum(class == "immune"), n_brain_tested = sum(class != "immune"),
                 max_PP4 = max(PP4), n_strong = sum(PP4 > 0.8),
                 n_strong_immune = sum(PP4 > 0.8 & class == "immune"),
                 n_strong_brain = sum(PP4 > 0.8 & class != "immune"),
                 best = group[which.max(PP4)], best_study = study[which.max(PP4)]),
             by = .(locus, gene_id, symbol, biotype)][order(-max_PP4)]
fwrite(by_gene, file.path(res_dir, "coloc_atlas_by_gene.tsv"), sep = "\t")
print(by_gene[max_PP4 > 0.5])

# ---- heatmap ------------------------------------------------------------------
keep <- by_gene[max_PP4 > 0.5, paste0(locus, ": ", symbol)]
ds_meta <- rbindlist(lapply(cfg$datasets, as.data.table))
ds_meta[, ds_label := paste0(study, " | ", group)]
ds_meta[, class := factor(class, levels = c("immune", "microglia", "brain"))]
ds_order <- ds_meta[order(class, study, group), ds_label]
h <- d[label %in% keep]
h <- merge(h, ds_meta[, .(dataset = id, ds_label)], by = "dataset")
grid <- CJ(label = keep, ds_label = ds_order)
h <- merge(grid, h[, .(label, ds_label, PP4, eqtl_beta_risk_allele)], by = c("label", "ds_label"), all.x = TRUE)
h[, ds_label := factor(ds_label, levels = ds_order)]
# flag locus-genes whose only strong colocalizations are with a non-lead signal
secondary <- setdiff(strong[tier == "Secondary / non-lead signal", label], strong[tags_lead_signal == TRUE, label])
h[label %in% secondary, label := paste0(label, " \u2020")]
keep <- ifelse(keep %in% secondary, paste0(keep, " \u2020"), keep)
h[, label := factor(label, levels = rev(sort(keep)))]
h[, mark := ifelse(!is.na(PP4) & PP4 > 0.8, ifelse(eqtl_beta_risk_allele > 0, "▲", "▼"), "")]
p <- ggplot(h, aes(ds_label, label, fill = PP4)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = mark), size = 2.4, colour = "white") +
  scale_fill_gradient(low = "#e8eef5", high = "#1f4e8c", na.value = "grey95", limits = c(0, 1), name = "PP4") +
  labs(x = NULL, y = NULL,
       title = "LBD risk-locus colocalization atlas (coloc.abf)",
       subtitle = "Grey = gene not tested in that dataset. ▲/▼: risk allele raises/lowers expression (PP4 > 0.8). \u2020 = secondary signal, not the GWAS lead (r\u00B2 < 0.5)") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6), panel.grid = element_blank())
ggsave(file.path(res_dir, "coloc_atlas_heatmap.pdf"), p, width = 13, height = 2 + 0.28 * length(keep))
ggsave(file.path(res_dir, "coloc_atlas_heatmap.png"), p, width = 13, height = 2 + 0.28 * length(keep), dpi = 200)
cat("Done.\n")
