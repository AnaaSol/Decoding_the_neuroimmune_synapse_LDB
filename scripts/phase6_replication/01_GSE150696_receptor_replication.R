#!/usr/bin/env Rscript
################################################################################
# Phase 6 / step 1: independent replication of the Phase 3b severity trend.
#
# Discovery (GSE178146, anterior cingulate snRNA-seq, neuron pseudobulk):
#   LIFR, LGALS9, CD86 rise monotonically with severity (padj < 0.05).
# Replication cohort: GSE150696 (Affymetrix HTA 2.0, prefrontal cortex BA9,
#   Control n=9, PDD n=12, DLB n=12, AD n=9) -- different donors, brain
#   region, platform and tissue resolution (bulk, not neuron-specific).
#
# Pre-specified test (same severity coding as discovery, no PD group here):
#   expr ~ severity (Control=0, PDD=2, DLB=3) + age + sex + PMI   [limma]
# Replication criterion (3 pre-specified genes): same sign as discovery AND
#   one-sided p < 0.05/3. AD vs Control is fit in the same model family as a
#   disease-specificity contrast (is any effect LBD-specific or generic
#   neurodegeneration?).
#
# Gene-level expression: median of exon-level PSR probesets per gene symbol
# (series matrix is Partek RMA log2, exon-level). Bulk-tissue signal mixes
# cell-type composition with per-cell expression -- see the neuron-marker
# adjustment sensitivity model below.
################################################################################
suppressPackageStartupMessages({ library(data.table); library(limma) })

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
raw_dir <- file.path(project_dir, "data/GSE150696_raw")
out_dir <- file.path(project_dir, "data/replication_GSE150696")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- sample metadata --------------------------------------------------------
hdr <- readLines(gzfile(file.path(raw_dir, "GSE150696_series_matrix.txt.gz")), n = 200)
field <- function(prefix) {
  l <- grep(paste0("^", prefix), hdr, value = TRUE)
  lapply(l, function(x) gsub('"', "", strsplit(x, "\t")[[1]][-1]))
}
title <- field("!Sample_title")[[1]]
gsm <- field("!Sample_geo_accession")[[1]]
ch <- field("!Sample_characteristics_ch1")
getch <- function(key) { v <- ch[[which(sapply(ch, function(x) startsWith(tolower(x[1]), key)))]]; sub("^[^:]*: ", "", v) }
meta <- data.table(gsm = gsm, title = title, group = sub(".*_", "", title),
                   age = as.numeric(getch("age")), sex = getch("sex"),
                   pmi = suppressWarnings(as.numeric(getch("pmi"))))
# one PDD donor (GSM4556883) has PMI recorded as "." in GEO -> median-impute
cat(sprintf("Median-imputing missing PMI for: %s\n", paste(meta[is.na(pmi), gsm], collapse = ", ")))
meta[is.na(pmi), pmi := median(meta$pmi, na.rm = TRUE)]
meta[, group := factor(group, levels = c("CTRL", "PDD", "DLB", "AD"))]
print(meta[, .N, by = group])

# ---- expression: PSR probesets -> gene median --------------------------------
ann <- fread(cmd = paste("zcat", file.path(raw_dir, "GPL17585_table.txt.gz"), "| grep -v '^[!#^]'"),
             select = c("ID", "gene_symbols"), quote = "")
ann <- ann[startsWith(ID, "PSR") & gene_symbols != "---" & !grepl(" ", gene_symbols)]
ex <- fread(cmd = paste("zcat", file.path(raw_dir, "GSE150696_series_matrix.txt.gz"), "| grep -v '^!' | grep -E '^\"(ID_REF|PSR)'"))
setnames(ex, "ID_REF", "ID")
ex <- merge(ann, ex, by = "ID")
gene_expr <- ex[, lapply(.SD, median), by = gene_symbols, .SDcols = gsm]
mat <- as.matrix(gene_expr[, ..gsm]); rownames(mat) <- gene_expr$gene_symbols
cat(sprintf("Gene-level matrix: %d genes x %d samples\n", nrow(mat), ncol(mat)))
# drop genes not expressed above array background in most samples
mat <- mat[rowMeans(mat > log2(32)) >= 0.5, ]
cat(sprintf("After expression filter: %d genes\n", nrow(mat)))
fwrite(data.table(gene = rownames(mat), mat), file.path(out_dir, "GSE150696_gene_level_log2.tsv.gz"), sep = "\t")

# ---- models -----------------------------------------------------------------
lbd <- meta$group %in% c("CTRL", "PDD", "DLB")
m <- meta[lbd]; m[, severity := c(CTRL = 0, PDD = 2, DLB = 3)[as.character(group)]]
fit_trend <- function(design_extra = NULL) {
  X <- model.matrix(~ severity + age + sex + pmi, data = m)
  if (!is.null(design_extra)) X <- cbind(X, design_extra)
  f <- eBayes(lmFit(mat[, m$gsm], X))
  tt <- topTable(f, coef = "severity", number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), logFC_per_step = tt$logFC, t = tt$t, p = tt$P.Value, padj = tt$adj.P.Val)
}
trend <- fit_trend()

# neuron-content sensitivity: bulk tissue confounds per-cell change with loss of
# neurons, so re-fit adjusting for a neuronal marker score (mean z of core
# pan-neuronal genes).
neuro_genes <- intersect(c("RBFOX3", "SNAP25", "SYT1", "STMN2", "GAP43", "NEFL"), rownames(mat))
neuro_score <- colMeans(t(scale(t(mat[neuro_genes, m$gsm]))))
trend_adj <- fit_trend(cbind(neuron_score = neuro_score))

# categorical contrasts incl. AD (all 42 samples)
Xc <- model.matrix(~ 0 + group + age + sex + pmi, data = meta)
colnames(Xc) <- sub("^group", "", colnames(Xc))
fc <- eBayes(contrasts.fit(lmFit(mat[, meta$gsm], Xc),
                           makeContrasts(PDD_vs_CTRL = PDD - CTRL, DLB_vs_CTRL = DLB - CTRL,
                                         AD_vs_CTRL = AD - CTRL, levels = Xc)))
cat_res <- rbindlist(lapply(c("PDD_vs_CTRL", "DLB_vs_CTRL", "AD_vs_CTRL"), function(k) {
  tt <- topTable(fc, coef = k, number = Inf, sort.by = "none")
  data.table(gene = rownames(tt), contrast = k, logFC = tt$logFC, p = tt$P.Value, padj = tt$adj.P.Val)
}))
cat_wide <- dcast(cat_res, gene ~ contrast, value.var = c("logFC", "p"))

# ---- replication table ------------------------------------------------------
disc <- fread(file.path(project_dir, "data/GSE178146_processed/receptor_doseresponse_pseudobulk.csv"))
panel <- disc[, .(gene, disc_log2FC = log2FC_per_severity_unit, disc_padj = padj_trend)]
rep <- merge(panel, trend, by = "gene", all.x = TRUE)
rep <- merge(rep, trend_adj[, .(gene, logFC_adj_neuron = logFC_per_step, p_adj_neuron = p)], by = "gene", all.x = TRUE)
rep <- merge(rep, cat_wide, by = "gene", all.x = TRUE)
rep[, prespecified := gene %in% c("LIFR", "LGALS9", "CD86")]
rep[, same_direction := sign(logFC_per_step) == sign(disc_log2FC)]
rep[, p_one_sided := ifelse(same_direction, p / 2, 1 - p / 2)]
rep[, replicated := prespecified & same_direction & p_one_sided < 0.05 / 3]
setorder(rep, -prespecified, p)
fwrite(rep, file.path(out_dir, "receptor_replication_GSE150696.csv"))
fwrite(trend[order(p)], file.path(out_dir, "genomewide_severity_trend_GSE150696.csv"))
print(rep[, .(gene, prespecified, disc_log2FC = round(disc_log2FC, 2), rep_logFC = round(logFC_per_step, 3),
              p = signif(p, 2), p_one_sided = signif(p_one_sided, 2), logFC_adj_neuron = round(logFC_adj_neuron, 3),
              p_adj_neuron = signif(p_adj_neuron, 2), AD_logFC = round(logFC_AD_vs_CTRL, 3), AD_p = signif(p_AD_vs_CTRL, 2),
              replicated)])

# ---- plot --------------------------------------------------------------------
pdf(file.path(out_dir, "replication_prespecified_receptors.pdf"), width = 10, height = 3.8)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
cols <- c(CTRL = "#3498db", PDD = "#e67e22", DLB = "#e74c3c", AD = "#7f8c8d")
for (g in c("LIFR", "LGALS9", "CD86")) {
  if (!g %in% rownames(mat)) { plot.new(); title(paste(g, "(not expressed)")); next }
  y <- mat[g, meta$gsm]
  boxplot(y ~ meta$group, col = adjustcolor(cols, 0.35), border = cols, outline = FALSE,
          ylab = "log2 expression (RMA)", xlab = "", main = g)
  stripchart(y ~ meta$group, vertical = TRUE, method = "jitter", pch = 16, col = cols, add = TRUE)
  r <- rep[gene == g]
  mtext(sprintf("trend p(1-sided)=%.3g", r$p_one_sided), side = 3, line = 0.1, cex = 0.7)
}
invisible(dev.off())
cat("Done.\n")
