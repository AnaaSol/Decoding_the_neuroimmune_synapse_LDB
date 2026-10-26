#!/usr/bin/env Rscript
################################################################################
# Phase 5 / step 5: where are the atlas effector genes expressed in LBD cortex,
# and does their expression track clinical severity?
#
# Effector genes = every locus-gene from the colocalization atlas that tags the
# GWAS lead signal (tiers 1-2 in coloc_atlas_strong_PP4gt0.8.tsv). Read from
# the atlas output, not hard-coded.
#
# Data: GSE178146 integrated object (ACC snRNA-seq, 28 donors, 2000 nuclei per
# donor after Phase 3 downsampling). Broad cell types are assigned with the
# SAME marker-argmax rule as Phase 3 (03_identify_neuronal_subtypes.R, Step 2),
# collapsing Neuron/Excitatory/Inhibitory to "Neuron".
#
# Test: per cell type, donor-level pseudobulk DESeq2 ~ severity_score with the
# same pre-specified coding as Phase 3b (Control=0, PD=1, PDD=2, DLB=3); BH
# across all expressed genes within the cell type, then subset to effectors.
################################################################################
suppressPackageStartupMessages({ library(Seurat); library(DESeq2); library(data.table); library(ggplot2) })

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
res_dir <- file.path(project_dir, "data/coloc_atlas/results")
out_dir <- file.path(project_dir, "data/coloc_atlas/effector_expression")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

strong <- fread(file.path(res_dir, "coloc_atlas_strong_PP4gt0.8.tsv"))
effectors <- unique(strong[tags_lead_signal == TRUE, symbol])
cat("Effector genes:", paste(effectors, collapse = ", "), "\n")

s <- readRDS(file.path(project_dir, "data/GSE178146_processed/GSE178146_integrated.rds"))
s[["RNA"]] <- JoinLayers(s[["RNA"]])
DefaultAssay(s) <- "RNA"
markers <- list(Neuron = c("RBFOX3", "SYT1", "SNAP25"), Excitatory = c("SLC17A7", "CAMK2A"),
                Inhibitory = c("GAD1", "GAD2"), Astrocyte = c("AQP4", "GFAP", "SLC1A3"),
                Oligodendrocyte = c("MBP", "MOG", "PLP1"), OPC = c("PDGFRA", "CSPG4"),
                Microglia = c("CX3CR1", "P2RY12", "AIF1"), Endothelial = c("CLDN5", "FLT1"))
Idents(s) <- s$seurat_clusters
avg <- AverageExpression(s, features = unlist(markers), layer = "data")$RNA
lab <- sapply(colnames(avg), function(cl) names(markers)[which.max(sapply(markers, function(g) mean(avg[intersect(g, rownames(avg)), cl])))])
names(lab) <- gsub("^g", "", colnames(avg))
s$cell_type <- unname(lab[as.character(s$seurat_clusters)])
s$cell_type[s$cell_type %in% c("Excitatory", "Inhibitory")] <- "Neuron"
print(table(s$cell_type))

counts <- GetAssayData(s, layer = "counts")
# append effector genes removed by Phase 3 QC's min-cells filter (from 05a)
rec_file <- file.path(out_dir, "raw_recovered_gene_counts.rds")
if (file.exists(rec_file)) {
  rec <- readRDS(rec_file)
  add <- setdiff(rownames(rec), rownames(counts))
  if (length(add)) {
    counts <- rbind(counts, rec[add, colnames(counts), drop = FALSE])
    cat("Appended from raw 10x matrices:", paste(add, collapse = ", "), "\n")
  }
}
eff_in <- intersect(effectors, rownames(counts))
cat("Effectors present in snRNA matrix:", paste(eff_in, collapse = ", "),
    " | absent:", paste(setdiff(effectors, eff_in), collapse = ", "), "\n")

# ---- expression by cell type -------------------------------------------------
expr_tab <- rbindlist(lapply(sort(unique(s$cell_type)), function(ct) {
  ix <- which(s$cell_type == ct)
  data.table(cell_type = ct, gene = eff_in, n_nuclei = length(ix),
             pct_expressing = 100 * Matrix::rowMeans(counts[eff_in, ix, drop = FALSE] > 0),
             cpm = 1e6 * Matrix::rowSums(counts[eff_in, ix, drop = FALSE]) / sum(counts[, ix]))
}))
fwrite(expr_tab, file.path(out_dir, "effector_expression_by_celltype.csv"))
print(dcast(expr_tab, gene ~ cell_type, value.var = "pct_expressing")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])

# ---- severity trend per cell type ------------------------------------------
sev <- c(Control = 0, PD = 1, PDD = 2, DLB = 3)
trend <- rbindlist(lapply(sort(unique(s$cell_type)), function(ct) {
  cells <- which(s$cell_type == ct)
  donors <- sort(unique(s$sample[cells]))
  pb <- sapply(donors, function(d) Matrix::rowSums(counts[, cells[s$sample[cells] == d], drop = FALSE]))
  ncell <- table(s$sample[cells])[donors]
  keep_d <- names(ncell)[ncell >= 20]                      # donors with too few nuclei are unreliable
  if (length(keep_d) < 12) return(NULL)
  pb <- pb[, keep_d, drop = FALSE]
  meta <- unique(data.frame(sample = s$sample, condition = s$condition))[, , drop = FALSE]
  meta <- meta[match(keep_d, meta$sample), ]; rownames(meta) <- meta$sample
  meta$severity_score <- sev[as.character(meta$condition)]
  pb <- pb[rowSums(pb >= 10) >= 0.5 * ncol(pb), , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(round(pb), meta, ~severity_score)
  dds <- suppressMessages(DESeq(dds, quiet = TRUE))
  r <- as.data.frame(results(dds, name = "severity_score"))
  r$gene <- rownames(r)
  r <- as.data.table(r)[gene %in% eff_in]
  if (nrow(r) == 0) return(NULL)
  r[, `:=`(cell_type = ct, n_donors = length(keep_d), n_genes_tested = nrow(pb))]
  r
}), fill = TRUE)
setnames(trend, c("log2FoldChange", "pvalue", "padj"), c("log2FC_per_severity_step", "p_trend", "padj_trend_celltype_genomewide"))
fwrite(trend[order(p_trend)], file.path(out_dir, "effector_severity_trend_by_celltype.csv"))
print(trend[order(p_trend), .(gene, cell_type, n_donors, baseMean = round(baseMean, 1),
                              log2FC = round(log2FC_per_severity_step, 3), p = signif(p_trend, 2),
                              padj = signif(padj_trend_celltype_genomewide, 2))])

p <- ggplot(expr_tab, aes(cell_type, gene, size = pct_expressing, colour = log1p(cpm))) +
  geom_point() + scale_colour_gradient(low = "#d6e2f0", high = "#1f4e8c", name = "log1p CPM") +
  scale_size_area(max_size = 9, name = "% nuclei") +
  labs(x = NULL, y = NULL, title = "Atlas effector genes in LBD anterior cingulate cortex (GSE178146)") +
  theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(out_dir, "effector_expression_by_celltype.pdf"), p, width = 7, height = 1.5 + 0.45 * length(eff_in))
cat("Done.\n")
