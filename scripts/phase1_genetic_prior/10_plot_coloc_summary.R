#!/usr/bin/env Rscript
################################################################################
# Combine coloc.abf results across all eQTL sources tested (brain ACC,
# BLUEPRINT CD4+ T cell) into one comparison table and plot, so the
# tissue-specific question ("does LBD risk converge with expression in the
# cell type actually implicated by the mechanistic hypothesis?") is visible
# directly, not just per-source.
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
results_dir <- file.path(project_dir, "data/coloc/results")

source_dirs <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
summary_files <- file.path(results_dir, source_dirs, "coloc_summary.csv")
summary_files <- summary_files[file.exists(summary_files)]

if (length(summary_files) == 0) stop("No per-source coloc_summary.csv files found under data/coloc/results/*/")

res <- do.call(rbind, lapply(summary_files, read.csv, stringsAsFactors = FALSE))

combined_file <- file.path(results_dir, "coloc_summary_all_tissues.csv")
write.csv(res, combined_file, row.names = FALSE)
cat(sprintf("Combined %d source(s) (%s) -> %s\n", length(summary_files),
            paste(unique(res$eqtl_source), collapse = ", "), combined_file))

gene_order <- res %>%
  group_by(symbol) %>%
  summarise(max_pp4 = max(PP4, na.rm = TRUE)) %>%
  arrange(desc(max_pp4)) %>%
  pull(symbol)
res$symbol <- factor(res$symbol, levels = gene_order)

res_long <- res %>%
  select(symbol, eqtl_source, PP3, PP4) %>%
  pivot_longer(cols = c(PP3, PP4), names_to = "hypothesis", values_to = "posterior_prob")

p <- ggplot(res_long, aes(x = symbol, y = posterior_prob, fill = hypothesis)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
  facet_wrap(~eqtl_source, ncol = 1) +
  scale_fill_manual(
    values = c(PP3 = "#e67e22", PP4 = "#3498db"),
    labels = c(PP3 = "PP3: distinct causal variants", PP4 = "PP4: shared causal variant")
  ) +
  labs(
    title = "GWAS-vs-eQTL colocalization across tissues (coloc.abf)",
    subtitle = "Dashed line = conventional 0.8 threshold for a confident call",
    x = NULL, y = "Posterior probability", fill = NULL
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_out <- file.path(results_dir, "coloc_PP3_PP4_summary_all_tissues.pdf")
ggsave(plot_out, p, width = 8, height = 4 + 3 * length(unique(res$eqtl_source)))
cat(sprintf("Saved: %s\n", plot_out))

if (length(unique(res$eqtl_source)) > 1) {
  cat("\n-- Cross-tissue comparison (PP4 per gene) --\n")
  wide <- res %>% select(symbol, eqtl_source, PP4) %>% pivot_wider(names_from = eqtl_source, values_from = PP4)
  print(wide)
}
