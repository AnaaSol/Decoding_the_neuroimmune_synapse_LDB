#!/usr/bin/env Rscript
################################################################################
# Permutation-based empirical null for the Phase 4 Prioritization Score.
#
# The manuscript already states, from inspection, that "no receptor carries
# significant genetic risk... this ranking is driven almost entirely by
# expression level, not genetic evidence." This script tests that claim
# formally instead of asserting it: permute the assignment of receptor GWAS
# Z-scores across the candidate receptors (breaking the true gene<->Z
# correspondence) many times, recompute Prioritization_Score and rank each
# time, and report how often OSM-LIFR's #1 rank (and each pair's exact score)
# would occur under random Z-score labeling alone.
#
# If permuting Z barely changes the ranking (expected here, since real Z's
# are all small/non-significant), that CONFIRMS the score is expression-
# dominated -- an empirical test, not just a stated caveat.
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

set.seed(42)
project_dir <- "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
net <- read.csv(file.path(project_dir, "results_final/weighted_neuroimmune_network.csv"), stringsAsFactors = FALSE)

cat("================================================================================\n")
cat("Phase 4b: Permutation null for the Prioritization Score\n")
cat("================================================================================\n\n")

n_pairs <- nrow(net)
observed_rank <- rank(-net$Prioritization_Score, ties.method = "first")
top_pair_idx <- which.min(observed_rank)
cat(sprintf("Observed top pair: %s-%s (score=%.4f, Z=%.3f, expr=%.4f)\n\n",
            net$ligand[top_pair_idx], net$receptor[top_pair_idx],
            net$Prioritization_Score[top_pair_idx], net$Receptor_Risk_Z[top_pair_idx],
            net$Max_Neuronal_Expression[top_pair_idx]))

n_perm <- 10000
top_pair_stays_top <- 0
rank_correlation <- numeric(n_perm)
score_at_top_pair <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  z_perm <- sample(net$Receptor_Risk_Z)
  score_perm <- net$Max_Neuronal_Expression * (1 + pmax(0, z_perm))
  rank_perm <- rank(-score_perm, ties.method = "first")
  if (rank_perm[top_pair_idx] == 1) top_pair_stays_top <- top_pair_stays_top + 1
  rank_correlation[i] <- cor(rank_perm, rank(-net$Max_Neuronal_Expression, ties.method = "first"),
                              method = "spearman")
  score_at_top_pair[i] <- score_perm[top_pair_idx]
}

p_empirical <- top_pair_stays_top / n_perm
cat(sprintf("De %d permutaciones de las etiquetas Z entre receptores:\n", n_perm))
cat(sprintf("  %s-%s sigue siendo el par #1 en %d permutaciones (%.1f%%)\n",
            net$ligand[top_pair_idx], net$receptor[top_pair_idx], top_pair_stays_top, 100 * p_empirical))
cat(sprintf("  Correlación de Spearman promedio (ranking con Z permutado vs. ranking por expresión sola): %.3f\n",
            mean(rank_correlation)))
cat(sprintf("  (1.0 significaría que el score con Z permutado es idéntico al ranking por expresión pura)\n\n"))

cat(sprintf("Interpretación: si %s-%s sigue #1 en %.0f%% de las permutaciones al azar de Z,\n",
            net$ligand[top_pair_idx], net$receptor[top_pair_idx], 100 * p_empirical))
cat("el ranking está dominado por la expresión neuronal, no por el riesgo genético --\n")
cat("confirma empíricamente, en vez de solo afirmar, la limitación ya declarada en el manuscrito.\n\n")

out <- data.frame(
  top_pair = paste0(net$ligand[top_pair_idx], "-", net$receptor[top_pair_idx]),
  n_permutations = n_perm,
  pct_permutations_top_pair_stays_1st = 100 * p_empirical,
  mean_spearman_vs_expression_only_rank = mean(rank_correlation),
  observed_score = net$Prioritization_Score[top_pair_idx],
  observed_Z = net$Receptor_Risk_Z[top_pair_idx]
)
out_file <- file.path(project_dir, "results_final/permutation_null_prioritization_score.csv")
write.csv(out, out_file, row.names = FALSE)
cat(sprintf("Guardado: %s\n\n", out_file))

cat("--- R Session Info ---\n")
print(sessionInfo())
