#!/usr/bin/env bash
################################################################################
# MAGMA competitive gene-set analysis: is LBD genetic risk enriched in
# curated immune pathways, as a pathway-level polygenic test? This has much
# more statistical power than the single-gene Bonferroni test (Phase 1) to
# detect convergence between genetic risk and immune function, since it
# aggregates signal across many genes rather than requiring any single gene
# to reach genome-wide significance.
#
# Input: LBD_gene_analysis.genes.raw (from the full genome-wide gene-based
# run, Phase 1) + a pre-specified immune-pathway gene-set file (KEGG, fetched
# by 11_fetch_kegg_genesets.py -- NOT selected post-hoc).
################################################################################
set -euo pipefail

PROJECT_DIR="/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
cd "$PROJECT_DIR"

MAGMA="tools/magma"
GENE_RESULTS="data/magma_results/LBD_gene_analysis.genes.raw"
SET_FILE="data/geneset_analysis/immune_pathways.set_annot.txt"
OUT_PREFIX="data/geneset_analysis/LBD_immune_geneset"

if [ ! -f "$GENE_RESULTS" ]; then
    echo "ERROR: $GENE_RESULTS not found. Run Phase 1 (parse_magma rule) first." >&2
    exit 1
fi

"$MAGMA" --gene-results "$GENE_RESULTS" \
         --set-annot "$SET_FILE" \
         --out "$OUT_PREFIX"

echo
echo "== Result: ${OUT_PREFIX}.gsa.out =="
cat "${OUT_PREFIX}.gsa.out"
