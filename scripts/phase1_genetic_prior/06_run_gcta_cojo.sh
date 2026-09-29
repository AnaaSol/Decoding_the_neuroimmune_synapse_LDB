#!/usr/bin/env bash
################################################################################
# GCTA-COJO conditional analysis of the chr19 APOE/BCAM locus
#
# Resolves whether BCAM (top MAGMA gene in the underpowered pre-filtered run,
# see REPORT.md D2/Sec 5.2) carries any signal independent of the true
# APOE-locus association. Reproduces, as a script, the manual analysis
# originally run interactively (REPORT.md Sec 11.3).
#
# Requires: plink (v1.9), gcta64 (v1.94+) on PATH; conda env with both,
# e.g. `conda activate base` if installed there via bioconda.
################################################################################
set -euo pipefail

PROJECT_DIR="/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
cd "$PROJECT_DIR"

OUT_DIR="data/gcta_cojo"
REF_PANEL="data/reference/g1000_eur"
MA_FILE="$OUT_DIR/chr19_APOE_locus.ma"
RSID_FILE="$OUT_DIR/chr19_rsids.txt"
BCAM_SNP="rs28399637"   # BCAM's most-associated SNP, per data/magma_results/LBD_gene_analysis.genes.annot

mkdir -p "$OUT_DIR"

echo "== Step 1: extract 1000G-EUR genotypes for the chr19 locus SNP set (by rsID) =="
plink --bfile "$REF_PANEL" \
      --chr 19 \
      --extract "$RSID_FILE" \
      --make-bed \
      --out "$OUT_DIR/chr19_locus_ref"

echo "== Step 2: stepwise selection of independent genome-wide-significant signals =="
gcta64 --bfile "$OUT_DIR/chr19_locus_ref" \
       --cojo-file "$MA_FILE" \
       --cojo-slct \
       --cojo-p 5e-08 \
       --out "$OUT_DIR/chr19_stepwise"

JMA_FILE="$OUT_DIR/chr19_stepwise.jma.cojo"
if [ ! -s "$JMA_FILE" ]; then
    echo "ERROR: stepwise selection produced no independent signals ($JMA_FILE empty/missing)"
    exit 1
fi

# All independently-selected signals (for joint conditioning)
tail -n +2 "$JMA_FILE" | cut -f2 > "$OUT_DIR/cond_snps_all3.txt"
N_SIGNALS=$(wc -l < "$OUT_DIR/cond_snps_all3.txt")
echo "  -> $N_SIGNALS independent signal(s): $(tr '\n' ' ' < "$OUT_DIR/cond_snps_all3.txt")"

# Single lead SNP = smallest joint p-value among the selected signals
tail -n +2 "$JMA_FILE" | sort -t$'\t' -k13,13g | head -1 | cut -f2 > "$OUT_DIR/cond_snp.txt"
echo "  -> lead SNP (smallest joint p): $(cat "$OUT_DIR/cond_snp.txt")"

echo "== Step 3: condition BCAM's top SNP ($BCAM_SNP) on the single lead signal =="
gcta64 --bfile "$OUT_DIR/chr19_locus_ref" \
       --cojo-file "$MA_FILE" \
       --cojo-cond "$OUT_DIR/cond_snp.txt" \
       --out "$OUT_DIR/chr19_cond_on_lead"

echo "== Step 4: condition BCAM's top SNP on ALL $N_SIGNALS jointly-selected signals =="
gcta64 --bfile "$OUT_DIR/chr19_locus_ref" \
       --cojo-file "$MA_FILE" \
       --cojo-cond "$OUT_DIR/cond_snps_all3.txt" \
       --out "$OUT_DIR/chr19_cond_on_all3"

echo
echo "== Summary: BCAM ($BCAM_SNP) conditional p-values =="
echo "Unconditional p:"
awk -F'\t' -v snp="$BCAM_SNP" 'NR==1 || $1==snp' "$MA_FILE"
echo
echo "Conditioned on single lead SNP ($(cat "$OUT_DIR/cond_snp.txt")):"
awk -F'\t' -v snp="$BCAM_SNP" 'NR==1 || $2==snp' "$OUT_DIR/chr19_cond_on_lead.cma.cojo"
echo
echo "Conditioned jointly on all $N_SIGNALS signals:"
awk -F'\t' -v snp="$BCAM_SNP" 'NR==1 || $2==snp' "$OUT_DIR/chr19_cond_on_all3.cma.cojo"
