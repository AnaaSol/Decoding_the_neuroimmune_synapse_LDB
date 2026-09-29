#!/bin/bash
################################################################################
# MAGMA Gene-Based Analysis Script
# Author: Ana
# Date: 2025-11-25
# Phase: 1 - Genetic Prior Construction
#
# Purpose: Run MAGMA to aggregate SNP-level p-values into gene-level Z-scores
#
# Requirements:
#   - MAGMA software (https://ctg.cptc.nl/software/magma)
#   - 1000 Genomes reference panel (European ancestry)
#   - Gene annotation file (NCBI37.3 or GRCh38)
#
# Output: Gene-level association scores ranked by Z-score
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Project paths
PROJECT_DIR="/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
DATA_DIR="${PROJECT_DIR}/data"
RESULTS_DIR="${DATA_DIR}/magma_results"
REFERENCE_DIR="${DATA_DIR}/reference"
TOOLS_DIR="${PROJECT_DIR}/tools"

# MAGMA executable
MAGMA="${TOOLS_DIR}/magma"

# Input files
MAGMA_INPUT="${DATA_DIR}/LBD_GWAS_magma_input.txt"

# MAGMA reference files (to be downloaded)
GENE_ANNOT="${REFERENCE_DIR}/NCBI38.gene.loc"  # Gene annotation file (GRCh38/hg38)
# NOTE: Most g1000_eur panels from CTG Lab are in hg19/GRCh37
# MAGMA will match by rsID which works for ~95% of SNPs
# For maximum accuracy, use a GRCh38 reference panel if available
REFERENCE_PANEL="${REFERENCE_DIR}/g1000_eur"    # 1000G EUR reference (prefix)

# Output files
OUTPUT_PREFIX="${RESULTS_DIR}/LBD_gene_analysis"

################################################################################
# Function: Print colored message
################################################################################
print_message() {
    local color=$1
    local message=$2
    echo -e "${color}${message}${NC}"
}

################################################################################
# Function: Check if file exists
################################################################################
check_file() {
    local file=$1
    local description=$2

    if [ ! -f "$file" ]; then
        print_message "$RED" "✗ Error: $description not found at: $file"
        return 1
    else
        print_message "$GREEN" "✓ Found: $description"
        return 0
    fi
}

################################################################################
# Function: Check MAGMA installation
################################################################################
check_magma() {
    if [ ! -f "$MAGMA" ]; then
        print_message "$RED" "✗ MAGMA executable not found at: $MAGMA"
        print_message "$YELLOW" "\nTo install MAGMA:"
        print_message "$YELLOW" "1. Download from: https://cncr.nl/research/magma/"
        print_message "$YELLOW" "2. Extract to: ${TOOLS_DIR}/"
        return 1
    else
        print_message "$GREEN" "✓ MAGMA installation found at: $MAGMA"
        "$MAGMA" --version
        return 0
    fi
}

################################################################################
# Function: Download reference files
################################################################################
download_references() {
    print_message "$BLUE" "\n=========================================="
    print_message "$BLUE" "Reference Files Setup"
    print_message "$BLUE" "=========================================="

    mkdir -p "$REFERENCE_DIR"

    print_message "$YELLOW" "\nMAGMA requires reference files:"
    print_message "$YELLOW" "1. Gene location file (NCBI build 38)"
    print_message "$YELLOW" "2. 1000 Genomes EUR reference panel"

    print_message "$YELLOW" "\nDownload instructions:"
    print_message "$YELLOW" "Visit: https://ctg.cptc.nl/software/magma"
    print_message "$YELLOW" "Download:"
    print_message "$YELLOW" "  - NCBI38.zip (gene annotation)"
    print_message "$YELLOW" "  - g1000_eur.zip (reference panel for European ancestry)"

    print_message "$YELLOW" "\nExtract to: ${REFERENCE_DIR}/"

    read -p "Have you downloaded the reference files? (y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_message "$RED" "Please download reference files first."
        exit 1
    fi
}

################################################################################
# Main Analysis Pipeline
################################################################################
main() {
    print_message "$BLUE" "\n=========================================="
    print_message "$BLUE" "MAGMA Gene-Based Analysis"
    print_message "$BLUE" "=========================================="

    # Check prerequisites
    print_message "$BLUE" "\n1. Checking prerequisites..."

    if ! check_magma; then
        exit 1
    fi

    if ! check_file "$MAGMA_INPUT" "MAGMA input file"; then
        print_message "$YELLOW" "Run: python3 02_prepare_magma_input.py"
        exit 1
    fi

    # Check reference files
    if [ ! -f "${GENE_ANNOT}" ] || [ ! -f "${REFERENCE_PANEL}.bim" ]; then
        download_references
    fi

    check_file "$GENE_ANNOT" "Gene annotation file" || exit 1
    check_file "${REFERENCE_PANEL}.bim" "Reference panel" || exit 1

    # Create results directory
    mkdir -p "$RESULTS_DIR"

    # Step 1: Annotate SNPs to genes
    print_message "$BLUE" "\n2. Annotating SNPs to genes..."
    print_message "$YELLOW" "This maps each SNP to genes based on genomic location"

    # MAGMA will auto-detect columns from header (SNP, CHR, BP)
    "$MAGMA" --annotate \
        --snp-loc "${MAGMA_INPUT}" \
        --gene-loc "${GENE_ANNOT}" \
        --out "${OUTPUT_PREFIX}"

    print_message "$GREEN" "✓ SNP annotation completed"

    # Step 2: Gene analysis
    print_message "$BLUE" "\n3. Running gene-based analysis..."
    print_message "$YELLOW" "Aggregating SNP p-values into gene-level statistics"

    "$MAGMA" --bfile "${REFERENCE_PANEL}" \
        --pval "${MAGMA_INPUT}" ncol=N \
        --gene-annot "${OUTPUT_PREFIX}.genes.annot" \
        --out "${OUTPUT_PREFIX}"

    print_message "$GREEN" "✓ Gene analysis completed"

    # Step 3: Process results
    print_message "$BLUE" "\n4. Processing results..."

    # Sort genes by Z-score
    RESULTS_FILE="${OUTPUT_PREFIX}.genes.out"
    SORTED_FILE="${OUTPUT_PREFIX}.genes.sorted.txt"

    if [ -f "$RESULTS_FILE" ]; then
        # Extract header and sort by Z-score (column 7) in descending order
        head -n 1 "$RESULTS_FILE" > "$SORTED_FILE"
        tail -n +2 "$RESULTS_FILE" | sort -k7 -rn >> "$SORTED_FILE"

        print_message "$GREEN" "✓ Results sorted by Z-score"

        # Display top 20 genes
        print_message "$BLUE" "\n=========================================="
        print_message "$BLUE" "Top 20 LBD Risk Genes (by Z-score)"
        print_message "$BLUE" "=========================================="
        head -n 21 "$SORTED_FILE" | column -t

        print_message "$GREEN" "\n✓ Full results saved to: ${SORTED_FILE}"

        # Generate summary statistics
        TOTAL_GENES=$(tail -n +2 "$SORTED_FILE" | wc -l)
        SIG_GENES=$(tail -n +2 "$SORTED_FILE" | awk '$8 < 0.05' | wc -l)
        BONF_GENES=$(tail -n +2 "$SORTED_FILE" | awk -v n="$TOTAL_GENES" '$8 < (0.05/n)' | wc -l)

        print_message "$BLUE" "\n=========================================="
        print_message "$BLUE" "Summary Statistics"
        print_message "$BLUE" "=========================================="
        print_message "$YELLOW" "Total genes analyzed: ${TOTAL_GENES}"
        print_message "$YELLOW" "Nominally significant (p < 0.05): ${SIG_GENES}"
        print_message "$YELLOW" "Bonferroni significant: ${BONF_GENES}"

    else
        print_message "$RED" "✗ Results file not found: $RESULTS_FILE"
        exit 1
    fi

    print_message "$GREEN" "\n=========================================="
    print_message "$GREEN" "Analysis Complete!"
    print_message "$GREEN" "=========================================="
    print_message "$YELLOW" "\nOutput files:"
    print_message "$YELLOW" "  - ${OUTPUT_PREFIX}.genes.annot (SNP-gene mapping)"
    print_message "$YELLOW" "  - ${OUTPUT_PREFIX}.genes.out (Full results)"
    print_message "$YELLOW" "  - ${SORTED_FILE} (Sorted by Z-score)"
}

# Run main pipeline
main "$@"
