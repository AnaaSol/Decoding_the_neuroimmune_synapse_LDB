#!/usr/bin/env python3
"""
MAGMA Results Parser and Visualization
Author: Ana
Date: 2025-11-25
Phase: 1 - Genetic Prior Construction

Purpose: Parse MAGMA gene-level results and create ranked gene table
  - Extract gene symbols, Z-scores, and p-values
  - Create ranked table for downstream analysis
  - Generate summary statistics and visualizations

Input: ../../data/magma_results/LBD_gene_analysis.genes.out
Output: ../../data/LBD_risk_genes_ranked.tsv
"""

import pandas as pd
import numpy as np
import sys
import os

def map_entrez_to_symbol(df, gene_loc_file):
    """
    Map Entrez Gene IDs to Gene Symbols using NCBI gene annotation file.

    This greatly improves usability of MAGMA results by showing gene names
    instead of numeric IDs.

    Parameters:
    -----------
    df : pd.DataFrame
        MAGMA results with 'GENE' column containing Entrez IDs
    gene_loc_file : str
        Path to NCBI gene location file (e.g., NCBI38.gene.loc)

    Returns:
    --------
    pd.DataFrame
        DataFrame with added 'SYMBOL' column
    """

    try:
        print(f"\nMapping Entrez IDs to Gene Symbols...")
        print(f"  Using annotation file: {gene_loc_file}")

        # NCBI gene.loc format: EntrezID Chr Start Stop Strand Symbol
        # No header, whitespace-separated
        annot = pd.read_csv(
            gene_loc_file,
            sep=r'\s+',
            header=None,
            names=['GENE', 'CHR_ANNOT', 'START_ANNOT', 'STOP_ANNOT', 'STRAND', 'SYMBOL']
        )

        # Keep only GENE and SYMBOL columns for mapping
        gene_map = annot[['GENE', 'SYMBOL']].drop_duplicates(subset='GENE')

        print(f"  Loaded {len(gene_map):,} gene annotations")

        # Merge with MAGMA results
        df_mapped = df.merge(gene_map, on='GENE', how='left')

        # Check how many genes were successfully mapped
        mapped_count = df_mapped['SYMBOL'].notna().sum()
        unmapped_count = df_mapped['SYMBOL'].isna().sum()

        print(f"  ✓ Mapped {mapped_count:,} genes ({mapped_count/len(df)*100:.1f}%)")

        if unmapped_count > 0:
            print(f"  ⚠ {unmapped_count} genes without symbol mapping")
            # For unmapped genes, use Entrez ID as fallback
            df_mapped.loc[df_mapped['SYMBOL'].isna(), 'SYMBOL'] = \
                'ENTREZ_' + df_mapped.loc[df_mapped['SYMBOL'].isna(), 'GENE'].astype(str)

        return df_mapped

    except FileNotFoundError:
        print(f"\n⚠ Warning: Gene annotation file not found: {gene_loc_file}")
        print("  Continuing with Entrez IDs only (no symbols)")
        df['SYMBOL'] = 'ENTREZ_' + df['GENE'].astype(str)
        return df

    except Exception as e:
        print(f"\n⚠ Warning: Could not map gene symbols ({e})")
        print("  Continuing with Entrez IDs only")
        df['SYMBOL'] = 'ENTREZ_' + df['GENE'].astype(str)
        return df

def parse_magma_results(magma_output, output_file, gene_loc_file=None, top_n=100):
    """
    Parse MAGMA gene-based analysis results.

    MAGMA output columns:
    - GENE: Gene ID (Entrez)
    - CHR: Chromosome
    - START: Start position
    - STOP: Stop position
    - NSNPS: Number of SNPs in gene
    - NPARAM: Number of parameters
    - N: Sample size
    - ZSTAT: Z-statistic (main metric of interest)
    - P: P-value

    Parameters:
    -----------
    magma_output : str
        Path to MAGMA .genes.out file
    output_file : str
        Path to output ranked gene table
    gene_loc_file : str, optional
        Path to NCBI gene location file for mapping Entrez IDs to symbols
    top_n : int
        Number of top genes to highlight
    """

    print("="*70)
    print("MAGMA Results Parser")
    print("="*70)

    # Check if input file exists
    if not os.path.exists(magma_output):
        print(f"\n✗ Error: MAGMA output file not found: {magma_output}")
        print("Please run 03_run_magma_analysis.sh first")
        sys.exit(1)

    # Read MAGMA results
    print(f"\nReading MAGMA results from: {magma_output}")
    df = pd.read_csv(magma_output, sep=r'\s+')
    print(f"Loaded {len(df):,} genes")

    # Map Entrez IDs to Gene Symbols (CRITICAL for usability)
    if gene_loc_file:
        df = map_entrez_to_symbol(df, gene_loc_file)
    else:
        print("\n⚠ Warning: No gene annotation file provided")
        print("  Results will show Entrez IDs only")
        df['SYMBOL'] = 'ENTREZ_' + df['GENE'].astype(str)

    # Sort by Z-score (descending)
    df_sorted = df.sort_values('ZSTAT', ascending=False).reset_index(drop=True)

    # Add rank column
    df_sorted['RANK'] = range(1, len(df_sorted) + 1)

    # Calculate -log10(p-value) for easier interpretation
    df_sorted['NEG_LOG10_P'] = -np.log10(df_sorted['P'])

    # Add significance flags
    bonferroni_threshold = 0.05 / len(df_sorted)
    df_sorted['BONFERRONI_SIG'] = df_sorted['P'] < bonferroni_threshold
    df_sorted['NOMINAL_SIG'] = df_sorted['P'] < 0.05

    # Reorder columns for output (SYMBOL first for readability)
    output_columns = [
        'RANK', 'SYMBOL', 'GENE', 'CHR', 'START', 'STOP',
        'NSNPS', 'ZSTAT', 'P', 'NEG_LOG10_P',
        'BONFERRONI_SIG', 'NOMINAL_SIG'
    ]
    df_output = df_sorted[output_columns].copy()

    # Save ranked gene table
    print(f"\nSaving ranked gene table to: {output_file}")
    df_output.to_csv(output_file, sep='\t', index=False)

    # Print summary statistics
    print("\n" + "="*70)
    print("Summary Statistics")
    print("="*70)

    total_genes = len(df_sorted)
    bonf_sig = df_sorted['BONFERRONI_SIG'].sum()
    nom_sig = df_sorted['NOMINAL_SIG'].sum()

    print(f"\nTotal genes analyzed: {total_genes:,}")
    print(f"Bonferroni significant (p < {bonferroni_threshold:.2e}): {bonf_sig}")
    print(f"Nominally significant (p < 0.05): {nom_sig}")

    print(f"\nZ-score statistics:")
    print(f"  Max: {df_sorted['ZSTAT'].max():.2f}")
    print(f"  Min: {df_sorted['ZSTAT'].min():.2f}")
    print(f"  Mean: {df_sorted['ZSTAT'].mean():.2f}")
    print(f"  Median: {df_sorted['ZSTAT'].median():.2f}")

    # Display top genes
    print(f"\n{'='*70}")
    print(f"Top {top_n} LBD Risk Genes (Ranked by Z-score)")
    print("="*70)
    print(df_output.head(top_n).to_string(index=False))

    # Display significant genes if any
    if bonf_sig > 0:
        print(f"\n{'='*70}")
        print(f"Bonferroni Significant Genes (n={bonf_sig})")
        print("="*70)
        sig_genes = df_output[df_output['BONFERRONI_SIG'] == True]
        print(sig_genes.to_string(index=False))

    print("\n" + "="*70)
    print("✓ MAGMA results parsed successfully!")
    print("="*70)

    return df_output

def create_gene_set_for_downstream(ranked_genes, output_dir, percentiles=[1, 5, 10, 25]):
    """
    Create gene sets at different percentile thresholds for downstream analysis.

    Parameters:
    -----------
    ranked_genes : pd.DataFrame
        Ranked gene table from MAGMA
    output_dir : str
        Directory to save gene sets
    percentiles : list
        Percentile cutoffs for gene sets
    """

    print(f"\n{'='*70}")
    print("Creating Gene Sets for Downstream Analysis")
    print("="*70)

    os.makedirs(output_dir, exist_ok=True)

    for pct in percentiles:
        n_genes = int(np.ceil(len(ranked_genes) * pct / 100))
        top_genes = ranked_genes.head(n_genes)

        output_file = os.path.join(output_dir, f"LBD_risk_genes_top{pct}pct.txt")

        # Save gene IDs only (for pathway analysis, etc.)
        top_genes['GENE'].to_csv(output_file, index=False, header=False)

        print(f"  Top {pct}% ({n_genes} genes) → {output_file}")

    print("\n✓ Gene sets created successfully!")

if __name__ == "__main__":
    # Define file paths (relative to project root)
    magma_results = "../../data/magma_results/LBD_gene_analysis.genes.out"
    output_file = "../../data/LBD_risk_genes_ranked.tsv"
    gene_sets_dir = "../../data/gene_sets"
    gene_annotation = "../../data/reference/NCBI38.gene.loc"

    try:
        # Parse MAGMA results with gene symbol mapping
        ranked_genes = parse_magma_results(
            magma_results,
            output_file,
            gene_loc_file=gene_annotation,
            top_n=50
        )

        # Create gene sets for downstream analysis
        create_gene_set_for_downstream(ranked_genes, gene_sets_dir)

        print("\n" + "="*70)
        print("PHASE 1 COMPLETE: Genetic Prior Constructed")
        print("="*70)
        print("\nNext steps:")
        print("  → Phase 2: Immune cell profiling (GSE161192)")
        print("  → Phase 3: Neuronal vulnerability analysis (GSE178146)")
        print("  → Phase 4: Risk-weighted interactomics")
        print("="*70)

    except Exception as e:
        print(f"\n✗ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
