#!/usr/bin/env python3
"""
MAGMA Input Preparation Script
Author: Ana
Date: 2025-11-25
Phase: 1 - Genetic Prior Construction

Purpose: Convert filtered GWAS data to MAGMA-compatible format
  - Extract SNP ID, P-value, and sample size
  - Format according to MAGMA requirements

Input: ../../data/LBD_GWAS_filtered.tsv.gz
Output: ../../data/LBD_GWAS_magma_input.txt
"""

import pandas as pd
import sys

def prepare_magma_input(gwas_file, output_file, metadata_file=None):
    """
    Prepare GWAS data for MAGMA gene-based analysis.

    MAGMA requires a file with columns:
    - SNP: SNP identifier (rsID)
    - CHR: Chromosome number (required for --snp-loc)
    - BP: Base pair position (required for --snp-loc)
    - P: P-value
    - N: Sample size (optional but recommended)

    Parameters:
    -----------
    gwas_file : str
        Path to filtered GWAS file
    output_file : str
        Path to output MAGMA input file
    metadata_file : str
        Path to metadata YAML file (to extract sample size)
    """

    print("="*70)
    print("MAGMA Input Preparation")
    print("="*70)

    # Read filtered GWAS data
    print(f"\nReading filtered GWAS data from: {gwas_file}")
    df = pd.read_csv(gwas_file, sep='\t', compression='gzip')
    print(f"Loaded {len(df):,} genome-wide significant variants")

    # Extract sample size from metadata if available
    sample_size = None
    if metadata_file:
        try:
            import yaml
            with open(metadata_file, 'r') as f:
                metadata = yaml.safe_load(f)
                # Try to find sample size in metadata
                # Adjust this based on actual metadata structure
                sample_size = metadata.get('sample_size', None)
                if sample_size:
                    print(f"Sample size from metadata: {sample_size:,}")
        except Exception as e:
            print(f"Could not read metadata: {e}")

    # If sample size not in metadata, try to infer or use a default
    if sample_size is None:
        # Check if there's a sample size column in the GWAS file
        if 'n' in df.columns or 'N' in df.columns or 'sample_size' in df.columns:
            n_col = [c for c in df.columns if c.lower() in ['n', 'sample_size']][0]
            sample_size = int(df[n_col].iloc[0])
            print(f"Sample size from GWAS file: {sample_size:,}")
        else:
            # For Chia et al. 2021 LBD GWAS (GCST90001390)
            # Total N = 6,618: 2,591 cases + 4,027 controls (European ancestry)
            # Reference: Chia R, et al. Nat Genet 2021;53:294-303
            sample_size = 6618  # Total N (cases + controls)
            print(f"Using known sample size for Chia et al. 2021: {sample_size:,}")
            print("  (2,591 cases + 4,027 controls, European ancestry)")

    # Prepare MAGMA input file
    print("\nPreparing MAGMA input format...")

    # Create output dataframe with required columns
    # CRITICAL: CHR and BP are required for MAGMA's --snp-loc annotation step
    # Convert CHR and BP to integers to avoid decimal notation (1.0 -> 1)
    magma_df = pd.DataFrame({
        'SNP': df['hm_rsid'],
        'CHR': df['hm_chrom'].astype(int),
        'BP': df['hm_pos'].astype(int),
        'P': df['p_value'],
        'N': sample_size
    })

    # Remove any variants without rsID
    initial_count = len(magma_df)
    magma_df = magma_df[magma_df['SNP'].notna()].copy()
    removed_no_rsid = initial_count - len(magma_df)

    if removed_no_rsid > 0:
        print(f"  - Removed {removed_no_rsid} variants without rsID")

    print(f"  - Final variant count for MAGMA: {len(magma_df):,}")

    # Sort by chromosome and position for better MAGMA performance
    magma_df = magma_df.sort_values(['CHR', 'BP']).reset_index(drop=True)

    # Save to file
    print(f"\nSaving MAGMA input to: {output_file}")
    magma_df.to_csv(output_file, sep='\t', index=False)

    # Display summary statistics
    print("\n" + "="*70)
    print("MAGMA Input Summary:")
    print("="*70)
    print(f"Total SNPs: {len(magma_df):,}")
    print(f"Sample size: {sample_size:,}")
    print(f"\nP-value range:")
    print(f"  Min: {magma_df['P'].min():.2e}")
    print(f"  Max: {magma_df['P'].max():.2e}")
    print(f"  Median: {magma_df['P'].median():.2e}")
    print("\nFirst 5 rows:")
    print(magma_df.head())
    print("="*70)

    return magma_df

if __name__ == "__main__":
    # Define file paths (relative to project root)
    gwas_input = "../../data/LBD_GWAS_filtered.tsv.gz"
    magma_output = "../../data/LBD_GWAS_magma_input.txt"
    metadata = "../../data/LBD_GWAS_metadata.yaml"

    try:
        # Run preparation
        magma_data = prepare_magma_input(gwas_input, magma_output, metadata)
        print("\n✓ MAGMA input preparation completed successfully!")

    except FileNotFoundError as e:
        print(f"\n✗ Error: Input file not found - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Error during preparation: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
