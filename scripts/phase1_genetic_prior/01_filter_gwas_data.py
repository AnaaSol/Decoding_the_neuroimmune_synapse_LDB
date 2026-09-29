#!/usr/bin/env python3
"""
GWAS Data Filtering Script
Author: Ana
Date: 2025-11-25

Purpose: Filter LBD GWAS harmonised data by:
  1. p_value <= 1e-5 
  2. hm_beta is not null

Input: data/LBD_GWAS_harmonised.tsv.gz
Output: data/LBD_GWAS_filtered.tsv.gz
"""

import argparse
import pandas as pd
import gzip
import sys

def filter_gwas_data(input_file, output_file, p_threshold=1.0):
    """
    Filter GWAS data based on p-value and hm_beta criteria.

    Parameters:
    -----------
    input_file : str
        Path to input GWAS file (gzipped)
    output_file : str
        Path to output filtered file (gzipped)
    p_threshold : float
        P-value threshold (default 1.0 = keep all SNPs, recommended for MAGMA).
        Use 5e-8 for genome-wide significant SNPs only, or 1e-5 for suggestive.
    """

    print(f"Reading GWAS data from: {input_file}")

    # Read the gzipped TSV file
    df = pd.read_csv(input_file, sep='\t', compression='gzip')

    # Get initial count
    initial_count = len(df)
    print(f"Initial variant count: {initial_count:,}")

    print("\nApplying filters...")
    print(f"  P-value threshold: {p_threshold} ({'no filter — all SNPs kept' if p_threshold >= 1.0 else 'filtered'})")
    df_filtered = df[df['p_value'] <= p_threshold].copy()
    p_value_filtered_count = initial_count - len(df_filtered)
    print(f"  - Removed {p_value_filtered_count:,} variants with p_value > {p_threshold}")

    # Filter 2: hm_beta is not null (not NA/NaN)
    df_filtered = df_filtered[df_filtered['hm_beta'].notna()].copy()
    beta_filtered_count = len(df[df['p_value'] <= p_threshold]) - len(df_filtered)
    print(f"  - Removed {beta_filtered_count:,} variants with null hm_beta")

    # Final count
    final_count = len(df_filtered)
    total_removed = initial_count - final_count
    retention_rate = (final_count / initial_count) * 100

    print(f"\nFiltering summary:")
    print(f"  - Total variants removed: {total_removed:,}")
    print(f"  - Final variant count: {final_count:,}")
    print(f"  - Retention rate: {retention_rate:.2f}%")

    # Save filtered data
    print(f"\nSaving filtered data to: {output_file}")
    df_filtered.to_csv(output_file, sep='\t', compression='gzip', index=False)
    print("Done!")

    return df_filtered

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter GWAS harmonised data for downstream analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--input",
        default="/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB/data/LBD_GWAS_harmonised.tsv.gz",
        help="Path to input harmonised GWAS file (gzipped TSV)"
    )
    parser.add_argument(
        "--output",
        default="/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB/data/LBD_GWAS_filtered.tsv.gz",
        help="Path to output file (gzipped TSV)"
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=1.0,
        help=(
            "P-value upper bound. Default 1.0 keeps all SNPs — required for a "
            "valid MAGMA genome-wide gene-based test. Use 5e-8 for genome-wide "
            "significant SNPs only, or 1e-5 for the (incorrect) suggestive filter "
            "used in the original analysis."
        )
    )
    args = parser.parse_args()

    try:
        filtered_data = filter_gwas_data(args.input, args.output,
                                         p_threshold=args.p_threshold)

        print("\n" + "="*60)
        print("P-value statistics of filtered data:")
        print(filtered_data['p_value'].describe())
        print("\nBeta coefficient statistics:")
        print(filtered_data['hm_beta'].describe())
        print("="*60)

    except FileNotFoundError:
        print("Error: Input file not found. Please check the path.")
        sys.exit(1)
    except Exception as e:
        print(f"Error during filtering: {str(e)}")
        sys.exit(1)
