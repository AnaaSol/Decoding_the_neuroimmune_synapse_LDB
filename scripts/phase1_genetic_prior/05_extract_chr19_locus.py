#!/usr/bin/env python3
"""
Extract chr19 APOE-locus summary statistics in GCTA .ma format, and the
matching rsID list for PLINK --extract.

This reproduces (as a script, not an ad hoc manual step) the extraction
originally done interactively for the BCAM/APOE GCTA-COJO conditional
analysis (see REPORT.md Sec 11.3). SNP matching is done by hm_rsid on the
GRCh38-harmonised GWAS file. Do NOT extract by physical position against
data/reference/g1000_eur directly -- that reference panel is GRCh37/hg19
(defect D13/D23): matching must happen via rsID, which PLINK does downstream
using this script's rsID list.
"""

import argparse
import gzip
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="LBD_GWAS_harmonised.tsv.gz")
    parser.add_argument("--chrom", type=int, required=True)
    parser.add_argument("--start", type=int, required=True, help="GRCh38 region start (bp)")
    parser.add_argument("--end", type=int, required=True, help="GRCh38 region end (bp)")
    parser.add_argument("--n-total", type=int, required=True, help="Total GWAS N (constant per-SNP)")
    parser.add_argument("--out-ma", required=True, help="Output .ma file (GCTA COJO format)")
    parser.add_argument("--out-rsids", required=True, help="Output rsID list for PLINK --extract")
    args = parser.parse_args()

    header_idx = {}
    n_written = 0
    n_seen_in_region = 0
    seen_rsids = set()

    with gzip.open(args.input, "rt") as fin, \
         open(args.out_ma, "w") as fma, \
         open(args.out_rsids, "w") as frs:

        header = fin.readline().rstrip("\n").split("\t")
        header_idx = {name: i for i, name in enumerate(header)}

        required = ["hm_rsid", "hm_chrom", "hm_pos", "hm_effect_allele",
                    "hm_other_allele", "hm_effect_allele_frequency",
                    "hm_beta", "standard_error", "p_value"]
        missing = [c for c in required if c not in header_idx]
        if missing:
            sys.exit(f"ERROR: missing expected columns in {args.input}: {missing}")

        fma.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n")

        for line in fin:
            f = line.rstrip("\n").split("\t")
            try:
                chrom = int(f[header_idx["hm_chrom"]])
                pos = int(f[header_idx["hm_pos"]])
            except (ValueError, IndexError):
                continue
            if chrom != args.chrom or not (args.start <= pos <= args.end):
                continue
            n_seen_in_region += 1

            rsid = f[header_idx["hm_rsid"]]
            a1 = f[header_idx["hm_effect_allele"]]
            a2 = f[header_idx["hm_other_allele"]]
            freq = f[header_idx["hm_effect_allele_frequency"]]
            b = f[header_idx["hm_beta"]]
            se = f[header_idx["standard_error"]]
            p = f[header_idx["p_value"]]

            if not rsid or rsid == "NA" or rsid in seen_rsids:
                continue
            if "NA" in (a1, a2, freq, b, se, p) or "" in (a1, a2, freq, b, se, p):
                continue

            seen_rsids.add(rsid)
            fma.write(f"{rsid}\t{a1}\t{a2}\t{freq}\t{b}\t{se}\t{p}\t{args.n_total}\n")
            frs.write(f"{rsid}\n")
            n_written += 1

    print(f"Region chr{args.chrom}:{args.start}-{args.end} — "
          f"{n_seen_in_region} SNPs in region, {n_written} with complete "
          f"rsID/effect data written to {args.out_ma} and {args.out_rsids}")


if __name__ == "__main__":
    main()
