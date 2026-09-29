#!/usr/bin/env python3
"""
Single-pass extraction of GWAS summary statistics for each coloc target gene's
cis-window (config: phase1c_coloc.target_genes, cis_window_bp), from the full
harmonised GWAS. One pass over the 374MB file regardless of gene count (some
windows overlap, e.g. the chr19 APOE-locus genes and chr4 SNCA/MMRN1).

Output: one file per gene, data/coloc/gwas/{symbol}_gwas.tsv, columns:
  rsid  effect_allele  other_allele  eaf  beta  se  pvalue
"""

import argparse
import gzip
import os
import sys

import yaml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", required=True, help="LBD_GWAS_harmonised.tsv.gz")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    window = cfg["phase1c_coloc"]["cis_window_bp"]
    genes = cfg["phase1c_coloc"]["target_genes"]

    regions = []
    for g in genes:
        regions.append({
            "symbol": g["symbol"],
            "chrom": int(g["chrom"]),
            "start": int(g["start"]) - window,
            "end": int(g["end"]) + window,
        })

    os.makedirs(args.out_dir, exist_ok=True)
    out_files = {}
    counts = {r["symbol"]: 0 for r in regions}
    for r in regions:
        path = os.path.join(args.out_dir, f"{r['symbol']}_gwas.tsv")
        fh = open(path, "w")
        fh.write("rsid\teffect_allele\tother_allele\teaf\tbeta\tse\tpvalue\n")
        out_files[r["symbol"]] = fh

    with gzip.open(args.input, "rt") as fin:
        header = fin.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        required = ["hm_rsid", "hm_chrom", "hm_pos", "hm_effect_allele",
                    "hm_other_allele", "hm_effect_allele_frequency",
                    "hm_beta", "standard_error", "p_value"]
        missing = [c for c in required if c not in idx]
        if missing:
            sys.exit(f"ERROR: missing columns {missing}")

        for line in fin:
            f = line.rstrip("\n").split("\t")
            try:
                chrom = int(f[idx["hm_chrom"]])
                pos = int(f[idx["hm_pos"]])
            except (ValueError, IndexError):
                continue

            for r in regions:
                if chrom == r["chrom"] and r["start"] <= pos <= r["end"]:
                    rsid = f[idx["hm_rsid"]]
                    a1 = f[idx["hm_effect_allele"]]
                    a2 = f[idx["hm_other_allele"]]
                    eaf = f[idx["hm_effect_allele_frequency"]]
                    b = f[idx["hm_beta"]]
                    se = f[idx["standard_error"]]
                    p = f[idx["p_value"]]
                    if not rsid or rsid == "NA" or "NA" in (a1, a2, eaf, b, se, p):
                        continue
                    out_files[r["symbol"]].write(f"{rsid}\t{a1}\t{a2}\t{eaf}\t{b}\t{se}\t{p}\n")
                    counts[r["symbol"]] += 1

    for fh in out_files.values():
        fh.close()

    for r in regions:
        print(f"{r['symbol']}\tchr{r['chrom']}:{r['start']}-{r['end']}\t{counts[r['symbol']]} SNPs")


if __name__ == "__main__":
    main()
