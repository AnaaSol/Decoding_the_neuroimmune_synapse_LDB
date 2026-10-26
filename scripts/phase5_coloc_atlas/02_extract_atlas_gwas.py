#!/usr/bin/env python3
"""
Phase 5 / step 2: extract LBD GWAS summary stats (GRCh38 harmonised,
GCST90001390) in +/- locus_window_bp around each lead SNP, in one pass.
Output: data/coloc_atlas/gwas/<locus>.tsv (pos, rsid, other, effect, beta, se, eaf, p)
"""
import gzip, os, yaml
P = "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
cfg = yaml.safe_load(open(f"{P}/config/pipeline_config.yaml"))["phase5_coloc_atlas"]
w = int(cfg["locus_window_bp"])
wins = [(str(l["chrom"]), l["lead_pos"] - w, l["lead_pos"] + w, l["locus"]) for l in cfg["loci"]]
outs = {l: open(f"{P}/data/coloc_atlas/gwas/{l}.tsv", "w") for *_, l in wins}
for fh in outs.values():
    fh.write("pos\trsid\tother_allele\teffect_allele\tbeta\tse\teaf\tpvalue\n")
with gzip.open(f"{P}/data/LBD_GWAS_harmonised.tsv.gz", "rt") as g:
    hdr = g.readline().rstrip("\n").split("\t")
    ix = {k: hdr.index(k) for k in ["hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele",
                                     "hm_beta", "standard_error", "hm_effect_allele_frequency", "p_value"]}
    for line in g:
        f = line.rstrip("\n").split("\t")
        c = f[ix["hm_chrom"]]
        for wc, s, e, l in wins:
            if c == wc and f[ix["hm_pos"]] not in ("", "NA") and s <= int(f[ix["hm_pos"]]) <= e:
                vals = [f[ix[k]] for k in ["hm_pos", "hm_rsid", "hm_other_allele", "hm_effect_allele",
                                           "hm_beta", "standard_error", "hm_effect_allele_frequency", "p_value"]]
                if "NA" not in vals and "" not in vals:
                    outs[l].write("\t".join(vals) + "\n")
for fh in outs.values():
    fh.close()
