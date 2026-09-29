#!/usr/bin/env python3
"""
Fetch immune-pathway gene sets from the KEGG REST API and write a MAGMA
row-based --set-annot file (SET_NAME gene1 gene2 ... per line, Entrez IDs --
KEGG's "hsa:<N>" gene IDs for Homo sapiens ARE NCBI Entrez Gene IDs, so no ID
mapping step is needed).

Pathways are pre-specified in config (phase1d_geneset.kegg_pathways), chosen
from mechanisms already discussed in the manuscript's own narrative (Th17/
IL-17A, TCR signaling/clonal expansion, cytokine receptors, CXCL12-CXCR4
chemokine signaling, OSM/JAK-STAT) -- not selected post-hoc after seeing
results.
"""

import argparse
import sys
import time
import urllib.request

import yaml


def fetch_pathway_genes(pathway_id):
    url = f"https://rest.kegg.jp/link/hsa/pathway:{pathway_id}"
    with urllib.request.urlopen(url, timeout=20) as resp:
        text = resp.read().decode("utf-8")
    genes = []
    for line in text.strip().split("\n"):
        if not line:
            continue
        _, gene_field = line.split("\t")
        entrez_id = gene_field.split(":")[1]
        genes.append(entrez_id)
    return genes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    pathways = cfg["phase1d_geneset"]["kegg_pathways"]

    with open(args.out, "w") as out:
        for p in pathways:
            pid, label = p["id"], p["label"]
            try:
                genes = fetch_pathway_genes(pid)
            except Exception as e:
                sys.exit(f"ERROR fetching {pid} ({label}): {e}")
            set_name = pid  # MAGMA set names: keep short/unique; label kept in a sidecar
            out.write(set_name + " " + " ".join(genes) + "\n")
            print(f"{pid}\t{label}\t{len(genes)} genes")
            time.sleep(1)  # be polite to the KEGG API


if __name__ == "__main__":
    main()
