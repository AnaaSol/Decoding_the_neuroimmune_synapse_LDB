#!/usr/bin/env python3
"""
Phase 5 / step 1: fetch eQTL summary statistics for every gene around each
genome-wide-significant LBD locus, for every dataset in
config['phase5_coloc_atlas']['datasets'], via remote tabix region queries
against the EBI eQTL Catalogue (no bulk download).

Output: data/coloc_atlas/eqtl/<dataset_id>/<locus>.tsv.gz with columns
  gene_id rsid pos ref alt maf pvalue beta se
Resumable: an existing, non-empty output file is skipped.

EBI asks for spaced tabix requests (https://www.ebi.ac.uk/eqtl/Data_access/),
hence the sleep between queries and strictly sequential execution.
"""
import gzip
import os
import subprocess
import sys
import time

import yaml

PROJECT_DIR = "/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
cfg = yaml.safe_load(open(os.path.join(PROJECT_DIR, "config/pipeline_config.yaml")))["phase5_coloc_atlas"]
out_root = os.path.join(PROJECT_DIR, "data/coloc_atlas/eqtl")
tbi_dir = os.path.join(PROJECT_DIR, "data/coloc_atlas/.tbi_cache")  # tabix drops the remote .tbi in cwd
os.makedirs(tbi_dir, exist_ok=True)

# eQTL Catalogue columns: 0 molecular_trait_id 1 chromosome 2 position 3 ref 4 alt
# 5 variant 6 ma_samples 7 maf 8 pvalue 9 beta 10 se ... 16 gene_id 18 rsid
KEEP = [16, 18, 2, 3, 4, 7, 8, 9, 10]
HEADER = "gene_id\trsid\tpos\tref\talt\tmaf\tpvalue\tbeta\tse\n"
win = int(cfg["locus_window_bp"])

only = set(sys.argv[1:])  # optional dataset-id filter
for ds in cfg["datasets"]:
    if only and ds["id"] not in only:
        continue
    ds_dir = os.path.join(out_root, ds["id"])
    os.makedirs(ds_dir, exist_ok=True)
    for loc in cfg["loci"]:
        out = os.path.join(ds_dir, f"{loc['locus']}.tsv.gz")
        if os.path.exists(out) and os.path.getsize(out) > 100:
            continue
        region = f"{loc['chrom']}:{max(1, loc['lead_pos'] - win)}-{loc['lead_pos'] + win}"
        t0 = time.time()
        for attempt in range(3):
            proc = subprocess.run(["tabix", ds["url"], region], cwd=tbi_dir,
                                  capture_output=True, text=True)
            if proc.returncode == 0:
                break
            time.sleep(20 * (attempt + 1))
        else:
            print(f"FAILED {ds['id']} {loc['locus']}: {proc.stderr.strip()[:200]}", flush=True)
            continue
        n = 0
        tmp = out + ".part"
        with gzip.open(tmp, "wt") as fh:
            fh.write(HEADER)
            for line in proc.stdout.splitlines():
                f = line.split("\t")
                if len(f) < 19:
                    continue
                fh.write("\t".join(f[i] for i in KEEP) + "\n")
                n += 1
        os.replace(tmp, out)
        print(f"{ds['id']} ({ds['study']}:{ds['group']}) {loc['locus']} {region}: "
              f"{n} rows, {time.time() - t0:.0f}s", flush=True)
        time.sleep(2)
print("DONE", flush=True)
