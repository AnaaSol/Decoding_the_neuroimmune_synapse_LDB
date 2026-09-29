#!/usr/bin/env bash
################################################################################
# Fetch per-gene cis-window eQTL summary stats from the EBI eQTL Catalogue via
# tabix (remote region query -- does NOT download the full multi-GB dataset
# file).
#
# Usage: 08_fetch_eqtl_regions.sh [SOURCE]
#   SOURCE = "brain_acc" (default) -- GTEx anterior cingulate cortex (QTD000151),
#            tissue-matched to the GSE178146 snRNA-seq dataset used in Phase 3.
#   SOURCE = "blueprint_tcell"     -- BLUEPRINT CD4+ naive T cell (QTD000031),
#            matched to the study's T-cell-mediated mechanistic hypothesis.
#
# EBI explicitly asks for tabix requests to be spaced out (frequent rapid
# requests can trigger their DoS-detection firewall) -- see
# https://www.ebi.ac.uk/eqtl/Data_access/ -- hence the `sleep` between genes.
################################################################################
set -euo pipefail

SOURCE="${1:-brain_acc}"

PROJECT_DIR="/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB"
cd "$PROJECT_DIR"

CONFIG="config/pipeline_config.yaml"
OUT_DIR="data/coloc/eqtl/${SOURCE}"
mkdir -p "$OUT_DIR"

FTP_URL=$(python3 -c "
import yaml
cfg = yaml.safe_load(open('$CONFIG'))['phase1c_coloc']
src = cfg['eqtl_source'] if '$SOURCE' == 'brain_acc' else next(s for s in cfg['eqtl_sources_extra'] if s['name'] == '$SOURCE')
print(src['ftp_url'])
")
WINDOW=$(python3 -c "import yaml; print(yaml.safe_load(open('$CONFIG'))['phase1c_coloc']['cis_window_bp'])")

echo "== eQTL source: $SOURCE =="
echo "== FTP: $FTP_URL =="

python3 -c "
import yaml
cfg = yaml.safe_load(open('$CONFIG'))
for g in cfg['phase1c_coloc']['target_genes']:
    start = max(1, int(g['start']) - $WINDOW)
    end = int(g['end']) + $WINDOW
    print(f\"{g['symbol']}\t{g['ensembl_id']}\t{g['chrom']}\t{start}\t{end}\")
" > "$OUT_DIR/.gene_regions.tsv"

while IFS=$'\t' read -r symbol ensembl chrom start end; do
    out="$OUT_DIR/${symbol}_eqtl_raw.tsv"
    echo "== $symbol ($ensembl) chr${chrom}:${start}-${end} =="
    {
        echo -e "molecular_trait_id\tchromosome\tposition\tref\talt\tvariant\tma_samples\tmaf\tpvalue\tbeta\tse\ttype\tac\tan\tr2\tmolecular_trait_object_id\tgene_id\tmedian_tpm\trsid"
        tabix "$FTP_URL" "${chrom}:${start}-${end}" | awk -F'\t' -v g="$ensembl" '$17==g'
    } > "$out"
    n=$(($(wc -l < "$out") - 1))
    echo "   -> $n eQTL rows for $ensembl written to $out"
    sleep 3
done < "$OUT_DIR/.gene_regions.tsv"

echo "Done ($SOURCE)."
