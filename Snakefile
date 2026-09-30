# ============================================================================
# Snakemake Pipeline — LBD Neuroimmune Synapse Study
# ============================================================================
# Run the full pipeline:
#   snakemake --cores 4
#
# Dry-run (show what would execute):
#   snakemake -n --reason
#
# Rebuild from a specific rule forward:
#   snakemake --forcerun extract_cd4_signatures --cores 4
#
# Visualize the DAG:
#   snakemake --dag | dot -Tsvg > dag.svg
# ============================================================================

configfile: "config/pipeline_config.yaml"

import os

PROJECT_DIR = config["project"]["dir"]
DATA_DIR    = os.path.join(PROJECT_DIR, "data")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results_final")
SCRIPTS_DIR = os.path.join(PROJECT_DIR, "scripts")
LOGS_DIR    = os.path.join(PROJECT_DIR, "logs")

# All 28 GSE178146 samples (verified from GEO, 2026-07-16)
P3_SAMPLES = (
    config["datasets"]["cortex_snrna"]["sample_groups"]["Control"] +
    config["datasets"]["cortex_snrna"]["sample_groups"]["DLB"]     +
    config["datasets"]["cortex_snrna"]["sample_groups"]["PDD"]     +
    config["datasets"]["cortex_snrna"]["sample_groups"]["PD"]
)

# Phase 2 sample IDs
P2_SAMPLES = list(config["datasets"]["csf_scrna"]["samples"].keys())

# ── ALL: Final outputs that trigger the full pipeline ───────────────────────

rule all:
    input:
        # Phase 1
        os.path.join(DATA_DIR, "LBD_risk_genes_ranked.tsv"),
        expand(os.path.join(DATA_DIR, "gene_sets", "LBD_risk_genes_top{pct}pct.txt"),
               pct=config["phase1"]["gene_set_percentiles"]),
        # Phase 1b (chr19 BCAM/APOE conditional analysis)
        os.path.join(DATA_DIR, "gcta_cojo", "chr19_cond_on_all3.cma.cojo"),
        # Phase 1c (GWAS-vs-eQTL colocalization, brain + T-cell)
        os.path.join(DATA_DIR, "coloc", "results", "coloc_PP3_PP4_summary_all_tissues.pdf"),
        # Phase 1d (MAGMA immune-pathway gene-set enrichment)
        os.path.join(DATA_DIR, "geneset_analysis", "LBD_immune_geneset.gsa.out"),
        # Phase 1e/1f (Mendelian randomization, SuSiE-coloc)
        os.path.join(DATA_DIR, "mr", "results", "mendelian_randomization_summary.csv"),
        os.path.join(DATA_DIR, "susie_coloc", "MMRN1_susie_coloc_blueprint_tcell.csv"),
        # Phase 2
        os.path.join(DATA_DIR, "GSE161192_processed", "DEG_expanded_vs_nonexpanded_CD4.csv"),
        os.path.join(DATA_DIR, "GSE161192_processed", "DEG_pseudobulk_expanded_vs_nonexpanded_CD4.csv"),
        os.path.join(DATA_DIR, "GSE161192_processed", "expanded_CD4_upregulated_genes.txt"),
        os.path.join(DATA_DIR, "GSE161192_processed", "clonotype_size_table.csv"),
        os.path.join(DATA_DIR, "GSE141578_processed", "pseudobulk_DEG_disease_vs_HC.csv"),
        # Phase 2c (batch sensitivity)
        os.path.join(DATA_DIR, "GSE141578_processed", "pseudobulk_DEG_batch_sensitivity.csv"),
        # Phase 3
        os.path.join(DATA_DIR, "GSE178146_processed", "active_ligand_receptor_pairs.csv"),
        os.path.join(DATA_DIR, "GSE178146_processed", "receptor_expression_DLB_vs_Control.csv"),
        os.path.join(DATA_DIR, "GSE178146_processed", "receptor_expression_by_neuron_subtype.csv"),
        os.path.join(DATA_DIR, "GSE178146_processed", "receptor_doseresponse_pseudobulk.csv"),
        # NOTE: Phase 3c (nichenet_ligand_activity) deliberately NOT in this
        # default target list -- it requires the isolated `nichenet_env`
        # conda environment (ggplot2>=4.0) and a manually-fetched 262MB prior
        # model (data/nichenet_prior/, not produced by any rule). Run it
        # explicitly: `snakemake nichenet_ligand_activity --cores 1`.
        # Phase 4
        os.path.join(RESULTS_DIR, "weighted_neuroimmune_network.csv"),
        os.path.join(RESULTS_DIR, "plots", "integrative_prioritization_network.pdf"),
        # Phase 4b (permutation null)
        os.path.join(RESULTS_DIR, "permutation_null_prioritization_score.csv"),


# ── PHASE 1: Genetic Prioritization ─────────────────────────────────────────

rule filter_gwas:
    """
    Filter (or pass-through) GWAS harmonised file for MAGMA.
    Default p_threshold=1.0 keeps all SNPs — required for a valid
    genome-wide gene-based test. The original 1e-5 filter left only
    454 SNPs → 34 genes (defect D6).
    """
    input:
        os.path.join(DATA_DIR, "LBD_GWAS_harmonised.tsv.gz")
    output:
        os.path.join(DATA_DIR, "LBD_GWAS_filtered.tsv.gz")
    params:
        p_threshold=config["phase1"]["magma_filter_p_threshold"]
    log:
        os.path.join(LOGS_DIR, "phase1_filter_gwas.log")
    shell:
        "python3 {SCRIPTS_DIR}/phase1_genetic_prior/01_filter_gwas_data.py "
        "--input {input} --output {output} --p-threshold {params.p_threshold} "
        "> {log} 2>&1"


rule prepare_magma:
    """
    Convert filtered GWAS to MAGMA-compatible format with correct N=6618.
    """
    input:
        gwas=os.path.join(DATA_DIR, "LBD_GWAS_filtered.tsv.gz"),
        meta=os.path.join(DATA_DIR, "LBD_GWAS_metadata.yaml")
    output:
        os.path.join(DATA_DIR, "LBD_GWAS_magma_input.txt")
    log:
        os.path.join(LOGS_DIR, "phase1_prepare_magma.log")
    shell:
        "python3 {SCRIPTS_DIR}/phase1_genetic_prior/02_prepare_magma_input.py "
        "> {log} 2>&1"


rule run_magma:
    """
    Run MAGMA SNPwise-mean gene-based test using GRCh38 coordinates and
    1000G-EUR LD reference. With all SNPs, expect ~18,000 genes analyzed.
    """
    input:
        magma_input=os.path.join(DATA_DIR, "LBD_GWAS_magma_input.txt"),
        gene_loc=config["phase1"]["magma_gene_loc"],
        bim=config["phase1"]["magma_reference"] + ".bim"
    output:
        genes_out=os.path.join(DATA_DIR, "magma_results", "LBD_gene_analysis.genes.out"),
        genes_annot=os.path.join(DATA_DIR, "magma_results", "LBD_gene_analysis.genes.annot")
    log:
        os.path.join(LOGS_DIR, "phase1_run_magma.log")
    shell:
        "bash {SCRIPTS_DIR}/phase1_genetic_prior/03_run_magma_analysis.sh "
        "> {log} 2>&1"


rule parse_magma:
    """
    Parse MAGMA output → ranked gene table + gene-set files.
    """
    input:
        os.path.join(DATA_DIR, "magma_results", "LBD_gene_analysis.genes.out")
    output:
        ranked=os.path.join(DATA_DIR, "LBD_risk_genes_ranked.tsv"),
        gene_sets=expand(
            os.path.join(DATA_DIR, "gene_sets", "LBD_risk_genes_top{pct}pct.txt"),
            pct=config["phase1"]["gene_set_percentiles"]
        )
    log:
        os.path.join(LOGS_DIR, "phase1_parse_magma.log")
    shell:
        "python3 {SCRIPTS_DIR}/phase1_genetic_prior/04_parse_magma_results.py "
        "> {log} 2>&1"


rule extract_chr19_locus:
    """
    Extract chr19 APOE/BCAM-locus summary stats (GCTA .ma format) + rsID list
    from the full harmonised GWAS, by rsID (never by BP -- reference panel is
    GRCh37/hg19, GWAS is GRCh38; see D13/D23).
    """
    input:
        os.path.join(DATA_DIR, "LBD_GWAS_harmonised.tsv.gz")
    output:
        ma=os.path.join(DATA_DIR, "gcta_cojo", "chr19_APOE_locus.ma"),
        rsids=os.path.join(DATA_DIR, "gcta_cojo", "chr19_rsids.txt")
    params:
        chrom=config["phase1b_gcta_cojo"]["region_chrom"],
        start=config["phase1b_gcta_cojo"]["region_start_bp"],
        end=config["phase1b_gcta_cojo"]["region_end_bp"],
        n_total=config["datasets"]["gwas"]["n_total"]
    log:
        os.path.join(LOGS_DIR, "phase1b_extract_chr19_locus.log")
    shell:
        "python3 {SCRIPTS_DIR}/phase1_genetic_prior/05_extract_chr19_locus.py "
        "--input {input} --chrom {params.chrom} --start {params.start} --end {params.end} "
        "--n-total {params.n_total} --out-ma {output.ma} --out-rsids {output.rsids} "
        "> {log} 2>&1"


rule gcta_cojo_bcam_apoe:
    """
    GCTA-COJO conditional analysis: is BCAM's top-MAGMA-gene signal
    independent of the APOE locus, or fully attributable to LD? (D2/Sec 5.2).
    Resolution (2026-08-19): BCAM's top SNP drops from p=6.6e-13 to pC=0.56
    once conditioned jointly on the locus's 3 true independent signals.
    """
    input:
        ma=os.path.join(DATA_DIR, "gcta_cojo", "chr19_APOE_locus.ma"),
        rsids=os.path.join(DATA_DIR, "gcta_cojo", "chr19_rsids.txt"),
        ref_bed=config["phase1b_gcta_cojo"]["ref_panel"] + ".bed"
    output:
        stepwise=os.path.join(DATA_DIR, "gcta_cojo", "chr19_stepwise.jma.cojo"),
        cond_all3=os.path.join(DATA_DIR, "gcta_cojo", "chr19_cond_on_all3.cma.cojo")
    log:
        os.path.join(LOGS_DIR, "phase1b_gcta_cojo.log")
    shell:
        "bash {SCRIPTS_DIR}/phase1_genetic_prior/06_run_gcta_cojo.sh "
        "> {log} 2>&1"


rule extract_coloc_gwas_loci:
    """
    Single-pass extraction of GWAS summary stats for each Phase 1c coloc
    target gene's cis-window (chr19 Bonferroni genes + SNCA/MSTO1/MMRN1 +
    Phase 3 candidate receptors TGFBR2/LGALS9/CD86).
    """
    input:
        os.path.join(DATA_DIR, "LBD_GWAS_harmonised.tsv.gz")
    output:
        expand(os.path.join(DATA_DIR, "coloc", "gwas", "{gene}_gwas.tsv"),
               gene=[g["symbol"] for g in config["phase1c_coloc"]["target_genes"]])
    log:
        os.path.join(LOGS_DIR, "phase1c_extract_gwas_loci.log")
    shell:
        "python3 {SCRIPTS_DIR}/phase1_genetic_prior/07_extract_coloc_gwas_loci.py "
        "--config config/pipeline_config.yaml --input {input} "
        "--out-dir " + os.path.join(DATA_DIR, "coloc", "gwas") +
        " > {log} 2>&1"


COLOC_GENES = [g["symbol"] for g in config["phase1c_coloc"]["target_genes"]]
COLOC_SOURCES = ["brain_acc"] + [s["name"] for s in config["phase1c_coloc"]["eqtl_sources_extra"]]

rule fetch_eqtl_regions:
    """
    Remote tabix region queries against the EBI eQTL Catalogue -- does NOT
    bulk-download the multi-GB dataset file. Requires network access and
    `tabix` on PATH. {source} = brain_acc (GTEx anterior cingulate cortex,
    tissue-matched to GSE178146) or blueprint_tcell (BLUEPRINT CD4+ naive
    T cell, matched to the study's T-cell-mediated mechanistic hypothesis).
    """
    output:
        expand(os.path.join(DATA_DIR, "coloc", "eqtl", "{{source}}", "{gene}_eqtl_raw.tsv"),
               gene=COLOC_GENES)
    log:
        os.path.join(LOGS_DIR, "phase1c_fetch_eqtl_regions_{source}.log")
    shell:
        "bash {SCRIPTS_DIR}/phase1_genetic_prior/08_fetch_eqtl_regions.sh {wildcards.source} "
        "> {log} 2>&1"


rule run_coloc:
    """
    coloc.abf colocalization: LBD GWAS vs eQTL, per target gene, per source.
    """
    input:
        gwas=expand(os.path.join(DATA_DIR, "coloc", "gwas", "{gene}_gwas.tsv"), gene=COLOC_GENES),
        eqtl=expand(os.path.join(DATA_DIR, "coloc", "eqtl", "{{source}}", "{gene}_eqtl_raw.tsv"), gene=COLOC_GENES)
    output:
        summary=os.path.join(DATA_DIR, "coloc", "results", "{source}", "coloc_summary.csv")
    log:
        os.path.join(LOGS_DIR, "phase1c_run_coloc_{source}.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase1_genetic_prior/09_run_coloc.R {wildcards.source} > {log} 2>&1"


rule plot_coloc_summary:
    """
    Combines all per-source coloc_summary.csv files into one cross-tissue
    comparison table + plot (e.g. MMRN1: PP4=0.97 in CD4+ T cells vs
    PP4=0.01/PP3=0.78 in brain ACC -- a tissue-specific colocalization).
    """
    input:
        expand(os.path.join(DATA_DIR, "coloc", "results", "{source}", "coloc_summary.csv"),
               source=COLOC_SOURCES)
    output:
        os.path.join(DATA_DIR, "coloc", "results", "coloc_PP3_PP4_summary_all_tissues.pdf"),
        os.path.join(DATA_DIR, "coloc", "results", "coloc_summary_all_tissues.csv")
    log:
        os.path.join(LOGS_DIR, "phase1c_plot_coloc_summary.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase1_genetic_prior/10_plot_coloc_summary.R > {log} 2>&1"


# ── PHASE 1d: MAGMA competitive gene-set analysis (immune pathways) ─────────

rule fetch_kegg_genesets:
    """
    Fetch pre-specified immune-pathway gene sets (KEGG) for MAGMA --set-annot.
    Not selected post-hoc -- see config phase1d_geneset for the fixed list.
    """
    output:
        os.path.join(DATA_DIR, "geneset_analysis", "immune_pathways.set_annot.txt")
    log:
        os.path.join(LOGS_DIR, "phase1d_fetch_kegg_genesets.log")
    shell:
        "python3 {SCRIPTS_DIR}/phase1_genetic_prior/11_fetch_kegg_genesets.py "
        "--config config/pipeline_config.yaml --out {output} > {log} 2>&1"


rule run_magma_geneset:
    """
    Competitive gene-set analysis: is LBD genetic risk enriched in curated
    immune pathways at the polygenic level (more power than single-gene
    Bonferroni tests)? Result (2026-09-29): none of the 5 pre-specified
    pathways reach significance (all P>0.18); an honest null at this level.
    """
    input:
        gene_results=os.path.join(DATA_DIR, "magma_results", "LBD_gene_analysis.genes.raw"),
        set_file=os.path.join(DATA_DIR, "geneset_analysis", "immune_pathways.set_annot.txt")
    output:
        os.path.join(DATA_DIR, "geneset_analysis", "LBD_immune_geneset.gsa.out")
    log:
        os.path.join(LOGS_DIR, "phase1d_run_magma_geneset.log")
    shell:
        "bash {SCRIPTS_DIR}/phase1_genetic_prior/12_run_magma_geneset.sh > {log} 2>&1"


# ── PHASE 1e/1f: Causal verification (MR, SuSiE-coloc) ──────────────────────

rule run_mendelian_randomization:
    """
    Two-sample MR (Wald ratio / IVW via GCTA-COJO instruments for chr19).
    Result: MMRN1 in T cells beta=-0.29, padj=2.5e-10; no instrument
    (eQTL p<1e-5) for the rest of the panel in either tissue.
    """
    input:
        expand(os.path.join(DATA_DIR, "coloc", "gwas", "{gene}_gwas.tsv"), gene=COLOC_GENES),
        expand(os.path.join(DATA_DIR, "coloc", "eqtl", "{source}", "{gene}_eqtl_raw.tsv"),
               source=COLOC_SOURCES, gene=COLOC_GENES)
    output:
        os.path.join(DATA_DIR, "mr", "results", "mendelian_randomization_summary.csv")
    log:
        os.path.join(LOGS_DIR, "phase1e_mendelian_randomization.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase1_genetic_prior/13_run_mendelian_randomization.R > {log} 2>&1"


rule run_susie_coloc_mmrn1:
    """
    SuSiE-coloc fine-mapping for MMRN1 (multi-causal-variant, vs. coloc.abf's
    single-variant assumption). Result: the credible set anchored on
    rs7680557 (same lead SNP as coloc.abf/MR) colocalizes at PP4=0.977,
    validating the simple coloc.abf PP4=0.97 as a specific signal, not a
    blend across the locus's multiple signals.
    """
    input:
        os.path.join(DATA_DIR, "coloc", "gwas", "MMRN1_gwas.tsv"),
        os.path.join(DATA_DIR, "coloc", "eqtl", "blueprint_tcell", "MMRN1_eqtl_raw.tsv")
    output:
        os.path.join(DATA_DIR, "susie_coloc", "MMRN1_susie_coloc_blueprint_tcell.csv")
    log:
        os.path.join(LOGS_DIR, "phase1f_susie_coloc_tcell.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase1_genetic_prior/14_run_susie_coloc_mmrn1.R blueprint_tcell > {log} 2>&1"


# ── PHASE 2: CSF Immune Profiling ──────────────────────────────────────────

rule process_scrna:
    """
    Load, QC, normalize, and cluster CSF T cells (GSE161192).
    nFeature_max raised to 6000 (was 2500) to retain activated clones.
    """
    input:
        expand(
            os.path.join(DATA_DIR, "GSE161192_raw",
                         "{gsm}_{sample}_filtered_feature_bc_matrix.zip"),
            zip,
            gsm=["GSM4905012","GSM4905013","GSM4905015","GSM4905014"],
            sample=["P1U","P1S","P2U","P2S"]
        )
    output:
        os.path.join(DATA_DIR, "GSE161192_processed", "CSF_Tcells_seurat_processed.rds")
    log:
        os.path.join(LOGS_DIR, "phase2_process_scrna.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase2_immune_profiling/01_process_scrna_seurat.R "
        "> {log} 2>&1"


rule integrate_tcr:
    """
    Integrate TCR contig annotations with Seurat object.
    Computes clonal expansion (aa-level CDR3, ≥2 cells primary; ≥3 cells strict).
    """
    input:
        seurat=os.path.join(DATA_DIR, "GSE161192_processed", "CSF_Tcells_seurat_processed.rds"),
        tcr=expand(
            os.path.join(DATA_DIR, "GSE161192_raw",
                         "{gsm}_{sample}_filtered_contig_annotations.csv.gz"),
            zip,
            gsm=["GSM4891327","GSM4891328","GSM4891329","GSM4891330"],
            sample=["P1U","P1S","P2U","P2S"]
        )
    output:
        seurat=os.path.join(DATA_DIR, "GSE161192_processed", "CSF_Tcells_with_TCR.rds"),
        all_cells=os.path.join(DATA_DIR, "GSE161192_processed", "all_cells_with_tcr_annotations.csv"),
        clone_sizes=os.path.join(DATA_DIR, "GSE161192_processed", "clonotype_size_table.csv")
    log:
        os.path.join(LOGS_DIR, "phase2_integrate_tcr.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase2_immune_profiling/02_integrate_tcr_data.R "
        "> {log} 2>&1"


rule extract_cd4_signatures:
    """
    Differential expression: expanded vs non-expanded CD4+ T cells.
    Wilcoxon (inflated significance, kept for comparison) +
    DESeq2 pseudobulk with ~stimulation + expansion design (preferred).
    CAVEAT: N=2 donors — results are directional only.
    """
    input:
        seurat=os.path.join(DATA_DIR, "GSE161192_processed", "CSF_Tcells_with_TCR.rds"),
        risk_genes=os.path.join(DATA_DIR, "LBD_risk_genes_ranked.tsv")
    output:
        deg_wilcox=os.path.join(DATA_DIR, "GSE161192_processed", "DEG_expanded_vs_nonexpanded_CD4.csv"),
        deg_pseudo=os.path.join(DATA_DIR, "GSE161192_processed", "DEG_pseudobulk_expanded_vs_nonexpanded_CD4.csv"),
        upregulated=os.path.join(DATA_DIR, "GSE161192_processed", "expanded_CD4_upregulated_genes.txt")
    log:
        os.path.join(LOGS_DIR, "phase2_extract_cd4_signatures.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase2_immune_profiling/03_extract_cd4_signatures.R "
        "> {log} 2>&1"


rule gse141578_pseudobulk:
    """
    GSE141578 — larger, unstimulated CSF cohort (n=11 HC, n=9 PD, n=2 DLB;
    same paper as GSE161192, Gate et al. 2021 Science). Independently tests
    the CXCR4/CXCL12 claim on the dataset the source paper actually used for
    that finding, rather than on the small stimulated GSE161192 cohort.
    CSF32 corrected from "excluded/unreported" to HC (verified from GSM4404059).
    """
    input:
        raw_dir=os.path.join(DATA_DIR, "GSE141578_raw", "csf_10x")
    output:
        deg=os.path.join(DATA_DIR, "GSE141578_processed", "pseudobulk_DEG_disease_vs_HC.csv"),
        upregulated=os.path.join(DATA_DIR, "GSE141578_processed", "cd4_disease_upregulated_genes.txt")
    log:
        os.path.join(LOGS_DIR, "phase2_gse141578_pseudobulk.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase2_immune_profiling/04_GSE141578_disease_vs_HC_pseudobulk.R "
        "> {log} 2>&1"


rule gse141578_batch_sensitivity:
    """
    Re-fits the GSE141578 pseudobulk matrix (a) without the batch covariate
    and (b) PD-only vs HC (drops the 2 DLB samples, both in batch2). Result:
    0 significant genes in both alternative models -- the CXCL12/CXCR4
    non-replication does not depend on the batch-adjustment choice.
    """
    input:
        os.path.join(DATA_DIR, "GSE141578_processed", "GSE141578_CD4_Tcells_seurat.rds")
    output:
        os.path.join(DATA_DIR, "GSE141578_processed", "pseudobulk_DEG_batch_sensitivity.csv")
    log:
        os.path.join(LOGS_DIR, "phase2c_batch_sensitivity.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase2_immune_profiling/05_GSE141578_batch_sensitivity.R "
        "> {log} 2>&1"


# ── PHASE 3: Neuronal Vulnerability ────────────────────────────────────────

rule load_qc_snnucseq:
    """
    Load and QC all 28 GSE178146 snRNA-seq samples.
    Applies verified sample→group mapping (Control/DLB/PDD/PD).
    Previous grepl("^C") bug misclassified PDC05-PDC91 as "PD" (defect D1).
    """
    input:
        raw_dir=os.path.join(DATA_DIR, "GSE178146_raw")
    output:
        per_sample=expand(
            os.path.join(DATA_DIR, "GSE178146_processed", "{sample}_seurat_qc.rds"),
            sample=P3_SAMPLES
        ),
        qc_meta=os.path.join(DATA_DIR, "GSE178146_processed", "qc_metadata_only.rds")
    log:
        os.path.join(LOGS_DIR, "phase3_load_qc_snnucseq.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase3_neuronal_vulnerability/01_load_and_qc_snnucseq.R "
        "> {log} 2>&1"


rule downsample_integrate_cluster:
    """
    Downsample (max {max_cells}/sample), Harmony batch correction, UMAP, clustering.
    """
    input:
        expand(
            os.path.join(DATA_DIR, "GSE178146_processed", "{sample}_seurat_qc.rds"),
            sample=P3_SAMPLES
        )
    output:
        os.path.join(DATA_DIR, "GSE178146_processed", "GSE178146_integrated.rds")
    params:
        max_cells=config["phase3"]["integration"]["max_cells_per_sample"]
    log:
        os.path.join(LOGS_DIR, "phase3_downsample_integrate_cluster.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase3_neuronal_vulnerability/02_downsampling_integrate_and_cluster.R "
        "> {log} 2>&1"


rule identify_neuronal_subtypes:
    """
    Classify neurons as excitatory (SLC17A7+) or inhibitory (GAD1/GAD2+),
    assign cortical layer identities.
    """
    input:
        os.path.join(DATA_DIR, "GSE178146_processed", "GSE178146_integrated.rds")
    output:
        os.path.join(DATA_DIR, "GSE178146_processed", "ACC_neurons_only.rds")
    log:
        os.path.join(LOGS_DIR, "phase3_identify_neuronal_subtypes.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase3_neuronal_vulnerability/03_identify_neuronal_subtypes.R "
        "> {log} 2>&1"


rule ligand_receptor:
    """
    Map CD4+ T cell ligands → neuronal receptors.
    Primary comparison: DLB vs Control (fixed from incorrect PD vs Control).
    CXCL12/CXCR4 and HLA/TCR NOT present in database — remove from manuscript claims.
    """
    input:
        neurons=os.path.join(DATA_DIR, "GSE178146_processed", "ACC_neurons_only.rds"),
        cd4_ligands=os.path.join(DATA_DIR, "GSE161192_processed", "expanded_CD4_upregulated_genes.txt")
    output:
        active_pairs=os.path.join(DATA_DIR, "GSE178146_processed", "active_ligand_receptor_pairs.csv"),
        exp_celltype=os.path.join(DATA_DIR, "GSE178146_processed", "receptor_expression_by_celltype.csv"),
        exp_subtype=os.path.join(DATA_DIR, "GSE178146_processed", "receptor_expression_by_neuron_subtype.csv"),
        exp_condition=os.path.join(DATA_DIR, "GSE178146_processed", "receptor_expression_DLB_vs_Control.csv")
    log:
        os.path.join(LOGS_DIR, "phase3_ligand_receptor.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase3_neuronal_vulnerability/04_ligand_receptor_analysis.R "
        "> {log} 2>&1"


rule receptor_doseresponse:
    """
    Donor-level pseudobulk DESeq2 trend test: does receptor expression scale
    monotonically with LBD clinical severity (Control=0 < PD=1 < PDD=2 <
    DLB=3, pre-specified ordering)? More powerful than 3 separate pairwise
    contrasts, and tests the actual "candidate hypothesis" framing in the
    manuscript directly rather than leaving it as a single DLB-vs-Control
    snapshot. Result (2026-09-29): LIFR, LGALS9, CD86 all show a significant
    (padj<0.05), monotonic increase with severity; TGFBR2 near-significant
    (padj=0.067) with the same direction.
    """
    input:
        neurons=os.path.join(DATA_DIR, "GSE178146_processed", "ACC_neurons_only.rds"),
        active_pairs=os.path.join(DATA_DIR, "GSE178146_processed", "active_ligand_receptor_pairs.csv")
    output:
        os.path.join(DATA_DIR, "GSE178146_processed", "receptor_doseresponse_pseudobulk.csv")
    log:
        os.path.join(LOGS_DIR, "phase3_receptor_doseresponse.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase3_neuronal_vulnerability/05_receptor_doseresponse_pseudobulk.R "
        "> {log} 2>&1"


rule nichenet_ligand_activity:
    """
    Directional ligand-activity analysis (NicheNet), run in an isolated
    conda env (nichenet_env) so its ggplot2>=4.0 dependency never touches
    this pipeline's main library (ggplot2 3.5.2). Result: OSM ranks 28/29
    candidate ligands by directional activity (AUROC=0.508), IL17A ranks
    29/29 -- reverses the mechanistic reading of OSM-LIFR from Phase 4.
    IMPORTANT: must `cd` away from the project directory before invoking
    this env's Rscript, or this project's own .Rprofile/renv activation
    silently contaminates .libPaths() with the main env's R 4.3.1 packages
    inside nichenet_env's R 4.3.3 -- causes confusing dyn.load failures
    that look like missing dependencies but are actually a version clash.
    """
    input:
        os.path.join(DATA_DIR, "nichenet_prior", "ligand_target_matrix_human.rds"),
        os.path.join(DATA_DIR, "nichenet_prior", "lr_network_human.rds"),
        os.path.join(PROJECT_DIR, "archive", "orphaned_nichenet_prep", "neuron_DLB_vs_Control_genomewide_DE.csv")
    output:
        os.path.join(DATA_DIR, "nichenet_results", "ligand_activities.csv")
    log:
        os.path.join(LOGS_DIR, "phase3c_nichenet.log")
    shell:
        "bash -c 'source $(conda info --base)/etc/profile.d/conda.sh && "
        "conda activate nichenet_env && cd /tmp && "
        "Rscript {SCRIPTS_DIR}/phase3_neuronal_vulnerability/06_nichenet_ligand_activity.R' "
        "> {log} 2>&1"


# ── PHASE 4: Weighted Interactomics ─────────────────────────────────────────

rule weighted_interactomics:
    """
    Integrative prioritization network: neuronal expression × GWAS Z.
    Score = Prioritization_Score (NOT Causality_Score; heuristic weighting only).
    """
    input:
        interactions=os.path.join(DATA_DIR, "GSE178146_processed", "active_ligand_receptor_pairs.csv"),
        expression=os.path.join(DATA_DIR, "GSE178146_processed", "receptor_expression_by_celltype.csv"),
        gwas=os.path.join(DATA_DIR, "LBD_risk_genes_ranked.tsv")
    output:
        network=os.path.join(RESULTS_DIR, "weighted_neuroimmune_network.csv"),
        scatter=os.path.join(RESULTS_DIR, "plots", "prioritization_scatter_plot.pdf"),
        plot=os.path.join(RESULTS_DIR, "plots", "integrative_prioritization_network.pdf")
    log:
        os.path.join(LOGS_DIR, "phase4_weighted_interactomics.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase4_weighted_interactomics/01_weighted_interactomics.R "
        "> {log} 2>&1"


rule permutation_null_prioritization:
    """
    Empirical null for the Prioritization Score: permutes receptor GWAS
    Z-scores across candidate pairs (10,000x). Result: OSM-LIFR stays #1 in
    83.7% of permutations, confirming the score is expression-dominated.
    """
    input:
        os.path.join(RESULTS_DIR, "weighted_neuroimmune_network.csv")
    output:
        os.path.join(RESULTS_DIR, "permutation_null_prioritization_score.csv")
    log:
        os.path.join(LOGS_DIR, "phase4b_permutation_null.log")
    shell:
        "Rscript {SCRIPTS_DIR}/phase4_weighted_interactomics/02_permutation_null_prioritization.R "
        "> {log} 2>&1"


# ── UTILITY: Validate pipeline outputs ──────────────────────────────────────

rule validate:
    """
    Run assertion checks against all outputs. Call after full pipeline:
      snakemake validate --cores 1
    """
    input:
        rules.all.input
    log:
        os.path.join(LOGS_DIR, "validation.log")
    shell:
        "Rscript {SCRIPTS_DIR}/utils/validate_pipeline.R > {log} 2>&1 && "
        "echo '✓ All validation checks passed'"
