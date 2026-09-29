# Phase A Audit Inventory
_Generated 2026-07-16. Source of truth: scripts and intermediate outputs, not manuscript prose._
_GSE178146 identity confirmed by user from GEO record._

---

## 1. Directory Map

```
Decoding_the_neuroimmune_synapse_LDB/
├── INVESTIGATION_REVIEW_SUMMARY.md   ← handoff document
├── MANUSCRIPT_Neuroimmune_Interactions_LBD.md
├── README.md
├── tools/
│   ├── magma                          ← MAGMA v1.10 binary (linux)
│   └── manual_v1.10.pdf
├── scripts/
│   ├── phase1_genetic_prior/
│   │   ├── 01_filter_gwas_data.py     ← Python; filters GWAS at p≤1e-5
│   │   ├── 02_prepare_magma_input.py  ← Python; formats SNP table for MAGMA
│   │   ├── 03_run_magma_analysis.sh   ← Bash; runs MAGMA annotate + gene-analysis
│   │   ├── 04_parse_magma_results.py  ← Python; maps Entrez→symbol, creates gene sets
│   │   └── magma.log                  ← stale (ran --help only; real log in magma_results/)
│   ├── phase2_immune_profiling/
│   │   ├── 01_process_scrna_seurat.R  ← R/Seurat v5; QC, normalization, UMAP
│   │   ├── 02_integrate_tcr_data.R    ← R/scRepertoire; TCR integration, clonality
│   │   └── 03_extract_cd4_signatures.R ← R; Wilcoxon DE expanded vs non-expanded
│   ├── phase3_neuronal_vulnerability/
│   │   ├── 01_load_and_qc_snnucseq.R  ← R; loads GSE178146; WRONG condition labels
│   │   ├── 02_downsampling_integrate_and_cluster.R ← R; downsample 2000/sample, merge
│   │   ├── 03_identify_neuronal_subtypes.R ← R; cell type annotation
│   │   └── 04_ligand_receptor_analysis.R   ← R; curated L-R database, "PD vs Control"
│   └── phase4_weighted_interactomics/
│       └── 01_weighted_interactomics.R ← R; igraph/ggraph network, "Causality_Score"
├── data/
│   ├── LBD_GWAS_harmonised.tsv.gz    ← raw GWAS (GRCh38, GCST90001390, N=6618)
│   ├── LBD_GWAS_filtered.tsv.gz      ← p≤1e-5 filtered (454 SNPs)
│   ├── LBD_GWAS_magma_input.txt      ← 453 data lines; N=2591 (wrong: should be 6618)
│   ├── LBD_GWAS_metadata.yaml        ← study metadata (N=6,618 confirmed)
│   ├── LBD_risk_genes_ranked.tsv     ← MAGMA output, 34 genes ranked by Z-score
│   ├── gene_sets/                     ← top 1/5/10/25% gene lists (Entrez IDs)
│   ├── magma_results/                 ← .genes.annot, .genes.out, .genes.sorted.txt, .log
│   ├── reference/                     ← NCBI38.gene.loc (GRCh38) + g1000_eur (hg19)
│   ├── GSE161192_raw/                 ← 10X matrices + TCR contigs (P1U,P1S,P2U,P2S)
│   ├── GSE161192_processed/           ← Seurat .rds, DEGs, clonotype CSVs, plots
│   ├── GSE178146_raw/                 ← 28 raw tarballs (GSM5380921–GSM5380948)
│   ├── GSE178146_extracted/           ← 28 per-sample feature-barcode matrices
│   │   ├── C36, C48                   ← Control samples
│   │   ├── PD060–PD747 (21 samples)  ← ALL disease groups (PD, PDD, DLB mixed)
│   │   └── PDC05–PDC91 (5 samples)   ← LIKELY Control group (see §5.D1)
│   └── GSE178146_processed/           ← per-sample QC .rds, BPCells store, integrated .rds
│                                        active_ligand_receptor_pairs.csv (16 pairs)
│                                        receptor_expression_PD_vs_Control.csv
│                                        qc_metadata_only.rds
└── results_final/
    ├── weighted_neuroimmune_network.csv
    └── plots/
        ├── causal_network_graph.pdf
        └── prioritization_scatter_plot.pdf
```

---

## 2. Pipeline Dependency Graph

```
Phase 1 ─────────────────────────────────────────────────────────────────
LBD_GWAS_harmonised.tsv.gz
  → [01_filter_gwas_data.py]     → LBD_GWAS_filtered.tsv.gz
  → [02_prepare_magma_input.py]  → LBD_GWAS_magma_input.txt
  → [03_run_magma_analysis.sh]   → magma_results/LBD_gene_analysis.genes.out
                                   magma_results/LBD_gene_analysis.genes.annot
  → [04_parse_magma_results.py]  → LBD_risk_genes_ranked.tsv
                                   gene_sets/LBD_risk_genes_top{1,5,10,25}pct.txt

Phase 2 ─────────────────────────────────────────────────────────────────
GSE161192_raw/ (10X matrices + TCR contigs)
  → [01_process_scrna_seurat.R]   → CSF_Tcells_seurat_processed.rds
  → [02_integrate_tcr_data.R]     → CSF_Tcells_with_TCR.rds
                                    CD4_expanded_clones.csv
                                    all_cells_with_tcr_annotations.csv
                                    expanded_CD4_clonotypes.csv
  → [03_extract_cd4_signatures.R] → DEG_expanded_vs_nonexpanded_CD4.csv
                                    expanded_CD4_upregulated_genes.txt (642 genes)
                                    expanded_CD4_downregulated_genes.txt

Phase 3 ─────────────────────────────────────────────────────────────────
GSE178146_extracted/ (per-sample 10X matrices)
  → [01_load_and_qc_snnucseq.R]  → {sample}_seurat_qc.rds ×27 (1 missing from expected 28)
  → [02_downsampling_integrate_and_cluster.R] → GSE178146_integrated.rds
  → [03_identify_neuronal_subtypes.R]         → ACC_neurons_only.rds
  → [04_ligand_receptor_analysis.R]
       reads: expanded_CD4_upregulated_genes.txt + ACC_neurons_only.rds
              + curated L-R database (hardcoded, 27 pairs defined)
       → active_ligand_receptor_pairs.csv (16 active pairs)
         receptor_expression_by_celltype.csv
         receptor_expression_PD_vs_Control.csv   ← MISLABELED (see D1)
         ligand_receptor_interaction_summary.csv

Phase 4 ─────────────────────────────────────────────────────────────────
active_ligand_receptor_pairs.csv
receptor_expression_by_celltype.csv
LBD_risk_genes_ranked.tsv
  → [01_weighted_interactomics.R] → weighted_neuroimmune_network.csv
                                    results_final/plots/causal_network_graph.pdf
                                    results_final/plots/prioritization_scatter_plot.pdf
```

---

## 3. Runnability Assessment

### Phase 1 – RUNNABLE (with caveats)

| Script | Status | Issue |
|--------|--------|-------|
| 01_filter_gwas_data.py | BROKEN (standalone) | Hardcoded path has **typo**: `Decoding_the_neuroinmune_synapse_LDB` (missing 'm') and missing leading `/`. Will `FileNotFoundError` from `__main__`. |
| 02_prepare_magma_input.py | RUNNABLE | Relative paths OK when run from `scripts/phase1_genetic_prior/`. Hardcodes N=2,591 (should be 6,618). |
| 03_run_magma_analysis.sh | RUNNABLE | Absolute paths correct. MAGMA binary present. Reference files present. |
| 04_parse_magma_results.py | RUNNABLE | Relative paths OK. Gene annotation present. |

**Critical upstream flaw:** GWAS pre-filtered to p≤1e-5 → only 454 SNPs → 34 genes. MAGMA requires all GWAS SNPs; pre-filtering invalidates the genome-wide gene-based test.

**Build note:** g1000_eur = hg19/GRCh37; NCBI38.gene.loc = GRCh38. Matched 426/454 SNPs by rsID (93.8%).

### Phase 2 – RUNNABLE (installed packages required)

| Script | Status | Issue |
|--------|--------|-------|
| 01_process_scrna_seurat.R | RUNNABLE | nFeature_max=2500 may exclude activated T cells. No `set.seed()`. |
| 02_integrate_tcr_data.R | RUNNABLE | Manual barcode-fix workaround present. `cloneCall="aa"` (not nucleotide). No seed. |
| 03_extract_cd4_signatures.R | RUNNABLE | Wilcoxon DE; no stimulation covariate; no seed. |

### Phase 3 – RUNNABLE but produces wrong results

| Script | Status | Issue |
|--------|--------|-------|
| 01_load_and_qc_snnucseq.R | RUNNABLE | **WRONG condition assignment** (see D1 below). |
| 02_downsampling_integrate_and_cluster.R | RUNNABLE | `set.seed(123)` for downsampling only. |
| 03_identify_neuronal_subtypes.R | RUNNABLE | Depends on Harmony, BPCells. |
| 04_ligand_receptor_analysis.R | RUNNABLE | Hardcoded L-R DB missing CXCL12/CXCR4 and HLA/TCR. Step 7 runs "PD vs Control" (wrong). |

### Phase 4 – RUNNABLE

| Script | Status | Issue |
|--------|--------|-------|
| 01_weighted_interactomics.R | RUNNABLE | "Causal" overclaimed. Simple heuristic score. Requires igraph, ggraph. |

### Missing infrastructure (all phases)

- No `requirements.txt`, `environment.yml`, `renv.lock`, or `sessionInfo()` output
- No seeds in Phase 2 DE or Phase 1 Python scripts
- No Makefile/Snakemake/{targets}
- No config file (disease names, thresholds, paths all hardcoded)

---

## 4. Confirmed Outputs (exist on disk, inspected)

| File | Phase | Key Content |
|------|-------|-------------|
| `LBD_GWAS_magma_input.txt` | 1 | 454 SNPs at p≤1e-5; N=2,591 throughout |
| `magma_results/LBD_gene_analysis.genes.out` | 1 | 34 genes, SNPwise-mean model, N=2,591 |
| `LBD_risk_genes_ranked.tsv` | 1 | 34 genes; BCAM #1 (Z=7.09); BCL3 absent |
| `gene_sets/LBD_risk_genes_top{1,5,10,25}pct.txt` | 1 | Entrez IDs only |
| `CSF_Tcells_seurat_processed.rds` | 2 | QC-filtered merged Seurat (4 samples) |
| `CSF_Tcells_with_TCR.rds` | 2 | + TCR + clonal status |
| `DEG_expanded_vs_nonexpanded_CD4.csv` | 2 | 3,671 genes; Wilcoxon; TPM4 #1 (padj=1.7e-45) |
| `expanded_CD4_upregulated_genes.txt` | 2 | 642 genes; CXCL12 absent |
| `expanded_CD4_downregulated_genes.txt` | 2 | Downregulated gene list |
| `CD4_expanded_clones.csv` | 2 | Expanded CD4+ cells table |
| `{sample}_seurat_qc.rds` × 27 | 3 | Per-sample QC objects (27/28 present) |
| `GSE178146_integrated.rds` | 3 | Integrated (downsampled) object |
| `ACC_neurons_only.rds` | 3 | Neuronal subset |
| `active_ligand_receptor_pairs.csv` | 3 | 16 pairs; NO CXCL12/CXCR4; NO HLA/TCR |
| `receptor_expression_PD_vs_Control.csv` | 3 | 16 receptors; labels: "Control" / "PD" (WRONG) |
| `receptor_expression_by_celltype.csv` | 3 | Expression per neuronal subtype |
| `weighted_neuroimmune_network.csv` | 4 | 16 interactions + Causality_Score |
| `results_final/plots/causal_network_graph.pdf` | 4 | igraph network |
| `results_final/plots/prioritization_scatter_plot.pdf` | 4 | Expression vs Z-score scatter |

---

## 5. Defects Identified

### BLOCKING / Critical

---

**D1 — GSE178146 condition labels are DOUBLY WRONG** ⚠️ TOP PRIORITY

**Dataset identity confirmed:** GSE178146 = Feleke et al. 2021, Acta Neuropathol, 142:449–474. Title: *"Cross-platform transcriptional profiling identifies common and distinct molecular pathologies in Lewy body diseases."* 28 anterior cingulate cortex samples; **n=7 per group: Control, PD, PDD, DLB.** This is the correct LBD-spectrum dataset.

**The 28 samples — fully verified from GEO individual sample pages:**

| Sample ID | GSM | Group | Current script label |
|-----------|-----|-------|----------------------|
| C36 | GSM5380921 | **Control** | Control ✓ |
| C48 | GSM5380922 | **Control** | Control ✓ |
| PD060 | GSM5380923 | **DLB** | "PD" ✗ |
| PD115 | GSM5380924 | **DLB** | "PD" ✗ |
| PD163 | GSM5380925 | **DLB** | "PD" ✗ |
| PD294 | GSM5380926 | **DLB** | "PD" ✗ |
| PD332 | GSM5380927 | **DLB** | "PD" ✗ |
| PD341 | GSM5380928 | **PDD** | "PD" ✗ |
| PD366 | GSM5380929 | **PDD** | "PD" ✗ |
| PD413 | GSM5380930 | **PD** | "PD" (accidentally correct) |
| PD415 | GSM5380931 | **PDD** | "PD" ✗ |
| PD416 | GSM5380932 | **PD** | "PD" (accidentally correct) |
| PD501 | GSM5380933 | **PDD** | "PD" ✗ |
| PD523 | GSM5380934 | **PD** | "PD" (accidentally correct) |
| PD531 | GSM5380935 | **PDD** | "PD" ✗ |
| PD563 | GSM5380936 | **PDD** | "PD" ✗ |
| PD566 | GSM5380937 | **DLB** | "PD" ✗ |
| PD666 | GSM5380938 | **PD** | "PD" (accidentally correct) |
| PD678 | GSM5380939 | **PDD** | "PD" ✗ |
| PD683 | GSM5380940 | **PD** | "PD" (accidentally correct) |
| PD706 | GSM5380941 | **DLB** | "PD" ✗ |
| PD732 | GSM5380942 | **PD** | "PD" (accidentally correct) |
| PD747 | GSM5380943 | **PD** | "PD" (accidentally correct) |
| PDC05 | GSM5380944 | **Control** | "PD" ✗ |
| PDC22 | GSM5380945 | **Control** | "PD" ✗ |
| PDC34 | GSM5380946 | **Control** | "PD" ✗ |
| PDC87 | GSM5380947 | **Control** | "PD" ✗ |
| PDC91 | GSM5380948 | **Control** | "PD" ✗ |

Summary: Control n=7 (C36, C48, PDC05, PDC22, PDC34, PDC87, PDC91) | PD n=7 | PDD n=7 | DLB n=7

**Script error in `01_load_and_qc_snnucseq.R`, line 80:**
```r
seurat_obj$condition <- ifelse(grepl("^C", sample_name), "Control", "PD")
```
This commits two simultaneous errors:
1. **PDC05–PDC91 (5 controls) are misclassified as "PD"** (they don't start with "C")
2. **All 21 PD*, PDD*, DLB* samples are collapsed into a single "PD" group** — the study's three disease subgroups are obliterated

**Downstream consequence:** `receptor_expression_PD_vs_Control.csv` compares [PD+PDD+DLB+5wrongly-labeled-Controls] vs [C36+C48 only]. This is:
- Wrong numerically (5 controls in wrong group)
- Scientifically meaningless (PD vs DLB have very different receptor profiles)
- The source of the Figure 4 "PD vs Control" label

**Fix (Phase C):** Replace line 80 of `01_load_and_qc_snnucseq.R` with a lookup table derived from the verified GEO metadata:

```r
# CORRECT — replace the broken ifelse
sample_group_map <- c(
  C36="Control", C48="Control",
  PDC05="Control", PDC22="Control", PDC34="Control", PDC87="Control", PDC91="Control",
  PD060="DLB", PD115="DLB", PD163="DLB", PD294="DLB", PD332="DLB", PD566="DLB", PD706="DLB",
  PD341="PDD", PD366="PDD", PD415="PDD", PD501="PDD", PD531="PDD", PD563="PDD", PD678="PDD",
  PD413="PD",  PD416="PD",  PD523="PD",  PD666="PD",  PD683="PD",  PD732="PD",  PD747="PD"
)
seurat_obj$condition <- sample_group_map[sample_name]
```

After fixing, the primary LBD comparison must be **DLB vs Control** (n=7 vs n=7). Secondary: PDD vs Control. PD vs Control is exploratory (shows disease specificity). Never collapse all three disease groups into one label. The current "PD vs Control" in all output files and figure titles is wrong for 14 of 21 disease samples.

---

**D2 — BCAM vs BCL3 contradiction (confirmed from data)**

MAGMA output (`LBD_risk_genes_ranked.tsv`): **BCAM = rank #1** (Z=7.09, p=6.6e-13). BCL3 (Entrez 602) is absent from all 34 analyzed genes. The Abstract/Results text (BCAM) is consistent with the actual MAGMA output. The Discussion text (BCL3) is factually wrong and must be corrected.

Caveats that must accompany BCAM claim:
- BCAM (1 SNP, chr19:44,809,059–44,821,421) shares the same LD block as APOE, APOC1, TOMM40, APOC4 — all at Z≈6.1 with BCAM at Z=7.09 from a single SNP
- The elevated BCAM score almost certainly reflects the chr19 APOE-locus signal
- BCL3 (chr19q13, Entrez 602) is not in the output because no SNPs near BCL3 passed p≤1e-5; this says nothing about its biology, only that the pre-filtering at 1e-5 before MAGMA is wrong (see D6)

---

**D3 — CXCL12–CXCR4 axis absent from pipeline output (confirmed)**

- `expanded_CD4_upregulated_genes.txt` (642 genes checked): **CXCL12 is absent**
- `active_ligand_receptor_pairs.csv` (16 pairs): **CXCR4 is absent**
- The Phase 3 L-R database includes CXCL12→CXCR4 as a defined pair (line ~120 of `04_ligand_receptor_analysis.R`) but it never activates because CXCL12 is not in the CD4 upregulated gene list
- Manuscript claims this axis is central to Figure 4: **FALSE** — Figure 4 shows OSM/LIFR, TGFB1/TGFBRs, IL17A/IL17RA/RC, etc.

Note: Gate et al. (Science 2021) DID report CXCR4 upregulation in CSF T cells. The absence here may reflect (a) the stimulated vs unstimulated confound, (b) the permissive ≥2-cell expansion threshold capturing too many non-specific clones, or (c) CXCL12 being expressed on neurons (ligand side) not T cells.

---

**D4 — MHC-II (HLA-DRA/B)–TCR axis absent from pipeline output (confirmed)**

Neither HLA-DRA, HLA-DRB1, nor any TCR chain gene is in the hardcoded `ligand_receptor_db` in `04_ligand_receptor_analysis.R`. This axis cannot appear in any output. The manuscript claim is entirely unsupported by the pipeline.

---

**D5 — "Causal" network is a heuristic prioritization formula**

`Causality_Score = Max_Neuronal_Expression × (1 + max(0, Z))`

This is a weighted sum with no causal inference component. No directionality test, no mediation analysis, no MR, no eQTL colocalization. The script title is "Weighted Interactomics (FIX: Graph Layout)" — even the author recognized it is not causal. Must be renamed throughout.

---

**D6 — MAGMA severely underpowered by wrong pre-filtering (blocking for Phase 1)**

The GWAS was pre-filtered to p≤1e-5 before MAGMA input → 454 SNPs → **only 34 genes analyzed**.

MAGMA's gene-based test is designed to operate on all GWAS SNPs simultaneously, aggregating statistics across SNPs within each gene while accounting for LD. Pre-filtering discards the majority of SNPs that collectively inform gene-level tests. The current analysis simply re-ranks the handful of genome-wide-suggestive loci — it is NOT a genome-wide gene-based test.

From the MAGMA log:
```
read 454 lines from file, containing valid SNP p-values for 426 SNPs in data
37 gene definitions read from file
found 34 genes containing valid SNPs in genotype data
```
A proper run on the full harmonized GWAS would analyze ~18,000 genes.

---

### Methodological / Statistical

**D7 — Pseudo-replication in DE (Phase 2)**
Wilcoxon treats cells as independent. N=2 donors only. padj values (1.7e-45 for TPM4) are unreliable. Pseudobulk per donor required.

**D8 — Stimulation confound in expansion signature**
`03_extract_cd4_signatures.R` models expanded vs non-expanded without `condition` (Stimulated/Unstimulated) as covariate. The cytoskeletal signature (TPM4, ACTB, PFN1, ACTG1, CFL1, VIM) is the classic naïve→effector activation module and will be enriched in Stimulated samples regardless of clonal status.

**D9 — QC ceiling too low (Phase 2)**
`nFeature_max = 2,500` for activated CSF T cells. Activated clones often exceed this. Most biologically relevant cells may be excluded.

**D10 — "Nuclei" vs "cells" terminology**
Phase 2 = scRNA-seq (cells). Phase 3 = snRNA-seq (nuclei). Scripts and manuscript use both terms interchangeably.

**D11 — CDR3 clonality: amino acid only, no nucleotide confirmation**
`cloneCall = "aa"` permits convergent clonotypes (different nucleotide sequences encoding the same amino acid CDR3) to be conflated. Stringent definition requires nucleotide identity.

**D12 — Wrong N in MAGMA**
Script hardcodes N=2,591 with comment "Cases + controls". Metadata confirms N=6,618. N=2,591 is the case count only. This underestimates effective sample size and deflates Z-scores.

**D13 — Reference panel build mismatch (Phase 1)**
g1000_eur = hg19; NCBI38.gene.loc = GRCh38. SNP→gene mapping partially affected.

**D14 — p-value threshold comment vs code mismatch**
Comment says "p_value <= 1e-8 (genome-wide significance)"; code uses `p_threshold = 1e-5`.

**D15 — No seeds**
Only downsampling in Phase 3 script 02 uses `set.seed(123)`. UMAP, clustering, DE, all unseeded.

**D16 — No reproducibility infrastructure**
No `renv.lock`, `requirements.txt`, `environment.yml`, or `sessionInfo()` anywhere.

**D17 — No pipeline orchestrator**
No Makefile/Snakemake/{targets}. Manual execution required in strict order.

**D18 — Hardcoded absolute path + typo (Phase 1)**
`01_filter_gwas_data.py` has `home/ana/Desktop/Decoding_the_neuroinmune_synapse_LDB/...` — missing leading `/`, and "neuroinmune" vs "neuroimmune".

**D19 — GWAS (Chia 2021, GCST90001390) not cited in manuscript**
The primary GWAS source confirmed from metadata is absent from the reference list.

---

## 6. Phase 1 MAGMA Output – Full 34-Gene List

| Rank | Symbol | Entrez | CHR | Z | P | Note |
|------|--------|--------|-----|---|---|------|
| 1 | BCAM | 4059 | 19 | 7.09 | 6.64e-13 | 1 SNP; chr19 LD block |
| 2 | GBA | 2629 | 1 | 7.00 | 1.28e-12 | Known LBD locus |
| 3 | APOC1 | 341 | 19 | 6.22 | 2.56e-10 | chr19 LD block |
| 4 | GON4L | 54856 | 1 | 6.14 | 4.02e-10 | Near GBA |
| 5 | APOE | 5819 | 19 | 6.11 | 5.00e-10 | chr19 LD block |
| 5= | TOMM40 | 10452 | 19 | 6.11 | 5.00e-10 | chr19 LD block |
| 5= | APOC4 | 348 | 19 | 6.11 | 5.00e-10 | chr19 LD block |
| … | … | … | … | … | … | |
| 34 | ? | 1012 | 16 | 4.44 | 4.48e-06 | |

BCL3 (Entrez 602, chr19) = **absent from all 34 genes**.
Ranks 5, 5=, 5= share identical Z/P — one LD block reported as three genes.

---

## 7. Phase 2 DEG Top Hits (Expanded vs Non-expanded CD4+, Wilcoxon, N=2 donors)

| Rank | Gene | log2FC | padj |
|------|------|--------|------|
| 1 | TPM4 | 1.16 | 1.7e-45 |
| 2 | ACTB | 0.85 | 8.3e-44 |
| 3 | PFN1 | 0.75 | 1.7e-40 |
| 4 | ACTG1 | 0.75 | 8.0e-40 |
| Top 30 | CFL1, VIM, CTSH, LGALS1, GAPDH, MYL6, ANXA2, S100A4, … | — | — |

3,671 total DEGs. All top hits are cytoskeletal/metabolic — T cell activation signature.
CXCL12: **absent** from upregulated genes. CXCR4: **absent** from neuronal L-R pairs.

---

## 8. Phase 3 Active L-R Pairs (all 16)

| Ligand (CD4+) | Receptor (Neuron) | Type |
|---------------|-------------------|------|
| IL17A | IL17RA, IL17RC | Inflammatory |
| IL26 | IL10RB | Inflammatory |
| TGFB1 | TGFBR1, TGFBR2, TGFBR3 | Anti-inf./Fibrotic |
| OSM | OSMR, LIFR | Inflam./Survival |
| TNFRSF1B | TNFRSF1A | Inflammatory |
| CD40LG | CD40 | Co-stimulation |
| CD70 | CD27 | Co-stimulation |
| LAG3 | CD86 | Immune checkpoint |
| HAVCR2 | LGALS9 | Immune checkpoint |
| ENTPD1 | ADORA1 | Purinergic |

CXCL12→CXCR4: **absent**. HLA-DRA/B→TCR: **absent** (not in L-R DB).

---

## 9. Dataset Provenance — CORRECTED

| Phase | GEO | Disease | Status |
|-------|-----|---------|--------|
| 1 – GWAS | GCST90001390 | LBD (Chia 2021) | ✓ Correct; uncited in manuscript |
| 2 – CSF scRNA | GSE161192 | PD (Gate 2021, Science 374:868) | ✓ Correct dataset |
| 3 – Cortex snRNA | GSE178146 | LBD spectrum: Control/PD/PDD/DLB (Feleke 2021, Acta Neuropathol 142:449) | ✓ **Correct dataset; WRONG condition labels in script** |

**GSE178146 is the right dataset.** The "PD vs Control" problem is a condition-labeling defect, not a wrong-dataset defect. The script collapses 4 groups into 2 with incorrect boundaries.

---

## 10. Next Steps for Phase B

1. **Obtain exact sample→group mapping** for GSE178146 (C36/C48/PDC05-91 = Controls; which PD* = PD/PDD/DLB). Source: Feleke 2021 supplementary Table 1 or GEO sample-level attributes.
2. Verify each manuscript headline claim against code + output tables.
3. Focus first on D1 (condition labels), D2 (BCAM/BCL3), D3 (CXCL12/CXCR4), D4 (HLA/TCR), D6 (MAGMA pre-filtering).

---

_End of Phase A Inventory. No pipeline files modified._

---

## 11. Addendum (2026-08-19) — GSE141578, GCTA-COJO, and updated runnability

This addendum covers items added to the repository after the original Phase A inventory (2026-07-16) and REPORT.md finalization (2026-07-17), discovered and resolved in a follow-up session. Full narrative and evidence: REPORT.md §11.

### 11.1 New files

```
scripts/phase2_immune_profiling/
  └── 04_GSE141578_disease_vs_HC_pseudobulk.R   ← R/Seurat v5 + DESeq2; disease vs HC pseudobulk, CD4+ T cells
data/
  ├── GSE141578_raw/csf_10x/{22 sample dirs}    ← standard 10x dirs, CSF32 staged in this session
  ├── GSE141578_raw/extracted/                   ← original tarballs incl. 13 unused PBMC samples
  ├── GSE141578_processed/
  │     ├── pseudobulk_DEG_disease_vs_HC.csv     ← 13,575 genes tested, DESeq2 ~batch+disease_group
  │     ├── cd4_disease_upregulated_genes.txt    ← 0 genes at padj<0.05
  │     ├── GSE141578_CSF_merged_seurat.rds
  │     ├── GSE141578_CD4_Tcells_seurat.rds
  │     └── plots/ (umap_overview.pdf, umap_cd4_cells.pdf, volcano_disease_vs_HC.pdf)
  └── gcta_cojo/                                  ← new; GCTA-COJO conditional analysis of chr19 locus
        ├── chr19_APOE_locus.ma                   ← 3,293-SNP summary stats, GCTA .ma format
        ├── chr19_locus_ref.{bed,bim,fam}          ← 1000G EUR, rsID-matched (3,074 SNPs, 503 indiv.)
        ├── chr19_stepwise.jma.cojo                ← 3 independent signals selected
        ├── chr19_cond_on_lead.cma.cojo            ← single-SNP conditional (superseded)
        └── chr19_cond_on_all3.cma.cojo            ← 3-SNP joint conditional (definitive)
tools (new, via conda, not previously in environment.yml):
  gcta64 v1.94.1, plink v1.90b6.21
```

### 11.2 Dependency graph addition

```
Phase 2b (new) ──────────────────────────────────────────────────────────
GSE141578_raw/csf_10x/ (22 samples)
  → [04_GSE141578_disease_vs_HC_pseudobulk.R] → pseudobulk_DEG_disease_vs_HC.csv
                                                  cd4_disease_upregulated_genes.txt
    (independent of Phase 2a/GSE161192; not currently consumed by Phase 3/4)

Phase 1 addendum ─────────────────────────────────────────────────────────
LBD_GWAS_harmonised.tsv.gz (chr19 subset) + data/reference/g1000_eur
  → [manual GCTA-COJO invocation, not yet a Snakemake rule]
  → data/gcta_cojo/chr19_stepwise.jma.cojo, chr19_cond_on_all3.cma.cojo
    (diagnostic; resolves BCAM/APOE independence question from §5, D2)
```

### 11.3 Runnability update

| Script | Status | Notes |
|--------|--------|-------|
| `04_GSE141578_disease_vs_HC_pseudobulk.R` | RUNNABLE (fixed, executed successfully) | Fixed: CSF32 metadata + staging, `JoinLayers()` before `GetAssayData()`. Ran to completion twice (with and without batch covariate). `ggrepel` added to `install_packages.R`. |
| GCTA-COJO conditional analysis | RUNNABLE (executed manually) | Not yet a Snakemake rule — one-off diagnostic. `gcta64`/`plink` now available via conda but not pinned in `environment.yml`. |
| NicheNet (`nichenetr`) install | **BLOCKED** | Missing system library `libcairo2-dev` (or conda-forge `cairo`), required by transitive dependency `gdtools`. Requires root (`sudo apt install libcairo2-dev`) or a longer conda-forge solve than attempted. Not installed; directional L-R inference (REPORT.md §5.7) remains unimplemented. |

### 11.4 New defects found and fixed this session

- **D20 — GSE141578 CSF32 double-fault.** Metadata table wrongly excluded CSF32 as "diagnosis not reported" (GEO explicitly lists `diagnosis: healthy control`); separately, CSF32's raw files were never staged into the `csf_10x/{sample}/` layout the script expects, so it was silently dropped even after the metadata fix. Both fixed. True cohort composition: 11 HC, 9 PD, 2 DLB (22 CSF samples), matching the GEO series-level design description.
- **D21 — `GetAssayData()` fails on multi-layer Seurat v5 object.** Same defect class as the Phase 3 script 4 fix in the original audit; `JoinLayers()` now called before pseudobulk count extraction in the new script.
- **D22 — `ggrepel` missing from `install_packages.R`.** Fixed.
- **D23 — Build-mismatch trap in ad hoc region extraction.** An initial attempt to extract the chr19 APOE-locus region from the `g1000_eur` reference panel by physical-position (BP) range, using GRCh38 coordinates, silently pulled the *wrong* genomic region because the reference panel is GRCh37/hg19 (previously flagged as D13, but not previously demonstrated to cause a concrete downstream error). Caught by checking that known APOE-defining SNPs (rs429358, rs7412) were absent from the BP-filtered extraction; fixed by matching on rsID instead of position. This is a cautionary, reusable finding for any future analysis touching `data/reference/g1000_eur`: **do not filter this file by GRCh38 physical position.**

### 11.5 Manuscript

`MANUSCRIPT_Neuroimmune_Interactions_LBD.md` was rewritten in full; see REPORT.md §11.5 for the complete rationale. All Results-section numbers were re-pulled directly from current output CSV/TSV files during this session (not copied from REPORT.md §3/§4 prose), which surfaced one additional, previously undocumented discrepancy: IL26 no longer reaches significance in the current GSE161192 DEG output (padj=0.077), contradicting REPORT.md §3/§4's "Confirmed" framing of the IL26 correction — REPORT.md has been updated accordingly.
