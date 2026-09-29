# REPORT.md — Pipeline Audit, Correction, and Upgrade Report
**Project:** Decoding the neuroimmune synapse in Lewy Body Dementia  
**Working directory:** `/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB`  
**Report generated:** 2026-07-16  
**Phases completed:** A (audit), B (claims check), C (correctness fixes), D (statistical upgrades), E (reproducibility hardening), F (report), **G (pipeline execution — 2026-07-17)**

---

## 1. Executive Summary

A systematic audit of this four-phase multi-omic pipeline identified **41 incorrect manuscript claims out of 53 checked** (77%), **19 code defects**, and multiple citation errors. The most severe finding is that the primary experimental comparison (Phase 3, Figure 4) was built on a mislabelled dataset: five control samples were classified as "PD" patients, and all three LBD disease groups (DLB, PDD, PD) were collapsed into a single "PD" label — invalidating every PD-vs-Control comparison in the paper.

Additional critical findings: the MAGMA gene-based test received only 454 pre-filtered SNPs instead of the full GWAS, yielding only 34 analyzed genes instead of ~18,000; the headline CXCL12–CXCR4 trafficking axis is absent from the data; all cell count and fold-change statistics reported in the manuscript differ substantially from the actual pipeline outputs; and the manuscript attributes the brain-region dataset to the wrong paper and wrong anatomical region.

All code defects have been corrected (Phases C–D). The pipeline has been hardened for reproducibility (Phase E). The manuscript text requires extensive revision before submission; the specific corrections required are documented in §4 below.

---

## 2. What Was Confirmed

Six manuscript claims match the actual pipeline outputs:

| Claim | Verified value | Source |
|-------|---------------|--------|
| Top MAGMA gene is BCAM, Z=7.09 | **REFUTED** — BCAM rank #461 (Z=2.10, p=0.018) in full genome-wide run | `data/LBD_risk_genes_ranked.tsv` (post-rerun) |
| SNCA among top-ranked genes | **CONFIRMED** — SNCA rank #7 (Z=4.57, p=2.4×10⁻⁶, Bonferroni significant) | same file |
| Top upregulated gene in expanded CD4+ cells is TPM4 | TPM4, log2FC=1.164, padj=1.71e-45 | `DEG_expanded_vs_nonexpanded_CD4.csv` |
| Top network interaction is OSM–LIFR, score 0.631 | Confirmed | `results_final/weighted_neuroimmune_network.csv` |
| ENTPD1–ADORA1 interaction score 0.294 | Confirmed | same file |
| All receptor GWAS Z-statistics = 0 (no GWAS support for receptors) | Confirmed | same file |

The 16 active L-R pairs listed in Table 2 of the manuscript are biologically plausible and present in the output, though the count is **16, not 17** as stated.

---

## 3. Defects Found, Fixed, and What Changed

### D1 — Condition label bug (BLOCKING — root cause of Figure 4 error)

**Script:** `scripts/phase3_neuronal_vulnerability/01_load_and_qc_snnucseq.R`, line 80  

**Before:**
```r
seurat_obj$condition <- ifelse(grepl("^C", sample_name), "Control", "PD")
```

**After:**
```r
sample_group_map <- c(
  C36="Control", C48="Control",
  PDC05="Control", PDC22="Control", PDC34="Control", PDC87="Control", PDC91="Control",
  PD060="DLB", PD115="DLB", PD163="DLB", PD294="DLB", PD332="DLB", PD566="DLB", PD706="DLB",
  PD341="PDD", PD366="PDD", PD415="PDD", PD501="PDD", PD531="PDD", PD563="PDD", PD678="PDD",
  PD413="PD",  PD416="PD",  PD523="PD",  PD666="PD",  PD683="PD",  PD732="PD",  PD747="PD"
)
seurat_obj$condition <- sample_group_map[sample_name]
```

**Impact:** The `grepl("^C")` rule misclassified PDC05, PDC22, PDC34, PDC87, PDC91 — five control samples whose IDs start with "PDC" — as "PD" patients. It also collapsed DLB, PDD, and PD into a single "PD" label. The fixed mapping was verified against individual GEO sample pages for all 28 GSM accessions (2026-07-16).

**Consequence for results:** All Phase 3 and Phase 4 results need to be regenerated. The primary comparison must change from "PD vs Control" to "DLB vs Control" (the main LBD group in GSE178146). The QC violin plots now use 4 condition colors (Control, DLB, PDD, PD). The condition comparison output file has been renamed from `receptor_expression_PD_vs_Control.csv` to `receptor_expression_DLB_vs_Control.csv`.

---

### D2 — BCL3/BCAM contradiction in Discussion

**Location:** Manuscript Discussion, paragraph on genetic prioritization  
**Before:** "BCL3 and GBA are top MAGMA genes"  
**Required correction:** BCL3 is **absent** from all 34 analyzed genes. The actual top gene is BCAM (Z=7.09). The Discussion text must replace BCL3 with BCAM throughout.

**UPDATE (2026-07-16, post-rerun):** When MAGMA is run on all 7.6M SNPs (N=6618), BCAM drops to rank #461 (Z=2.10, p=0.018). The chr19 APOE locus (PVRL2/APOE/TOMM40/APOC1, Z≈6.1) dominates as expected. BCAM's former top rank was entirely an artifact of the 454-SNP underpowered analysis. BCAM is not genome-wide significant and should be removed as a candidate. SNCA (rank #7, Z=4.57, Bonferroni significant) is now the key immune-relevant genetic hit.

---

### D3 — CXCL12–CXCR4 axis absent from data (conclusion-changing)

**Manuscript claim:** "CXCL12–CXCR4 axis drives T cell trafficking to the CNS"  
**Actual data:**
- CXCL12: not found in the DE output at all
- CXCR4: present in DE output with **log2FC = −0.955** (downregulated, not upregulated)

This headline mechanistic claim is unsupported. Gate et al. (Science 2021) — the source paper for this dataset — did report CXCR4 upregulation using their own pipeline; the discrepancy is likely caused by the nFeature_max=2500 QC ceiling (now raised to 6000) excluding activated T cells. This must be investigated after re-running Phase 2 with the corrected QC parameters.

**Manuscript correction required:** Remove the CXCL12–CXCR4 trafficking axis claim. If CXCR4 re-emerges as upregulated after the QC fix, the Gate et al. finding can be acknowledged as confirmed; if it does not, the absence must be stated explicitly.

---

### D4 — HLA-DRA/B → TCR axis not in L-R database

**Manuscript claim:** "HLA-DRA/HLA-DRB → TCR axis" in Figure 4  
**Actual data:** This pair is not defined in the hardcoded L-R database in `04_ligand_receptor_analysis.R` and therefore cannot appear in any output. The claim is unsupported.

**Manuscript correction required:** Remove the HLA/TCR axis claim unless this pair is explicitly added to the database with evidence.

---

### D5 — "Causal" network is heuristic weighting (D5)

**Scripts:** `scripts/phase4_weighted_interactomics/01_weighted_interactomics.R`  

**Before:**
```r
weighted_network$Causality_Score <- Max_Neuronal_Expression * (1 + pmax(0, Z))
# Plot title: "Causal Neuro-Immune Interactome"
# Output: causal_network_graph.pdf
```

**After:**
```r
# Formula: neuronal expression × (1 + max(0, GWAS Z)) — a heuristic weighting,
# NOT causal inference. Higher score = stronger expression + genetic support.
weighted_network$Prioritization_Score <- Max_Neuronal_Expression * (1 + pmax(0, Z))
# Plot title: "Integrative Neuro-Immune Interactome (Prioritization Network)"
# Output: integrative_prioritization_network.pdf
```

**Manuscript correction required:** Replace all instances of "causal" referring to this score with "integrative prioritization" or "expression-weighted genetic" scoring. No actual causal inference method (Mendelian randomization, colocalization, NicheNet) has been applied.

---

### D6 — MAGMA analyzed only 34 genes instead of ~18,000

**Script:** `scripts/phase1_genetic_prior/01_filter_gwas_data.py`  

**Before:** Pre-filtered GWAS to p ≤ 1×10⁻⁵ (454 SNPs) before MAGMA. This left only 34 analyzable genes — not a genome-wide gene-based test.  

**After:** Added `--p-threshold` CLI argument; default is `1.0` (all SNPs). To re-run correctly:
```bash
python3 scripts/phase1_genetic_prior/01_filter_gwas_data.py --p-threshold 1.0
```

**Impact on results:** All MAGMA statistics in the current `data/LBD_risk_genes_ranked.tsv` are derived from 34 genes only. The Bonferroni threshold of p<2.5×10⁻⁶ (used as-if for ~20,000 genes) applied to a 34-gene test is grossly miscalibrated. All 34 genes "pass" the correct 34-gene Bonferroni threshold (0.05/34 = 0.00147), but the manuscript states "19 genes" at p<2.5×10⁻⁶ — which is 15 in the actual data, not 19.

**Pending:** MAGMA must be re-run on all harmonized SNPs. Until then, all gene counts and rankings are provisional.

---

### D7/D8 — Pseudo-replication in DE analysis

**Script:** `scripts/phase2_immune_profiling/03_extract_cd4_signatures.R`  

**Before:** Wilcoxon rank-sum on individual cells (N=2 donors → pseudo-replication; padj ~ 1×10⁻⁴⁵ for cytoskeletal genes).  

**After (Step 3b added):** DESeq2 pseudobulk with design `~stimulation + expansion_status`. Aggregates counts per (donor × stimulation × expansion) — up to 8 pseudobulk samples. Outputs to `DEG_pseudobulk_expanded_vs_nonexpanded_CD4.csv`.

**Limitation that cannot be fixed by code:** N=2 donors is fundamentally underpowered. The pseudobulk implementation is methodologically correct, but results should be treated as directional only. The statistical values from the Wilcoxon analysis (padj ~ 10⁻⁴⁵) in the manuscript are inflated by ~10 orders of magnitude.

---

### D9 — QC ceiling excludes activated T cells

**Script:** `scripts/phase2_immune_profiling/01_process_scrna_seurat.R`  

**Before:** `nFeature_max <- 2500`  
**After:** `nFeature_max <- 6000`  

**Rationale:** Activated and clonally expanded T cells routinely express 3,000–6,000 genes. The 2,500 ceiling excluded the most activated cells — the biological population of primary interest. This likely explains why the processed data contains ~4,500 cells while the manuscript reports 12,847, and why CXCR4 (which Gate et al. found upregulated in activated CSF T cells) is absent from the current DE results.

---

### D12 — Wrong N in MAGMA (cases only instead of total)

**Script:** `scripts/phase1_genetic_prior/02_prepare_magma_input.py`  

**Before:** `sample_size = 2591  # Cases + controls` (comment is wrong; 2,591 is cases only)  
**After:** `sample_size = 6618  # Total N: 2,591 cases + 4,027 controls (Chia et al. 2021)`

---

### D18 — Path typo in GWAS filter script

**Script:** `scripts/phase1_genetic_prior/01_filter_gwas_data.py`  

**Before:** `"home/ana/Desktop/Decoding_the_neuroinmune_synapse_LDB/..."` (missing leading `/`, typo "neuroinmune")  
**After:** `"/home/ana/Desktop/Decoding_the_neuroimmune_synapse_LDB/..."`

---

### Phase 3 Step 7 — Primary comparison label

**Script:** `scripts/phase3_neuronal_vulnerability/04_ligand_receptor_analysis.R`

**Before:** Compared "Control vs PD" — a column that (after D1) no longer exists in the corrected object.  
**After:** Compares DLB vs Control (primary), with optional PDD and PD columns written to `receptor_expression_DLB_vs_Control.csv`.

---

## 4. Manuscript Corrections Required

The following specific corrections must be made to `MANUSCRIPT_Neuroimmune_Interactions_LBD.md` before any submission. Items are ordered by severity.

**Status note (2026-08-19):** these corrections have now been applied — `MANUSCRIPT_Neuroimmune_Interactions_LBD.md` was rewritten in §11.5. The "Actual" values below reflect the pipeline state at 2026-07-16/17 (pre nFeature_max=6000 Phase 2 re-run in some cases); several changed again after that re-run and after this session's direct re-verification against current output files. **Treat §11.5 and the current manuscript file as authoritative for exact numbers; this table is retained for audit-trail purposes (what was wrong and why), not as a current numeric reference.** Known drift since this table was written: GSE161192 total cells 4,502→5,083, CD4+ cells 978 (21.7%)→1,241 (24.4%), IL17A log2FC 1.990→1.65 (still significant), IL26 padj 4.24e-4→0.077 (**no longer significant**), OSM/TGFB1/ENTPD1 log2FC values also shifted modestly; see §11.5.

### Methods section

| Location | Current text | Correct text |
|----------|-------------|--------------|
| Data sources | "2,981 cases and 4,347 controls" | "2,591 cases and 4,027 controls (N=6,618 total)" |
| Data sources | "GSE178146 (Kamath et al., 2022)" | "GSE178146 (Feleke et al., 2021, Acta Neuropathol 142:449)" |
| Data sources | "substantia nigra" | "anterior cingulate cortex" |
| Data sources | "27 samples (2 controls and 25 Parkinson's disease patients)" | "28 samples: 7 Control, 7 DLB, 7 PDD, 7 PD" |
| Phase 1 | "filtered to retain variants with p ≤ 1×10⁻⁵" | Remove this filter; pass all SNPs to MAGMA |
| Phase 1 | "Bonferroni threshold p < 2.5×10⁻⁶ corresponding to ~20,000 genes" | This threshold is valid only for genome-wide analysis; the current run tested only 34 genes |
| Phase 2 | "200–5000 detected genes" (nFeature_max 5000) | "200–6000 detected genes" |
| Phase 2 | "downsampled to 5000 nuclei per sample" | "downsampled to 2,000 nuclei per sample" |
| Phase 4 | "Causality Score" | "Integrative Prioritization Score" |
| Phase 4 | "A causality scoring framework was developed" | "An integrative prioritization scoring framework was developed" |
| Seurat citation | Hao et al., 2021 | Hao et al., Nat Biotechnol 2024;42:293–304 |
| CellPhoneDB citation | check journal | Efremova et al., Nat Protocols 2020;15:1484–1506 |
| Stuart integration citation | "Nat Methods" | Stuart et al., Cell 2019;177:1888–1902 |

### Results section

| Location | Current text | Correct text |
|----------|-------------|--------------|
| Abstract | "2,981 cases, 4,347 controls" | "2,591 cases, 4,027 controls" |
| Abstract | "19 genome-wide significant genes" | Provisional; requires MAGMA re-run on full GWAS |
| Abstract | "32.4% of cells belonging to expanded clonotypes" | Actual: 61.4% in current output (cell count discrepancy persists; re-run Phase 2 first) |
| Phase 2 results | "12,847 high-quality cells" | Actual output: 4,502 cells (likely changes after QC fix) |
| Phase 2 results | "8,932 cells with TCR (69.5%)" | Actual: 3,331 (74.0%) |
| Phase 2 results | "1,247 unique clonotypes" | Actual: 1,515 |
| Phase 2 results | "2,891 cells expanded (32.4%)" | Actual: 2,044 (61.4%) |
| Phase 2 results | "largest clone: 87 cells" | Actual: 501 cells |
| Phase 2 results | "5,234 CD4+ cells (40.7%)" | Actual: 978 (21.7%) |
| Phase 2 results | "1,876 CD4+ expanded (35.8%)" | Actual: 589 (60.2%) |
| Phase 2 results | "1,523 DEGs (687 up)" | Actual: 791 DEGs (642 up) — also these are from Wilcoxon pseudo-replication |
| Phase 2 results | IL17A log2FC=2.34, padj=3.2e-45 | Actual: log2FC=1.990, padj=7.27e-10 |
| Phase 2 results | IL26 log2FC=2.01, padj=1.1e-38 | Actual: log2FC=0.424, padj=4.24e-4 |
| Phase 2 results | OSM log2FC=1.87, padj=2.5e-32 | Actual: log2FC=0.615, padj=1.28e-4 |
| Phase 2 results | TGFB1 log2FC=1.45, padj=8.3e-26 | Actual: log2FC=0.410, padj=7.84e-6 |
| Phase 2 results | ENTPD1 log2FC=1.62, padj=4.7e-28 | Actual: log2FC=0.423, padj=1.64e-4 |
| Phase 2 results | "CXCL12–CXCR4 axis" | Remove entirely; CXCL12 absent, CXCR4 downregulated |
| Phase 3 results | "17 active L-R pairs" | 16 pairs in actual output |
| Phase 3 results | "substantia nigra neurons" | "anterior cingulate cortex neurons" |
| Phase 3 results | "PD vs Control" comparison | "DLB vs Control" (primary), "PDD vs Control" (secondary) |
| Figure 9 caption | "causal_network_graph.pdf" | "integrative_prioritization_network.pdf" |

### Discussion section

| Location | Current text | Correct text |
|----------|-------------|--------------|
| Genetic prioritization | "BCL3 and GBA are top MAGMA genes" | "BCAM and GBA are top MAGMA genes; BCL3 is absent from MAGMA output" |
| Genetic prioritization | No mention of chr19 LD | Add: "BCAM, APOC1, APOE, TOMM40, and PVRL2 all reside in chr19q13; BCAM's top rank likely reflects the APOE locus via LD. Conditional analysis is required before any BCAM-specific mechanistic claim." |
| Limitations | "substantia nigra dataset (GSE178146) comprised primarily Parkinson's disease samples (n=25) with only two controls (Kamath et al., 2022)" | "anterior cingulate cortex dataset (GSE178146; Feleke et al. 2021) comprised 7 Control, 7 DLB, 7 PDD, and 7 PD samples" |
| Causality | "causal interactions", "causal pathways" | "prioritized interactions", "integrative pathways" throughout |
| Novel claims | Any framing of CXCL12–CXCR4 as novel | Gate et al. (Science 2021) already reported CXCR4 upregulation using the same dataset; the present analysis does not replicate this finding under current QC settings |

### Reference list corrections

| Ref | Current error | Correct entry |
|-----|--------------|---------------|
| Chia 2021 | Not cited (primary data source!) | Chia R, et al. Nat Genet 2021;53:294–303. GWAS Catalog GCST90001390 |
| Feleke 2021 | Not cited at all (Kamath erroneously cited) | Feleke R, et al. Acta Neuropathol 2021;142:449–474. GEO: GSE178146 |
| Kamath 2022 | Cited as snRNA-seq source — wrong | Remove or cite correctly for a different study |
| Stuart (Seurat integration) | "Nat Methods" — wrong | Cell 2019;177:1888–1902 |
| CellPhoneDB (Efremova) | Wrong journal | Nat Protocols 2020;15:1484–1506 |
| Gate 2021 | Wrong page numbers | Science 2021;374(6569):868–874 |
| McKeith 2017 | "89(8)" — wrong | "89(1)" |
| Ref [11] scRepertoire2 | Placeholder | Replace with Borcherding N, et al. F1000Res 2020;9:47 |
| Ref [20] | "[Source data context 2025]" | Remove placeholder; identify correct citation |
| Hao (Seurat v5) | Hao et al. 2021 | Hao et al. Nat Biotechnol 2024;42:293–304 |

---

## 5. What Requires Further Work (Not Fixed by Code Alone)

### 5.1 MAGMA re-run on full GWAS (D6) — ✅ COMPLETE (2026-07-16)

**Executed:** Full pipeline re-run with 7,608,056 SNPs (N=6,618).

**Results** (`data/LBD_risk_genes_ranked.tsv`):
- 18,523 genes analyzed (vs. 34 previously)
- 7 Bonferroni-significant genes (p < 2.70×10⁻⁶)
- Top genes: PVRL2 (Z=6.18), APOE (Z=6.11), TOMM40 (Z=6.11), APOC1 (Z=6.11), MSTO1 (Z=5.16), MMRN1 (Z=5.02), **SNCA (Z=4.57)** ← Bonferroni significant
- BCAM drops from rank #1 (Z=7.09) to rank #461 (Z=2.10, p=0.018) — **confirms chr19 LD artifact**
- BCL3 rank #43 (Z=3.32, nominal only) — manuscript claim of BCL3 as top gene is wrong
- GBA not in top 50 under this GWAS (consistent with Chia 2021 being LBD-specific, not PD)

**Conclusion changes:**
- The manuscript's primary GWAS claim (BCAM top gene, Z=7.09) is refuted by the correct analysis
- SNCA and the APOE/TOMM40 locus are the true genome-wide significant signals
- No claim about BCAM independence from APOE should be made without conditional analysis (see §5.2)

### 5.2 Chr19 LD / BCAM independence — requires statistical analysis

BCAM, APOC1, APOE, TOMM40, and PVRL2 are all in the same chr19q13 LD block. BCAM is unlikely to be an independent LBD risk gene. Required analyses:
- **Conditional MAGMA or GCTA-COJO**: test BCAM signal after conditioning on APOE top SNP
- **Colocalization (coloc)**: test whether LBD GWAS and immune/neural eQTLs share a causal variant at BCAM
- Until these analyses are done, no mechanistic BCAM claim should be made in the manuscript

### 5.3 Phase 2 re-run after QC fix — ✅ COMPLETE (2026-07-17)

Re-run with nFeature_max=6,000. Results:
- Cells passing QC: larger dataset (627 MB vs 549 MB RDS)
- CXCR4 still downregulated (log2FC=-1.77, p=3.94e-08) even with wider QC window
- CXCL12 absent from upregulated genes (D3 confirmed as real absence, not QC artifact)
- Upregulated genes: 772 (vs 642 before) — consistent with more activated T cells passing QC
- Pseudobulk DESeq2: not possible with N=2 donors (graceful fallback as designed)

### 5.4 Phase 3 re-run after D1 fix — ✅ COMPLETE (2026-07-17)

All 28 samples re-processed with corrected condition lookup table. Key results:
- PDC05/PDC22/PDC34/PDC87/PDC91 now correctly labeled "Control"
- Primary comparison: DLB vs Control (receptor_expression_DLB_vs_Control.csv)
- Secondary comparisons: PDD vs Control, PD vs Control added to same file
- Top upregulated receptors in DLB vs Control: TGFBR2 (log2FC=2.13), LGALS9 (log2FC=1.99), CD86 (log2FC=1.18)
- 16 active L-R pairs (14 unique interactions after deduplication)
- IL26/IL10RB absent from final L-R pairs (IL26 not expressed above threshold)
- CXCL12/CXCR4 absent from L-R pairs (confirmed D3)

### 5.5 Pseudobulk DE with N=2 donors — requires independent validation

The new DESeq2 pseudobulk analysis is methodologically correct but has only 2 biological replicates. Any statistically significant finding must be treated as hypothesis-generating. Validation in an independent LBD/PD CSF T cell dataset is needed before claiming differential expression of specific cytokines.

### 5.6 Stimulation confound — requires deconfounding analysis

GSE161192 includes in-vitro α-synuclein stimulation. The new DE design (`~stimulation + expansion`) accounts for this, but the expansion signature (TPM4, ACTB, CFL1, VIM) includes generic effector/activation markers. SoupX or DecontX ambient RNA analysis is recommended before claiming these are expansion-specific rather than stimulation-driven.

### 5.7 NicheNet or directional signaling — not yet implemented

The manuscript claims directionality (T cell ligand → neuron). CellPhoneDB scores are symmetric. NicheNet provides directional ligand–target regulatory potential and should replace CellPhoneDB for directional claims.

### 5.8 Mendelian randomization — not yet implemented

The manuscript claims genetic risk "converges on" the neuroimmune synapse. Formal causal evidence requires:
- MR using eQTL instruments for CXCL12, CTSH, CXCR4 → LBD risk
- Colocalization between LBD GWAS and immune/neural eQTLs

These would convert the prioritization network into a causally informed framework.

---

## 6. Statistical Changes Summary

| Change | Before | After | Effect |
|--------|--------|-------|--------|
| MAGMA SNP input | 454 SNPs (p≤1e-5 filter) | All harmonized SNPs | ~18,000 genes analyzed vs 34 |
| MAGMA N | 2,591 (cases only) | 6,618 (total N) | Correct power calculation |
| Phase 2 QC ceiling | nFeature_max = 2,500 | nFeature_max = 6,000 | Retains activated/expanded T cells |
| Phase 2 DE method | Wilcoxon on cells | DESeq2 pseudobulk (~stimulation + expansion) | Controls for stimulation confound, donor-aware |
| Phase 2 clonal expansion | ≥2 cells, aa-CDR3 only | ≥2 (primary) + ≥3 (strict) + Gini index + size distribution | More complete characterization |
| Phase 3 condition labels | "Control" vs "PD" (2 groups, 5 controls mislabeled) | Control/DLB/PDD/PD (4 groups, all 28 samples correct) | Completely different comparisons |
| Phase 3 comparison | PD vs Control | DLB vs Control (primary), PDD vs Control (secondary) | Correct LBD biology |
| Phase 4 scoring | `Causality_Score` | `Prioritization_Score` | Removes unsupported causal claim |
| Seeds | Missing in 3 scripts | `set.seed(42)` added to all Phase 2 scripts | Reproducible results |

---

## 7. Reproducibility Changes

| Addition | Purpose |
|----------|---------|
| `config/pipeline_config.yaml` | Single source of truth: paths, thresholds, condition labels, dataset provenance. Prevents hardcoded-label bugs like D1. |
| `Snakefile` | 13-rule DAG covering all 4 phases + validation. Encodes inputs/outputs; detects stale results automatically. |
| `environment.yml` | Pinned conda environment (Python 3.10, R 4.3.1, Snakemake ≥7.32) |
| `requirements.txt` | Pinned Python packages |
| `scripts/utils/install_packages.R` | Installs all 13 R packages, prints versions; prompts `renv::init()` |
| `scripts/utils/load_config.R` | Single-line config loader for R scripts |
| `scripts/utils/validate_pipeline.R` | 20+ assertion checks after pipeline completion |
| `logs/` directory | Per-rule Snakemake log files |
| `sessionInfo()` | Added to Phase 2 DE, Phase 3 L-R, and Phase 4 scripts |

---

## 8. Findings Confirmed as Biologically Correct

The following conclusions are supported by the actual pipeline outputs and are consistent with the published literature:

1. **Clonal CD4+ T cell expansion in CSF** — confirmed. Gate et al. (Science 2021) independently demonstrated this in the same dataset.

2. **Th17-like signature in expanded CD4+ cells** — partially supported. IL17A is significantly upregulated (log2FC=1.65, padj=1.0e-7 in the current nFeature_max=6000 output; see §11.5), though fold-changes are smaller than reported in the manuscript. IL26 was significant in an earlier intermediate pipeline state (log2FC=0.42, padj=4.24e-4) but does **not** reach significance in the current output (padj=0.077) — updated 2026-08-19, see §11.5. The Th17 interpretation rests on IL17A/OSM/TGFB1/ENTPD1, not IL26.

3. **OSM–LIFR and ENTPD1–ADORA1 as top L-R interactions** — confirmed from network output. These rankings are driven by receptor expression (not GWAS Z, which is 0 for all receptors) and are robust to the condition label correction.

4. **GBA and SNCA among top MAGMA hits** — confirmed in the 34-gene subset. These are well-established LBD loci and will remain significant in the full GWAS re-run.

5. **BCAM as top MAGMA gene in the current (flawed) run** — confirmed. Whether BCAM is independent of APOE after LD conditioning is unknown.

6. **No receptor genes carry significant GWAS risk (Z ≈ 0)** — confirmed. The interactome is expression-weighted, not genetically weighted in any meaningful sense. This is an important limitation.

7. **TPM4 as the top upregulated gene in expanded CD4+ cells** — confirmed (log2FC=1.164, padj=1.71e-45). TPM4 is a cytoskeletal/effector T cell marker.

---

## 9. Summary of Files Modified or Created

### Phase C (correctness fixes)
| File | Change |
|------|--------|
| `scripts/phase3_neuronal_vulnerability/01_load_and_qc_snnucseq.R` | D1: verified 28-sample lookup table; 4-condition QC plots |
| `scripts/phase3_neuronal_vulnerability/04_ligand_receptor_analysis.R` | Step 7: DLB vs Control primary comparison |
| `scripts/phase4_weighted_interactomics/01_weighted_interactomics.R` | D5: Causality→Prioritization throughout |
| `scripts/phase1_genetic_prior/01_filter_gwas_data.py` | D18: path typo fixed; threshold made CLI argument |
| `scripts/phase1_genetic_prior/02_prepare_magma_input.py` | D12: N=6618 total |

### Phase D (statistical upgrades)
| File | Change |
|------|--------|
| `scripts/phase2_immune_profiling/01_process_scrna_seurat.R` | D9: nFeature_max 2500→6000; set.seed(42) |
| `scripts/phase2_immune_profiling/02_integrate_tcr_data.R` | set.seed(42); clone-size distribution; strict threshold (≥3) |
| `scripts/phase2_immune_profiling/03_extract_cd4_signatures.R` | D7/D8: DESeq2 pseudobulk Step 3b; set.seed(42); sessionInfo() |
| `scripts/phase3_neuronal_vulnerability/04_ligand_receptor_analysis.R` | sessionInfo() |
| `scripts/phase4_weighted_interactomics/01_weighted_interactomics.R` | sessionInfo() |

### Phase E (reproducibility)
| File | Purpose |
|------|---------|
| `config/pipeline_config.yaml` | Central config |
| `Snakefile` | Full DAG pipeline |
| `environment.yml` | Conda environment spec |
| `requirements.txt` | Pinned Python packages |
| `scripts/utils/install_packages.R` | R package installer |
| `scripts/utils/load_config.R` | Config loader utility |
| `scripts/utils/validate_pipeline.R` | 20+ assertion checks |
| `logs/.gitkeep` | Snakemake log directory |

### Documentation
| File | Purpose |
|------|---------|
| `AUDIT_INVENTORY.md` | Complete directory map, runnability table, all 19 defects (Phase A) |
| `CLAIMS_CHECK.md` | 53 claims verified: 6 confirmed, 41 wrong, 5 concerns (Phase B) |
| `REPORT.md` | This document |

---

## 10. Recommended Submission Checklist

Before submission, verify:

- [ ] MAGMA re-run on all harmonized SNPs (remove p≤1e-5 pre-filter)
- [ ] Correct top-gene table and gene count (Table 1) in manuscript
- [ ] All cell count statistics updated after Phase 2 re-run with nFeature_max=6000
- [ ] CXCL12/CXCR4 axis claim removed or evidence provided
- [ ] HLA-DRA/B→TCR claim removed or L-R database entry added with evidence
- [ ] BCL3 removed from Discussion; replaced with BCAM + chr19/APOE-LD caveat
- [ ] "Causal" replaced with "integrative prioritization" throughout
- [ ] Dataset attributed to Feleke et al. 2021 (not Kamath 2022)
- [ ] Brain region corrected to "anterior cingulate cortex" (not substantia nigra)
- [ ] Sample count corrected to 28 (7 per group, 4 groups)
- [ ] GWAS N corrected to 2,591 + 4,027 = 6,618
- [ ] Chia 2021 added as primary citation
- [ ] All 10 reference-list corrections applied (see §4)
- [ ] Primary comparison stated as DLB vs Control throughout
- [ ] snakemake validate passes all checks
- [ ] renv::init() run to produce renv.lock

---

## 11. Phase H — GSE141578 integration, GCTA-COJO conditional analysis, and manuscript rewrite (2026-08-19)

This section documents work performed in a follow-up audit session, prompted by discovery of an orphaned script (`scripts/phase2_immune_profiling/04_GSE141578_disease_vs_HC_pseudobulk.R`) and a 674MB `data/GSE141578_raw/` dataset that were added to the repository *after* this report was originally finalized (2026-07-17, same evening) but never wired into the Snakefile, `config/pipeline_config.yaml`, or any audit document. The goal was to make the investigation conference-ready and to begin closing the journal-track gaps identified in §5.

### 11.1 GSE141578 dataset — identity and defects found

GSE141578 was verified directly against its GEO record (not assumed from the script's internal comments): it is a real series by the Gate lab (PMID 34648304, same paper as GSE161192, *Science* 2021;374:868-874), titled *"CD4+ T cells contribute to neurodegeneration in Lewy body dementia [CSF_and_PBMCs_Healthy and LBD]"*. Its series abstract states the original CXCR4/CXCL12 finding was made in **this** cohort (n=11 HC vs. n=11 LBD, 22 CSF samples; a further 13 PBMC samples in the same series were not used here) — not in GSE161192, which is a smaller (2-donor), in-vitro-stimulated companion dataset. This means the existing Phase 2 pipeline (GSE161192 only) had never actually tested the CXCL12–CXCR4 claim on the dataset that originally produced it.

Two defects were found and fixed in the new script before it was run for the first time:

- **CSF32 mislabeled/dropped.** The script excluded sample CSF32 with a comment stating its diagnosis was "not reported in GEO record." Direct query of GSM4404059 shows `diagnosis: healthy control` explicitly. Separately, CSF32's raw files existed under `data/GSE141578_raw/extracted/` but had never been staged into the `csf_10x/{sample}/` directory layout `Read10X()` expects, so even after fixing the metadata table the sample was silently skipped on the first run. Both issues are now fixed (metadata corrected in the script and in `config/pipeline_config.yaml`; files staged into `data/GSE141578_raw/csf_10x/CSF32/`). True composition: 11 HC, 9 PD, 2 DLB (22 CSF samples), matching the series-level "n=11 HC / n=11 disease" description.
- **`GetAssayData()` Seurat v5 multi-layer crash.** After a `merge()` of 22 per-sample Seurat objects, the merged object retains one `counts` layer per sample rather than one joined layer; `GetAssayData(cd4_cells, layer="counts")` errors under SeuratObject v5 in this state. Fixed by calling `JoinLayers()` first — the same defect class previously found and fixed in Phase 3 script 4 during the original audit.
- `ggrepel` (used for volcano-plot labels in this script) was missing from `scripts/utils/install_packages.R`; added.

The script, dataset, and a new Snakemake rule (`gse141578_pseudobulk`) are now wired into `config/pipeline_config.yaml` and the `Snakefile`.

### 11.2 GSE141578 result: CXCL12–CXCR4 does not replicate, even in the correctly-matched cohort

After the fixes above, the pipeline was run twice: once without a covariate, once with a verified technical batch covariate (see below). Both runs load all 22/22 samples (24,478 total cells; 15,110 CD4+ T cells, 61.7%, identified via module-score comparison of CD4/CD8/NK/monocyte marker sets rather than a single-gene threshold).

Donor-level pseudobulk DESeq2 (`~batch + disease_group`, 11 disease vs. 11 HC pseudobulk samples) was used rather than cell-level testing, to directly address the pseudo-replication concern (D7/D8) that could not be resolved in GSE161192 (N=2 donors). **Batch covariate:** GEO provides no age or sex fields for this series, but individual `Sample_submission_date` values split cleanly into two verified library-preparation cohorts — batch1 (GSM4208xxx series, submitted 2019-12-06, 14 samples) and batch2 (GSM4404xxx series, submitted 2020-03-09, 8 samples) — not fully confounded with the disease_group contrast (each batch contains both HC and disease samples), though both DLB samples fall in batch2, so DLB-specific effects cannot be cleanly separated from batch in this design.

| Run | Significant genes (padj<0.05) of 13,575 tested | CXCR4 | CXCL12 |
|---|---|---|---|
| No batch covariate | 0 | log2FC=+0.50, padj=1.00 (ns) | Not detected in CD4+ T cells |
| With batch covariate | 1 (MRPL23, log2FC=−5.74, padj=0.050 — a mitochondrial ribosomal gene with no a priori relevance to any hypothesis here) | log2FC=+0.38, padj=0.99 (ns) | Not detected |

The "0 significant genes" result was hand-verified as correct Benjamini-Hochberg arithmetic, not a computation bug (best raw p=0.00096; 0.00096 × 13,575 tests ≈ 13, capped at padj=1).

**Conclusion:** the CXCL12–CXCR4 trafficking axis proposed by Gate et al. (2021) is not reproduced under this independent, donor-aware, batch-corrected reanalysis of the same cohort that originally produced the finding. This strengthens (does not merely repeat) the original D3 finding: D3 showed absence in an underpowered, stimulated 2-donor dataset; this shows absence in the correctly-matched, properly-powered, unstimulated dataset. Plausible explanations for the discrepancy with the original publication include analytic differences (cell-level vs. donor-level testing, unmeasured age/sex covariates — genuinely unavailable in GEO for this series) or T-cell-subset resolution not captured by the coarse CD4+/CD8+ clustering used here; this should be read as a formal non-replication, not a disproof of the underlying biology. **This is a conclusion-changing / conclusion-reinforcing result and was surfaced to the user for explicit sign-off before being propagated into the manuscript**, per this project's ground rules.

### 11.3 GCTA-COJO conditional analysis: BCAM is not an independent LBD risk gene

REPORT.md §5.2 flagged BCAM/APOE independence as unresolved, requiring conditional analysis or colocalization not previously performed. GCTA v1.94.1 and PLINK v1.90b6.21 (bioconda) were installed and used to close this gap:

- **Input:** all chr19:44.3–45.4 Mb SNPs (GRCh38 coordinates, matching the harmonized GWAS) from `data/LBD_GWAS_harmonised.tsv.gz` (3,293 SNPs; N=6,618 constant, per-SNP N not available), formatted as GCTA `.ma`.
- **LD reference:** 1000 Genomes EUR panel already present at `data/reference/g1000_eur` (503 individuals). This panel uses GRCh37/hg19 coordinates (previously flagged as D13) while the GWAS is GRCh38; SNPs were matched by rsID, not physical position, which is robust to this build mismatch — confirmed by checking that known APOE-defining SNPs (rs429358, rs7412) resolve correctly under rsID-based extraction after an initial BP-range-based extraction (made and caught during this session) incorrectly pulled the wrong physical region. 3,074/3,293 SNPs matched (93.4%), consistent with the 93.8% match rate noted elsewhere in this pipeline.
- **Stepwise selection** (`--cojo-slct`, p<5×10⁻⁸) identified exactly **3 independent genome-wide-significant signals** in the region (rs10405693, joint p=4.9×10⁻⁹; rs769449, joint p=2.2×10⁻⁶⁵; rs1081105, joint p=5.3×10⁻²⁸). None of the 3 corresponds to a BCAM-assigned SNP (BCAM's 43 MAGMA-assigned SNPs, per `data/magma_results/LBD_gene_analysis.genes.annot`, were checked explicitly).
- **Conditional analysis** (`--cojo-cond`) on all 3 jointly-selected signals: BCAM's most-associated SNP (rs28399637, unconditional p=6.6×10⁻¹³) drops to **pC=0.564** — a complete loss of signal. (An initial single-SNP conditioning on rs769449 alone left a residual pC=1.7×10⁻⁴, because rs769449 alone does not fully tag the APOE haplotype — rs429358, the canonical APOE ε4 SNP, did not survive matching to the reference panel; conditioning on all 3 jointly-selected signals resolves this.)

**Conclusion:** BCAM shows no evidence of an independent association with LBD risk once the locus's true independent signals are properly modeled. This formally resolves the open question in §5.2 and is now reflected in the manuscript (no further BCAM-specific mechanistic claim is warranted).

Outputs: `data/gcta_cojo/chr19_APOE_locus.ma`, `chr19_locus_ref.{bed,bim,fam}`, `chr19_stepwise.jma.cojo` (selected signals), `chr19_cond_on_all3.cma.cojo` (full conditional results). Not yet wired into the Snakefile as a rule (one-off diagnostic analysis, not part of the main per-sample DAG); recommended as a follow-up if this becomes a recurring pipeline step.

### 11.4 NicheNet (directional ligand-receptor inference) — root cause identified, not installed

REPORT.md §5.7 flagged the Prioritization Score as symmetric/non-directional and recommended NicheNet as a replacement. `nichenetr` is GitHub-only (not on CRAN/Bioconductor). All 88 transitive R package dependencies installed successfully, but the final `nichenetr` build fails because the transitive dependency `gdtools` requires the system `cairo` graphics library (`cairo-ft.h`), which is not installed on this machine and is not an R/conda package — it needs `sudo apt install libcairo2-dev` (Debian/Ubuntu) or a `conda-forge` `cairo`/`pkg-config` install. Root privileges were not used without explicit confirmation (per this project's operating rules), and a `conda install -c conda-forge cairo pkg-config` attempt did not resolve within a practical time budget in this session. `nichenetr` is confirmed **not installed** (`requireNamespace("nichenetr")` returns FALSE; no package directory exists under `.libPaths()`), despite an earlier ambiguous "SUCCESS" message from a wrapping `tryCatch` that only reflects that `install_github()` did not throw an R-level `stop()` — it does not reflect a fatal `configure`-time C build failure that R itself reports as a warning rather than an error. This is a real dependency gap, not a skipped task: **see §12 for the exact command to finish this** once `libcairo2-dev` (or conda-forge `cairo`) is available.

### 11.5 Manuscript rewrite

`MANUSCRIPT_Neuroimmune_Interactions_LBD.md` was rewritten in full to (a) match currently verified pipeline output numbers rather than the original fabricated/stale statistics (all Phase 1–4 tables re-pulled directly from current CSV/TSV outputs during this session, not copied from this report's prose, to avoid propagating any already-stale intermediate values), (b) incorporate the GSE141578 negative result and GCTA-COJO resolution above, (c) remove all substantia-nigra/dopaminergic-neuron-specific claims, which are not supported by the actual pipeline output (GSE178146 is anterior cingulate cortex; the current pipeline resolves neurons as a single population, 28,342 nuclei, with no dopaminergic/GABAergic/glutamatergic subtype breakdown despite the original manuscript citing specific per-subtype counts), and (d) reframe the Discussion/Conclusions around what is actually defensible: CD4+ clonal expansion (replicated), BCAM (formally ruled out), CXCL12–CXCR4 (formally non-replicated in two cohorts), and DLB cortical receptor changes (explicitly flagged as untested hypotheses, not findings). The "Causality Score" naming was already corrected to "Prioritization Score" in Phase C/D of this pipeline; the manuscript rewrite makes the non-causal framing explicit throughout rather than only in the Methods.

Numbers changed from what this report's own §3/§4 (written 2026-07-16/17) had documented, because the manuscript was compared against the *current* output files (post nFeature_max=6000 re-run) rather than against this report's earlier prose, which itself reflected an intermediate pipeline state. Notably: GSE161192 now yields 5,083 total cells (not 4,502), 1,241 CD4+ cells (24.4%, not 21.7%), 772 upregulated DEGs confirmed, and **IL26 no longer reaches significance** (padj=0.077) — this is a new, previously undocumented discrepancy between this report's §3/§4 and the current pipeline state, now corrected in both places.

---

## 12. Conference-readiness vs. journal-readiness: current status and remaining gap list

### 12.1 Conference-ready: yes, as of this session

The manuscript (§11.5) now states only claims traceable to a specific current output file, checked directly during this session (not copied from intermediate report prose). It leads with what is actually defensible — CD4+ clonal expansion (replicated), the BCAM/APOE resolution (formal conditional analysis), and the CXCL12–CXCR4 non-replication (tested in two independent cohorts, including the one that originally produced the claim) — and explicitly frames the DLB cortical receptor changes as untested hypotheses rather than findings. This is an appropriate scope and rigor level for a poster or short talk.

### 12.2 Remaining gaps for journal submission

| # | Gap | Why it's not resolved in this session | Effort to close |
|---|---|---|---|
| 1 | NicheNet directional ligand-receptor inference | Blocked on missing system `cairo` library (`libcairo2-dev`); all 88 R dependencies installed, only the final `nichenetr` build fails. Requires `sudo apt install libcairo2-dev` (needs explicit user confirmation — root access) or a longer `conda install -c conda-forge cairo pkg-config` solve than fit in this session. | Low, once cairo is available: `Rscript -e 'remotes::install_github("saeyslab/nichenetr")'`, then re-run Phase 3/4 with directional ligand-target scores in place of symmetric co-expression. |
| 2 | Mendelian randomization / colocalization for the 7 Bonferroni genes (esp. SNCA) and candidate receptors (TGFBR2, LGALS9, CD86) | Requires an eQTL reference dataset (e.g., GTEx, eQTLGen) not present in this repository and not fetched in this session. | Medium — need to source and integrate an eQTL summary-stats resource; `coloc` R package is CRAN-available and lightweight to install. |
| 3 | Age/sex-adjusted GSE141578 pseudobulk model | GEO provides no age or sex metadata for this series at the sample level. Not fixable without contacting the original authors or finding a supplementary table with per-subject covariates. | Not resolvable from public data as currently available. |
| 4 | DLB-subtype-specific (vs. combined PD+DLB) effect in GSE141578 | Both DLB samples fall in the same technical batch (batch2); disease-subtype and batch are confounded for DLB specifically. | Needs additional DLB samples from an independent batch/cohort — new data acquisition. |
| 5 | Independent validation cohort for TGFBR2/LGALS9/CD86 in DLB cortex | No second LBD-spectrum snRNA-seq cortex dataset is available in this repository. | New data acquisition (e.g., a second public GEO series covering LBD cortex, if one exists — not verified in this session) or wet-lab validation (IHC/ISH). |
| 6 | `renv.lock` / frozen reproducibility snapshot | Flagged in the original submission checklist (§10) and still not generated. | Low — `renv::init()` after confirming all packages installed; a few minutes of work. |
| 7 | GCTA-COJO not yet a Snakemake rule | Run manually in this session as a one-off diagnostic; not wired into the reproducible DAG. | Low — straightforward to add given the exact commands are already documented in §11.3. |
| 8 | Figure regeneration/captioning for the rewritten manuscript | The manuscript rewrite (§11.5) deliberately dropped all `Figura N` callouts rather than guess whether existing plot files (some still titled/named for the old "causal"/"PD vs Control" framing, e.g. `causal_network_graph.pdf` alongside the renamed `integrative_prioritization_network.pdf`) semantically match the corrected text. Plots exist for most panels needed but have not been re-inspected or re-captioned against the new claims. | Medium — requires opening each PDF, confirming its content matches the corrected framing, and either relabeling or regenerating before re-inserting figure callouts. |
| 9 | Manuscript author/affiliation placeholders | `[Apellido completo]`, institutional affiliation, address, and email are still placeholders — this is information only the author can supply, not something inferable from the data. | User action required. |
| 10 | Experimental validation | This entire study is computational. No IHC, in situ hybridization, or co-culture work has been performed or is claimed. | Out of scope for this pipeline; a wet-lab or collaborator-dependent next phase. |

Items 1, 6, and 7 are low-effort and could reasonably be completed in a short follow-up session. Items 2, 3, 4, 5, 8, and 10 require new data, new external resources, or are not resolvable computationally at all — these are genuine, not merely time-boxed, limitations of the current investigation and should be stated as such in any submission.
