# CLAIMS_CHECK.md — Phase B Verification Report
_Generated 2026-07-16. Every claim checked against actual script code and output files._
_Status codes: ✓ CONFIRMED | ✗ WRONG | ~ PARTIALLY CORRECT | ? UNVERIFIABLE from current data_

---

## Phase 1 — Genetic Prioritization (MAGMA)

| # | Manuscript Claim | Actual Value | Source | Status |
|---|-----------------|--------------|--------|--------|
| 1.1 | GWAS N = "2,981 cases and 4,347 controls" | N=2,591 cases + 4,027 controls = 6,618 total | `data/LBD_GWAS_metadata.yaml` | ✗ WRONG |
| 1.2 | MAGMA run on genome-wide SNPs (~whole GWAS) | Only 454 pre-filtered SNPs (p≤1e-5) input | `data/magma_results/LBD_gene_analysis.log` | ✗ WRONG |
| 1.3 | "19 genes reach genome-wide significance (p<2.5×10⁻⁶)" | 15 genes at that threshold (34 total analyzed) | `data/LBD_risk_genes_ranked.tsv` | ✗ WRONG |
| 1.4 | Top MAGMA gene = BCAM (Z=7.09) | BCAM rank #1, Z=7.09 ✓ | `data/LBD_risk_genes_ranked.tsv` | ✓ CONFIRMED |
| 1.5 | Discussion: "BCL3 and GBA are top MAGMA genes" | BCL3 **absent** from all 34 analyzed genes | `data/LBD_risk_genes_ranked.tsv` | ✗ WRONG |
| 1.6 | SNCA among top-ranked genes | SNCA rank #9, Z=5.81 | `data/LBD_risk_genes_ranked.tsv` | ✓ CONFIRMED |
| 1.7 | MAGMA uses Bonferroni correction | Applied at p<2.5×10⁻⁶ (correct for GWAS-scale), but gene-level Bonferroni would be 0.05/34=0.00147 since only 34 genes were tested | `scripts/phase1_genetic_prior/04_parse_magma_results.py` | ~ PARTIALLY CORRECT |
| 1.8 | GWAS citation: Chia et al. (implied) | GWAS **not cited** in manuscript despite being the primary data source | `INVESTIGATION_REVIEW_SUMMARY.md` | ✗ WRONG |
| 1.9 | `02_prepare_magma_input.py`: sample_size labeled "Cases + controls" = 2,591 | 2,591 is cases only; 6,618 is correct total | `scripts/phase1_genetic_prior/02_prepare_magma_input.py` | ✗ WRONG |
| 1.10 | `01_filter_gwas_data.py` file path | Hardcoded path has typo ("neuroinmune" missing 'm') and missing leading `/`; script cannot run standalone | `scripts/phase1_genetic_prior/01_filter_gwas_data.py` | ✗ WRONG |
| 1.11 | Comment in filter script says p≤1e-8 | Actual filter threshold is p≤1e-5 | `scripts/phase1_genetic_prior/01_filter_gwas_data.py` | ✗ WRONG |

**Critical note on 1.2–1.3:** Because MAGMA received only 454 SNPs (pre-filtered at p≤1e-5), it mapped to only 34 genes — not the ~18,000 expected from a genome-wide gene-based test. This invalidates the MAGMA analysis as currently run. All gene counts, significance thresholds, and gene rankings must be treated as provisional until MAGMA is re-run on all harmonized SNPs.

**Critical note on 1.4 + 1.5 (BCAM/BCL3 contradiction):** BCAM is the true top gene in the current output. BCL3 is absent. However, BCAM, APOC1, APOE, TOMM40, and PVRL2 all reside in the same chr19q13 LD block as APOE. BCAM's top rank likely reflects the APOE locus via LD contamination, not an independent BCAM signal. This is not disclosed anywhere in the manuscript.

---

## Phase 2 — CSF Immune Profiling (scRNA-seq + TCR-seq)

| # | Manuscript Claim | Actual Value | Source | Status |
|---|-----------------|--------------|--------|--------|
| 2.1 | Total cells after QC: 12,847 | 4,502 cells | `data/GSE161192_processed/all_cells_with_tcr_annotations.csv` | ✗ WRONG |
| 2.2 | Cells with TCR: 8,932 (69.5%) | 3,331 (74.0%) | same file | ✗ WRONG |
| 2.3 | Unique clonotypes: 1,247 | 1,515 | same file | ✗ WRONG |
| 2.4 | Expanded cells: 2,891 (32.4% of TCR) | 2,044 (61.4% of TCR) | same file | ✗ WRONG |
| 2.5 | Largest clone: 87 cells | 501 cells | same file | ✗ WRONG |
| 2.6 | CD4+ cells: 5,234 (40.7% of total) | 978 (21.7% of total) | same file | ✗ WRONG |
| 2.7 | CD4+ expanded: 1,876 (35.8% of CD4+) | 589 (60.2% of CD4+) | same file | ✗ WRONG |
| 2.8 | DE: 1,523 significant genes (padj<0.05) | 791 significant genes | `data/GSE161192_processed/DEG_expanded_vs_nonexpanded_CD4.csv` | ✗ WRONG |
| 2.9 | DE: 687 upregulated genes | 642 upregulated | same file | ✗ WRONG |
| 2.10 | IL17A: log2FC=2.34, padj=3.2e-45 | log2FC=1.990, padj=7.27e-10 | same file | ✗ WRONG |
| 2.11 | IL26: log2FC=2.01, padj=1.1e-38 | log2FC=0.424, padj=4.24e-4 | same file | ✗ WRONG |
| 2.12 | OSM: log2FC=1.87, padj=2.5e-32 | log2FC=0.615, padj=1.28e-4 | same file | ✗ WRONG |
| 2.13 | TGFB1: log2FC=1.45, padj=8.3e-26 | log2FC=0.410, padj=7.84e-6 | same file | ✗ WRONG |
| 2.14 | ENTPD1: log2FC=1.62, padj=4.7e-28 | log2FC=0.423, padj=1.64e-4 | same file | ✗ WRONG |
| 2.15 | CXCL12 upregulated; CXCL12–CXCR4 axis drives T cell trafficking | CXCL12 **not found** in DEG output; CXCR4 log2FC=**−0.955** (DOWNREGULATED) | same file | ✗ WRONG |
| 2.16 | TPM4 is top upregulated gene | TPM4 log2FC=1.164, padj=1.71e-45 ✓ | same file | ✓ CONFIRMED |
| 2.17 | QC wording: "nuclei filtered" | These are CSF **cells**, not nuclei | `scripts/phase2_immune_profiling/01_process_scrna_seurat.R` | ✗ WRONG |
| 2.18 | QC feature ceiling not stated explicitly | Script uses nFeature_max=2,500 — may exclude most activated T cells | `scripts/phase2_immune_profiling/01_process_scrna_seurat.R` | ~ CONCERN |
| 2.19 | Expansion defined by clonal TCR | cloneCall="aa" (amino acid CDR3), not nucleotide-level — permissive | `scripts/phase2_immune_profiling/02_integrate_tcr_data.R` | ~ CONCERN |
| 2.20 | DE controls for stimulation condition | DE in `03_extract_cd4_signatures.R` does NOT model Stimulated vs Unstimulated | `scripts/phase2_immune_profiling/03_extract_cd4_signatures.R` | ✗ WRONG |
| 2.21 | Pseudobulk / donor-aware DE | Wilcoxon on individual cells; N=2 donors — pseudo-replication | `scripts/phase2_immune_profiling/03_extract_cd4_signatures.R` | ✗ WRONG |

**Critical note on 2.1–2.7 (cell count discrepancies):** The manuscript numbers are approximately 3× higher than observed in the processed output files for total cells and CD4+ counts. This suggests either (a) the QC script was re-run with different parameters after the reported analysis, or (b) the manuscript numbers were not derived from the current pipeline outputs. The expansion rate discrepancy is also qualitatively reversed: manuscript reports 32.4% expansion rate, but actual data shows 61.4% — a fundamentally different biological conclusion.

**Critical note on 2.15 (CXCL12/CXCR4):** This is a conclusion-changing defect. The manuscript's central mechanistic claim about T cell trafficking via CXCL12–CXCR4 is unsupported. CXCL12 is absent from DE output; CXCR4 is downregulated. Note: Gate et al. (Science 2021) — the source of this dataset — already reported CXCR4 upregulation in CSF T cells. If the current analysis contradicts Gate et al. using the same dataset, the QC or analysis pipeline has a problem (likely the QC ceiling at 2,500 features removing activated T cells).

---

## Phase 3 — Cortical snRNA-seq (GSE178146)

| # | Manuscript Claim | Actual Value | Source | Status |
|---|-----------------|--------------|--------|--------|
| 3.1 | Dataset: "Kamath et al., 2022" | Dataset is **Feleke et al., 2021** (Acta Neuropathol 142:449); PMID 34383113 | GEO page GSE178146 | ✗ WRONG |
| 3.2 | Brain region: "substantia nigra" | Region is **anterior cingulate cortex** | GSE178146 GEO metadata | ✗ WRONG |
| 3.3 | "27 samples (2 controls and 25 PD patients)" | **28 samples**: 7 Control, 7 DLB, 7 PDD, 7 PD | GSE178146 individual sample pages | ✗ WRONG |
| 3.4 | Comparison: "PD vs Control" | Samples labeled "PD" include 5 true controls (PDC05, PDC22, PDC34, PDC87, PDC91) due to `ifelse(grepl("^C"))` bug | `scripts/phase3_neuronal_vulnerability/01_load_and_qc_snnucseq.R:80` | ✗ WRONG |
| 3.5 | Study is about LBD | The correct primary comparison is **DLB vs Control** (and PDD vs Control); the current analysis collapses all 3 disease groups into "PD" | same script | ✗ WRONG |
| 3.6 | Downsampled to 5,000 nuclei per sample | Script uses `MAX_CELLS=2000` | `scripts/phase3_neuronal_vulnerability/02_downsampling_integrate_and_cluster.R` | ✗ WRONG |
| 3.7 | 17 active ligand-receptor pairs | **16 pairs** in actual output | `data/GSE178146_processed/active_ligand_receptor_pairs.csv` | ✗ WRONG |
| 3.8 | CXCL12–CXCR4 axis in Figure 4 | CXCL12 not in CD4+ upregulated genes; pair not activating in L-R analysis | `data/GSE161192_processed/expanded_CD4_upregulated_genes.txt`; `scripts/phase3_neuronal_vulnerability/04_ligand_receptor_analysis.R` | ✗ WRONG |
| 3.9 | HLA-DRA/B → TCR axis in Figure 4 | HLA-DRA/B–TCR pair **not in hardcoded L-R database** | `scripts/phase3_neuronal_vulnerability/04_ligand_receptor_analysis.R` | ✗ WRONG |

**Critical note on 3.1–3.5 (dataset identity + condition labeling):**
These are the most severe defects in the pipeline. The correct dataset is being used (Feleke 2021 LBD spectrum study) but:
- It is wrongly attributed (Kamath 2022, which is a different dataset)
- The brain region is misreported (substantia nigra vs anterior cingulate cortex)
- The sample count is wrong (27 vs 28)
- Five control samples (PDC05, PDC22, PDC34, PDC87, PDC91) that begin with "PD" are misclassified as disease
- All disease groups (DLB, PDD, PD) are collapsed into a single "PD" label

This means every Figure 3 and Figure 4 result comparing "PD vs Control" is based on:
- A mis-attributed brain region
- A wrong disease label (should be "DLB" for the primary comparison)
- Contaminated control group (5 controls mis-assigned to disease)
- Loss of the DLB/PDD/PD distinction that is central to the LBD spectrum

**Verified complete sample→group mapping (from GEO individual sample pages, 2026-07-16):**
```
Control (n=7): C36, C48, PDC05, PDC22, PDC34, PDC87, PDC91
DLB     (n=7): PD060, PD115, PD163, PD294, PD332, PD566, PD706
PDD     (n=7): PD341, PD366, PD415, PD501, PD531, PD563, PD678
PD      (n=7): PD413, PD416, PD523, PD666, PD683, PD732, PD747
```

---

## Phase 4 — Weighted Interactomics / Network

| # | Manuscript Claim | Actual Value | Source | Status |
|---|-----------------|--------------|--------|--------|
| 4.1 | Top interaction: OSM–LIFR (score 0.631) | Confirmed ✓ | `results_final/weighted_neuroimmune_network.csv` | ✓ CONFIRMED |
| 4.2 | ENTPD1–ADORA1 score 0.294 | Confirmed ✓ | same file | ✓ CONFIRMED |
| 4.3 | Network represents "causal" neuro-immune interactions | Formula is `Expression × (1 + max(0, Z))` — a heuristic prioritization score, not causal inference | `scripts/phase4_weighted_interactomics/01_weighted_interactomics.R` | ✗ WRONG |
| 4.4 | Receptor Z-scores reflect GWAS risk | All receptor Z-stats in network = 0 (no GWAS support for any receptor gene) | `results_final/weighted_neuroimmune_network.csv` | ~ CONCERN |

---

## Methods Section — Additional Errors

| # | Claim | Actual Value | Status |
|---|-------|--------------|--------|
| M1 | Seurat integration ref: Stuart et al., Nat Methods | Should be Stuart et al., Cell 2019;177:1888–1902 | ✗ WRONG |
| M2 | CellPhoneDB ref: wrong journal | Should be Efremova et al., Nat Protocols 2020;15:1484–1506 | ✗ WRONG |
| M3 | Gate et al. page numbers: wrong | Should be Science 2021;374(6569):868–874 | ✗ WRONG |
| M4 | McKeith 2017: vol 89(8) | Should be vol 89(1) | ✗ WRONG |
| M5 | Ref [11] (scRepertoire2): placeholder | Incomplete citation | ✗ WRONG |
| M6 | Ref [20]: "[Source data context 2025]" | Placeholder — not a real citation | ✗ WRONG |
| M7 | Chia 2021 GWAS (primary data source) | **Not cited anywhere in manuscript** | ✗ WRONG |
| M8 | Seurat v5 ref: Hao et al. | Should specify Nat Biotechnol 2024;42:293–304 | ~ |

---

## Summary Statistics

| Category | Total claims checked | Confirmed ✓ | Wrong ✗ | Concern ~ |
|----------|---------------------|-------------|---------|-----------|
| Phase 1 (MAGMA) | 11 | 2 | 8 | 1 |
| Phase 2 (CSF scRNA) | 21 | 2 | 16 | 3 |
| Phase 3 (snRNA/L-R) | 9 | 0 | 9 | 0 |
| Phase 4 (Network) | 4 | 2 | 1 | 1 |
| Methods/References | 8 | 0 | 7 | 1 |
| **TOTAL** | **53** | **6** | **41** | **5** |

---

## Prioritized Fix List (for Phase C)

### BLOCKING — Must fix before any re-submission
1. **D1** — `01_load_and_qc_snnucseq.R:80`: Replace `ifelse(grepl("^C"))` with verified lookup table. All Phase 3 results must be regenerated with DLB/PDD/PD/Control labels.
2. **D2** — Manuscript Discussion: Remove BCL3 claim; replace with BCAM + chr19/APOE-LD caveat.
3. **D3** — Remove CXCL12–CXCR4 axis claim. State that CXCL12 was not DE and CXCR4 was downregulated. Investigate whether nFeature_max=2500 QC ceiling explains the discrepancy with Gate et al. findings.
4. **D4** — Remove HLA-DRA/B→TCR axis claim or add the pair to the L-R database with evidence.
5. **D5** — Rename "causal" → "integrative prioritization" throughout manuscript and script.
6. **D6** — Re-run MAGMA on all harmonized SNPs (not pre-filtered 454 SNPs).
7. **D7/D8** — Re-run DE as pseudobulk (DESeq2/edgeR) with stimulation condition as covariate.
8. **Brain region** — Correct Methods: "anterior cingulate cortex" not "substantia nigra".
9. **Dataset citation** — Correct Methods: Feleke et al. 2021, not Kamath et al. 2022.
10. **Sample count** — Correct Methods: 28 samples (7 per group), not "27 samples (2 controls and 25 PD)".
11. **GWAS sample size** — Correct to 2,591 cases + 4,027 controls = 6,618 total.
12. **GWAS citation** — Add Chia et al. Nat Genet 2021 as primary data source.
13. **All DEG stats** (IL17A, IL26, OSM, TGFB1, ENTPD1) — Update fold-changes and p-values to match actual output.
14. **All cell counts** — Reconcile manuscript numbers (3× inflated) with actual output.

### HIGH PRIORITY — Fix before final submission
15. **D9** — Re-examine nFeature_max=2500 QC ceiling; increase and re-run if activated T cells are excluded.
16. **D12** — Update MAGMA N to 6,618.
17. **D18** — Fix path typo in `01_filter_gwas_data.py`.
18. **chr19 LD** — Add APOE conditional analysis / colocalization before any BCAM mechanistic claim.
19. **Reference list** — Fix all 8 citation errors (Stuart, CellPhoneDB, Gate, McKeith, refs 11/20, add Chia).
20. **Downsampling** — Correct Methods to state 2,000 nuclei per sample (or change script to 5,000 and re-run).
