# Investigation Review & Upgrade Summary
_Handoff document for Claude Code. Read this fully before touching any script._

## 0. Purpose of this document
This summary captures the current state of a multi-omic bioinformatics investigation on
**Lewy Body Dementia (LBD)** and the neuroimmune synapse, plus every defect and
optimization target found during an expert manuscript review. Use it as ground truth and
context. The manuscript is only one output — **the scripts and their intermediate/final
outputs contain the real evidence**, so verify claims against code and data, not prose.

---

## 1. What the investigation is

**Title (working):** Genetic risk convergence on CD4⁺ T cell–neuron interactions in the
LBD neuroimmune synapse.

**Central hypothesis:** LBD pathogenesis is driven by a genetically primed neuroimmune
synapse where clonally expanded, antigen-experienced CD4⁺ T cells interact with vulnerable
cortical neurons (Layer 5/6 excitatory).

**Design — a four-phase computational pipeline:**
1. **Phase 1 – Genetic prioritization.** MAGMA v1.10 gene-based test on LBD GWAS summary
   stats. SNPwise-mean model, 1000G-EUR LD reference, GRCh38, Bonferroni correction.
   Produces ranked risk genes + gene sets at top 1/5/10%.
2. **Phase 2 – CSF immune profiling.** scRNA-seq + TCR-seq of CSF T cells (Seurat v5 +
   scRepertoire). Clonal expansion defined by identical CDR3 aa in ≥2 cells. Wilcoxon DE:
   expanded vs non-expanded CD4⁺ → "expansion signature".
3. **Phase 3 – Cortical snRNA-seq.** Neuronal subtyping (Harmony integration, scDblFinder
   doublets). Excitatory (SLC17A7+) / inhibitory (GAD1/GAD2+), layer classification.
4. **Phase 4 – Interactome + network.** Ligand-receptor scoring (CellPhoneDB, OmniPath)
   between Phase-2 clones and Phase-3 neurons; weighted "causal" network integrating
   Phase-1 Z-scores, Phase-2 log₂FC, Phase-3 interaction weights.

**Data sources (verify provenance in code):**
- GWAS: **Chia R, et al. Nat Genet 2021;53:294–303** (2,591 cases / 4,027 controls; GWAS
  Catalog **GCST90001390**). NOTE: this is currently **used but not cited** in the manuscript.
- CSF scRNA/TCR: **GSE161192** (originates from Gate et al., Science 2021).
- Cortex snRNA: **GSE178146**.

---

## 2. Verified ground-truth facts (do not reintroduce the old errors)

- Real Chia 2021 genome-wide loci = **GBA, BIN1, TMEM175, SNCA-AS1, APOE**. BCAM was **not**
  an independent genome-wide hit.
- BCAM, APOC1, APOE (and BCL3) all lie in the **same chr19q13 LD block**. A high MAGMA
  gene score for BCAM most likely reflects the **APOE locus via LD**, not an independent
  adhesion-molecule mechanism.
- Gate et al. (Science 2021;374(6569):868–874) **already reported** clonal CD4⁺ expansion,
  CXCR4 upregulation, the CXCL12–CXCR4 trafficking axis and Th17/IL-17A in this CSF data.
  Novelty of the present work must be framed around Phase-1 integration, the TPM4/CTSH
  signature, and the L5/6 ligand-receptor mapping — not around re-finding CXCR4.
- Correct method citations: Seurat v5 = **Hao et al., Nat Biotechnol 2024;42:293–304**;
  CellPhoneDB = **Efremova et al., Nat Protocols 2020;15:1484–1506**; Seurat integration
  (Stuart) = **Cell 2019;177:1888–1902** (not Nature Methods).

---

## 3. CRITICAL defects to fix (blocking)

1. **BCAM vs BCL3 contradiction.** Abstract/Results say BCAM is top gene; Discussion says
   "BCL3 and GBA". Only one can be correct — resolve against the actual MAGMA output table.
2. **Figure 4 mislabeled "Parkinson's Disease (PD)"** ("PD vs Control"). The whole study is
   LBD. Determine from the scripts whether this is (a) a wrong hardcoded label/title, or
   (b) the wrong dataset was loaded. This is the single most important thing to verify.
3. **Figure 4A network does not contain the headline axes.** Text centers on CXCL12–CXCR4
   and MHC-II (HLA-DRA/B)–TCR, but the plotted network shows cytokine receptors
   (IL17RA/RC, IL10RB, TGFBR1/2/3, OSMR, LIFR, TNFRSF1A, LGALS9, CD86, CD40, IL2RA/RB).
   Reconcile: is the figure filtering out the reported axes, or is the text unsupported?
4. **"Causal" is overstated.** The network is a weighted integration, not causal inference.
   Either rename (integrative/prioritization) or add real causal methods (see §5).

---

## 4. Statistical / methodological issues to optimize

- **Pseudo-replication.** DE p-values (e.g. padj ≈ 1e-45) treat cells as independent.
  Re-run DE as **pseudobulk per donor** (DESeq2/edgeR on aggregated counts). Report n donors
  driving each gene.
- **Clonal-expansion threshold too permissive** ("≥2 cells, identical CDR3 aa"). Report
  clone-size distributions, confirm nucleotide-level identity, test a stricter cutoff.
- **Confound: "stimulated vs unstimulated" vs "expanded vs non-expanded."** GSE161192 has an
  in-vitro α-synuclein stimulation arm. Ensure the expansion signature (TPM4/ACTB/CTSH) is
  not just an activation/stimulation effect. Model both factors.
- **Cytoskeletal signature caution.** TPM4/ACTB/PFN1/ACTG1/CFL1/VIM + loss of IL7R/TCF7 is
  the generic naïve→effector transition and an ambient-RNA-prone module. Check for
  ambient/contamination (e.g. SoupX/DecontX) before claiming specificity.
- **QC wording/params.** CSF section says "nuclei" but these are cells; feature ceiling of
  2,500 is low for activated T cells and may exclude the most active clones. Re-examine QC
  thresholds and justify them.
- **chr19 / BCAM.** Add conditional/joint analysis on APOE, and colocalization, before any
  BCAM mechanistic claim.

---

## 5. Steps to move toward causality (implement where feasible)

- Conditional/joint MAGMA or GCTA-COJO on chr19 (separate BCAM from APOE).
- Colocalization (coloc/SMR) between LBD GWAS and immune/neural eQTLs for BCAM, CTSH, CXCR4.
- Mendelian randomization using eQTL/pQTL instruments (e.g. genetically predicted CXCL12 /
  CTSH → LBD risk).
- Directional signaling with **NicheNet** (ligand → neuronal target-gene regulatory
  potential) to replace symmetric CellPhoneDB scores where directionality is claimed.
- Independent replication on a second LBD cohort not derived from Gate et al.

---

## 6. Reference-list cleanup (regenerate from DOIs)
Known-wrong entries: Stuart/Seurat (wrong journal), CellPhoneDB (wrong journal), Gate
(wrong pages), McKeith 2017 (89(1) not 89(8)), ref [11] scRepertoire2 (placeholder), ref
[20] (placeholder "[Source data context 2025]"). Add missing **Chia 2021** GWAS citation and
proper GEO dataset attributions. Verify [3],[4],[19],[21],[26] individually.

---

## 7. Upgrade / optimization goals (for the scripts themselves)
- Reproducibility: pin package versions, set/record seeds everywhere, emit `sessionInfo()` /
  `renv.lock` (R) and `requirements.txt`/`environment.yml` (Python); add a config file so
  dataset paths, thresholds, and labels are not hardcoded (this is how the "PD" label likely
  crept in).
- Turn the four phases into a reproducible, re-runnable pipeline (Makefile / Snakemake /
  targets) with clear inputs→outputs and cached intermediates.
- Every figure should be regenerated from a script that reads the actual result tables, with
  titles/labels derived from a config (no hardcoded disease names).
- Add assertions/tests: dataset identity checks, cell/gene count sanity checks, "no NA in Z
  columns", CDR3 matching unit test, and a check that figure inputs match the reported axes.
- Optimize performance where relevant (vectorize, avoid dense matrix blow-ups, parallelize
  MAGMA/permutations, use sparse ops in Seurat).

---

## 8. Definition of done
- Pipeline re-runs end-to-end from raw inputs to every figure and table.
- No LBD/PD label mismatches; Phase-1 top-gene claim is internally consistent and
  LD-aware.
- Figures match the text (or the text is corrected to match the real results).
- DE is donor-aware (pseudobulk); expansion vs stimulation deconfounded.
- Reference list and data-availability regenerated and verified.
- A written REPORT.md summarizes what changed, what was confirmed, and what still needs wet-lab/statistical validation.
