# Decoding the Neuroimmune Synapse in Lewy Body Dementia

**Integración multi-ómica de riesgo genético, perfilado inmune y vulnerabilidad neuronal**

[![License](https://img.shields.io/badge/License-See%20LICENSE-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-4.3.1-276DC3.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](https://www.python.org/)
[![Snakemake](https://img.shields.io/badge/Snakemake-DAG%20validado-039475.svg)](Snakefile)

> **Nota sobre este README:** esta versión reemplaza por completo una versión anterior que contenía cifras fabricadas/no verificadas (top gen incorrecto, dataset mal atribuido, conteos de células ~3× inflados). Todo lo que sigue está verificado directamente contra los archivos de salida del pipeline. El proceso de auditoría que llevó a esta reescritura está documentado en `archive/` (histórico, no vigente).

## Qué es esto

Pipeline computacional de cuatro fases (más tres extensiones) que integra GWAS, scRNA-seq con perfilado de repertorio TCR, y snRNA-seq para poner a prueba, con verificación cruzada entre cohortes y métodos, si el riesgo genético de la demencia por cuerpos de Lewy (DCL) converge con señales inmunes mediadas por linfocitos T CD4+ hacia neuronas corticales vulnerables.

Para el flujo completo fase por fase, ver **[`PIPELINE_DIAGRAM.md`](PIPELINE_DIAGRAM.md)**. Para la redacción científica completa, ver **[`MANUSCRIPT_Neuroimmune_Interactions_LBD.md`](MANUSCRIPT_Neuroimmune_Interactions_LBD.md)**.

## Hallazgos principales (verificados contra archivos de salida, septiembre 2026)

1. **Expansión clonal de linfocitos T CD4+ en LCR** — confirmada en dos cohortes independientes (GSE161192, GSE141578). El hallazgo más sólido del estudio.
2. **BCAM descartado formalmente como gen de riesgo independiente** — MAGMA genome-wide (18.523 genes) + GCTA-COJO condicional muestran que su asociación es enteramente atribuible al desequilibrio de ligamiento con el locus APOE (p condicional: 6,6×10⁻¹³ → 0,56).
3. **El eje CXCL12–CXCR4 no replica** en ninguna de las dos cohortes de LCR, incluyendo la cohorte de origen de esa hipótesis publicada.
4. **MMRN1 colocaliza fuertemente con expresión en linfocitos T CD4+** (PP4=0,97) pero no en corteza cerebral (variantes distintas) — evidencia de convergencia genética-inmune tejido-específica.
5. **LIFR, LGALS9 y CD86 muestran un gradiente de expresión significativo y monótono** a lo largo del espectro clínico Control→PD→PDD→DLB (pseudobulk DESeq2, padj<0,05).
6. Sin evidencia de enriquecimiento poligénico en pathways inmunes (Th17, TCR, citoquinas, quimiocinas, JAK-STAT) — resultado nulo, reportado como tal.
7. **Mendelian randomization y SuSiE-coloc confirman independientemente el hallazgo de MMRN1** (beta=-0,29, padj=2,5×10⁻¹⁰ en linfocitos T; señal específica validada a nivel de credible set).
8. **NicheNet (análisis direccional) revierte la interpretación mecanística del eje OSM–LIFR**: OSM rankea 28 de 29 ligandos candidatos por actividad dirigida sobre la respuesta transcripcional de DLB en neuronas ACC; IL17A rankea último. El gradiente de severidad de LIFR (hallazgo #5) se mantiene, pero sin evidencia de que OSM lo explique.

## Arquitectura del pipeline

| Fase | Qué hace | Scripts |
|---|---|---|
| **1 — Priorización genética** | MAGMA gen-based test sobre el GWAS completo de DCL | `scripts/phase1_genetic_prior/01-04` |
| **1b — Resolución BCAM** | GCTA-COJO condicional en el locus chr19 | `scripts/phase1_genetic_prior/05-06` |
| **1c — Colocalización GWAS-eQTL** | coloc.abf vs. GTEx corteza y BLUEPRINT células T (tabix remoto, sin descarga masiva) | `scripts/phase1_genetic_prior/07-10` |
| **1d — Enrichment de pathways** | MAGMA gene-set sobre pathways inmunes de KEGG (pre-especificados) | `scripts/phase1_genetic_prior/11-12` |
| **1e — Mendelian randomization** | Wald ratio / IVW (instrumentos GCTA-COJO para chr19) | `scripts/phase1_genetic_prior/13` |
| **1f — SuSiE-coloc (fine-mapping)** | Colocalización multi-variante para MMRN1 | `scripts/phase1_genetic_prior/14` |
| **2a — Perfilado inmune LCR** | scRNA-seq+TCR-seq, clonalidad, DE Wilcoxon+pseudobulk | `scripts/phase2_immune_profiling/01-03` |
| **2b — Réplica independiente LCR** | Pseudobulk DESeq2 a nivel de donante en cohorte más grande | `scripts/phase2_immune_profiling/04` |
| **2c — Sensibilidad a batch** | DESeq2 sin covariable de lote + subset PD-only | `scripts/phase2_immune_profiling/05` |
| **3 — Vulnerabilidad neuronal cortical** | snRNA-seq, subtipificación neuronal por capa, mapeo ligando-receptor | `scripts/phase3_neuronal_vulnerability/01-04` |
| **3b — Gradiente de severidad clínica** | DESeq2 pseudobulk, severidad ordinal pre-especificada | `scripts/phase3_neuronal_vulnerability/05` |
| **3c — Actividad ligando dirigida** | NicheNet (entorno conda aislado, ver Limitaciones) | `scripts/phase3_neuronal_vulnerability/06` |
| **4 — Red de priorización integrada** | Score heurístico expresión × riesgo genético (explícitamente no causal) | `scripts/phase4_weighted_interactomics/01` |
| **4b — Verificación por permutación** | Null empírico del ranking (10.000 permutaciones) | `scripts/phase4_weighted_interactomics/02` |

Todo el pipeline está orquestado por `Snakefile` + `config/pipeline_config.yaml` (única fuente de verdad para rutas, umbrales y etiquetas — ningún valor debería estar hardcodeado en los scripts en producción).

```bash
# Pipeline completo
snakemake --cores 4

# Dry-run (ver qué se ejecutaría)
snakemake -n

# Validar outputs
snakemake validate --cores 1
```

> Nota de entorno: Snakemake no viene preinstalado en el `PATH` por defecto de este proyecto. Se usó un entorno conda dedicado (`conda create -n snakecheck -c bioconda -c conda-forge snakemake-minimal`) para validar el DAG completo.

## Fuentes de datos

| Dataset | Accesión | Referencia | Uso |
|---|---|---|---|
| GWAS DCL | GCST90001390 | Chia et al. 2021, *Nat Genet* 53:294-303 | Fase 1/1b/1c/1d — N=6.618 (2.591 casos, 4.027 controles) |
| CSF scRNA+TCR (estimulado) | GSE161192 | Gate et al. 2021, *Science* 374:868-874 | Fase 2a — 2 donantes |
| CSF scRNA (sin estimular) | GSE141578 | Gate et al. 2021 (mismo estudio) | Fase 2b — 22 donantes, cohorte de origen de la hipótesis CXCR4 |
| snRNA-seq corteza | GSE178146 | Feleke et al. 2021, *Acta Neuropathol* 142:449-474 | Fase 3/3b — 28 muestras, corteza cingulada anterior (no sustancia nigra) |
| eQTL corteza | GTEx v8 (vía eQTL Catalogue, QTD000151) | — | Fase 1c |
| eQTL células T CD4+ | BLUEPRINT (vía eQTL Catalogue, QTD000031) | — | Fase 1c |
| Gene sets inmunes | KEGG (hsa04659, hsa04660, hsa04060, hsa04062, hsa04630) | — | Fase 1d |
| Referencia LD | 1000 Genomes Phase 3, EUR | — | Fase 1/1b (GRCh37/hg19 — atención al desajuste de build, ver Métodos) |

## Dependencias

- **R 4.3.1** con `renv.lock` (201 paquetes fijados) — `renv::restore()` para reproducir el entorno exacto.
- **Python 3.10+** — `requirements.txt`.
- **MAGMA v1.10** (`tools/magma`), **GCTA v1.94+**, **PLINK v1.9**, **tabix/bgzip** (htslib).
- **Snakemake ≥7.32** (no incluido en `environment.yml` por defecto del sistema; ver nota arriba).
- **Docker** (opcional) — ver `Dockerfile`: imagen autocontenida (2,1GB) para el critical path de Fase 1/1b/1c (MAGMA+GCTA+PLINK+coloc+MR vía conda-forge/bioconda, sin depender del entorno de esta máquina). Validada: reproduce `coloc_summary.csv` byte a byte contra el resultado del host. Nota: con Docker instalado vía snap, el bind-mount de volúmenes fuera de `$HOME` puede fallar silenciosamente (confinamiento de snap) — montar desde una ruta dentro de `$HOME` o extraer resultados con `docker cp`/redirección de stdout en vez de `-v /tmp/...`.

```bash
docker build -t lbd-neuroimmune-phase1 .
docker run --rm \
  -v "$(pwd)/data/coloc/gwas:/pipeline/data/coloc/gwas:ro" \
  -v "$(pwd)/data/coloc/eqtl:/pipeline/data/coloc/eqtl:ro" \
  lbd-neuroimmune-phase1 \
  bash -c "cd /pipeline && Rscript scripts/phase1_genetic_prior/09_run_coloc.R brain_acc"
```

## Estructura del repositorio

```
.
├── PIPELINE_DIAGRAM.md          ← diagrama + tabla fase-por-fase
├── MANUSCRIPT_Neuroimmune_Interactions_LBD.md  ← redacción científica completa
├── config/pipeline_config.yaml  ← única fuente de verdad (rutas, umbrales, etiquetas)
├── Snakefile                    ← DAG completo, 30 reglas
├── renv.lock                    ← entorno R fijado
├── scripts/
│   ├── phase1_genetic_prior/    ← 1, 1b, 1c, 1d (12 scripts)
│   ├── phase2_immune_profiling/ ← 2a, 2b (4 scripts)
│   ├── phase3_neuronal_vulnerability/  ← 3, 3b (5 scripts)
│   ├── phase4_weighted_interactomics/  ← 4 (1 script)
│   └── utils/                   ← instalación de paquetes, validación
├── data/                        ← datos crudos y procesados (no versionados)
├── results_final/               ← red integrada + figuras
└── archive/                     ← documentos de auditoría histórica (no vigentes) + código huérfano
```

## Limitaciones conocidas (ver Discusión del manuscrito para el detalle completo)

- N=2 donantes en GSE161192 — resultados de esa cohorte específicamente son exploratorios, no confirmatorios.
- El score de Fase 4 es una heurística de expresión ponderada por riesgo genético, **no un método de inferencia causal** (no hay Mendelian randomization, mediación, ni fine-mapping tipo SuSiE).
- La colocalización (Fase 1c) usa coloc.abf estándar (un solo variante causal asumido por locus); el locus chr19 tiene 3 señales independientes reales (Fase 1b) — una limitación metodológica conocida, no oculta.
- `scripts/utils/load_config.R` existe pero **no está actualmente invocado** por los scripts de fase — cada script todavía hardcodea su propia ruta de proyecto. Queda como trabajo pendiente si se refactoriza el pipeline.
- NicheNet (inferencia direccional ligando-receptor) está instalado en un **entorno conda aislado** (`nichenet_env`), no en el entorno principal del pipeline, específicamente porque su dependencia final requiere ggplot2≥4.0 (incompatible con los gráficos existentes, que usan ggplot2 3.5.2). El script `scripts/phase3_neuronal_vulnerability/06_nichenet_ligand_activity.R` debe ejecutarse con ese entorno activado (`conda activate nichenet_env`, y siempre desde un directorio de trabajo distinto de la raíz del proyecto para evitar que `.Rprofile`/renv contaminen `.libPaths()` con paquetes de R 4.3.1 dentro de un entorno de R 4.3.3 — ver comentarios en el script). Resultado: OSM e IL17A rankean entre los últimos de 29 ligandos candidatos por actividad dirigida, revirtiendo la interpretación mecanística previa del eje OSM–LIFR (ver manuscrito, Fase 3c).

## Autoría

**Ana [Apellido]** — [Afiliación institucional completa]

Ver la sección de declaración de IA en el manuscrito para el rol de herramientas asistidas por IA en el proceso de auditoría y análisis.

## Licencia

Ver [`LICENSE`](LICENSE).
