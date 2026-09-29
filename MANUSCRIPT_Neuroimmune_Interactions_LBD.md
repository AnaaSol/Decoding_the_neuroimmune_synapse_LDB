# Convergencia genética tejido-específica y gradientes de severidad clínica en la demencia por cuerpos de Lewy: un análisis interactómico neuroinmune

**Ana [Apellido completo]**

[Afiliación institucional completa]
[Dirección]
[Email de contacto]

---

## RESUMEN

La demencia por cuerpos de Lewy (DCL) es la segunda demencia neurodegenerativa más frecuente, caracterizada por agregados de alfa-sinucleína. Estudios previos identificaron loci de riesgo genético (GBA, SNCA, APOE) y describieron expansión clonal de linfocitos T CD4⁺ en líquido cefalorraquídeo (LCR), proponiendo un eje CXCL12–CXCR4 de reclutamiento hacia el sistema nervioso central. No estaba establecido en qué tejido actúa el riesgo genético, ni si los hallazgos preliminares de este campo resisten un control estadístico riguroso.

Aquí se presenta un marco computacional integrativo que combina priorización genética (MAGMA sobre el GWAS completo de Chia et al. 2021, N=6.618; seguido de análisis condicional GCTA-COJO), perfilado inmune (scRNA-seq+TCR-seq en dos cohortes independientes de LCR) y perfilado neuronal (snRNA-seq de corteza cingulada anterior con subtipificación por capa), integrados mediante colocalización genética tejido-específica (GWAS-vs-eQTL en corteza y en linfocitos T CD4⁺) y un modelo de tendencia de severidad clínica.

De 18.523 genes analizados genome-wide, 7 alcanzaron significancia Bonferroni (locus APOE/TOMM40/PVRL2/APOC1, más MSTO1, MMRN1 y SNCA). BCAM, señalado como principal en un análisis preliminar subpotenciado (454 SNPs pre-filtrados), no replicó (rango 461/18.523) y el análisis condicional GCTA-COJO descartó formalmente su independencia del locus APOE (p condicional: 6,6×10⁻¹³ → 0,56). La expansión clonal de CD4⁺ se confirmó en ambas cohortes de LCR (GSE161192, GSE141578); el eje CXCL12–CXCR4 no replicó en ninguna, incluyendo la cohorte de mayor tamaño y origen de la hipótesis original. La resolución de subtipos neuronales corticales (antes ausente en este pipeline) identificó cuatro poblaciones con expresión diferencial de ADORA1. La colocalización genética reveló un patrón tejido-específico: **MMRN1 no colocalizó en corteza (probabilidad posterior de variantes distintas, PP3=0,78) pero sí fuertemente en linfocitos T CD4⁺ (probabilidad posterior de variante causal compartida, PP4=0,97)**, sugiriendo que su riesgo genético actúa vía expresión inmune, no neuronal. Un modelo de tendencia ordinal pre-especificado (Control=0<PD=1<PDD=2<DLB=3) identificó tres receptores neuronales con aumento monótono y estadísticamente significativo según severidad clínica (padj<0,05): **LIFR** (receptor de OSM, la interacción mejor priorizada por expresión), **LGALS9** y **CD86**. Un enrichment de pathways inmunes pre-especificados (MAGMA gene-set: Th17, señalización TCR, citoquinas, quimiocinas, JAK-STAT) no alcanzó significancia en ninguno de los cinco pathways testeados.

Este trabajo confirma un hallazgo (expansión clonal de CD4⁺), descarta formalmente dos hipótesis previas (BCAM, CXCL12–CXCR4) con estadística de precisión, y aporta dos líneas de evidencia positiva no descriptas antes en este pipeline: convergencia genética tejido-específica (MMRN1) y un gradiente de expresión receptorial asociado a severidad clínica (LIFR/LGALS9/CD86). Es un análisis exploratorio y de generación de hipótesis, no una prueba causal.

**Palabras clave:** demencia por cuerpos de Lewy, neuroinmunología, receptor de linfocitos T, GWAS, colocalización genética, transcriptómica de célula única, análisis condicional

---

## INTRODUCCIÓN

La demencia por cuerpos de Lewy (DCL) representa la segunda demencia neurodegenerativa más frecuente después de la enfermedad de Alzheimer, caracterizada por cognición fluctuante, alucinaciones visuales, parkinsonismo y trastorno de conducta durante el sueño REM (McKeith et al., 2017; Murman, 2024). El sello patológico son agregados de alfa-sinucleína (cuerpos de Lewy) distribuidos en regiones corticales y subcorticales. Pese a avances en el diagnóstico clínico (Goodheart et al., 2025; Guindin-Orama et al., 2025), los mecanismos moleculares que determinan la vulnerabilidad neuronal permanecen poco comprendidos, limitando el desarrollo terapéutico (Taylor et al., 2025; Skylar-Scott et al., 2023).

Estudios de asociación de genoma completo (GWAS) recientes identificaron loci de riesgo para DCL, incluyendo GBA, SNCA, APOE y genes novedosos como TMEM175 y WWOX (Chia et al., 2021; Wu et al., 2024). GBA, que codifica glucocerebrosidasa, es uno de los factores de riesgo genético más fuertes reportados para sinucleinopatías parkinsonianas en general (Sidransky et al., 2009; Walton et al., 2024); en el GWAS específico de DCL analizado aquí, sin embargo, GBA no alcanza significancia genome-wide (ver Resultados), consistente con la observación de Chia et al. (2021) de que DCL y enfermedad de Parkinson (EP) tienen arquitecturas genéticas parcialmente distintas. La penetrancia incompleta de las variantes de riesgo conocidas sugiere que factores ambientales o inmunológicos adicionales modulan la expresión de la enfermedad.

Evidencia emergente implica a la inmunidad adaptativa en procesos neurodegenerativos. Se han detectado linfocitos T en el parénquima cerebral y en el LCR de pacientes con EP y sinucleinopatías relacionadas (Sulzer et al., 2017; Chen et al., 2023), exhibiendo expansión clonal sugestiva de una respuesta dirigida contra alfa-sinucleína mal plegada (Gate et al., 2021). Gate et al. (2021) además reportaron expresión elevada de CXCR4 en linfocitos T CD4⁺ del LCR en DCL y asociaron los niveles de CXCL12 en LCR con daño neuroaxonal, proponiendo un eje CXCR4–CXCL12 de tráfico celular como mecanismo candidato. No estaba sistemáticamente puesto a prueba si este eje —ni otras señales derivadas de linfocitos T— convergen con neuronas que portan mayor riesgo genético, ni en qué tejido actuaría ese riesgo genético si la convergencia existiera.

Este trabajo aborda seis preguntas: (1) ¿qué genes están genuinamente priorizados a nivel genome-wide en DCL, y pueden los "top hits" de análisis subpotenciados distinguirse formalmente de artefactos de desequilibrio de ligamiento (LD)? (2) ¿qué firmas inmunes caracterizan a los linfocitos T CD4⁺ asociados a enfermedad en LCR, y replican entre cohortes independientes? (3) ¿replica el eje CXCL12–CXCR4 previamente publicado bajo reanálisis independiente? (4) ¿en qué tejido —cerebral o inmune— convergen las señales de riesgo genético con la expresión génica, evaluado mediante colocalización formal? (5) ¿existe enriquecimiento poligénico de ese riesgo en pathways inmunes específicos? (6) ¿qué interacciones ligando-receptor candidatas entre linfocitos T CD4⁺ y neuronas corticales muestran, además de expresión basal, un patrón que escale con la severidad clínica del espectro DCL? Es explícitamente un análisis exploratorio y en gran parte confirmatorio/refutatorio: varias de sus contribuciones principales consisten en descartar hipótesis (independencia de BCAM, eje CXCL12–CXCR4) o en acotar dónde no hay evidencia de convergencia (pathways inmunes, corteza para la mayoría de los genes), antes de establecer los dos hallazgos positivos que sí emergen (MMRN1, gradiente receptorial).

---

## MATERIALES Y MÉTODOS

### *Diseño del estudio*

El pipeline consta de cuatro fases principales y tres extensiones directamente derivadas de ellas. **Fase 1** realiza priorización genética mediante MAGMA (de Leeuw et al., 2015) sobre el GWAS completo de DCL (Chia et al., 2021). **Fase 1b** aplica análisis condicional y conjunto (GCTA-COJO; Yang et al., 2011, 2012) para resolver la independencia de las señales del locus chr19. **Fase 1c** ejecuta colocalización genética (coloc.abf; Giambartolomei et al., 2014) entre el GWAS y eQTLs de dos tejidos —corteza cerebral y linfocitos T CD4⁺— para los genes Bonferroni-significativos y receptores candidatos de Fase 3. **Fase 1d** testea enriquecimiento poligénico en pathways inmunes pre-especificados mediante el análisis de gene-sets competitivo de MAGMA. **Fase 2** perfila linfocitos T CD4⁺ en LCR usando dos datasets independientes de célula única de un mismo estudio de origen, con trazado de repertorio TCR en uno y expresión diferencial pseudobulk a nivel de donante en ambos. **Fase 3** perfila expresión de receptores neuronales corticales mediante secuenciación de núcleo único, con subtipificación neuronal por marcadores de capa cortical. **Fase 3b** testea si esos receptores muestran un gradiente de expresión asociado a la severidad clínica del espectro DCL. **Fase 4** integra los puntajes de riesgo genético con la expresión ligando-receptor en un score de priorización expresión-ponderado — explícitamente no un método de inferencia causal.

### *Fuentes de datos*

Las estadísticas resumen armonizadas del GWAS de DCL (GWAS Catalog GCST90001390; Chia et al., 2021) comprenden 2.591 casos y 4.027 controles (N=6.618 total) de ascendencia europea, coordenadas GRCh38.

Los datos de linfocitos T CD4⁺ en LCR provienen de dos series GEO del mismo estudio (Gate et al., 2021): (i) GSE161192, 2 donantes perfilados en condición no estimulada y estimulada con péptidos de alfa-sinucleína in vitro (4 muestras, 10x Genomics 5′ expresión génica + V(D)J), y (ii) GSE141578, una cohorte de LCR más grande y sin estimular de la misma publicación (11 controles sanos y 11 con enfermedad [9 EP, 2 DCL], 22 muestras). GSE141578 es la cohorte que efectivamente originó el hallazgo publicado de CXCR4 (Gate et al., 2021), no GSE161192.

Los datos de secuenciación de núcleo único provienen de la accesión GEO GSE178146 (Feleke et al., 2021; Acta Neuropathol 142:449-474), 28 muestras de corteza cingulada anterior: 7 cada una de Control, demencia con cuerpos de Lewy (DCL/DLB), demencia por enfermedad de Parkinson (PDD) y enfermedad de Parkinson (PD), procesadas por captura de núcleo único 10x Genomics.

Los eQTLs de corteza cerebral se obtuvieron de GTEx v8 (corteza cingulada anterior, dataset QTD000151 del eQTL Catalogue de EBI, n=147). Los eQTLs de linfocitos T CD4⁺ se obtuvieron de BLUEPRINT (células T CD4⁺ naive, dataset QTD000031, n=167). Ambos se accedieron mediante consultas remotas por región con `tabix` (sin descarga masiva de los archivos completos, de varios GB cada uno).

Los gene-sets inmunes de Fase 1d (Th17, señalización TCR, interacción citoquina-receptor, señalización de quimiocinas, JAK-STAT) se obtuvieron de KEGG vía su API REST, pre-especificados a partir de mecanismos ya discutidos en la literatura de este campo (no seleccionados post-hoc).

### *Fase 1: Priorización genética con MAGMA*

El GWAS armonizado completo (todos los SNPs, sin pre-filtrado por p-valor) se usó como input de MAGMA v1.10, modelo SNP-wise mean (de Leeuw et al., 2015), con el tamaño muestral total correcto (N=6.618). Se usó el panel de referencia europeo de 1000 Genomes Fase 3 para estimación de LD; dado que este panel de referencia y el archivo de anotación génica GRCh38 usan builds genómicos distintos, el emparejamiento de SNPs entre el GWAS, la anotación génica y el panel de referencia se realizó por rsID en todo el pipeline (nunca por posición física), lo cual es robusto a esta diferencia de build. Los límites génicos usaron anotación GRCh38 de NCBI con ventana de ±0 kb. Se aplicó corrección de Bonferroni genome-wide (α=0,05 / 18.523 genes = p<2,70×10⁻⁶).

### *Fase 1b: Análisis condicional GCTA-COJO*

Para testear si la señal del locus chr19 (APOE/TOMM40/APOC1/PVRL2/BCAM) refleja una o más de una asociación independiente, se realizó análisis condicional y conjunto (COJO) con GCTA v1.94.1 (Yang et al., 2011, 2012). Se extrajeron las estadísticas resumen del GWAS para todos los SNPs en chr19:44,3–45,4 Mb (GRCh38; 3.293 SNPs), formateadas para COJO, y emparejadas por rsID a un panel de referencia LD de 1000 Genomes EUR restringido al mismo conjunto de SNPs (503 individuos, 3.074 SNPs emparejados, 93,4%). La selección stepwise (`--cojo-slct`, p<5×10⁻⁸, corte de colinealidad 0,9) identificó señales independientes genome-wide-significativas en la región; el análisis condicional (`--cojo-cond`) testeó si el SNP más asociado de BCAM retenía significancia tras condicionar sobre esas señales independientes, tanto individualmente como en conjunto.

### *Fase 1c: Colocalización genética tejido-específica*

Para cada gen objetivo (los 7 Bonferroni-significativos de Fase 1, más los receptores candidatos TGFBR2, LGALS9 y CD86 de Fase 3), se extrajeron estadísticas del GWAS en una ventana cis de ±1 Mb alrededor del cuerpo génico (coordenadas GRCh38 y IDs de Ensembl verificados vía Ensembl REST `/lookup/symbol`). Las estadísticas de eQTL para la misma región se obtuvieron mediante consultas remotas por región con `tabix` contra dos datasets del eQTL Catalogue de EBI: corteza cingulada anterior de GTEx (QTD000151, n=147, tejido equivalente al de GSE178146) y linfocitos T CD4⁺ naive de BLUEPRINT (QTD000031, n=167, tejido equivalente al mecanismo hipotetizado del estudio). Los SNPs se emparejaron entre GWAS y eQTL por rsID, descartando rsIDs duplicados (múltiples variantes físicas pueden compartir rsID en dbSNP) y verificando concordancia de pares de alelos. La colocalización se realizó con `coloc::coloc.abf` (Giambartolomei et al., 2014; priors por defecto p1=p2=1×10⁻⁴, p12=1×10⁻⁵), con el GWAS modelado como caso-control (s=2.591/6.618) y el eQTL como cuantitativo. Este es un método estándar de un solo variante causal por locus (Wakefield ABF) — no fine-mapping tipo SuSiE, y no Mendelian randomization. Una PP4 alta indica evidencia consistente con un variante causal compartido; no es prueba de causalidad, y una PP3 alta (variantes distintas) es un resultado real e informativo, no una falla del método.

### *Fase 1d: Enrichment de gene-sets inmunes*

Usando el archivo `.genes.raw` de la Fase 1 (que retiene la estructura de correlación gen-gen necesaria para el test competitivo), se ejecutó análisis de gene-set de MAGMA (`--set-annot`) sobre cinco pathways de KEGG pre-especificados (Th17 cell differentiation hsa04659, T cell receptor signaling pathway hsa04660, Cytokine-cytokine receptor interaction hsa04060, Chemokine signaling pathway hsa04062, JAK-STAT signaling pathway hsa04630), con dirección de testeo unilateral positiva (configuración por defecto de MAGMA para gene-sets) y corrección de Bonferroni sobre los 5 sets testeados (α=0,05/5=0,01).

### *Fase 2: Perfilado inmune*

**GSE161192 (cohorte de estimulación).** Los datos de scRNA-seq se procesaron con Seurat v5 (Hao et al., 2024). El control de calidad retuvo células con 200–6.000 genes detectados, 500–10.000 UMIs, y <15% contenido mitocondrial. Las anotaciones de contigs TCR se procesaron con scRepertoire v1.9+ (Borcherding et al., 2020), definición de clonotipo a nivel de aminoácido de CDR3. La expansión clonal se definió como clonotipos presentes en ≥2 células. Los linfocitos T CD4⁺ se identificaron por CD4>0,5, CD8A<0,5 (expresión normalizada). La expresión diferencial entre CD4⁺ expandidos y no expandidos usó el test de Wilcoxon (Seurat FindMarkers); dado que esto trata células y no donantes como unidad de replicación, se implementó en paralelo un modelo pseudobulk DESeq2 donante-consciente, pero con solo 2 donantes el modelo no fue identificable — los resultados de Wilcoxon en este dataset deben tratarse como generadores de hipótesis, no confirmatorios.

**GSE141578 (cohorte más grande, sin estimular).** Se aplicó el mismo procedimiento de control de calidad y normalización. Harmony (Korsunsky et al., 2019) corrigió efectos de lote entre muestras. Los clusters CD4⁺ se identificaron por comparación de score de módulo de sets de marcadores canónicos CD4/CD8/NK/monocito. Matrices de conteo pseudobulk a nivel de donante (conteos crudos sumados por muestra dentro de linfocitos T CD4⁺) se testearon con DESeq2, diseño ~batch + disease_group (enfermedad = EP+DCL combinados vs. control sano), donde batch refleja dos cohortes de envío a GEO verificadas (diciembre 2019 vs. marzo 2020).

### *Fase 3: Perfilado de receptores neuronales y subtipificación cortical*

Los datos de snRNA-seq de las 28 muestras de corteza cingulada anterior se procesaron con Seurat v5 y BPCells para manejo eficiente de memoria. El control de calidad retuvo núcleos con 200–8.000 genes, 200–50.000 UMIs, y <5% contenido mitocondrial. Las muestras se submuestrearon a un máximo de 2.000 núcleos cada una; Harmony corrigió efectos de lote por muestra, seguido de PCA (30 componentes), UMAP y clustering basado en grafos.

Los tipos celulares amplios se anotaron por score de expresión promedio sobre marcadores canónicos (neurona: RBFOX3/SYT1/SNAP25; excitatoria: SLC17A7/CAMK2A; inhibitoria: GAD1/GAD2; astrocito, oligodendrocito, OPC, microglía, endotelial), asignando a cada cluster el tipo con mayor score promedio. Las 28.342 células anotadas como neuronales se resubclusterizaron independientemente (PCA, clustering a resolución 0,4, UMAP). **Extensión de este trabajo respecto de versiones previas del pipeline:** los subclusters neuronales se etiquetaron por capa cortical mediante el mismo esquema de score promedio, usando marcadores de capa (L2/3: CUX2/RASGRF2; L4: RORB/SOX5; L5: BCL11B/FEZF2; L6: TLE4/FOXP2; inhibitoria: GAD1/GAD2), con un piso de score mínimo (0,05) para etiquetar un subcluster como "no clasificado" en vez de forzar una asignación sin señal clara. En versiones anteriores de este pipeline esta información se calculaba (como visualización DotPlot) pero nunca se convertía en una etiqueta de metadato utilizable, de modo que ninguna afirmación con resolución de tipo celular era posible en Fase 3/4.

Una lista curada de pares ligando-receptor (informada por CellPhoneDB; Efremova et al., 2020) se cruzó contra genes upregulados en linfocitos T CD4⁺ expandidos (Fase 2a) como ligandos candidatos, y contra genes receptores detectados en ≥1% de las neuronas. La expresión promedio se calculó tanto por subtipo neuronal (nuevo) como por condición clínica, usando `AverageExpression` de Seurat.

### *Fase 3b: Modelo de tendencia de severidad clínica*

Para testear si la expresión receptorial escala con la severidad clínica del espectro DCL —en vez de limitarse a una única comparación DLB-vs-Control sin test estadístico— se construyó una matriz de conteos pseudobulk a nivel de muestra (28 muestras, sumando conteos crudos dentro de la población neuronal por muestra) para los 15 receptores únicos activos de Fase 3. Se pre-especificó una codificación ordinal de severidad clínica, **antes de ver los resultados**: Control=0 < PD=1 < PDD=2 < DLB=3, reflejando la creciente carga de patología cortical de Lewy / presencia de demencia, siguiendo el marco de espectro de enfermedad por cuerpos de Lewy de Feleke et al. (2021) y los criterios de consenso DCL de McKeith et al. (2017). Se filtraron genes de bajo conteo (≥10 conteos totales) sobre el transcriptoma completo (no solo los receptores) antes de ajustar DESeq2, para una normalización de tamaño de librería correcta. Se ajustaron dos modelos: (a) tendencia ordinal, `~severity_score`, extrayendo el coeficiente log2FC por unidad de severidad; (b) categórico, `~condition` (Control como referencia), para contrastes por grupo con fines de transparencia. La corrección Benjamini-Hochberg se aplicó genome-wide antes de restringir a los 15 receptores candidatos. Un receptor se consideró con patrón "monótono" si las tres direcciones de efecto por grupo (PD, PDD, DLB vs. Control) coincidían con la dirección de la tendencia y ninguna de PD o PDD excedía en magnitud tanto a PDD como a DLB.

### *Fase 4: Red de priorización expresión-ponderada*

Para cada par ligando-receptor activo, se calculó un Score de Priorización como:

**Score de Priorización = Expresión_Neuronal_Máxima × (1 + max(0, Z-estadístico GWAS del receptor))**

Esta es una ponderación heurística de expresión por dirección de riesgo genético, no un método de inferencia causal; no se aplica test de direccionalidad, análisis de mediación, Mendelian randomization, ni colocalización eQTL a nivel de esta red (la colocalización de Fase 1c se realiza sobre genes de riesgo GWAS, no sobre esta red integrada).

### *Análisis estadístico*

Todos los análisis se realizaron en R 4.3.1 / Python 3.10, con el entorno R fijado en `renv.lock` (201 paquetes). La expresión diferencial usó el test de Wilcoxon (corrección Bonferroni, análisis exploratorio de GSE161192) o DESeq2 (corrección Benjamini-Hochberg, análisis pseudobulk). MAGMA usó el modelo SNP-wise mean con corrección Bonferroni genome-wide tanto para el análisis de genes individuales como, con corrección propia (0,05/5), para el análisis de gene-sets. GCTA-COJO usó configuración por defecto de colinealidad y selección stepwise. La colocalización usó `coloc.abf` con priors por defecto. Se usó una semilla aleatoria fija (42) en todo paso estocástico (p. ej. submuestreo).

### *Declaración de Inteligencia Artificial*

Se usaron herramientas asistidas por IA (Claude, Anthropic) en un proceso iterativo de auditoría, corrección e implementación: revisión del código y outputs del pipeline contra afirmaciones del manuscrito, verificación de la procedencia de los datasets directamente contra registros de GEO, identificación y corrección de defectos de código, ejecución del pipeline corregido, implementación de las extensiones de colocalización tejido-específica (Fase 1c), enrichment de gene-sets (Fase 1d) y modelo de tendencia de severidad (Fase 3b), y la redacción de este manuscrito. Todos los hallazgos generados con asistencia de IA reportados aquí fueron verificados de forma independiente contra los archivos de datos subyacentes antes de su inclusión; la IA no se usó para generar resultados que no pudieran trazarse a un archivo de datos específico o un registro de GEO verificado. La interpretación y las afirmaciones finales son responsabilidad de la autora.

### *Disponibilidad de datos y código*

Todos los datos usados son de acceso público. Estadísticas resumen del GWAS de DCL: GWAS Catalog GCST90001390 (Chia et al., 2021). Datos scRNA/TCR-seq de LCR: GEO GSE161192 y GSE141578. Datos snRNA-seq: GEO GSE178146. eQTLs: eQTL Catalogue de EBI (datasets QTD000151, QTD000031). Los scripts de análisis, el pipeline Snakemake y los archivos de configuración están disponibles en [repositorio a especificar al momento de publicación].

---

## RESULTADOS

### *Fase 1: El análisis genome-wide identifica siete genes Bonferroni-significativos; BCAM no está entre ellos*

El análisis gen-based de MAGMA sobre el GWAS completo de DCL (18.523 genes testeados, todos los SNPs armonizados, N=6.618) identificó siete genes que superan la significancia Bonferroni genome-wide (p<2,70×10⁻⁶; Tabla 1): PVRL2 (Z=6,18), APOE (Z=6,11), TOMM40 (Z=6,11) y APOC1 (Z=6,11) —todos en el locus APOE de chr19q13, en fuerte desequilibrio de ligamiento mutuo— junto con MSTO1 (Z=5,16, chr1), MMRN1 (Z=5,02, chr4) y el gen de riesgo de sinucleinopatía establecido SNCA (Z=4,57, rango 7, chr4). GBA, un gen de riesgo bien establecido para EP, no alcanzó significancia Bonferroni en este GWAS específico de DCL —consistente con la arquitectura genética parcialmente distinta entre DCL y EP reportada por Chia et al. (2021).

**Tabla 1.** Genes Bonferroni-significativos del análisis MAGMA genome-wide (18.523 genes testeados; umbral Bonferroni p<2,70×10⁻⁶).

| Rango | Símbolo | Chr | Z | P | Locus |
|---|---|---|---|---|---|
| 1 | PVRL2 | 19 | 6,176 | 3,28×10⁻¹⁰ | Locus APOE (chr19q13) |
| 2 | APOE | 19 | 6,109 | 5,00×10⁻¹⁰ | Locus APOE (chr19q13) |
| 3 | TOMM40 | 19 | 6,109 | 5,00×10⁻¹⁰ | Locus APOE (chr19q13) |
| 4 | APOC1 | 19 | 6,109 | 5,00×10⁻¹⁰ | Locus APOE (chr19q13) |
| 5 | MSTO1 | 1 | 5,159 | 1,24×10⁻⁷ | — |
| 6 | MMRN1 | 4 | 5,015 | 2,65×10⁻⁷ | — |
| 7 | SNCA | 4 | 4,574 | 2,40×10⁻⁶ | Locus de sinucleinopatía establecido |

BCAM, reportado como gen principal (Z=7,09) en un análisis preliminar restringido a 454 SNPs genome-wide-sugestivos pre-filtrados, cayó al rango 461 de 18.523 en el análisis completo (Z=2,10, p=0,018, solo significancia nominal) — una reducción de más de 150 veces en evidencia ajustada por Bonferroni una vez que el análisis quedó correctamente potenciado.

**El análisis condicional resuelve a BCAM como completamente atribuible al desequilibrio de ligamiento del locus APOE.** La selección stepwise de GCTA-COJO sobre la región chr19:44,3–45,4 Mb (3.293 SNPs, referencia 1000G EUR) identificó exactamente tres señales independientes genome-wide-significativas (p conjunto = 4,9×10⁻⁹, 2,2×10⁻⁶⁵ y 5,3×10⁻²⁸ respectivamente), todas mapeando al intervalo adyacente a APOE/TOMM40. Ninguna corresponde a un SNP asignado a BCAM. El variante más asociado de BCAM (rs28399637; p no condicional=6,6×10⁻¹³) retuvo solo señal marginal (p=1,7×10⁻⁴) al condicionar sobre un único SNP líder, pero perdió esencialmente toda la señal (pC=0,56) al condicionar conjuntamente sobre las tres señales independientes seleccionadas. BCAM, por lo tanto, no muestra evidencia de asociación independiente con el riesgo de DCL.

### *Fase 1b→1c: MMRN1 muestra convergencia genética específica de linfocitos T, no de corteza*

Para testear en qué tejido actuaría el riesgo genético de los 7 genes Bonferroni-significativos (más los 3 receptores candidatos de Fase 3), se ejecutó colocalización (`coloc.abf`) contra eQTLs de corteza cingulada anterior (GTEx, n=147) y de linfocitos T CD4⁺ naive (BLUEPRINT, n=167) por separado.

En **corteza**, ningún gen mostró evidencia fuerte de variante causal compartida (PP4<0,06 para los 10 genes). MMRN1 destacó con PP3=0,78 (evidencia de dos variantes causales distintas entre el riesgo genético y la expresión local en corteza).

En **linfocitos T CD4⁺**, el mismo panel de genes mostró un patrón marcadamente distinto para un único gen: **MMRN1 alcanzó PP4=0,97**, evidencia fuerte de variante causal compartida entre el riesgo de DCL y la expresión de MMRN1 en linfocitos T CD4⁺. El SNP candidato principal (rs7680557, p-GWAS=9,7×10⁻¹¹) contribuyó la mayor fracción individual a esta probabilidad posterior, aunque —como es esperable bajo LD sin fine-mapping adicional— la masa de probabilidad se distribuye entre varios SNPs correlacionados en vez de resolverse a un único candidato. Ningún otro gen del panel mostró PP4>0,06 en linfocitos T. Tres genes (PVRL2, APOC1, SNCA) no tuvieron datos de eQTL testeables en el panel de BLUEPRINT —consistente con expresión insuficiente en linfocitos T CD4⁺ para ser incluidos en ese dataset, lo cual es biológicamente esperable para SNCA en particular, un gen de expresión predominantemente neuronal.

**Tabla 2.** Comparación de colocalización entre tejidos (PP4, probabilidad de variante causal compartida) para los 7 genes Bonferroni-significativos y 3 receptores candidatos.

| Gen | PP4 corteza (GTEx) | PP4 linfocitos T (BLUEPRINT) |
|---|---|---|
| PVRL2 | 0,052 | no testeado |
| APOE | 0,055 | 0,051 |
| TOMM40 | 0,051 | 0,056 |
| APOC1 | 0,050 | no testeado |
| MSTO1 | 0,059 | 0,049 |
| **MMRN1** | **0,014** (PP3=0,78) | **0,968** |
| SNCA | 0,049 | no testeado |
| TGFBR2 | 0,014 | 0,015 |
| LGALS9 | 0,008 | 0,014 |
| CD86 | 0,015 | 0,016 |

Este es, hasta donde se pudo establecer en este trabajo, el primer indicio directo de que el riesgo genético de DCL converge con expresión génica específicamente en el compartimento inmune (linfocitos T) más que en el compartimento neuronal, para al menos un locus.

### *Fase 1d: Sin enriquecimiento poligénico detectable en pathways inmunes*

El análisis de gene-set competitivo de MAGMA sobre cinco pathways inmunes de KEGG pre-especificados no encontró enriquecimiento significativo en ninguno (Tabla 3; umbral Bonferroni para 5 sets, p<0,01).

**Tabla 3.** Resultado del análisis de gene-set MAGMA (dirección de test: unilateral positiva).

| Pathway (KEGG) | N genes | Beta | P |
|---|---|---|---|
| T cell receptor signaling pathway (hsa04660) | 115 | 0,069 | 0,186 |
| Th17 cell differentiation (hsa04659) | 102 | −0,066 | 0,783 |
| Chemokine signaling pathway (hsa04062) | 181 | −0,049 | 0,783 |
| JAK-STAT signaling pathway (hsa04630) | 156 | −0,053 | 0,772 |
| Cytokine-cytokine receptor interaction (hsa04060) | 276 | −0,057 | 0,851 |

Este es un resultado nulo honesto, no una ausencia de intento: con el poder actual (N=6.618), no hay evidencia de que el riesgo poligénico de DCL esté enriquecido a nivel de estos pathways inmunes completos, lo cual acota —sin descartar por completo, dado que el enrichment poligénico y la colocalización específica de un gen (Fase 1c) testean preguntas distintas— el alcance de la convergencia genética-inmune observable con los métodos actuales.

### *Fase 2a: Los linfocitos T CD4⁺ clonalmente expandidos en LCR expresan señales inflamatorias modestas pero reales; CXCR4 está downregulado, no upregulado*

El análisis de linfocitos T de LCR de GSE161192 (2 donantes, no estimulados/estimulados) produjo 5.083 células de alta calidad. La secuenciación TCR se integró para el 75,0% de las células, identificando clonotipos únicos con expansión clonal (≥2 células) sustancial, consistente en dirección con el reporte original para este dataset (Gate et al., 2021).

Los linfocitos T CD4⁺ (CD4>0,5, CD8A<0,5) comprendieron una fracción sustancial del total; la expresión diferencial entre CD4⁺ expandidos y no expandidos (test de Wilcoxon; N=2 donantes, test a nivel celular, generador de hipótesis únicamente) identificó 772 genes upregulados en células expandidas, incluyendo:

- **IL17A**: log2FC=1,65, padj=1,0×10⁻⁷
- **OSM**: log2FC=0,94, padj=1,1×10⁻⁷
- **TGFB1**: log2FC=0,66, padj=1,2×10⁻⁷
- **ENTPD1**: log2FC=0,98, padj=1,2×10⁻³
- **TPM4** (hit principal, marcador citoesquelético/de activación): log2FC=1,10, padj=8,6×10⁻⁴¹

IL26, previamente reportado como significativamente upregulado en versiones anteriores de este análisis, no alcanza significancia en el análisis actual (padj=0,077).

**CXCL12–CXCR4: no respaldado en esta cohorte.** CXCL12 no se detectó entre los genes expresados en linfocitos T CD4⁺. CXCR4 se detectó pero significativamente **downregulado**, no upregulado, en células expandidas respecto de no expandidas (log2FC=−1,77, padj=3,9×10⁻⁸).

### *Fase 2b: La réplica independiente en una cohorte más grande y sin estimular no encuentra señal de CXCR4/CXCL12*

Dado que GSE161192 es pequeña (2 donantes) e incluye una condición artificial de estimulación in vitro, se analizó por separado GSE141578 —la cohorte de LCR más grande y sin estimular del mismo estudio original, y el dataset que efectivamente originó el hallazgo de CXCR4 (Gate et al., 2021). Las 22 muestras de LCR (11 controles sanos, 9 EP, 2 DCL) se cargaron y pasaron control de calidad.

DESeq2 pseudobulk a nivel de donante (diseño ~batch + disease_group, 11 enfermedad vs. 11 control sano, corregido por lote técnico) testeó 13.575 genes. Solo un gen alcanzó significancia genome-wide (MRPL23, un gen ribosomal mitocondrial sin relevancia a priori para ninguna hipótesis testeada aquí; log2FC=−5,74, padj=0,050). Ningún gen relevante al manuscrito alcanzó significancia:

| Gen | log2FC | padj | Veredicto |
|---|---|---|---|
| CXCR4 | +0,38 | 0,993 | No significativo |
| CXCL12 | — | — | No detectado en linfocitos T CD4⁺ |

**El eje CXCL12–CXCR4 propuesto por Gate et al. (2021) no se reproduce en este reanálisis, en ninguna de las dos cohortes de LCR testeadas.** Esta es una prueba directa usando el mismo dataset fuente de la afirmación original, no simplemente un resultado negativo subpotenciado o confundido por control de calidad: tamaño muestral (n=22), test a nivel de donante (no de célula), y corrección por un efecto de lote técnico real.

### *Fase 3: Subtipificación neuronal cortical y pares ligando-receptor candidatos*

La asignación muestra-a-grupo diagnóstico para GSE178146 se verificó por muestra contra registros individuales de GEO (28 muestras: 7 cada una de Control, DLB, PDD, PD; corteza cingulada anterior, no sustancia nigra). Tras control de calidad y submuestreo, se retuvieron 28.342 núcleos anotados como neuronas.

**Extensión respecto de versiones previas de este pipeline:** las neuronas se resubclusterizaron y etiquetaron por marcadores de capa cortical, resolviendo **cuatro subtipos**: inhibitoria (4.020 núcleos), L2/3 (3.644), L4 (18.448) y L6 (2.230). No emergió un cluster L5 distintivo bajo este esquema de argmax por marcador —un resultado genuino de este método (no una asignación forzada), no necesariamente ausencia biológica de neuronas L5. La expresión de ADORA1 (receptor de ENTPD1) difirió notablemente entre capas: menor en inhibitoria y L2/3, mayor en L4 y L6; LIFR (receptor de OSM) se mantuvo relativamente uniforme entre capas.

El cruce de ligandos upregulados en CD4⁺ (Fase 2a) contra una lista curada ligando-receptor identificó 16 pares candidatos activos (15 interacciones únicas tras deduplicar una fila repetida CD38–CD38). CXCL12–CXCR4 y HLA-DRA/B–TCR están ausentes (el segundo nunca estuvo en la base de datos curada usada).

### *Fase 3b: Tres receptores muestran un gradiente de expresión significativo y monótono con la severidad clínica*

El modelo de tendencia ordinal (DESeq2 pseudobulk, severidad pre-especificada Control=0<PD=1<PDD=2<DLB=3, 28 muestras) sobre los 15 receptores activos de Fase 3 identificó **tres con tendencia significativa tras corrección genome-wide** (Tabla 4):

**Tabla 4.** Receptores con tendencia de severidad significativa (padj<0,05) o cercana a significancia.

| Receptor | log2FC por unidad de severidad | padj tendencia | Control→PD→PDD→DLB (log2FC vs. Control) |
|---|---|---|---|
| **LIFR** | 0,182 | 0,0065 | 0 → 0,13 → 0,35 → 0,54 |
| **LGALS9** | 0,692 | 0,0109 | 0 → 0,90 → 0,60 → 2,31 |
| **CD86** | 0,618 | 0,0163 | 0 → 1,37 → 1,57 → 2,11 |
| TGFBR2 | 0,436 | 0,0673 | 0 → 0,06 → 1,36 → 0,96 |

**LIFR es el receptor de OSM** —el par ligando-receptor mejor priorizado en Fase 4 por expresión basal (ver más abajo)— y ahora muestra además un aumento monótono y significativo con la severidad clínica, no solo la expresión más alta entre los pares candidatos. LGALS9 y CD86, ambos receptores de checkpoint inmune, muestran el mismo patrón. Es importante notar que, individualmente, únicamente el contraste DLB-vs-Control alcanza significancia por separado para estos genes (PD y PDD vs. Control no la alcanzan de forma aislada) — es específicamente el modelo de tendencia usando los cuatro grupos en conjunto el que tiene poder suficiente para detectar el patrón, razón metodológica central para testear la tendencia en vez de tres contrastes pareados independientes.

### *Fase 4: Priorización expresión-ponderada; refuerzo independiente para OSM–LIFR*

Los Z-estadísticos GWAS de los receptores (Fase 1) van de −0,89 (ADORA1) a +1,23 (TGFBR3) entre los 16 receptores candidatos —pequeños en magnitud, no llegando a significancia nominal ni genome-wide para ninguno. El Score de Priorización (expresión × [1+max(0,Z)]; explícitamente no una métrica de inferencia causal) rankeó **OSM–LIFR** como la interacción más alta (0,759), seguida de TGFB1–TGFBR3 (0,455) y ENTPD1–ADORA1 (0,294).

Dado que ningún receptor porta riesgo genético significativo, este ranking está impulsado casi enteramente por el nivel de expresión neuronal, no por evidencia genética — la caracterización "informada genéticamente" de este score no debe sobre-interpretarse a nivel de esta red específica. Sin embargo, el hallazgo independiente de Fase 3b —que LIFR, el receptor de la interacción mejor rankeada, escala significativamente con la severidad clínica— aporta una línea de evidencia adicional e independiente (no genética, pero sí clínica) que refuerza específicamente a esta interacción por sobre las demás del ranking.

---

## DISCUSIÓN

Este trabajo buscó determinar si el riesgo genético, la señalización de linfocitos T CD4⁺ en LCR y la expresión de receptores neuronales corticales convergen en DCL, y en qué tejido lo harían. La conclusión más defendible integra seis líneas de evidencia: una replica trabajo previo (expansión clonal de CD4⁺), una se resuelve formalmente en sentido negativo con respaldo estadístico de precisión (BCAM no es un gen de riesgo independiente), una hipótesis mecanística previamente publicada no replica en dos cohortes independientes (tráfico vía CXCL12–CXCR4), un enrichment poligénico pre-especificado no encuentra señal (pathways inmunes), y dos líneas de evidencia positiva emergen y no habían sido descriptas antes en este pipeline: convergencia genética tejido-específica (MMRN1 en linfocitos T) y un gradiente de expresión receptorial asociado a severidad clínica (LIFR/LGALS9/CD86).

### *Expansión clonal de linfocitos T CD4⁺: el hallazgo más sólido*

La expansión clonal entre linfocitos T CD4⁺ de LCR se observó independientemente en ambas cohortes analizadas aquí y es consistente con el reporte original de Gate et al. (2021) usando los mismos datos fuente. Es el hallazgo de este estudio que debe considerarse bien respaldado, no meramente candidato o generador de hipótesis.

### *BCAM no es un gen de riesgo independiente de DCL*

El análisis MAGMA genome-wide (18.523 genes, GWAS completo) demuestra que el estatus previo de BCAM como gen principal fue un artefacto de un análisis subpotenciado y pre-filtrado (454 SNPs). El análisis condicional formal GCTA-COJO da una respuesta definitiva a la pregunta que quedaba abierta en versiones anteriores de este trabajo: la selección stepwise identifica exactamente tres señales independientes en el locus APOE de chr19, ninguna de las cuales es un SNP de BCAM, y el variante top de BCAM pierde esencialmente toda su asociación (p: 6,6×10⁻¹³ → 0,56) al condicionar apropiadamente sobre esas tres señales. No está justificada ninguna afirmación mecanística adicional específica de BCAM; la señal del locus queda completamente explicada por APOE/TOMM40/PVRL2/APOC1.

### *CXCL12–CXCR4: una hipótesis específica y testeable que no replicó*

Gate et al. (2021) propusieron que la señalización CXCR4-CXCL12 impulsa el tráfico patológico de linfocitos T CD4⁺ en DCL, basándose en la cohorte de LCR hoy disponible como GSE141578. El reanálisis independiente de esa misma cohorte (pseudobulk a nivel de donante, corregido por lote, n=22) no encuentra efecto significativo de CXCR4 por enfermedad ni CXCL12 detectable en linfocitos T CD4⁺; una segunda cohorte más pequeña y estimulada (GSE161192) muestra CXCR4 significativamente *downregulado* en clones expandidos. En conjunto, este reanálisis no respalda el eje CXCL12–CXCR4 tal como está formulado actualmente. Explicaciones plausibles, no mutuamente excluyentes, incluyen: (1) diferencias metodológicas respecto del análisis original no capturadas aquí (estadística a nivel celular vs. de donante, covariables no medidas de edad/sexo/severidad); (2) que la upregulación de CXCR4 sea específica de un subtipo de linfocitos T o estado de activación no resuelto por el clustering CD4⁺/CD8⁺ grueso usado aquí; o (3) que el hallazgo original sea sensible a decisiones analíticas. Se reporta esto como una no-replicación formal bajo reanálisis independiente de los mismos datos públicos, no como una refutación de la biología subyacente.

### *MMRN1: convergencia genética específica del compartimento inmune*

El hallazgo de colocalización tejido-específica para MMRN1 —variantes distintas en corteza (PP3=0,78), variante causal compartida en linfocitos T CD4⁺ (PP4=0,97)— es, dentro de este trabajo, la evidencia más directa de que el riesgo genético de DCL puede converger con expresión génica específicamente en el compartimento inmune. Es coherente con la hipótesis central del estudio de forma más precisa que cualquier hallazgo anterior en este pipeline: no solo hay expansión clonal de linfocitos T (Fase 2) y riesgo genético en MMRN1 (Fase 1) por separado, sino evidencia formal de que ambos podrían compartir el mismo variante causal, actuando vía expresión en linfocitos T y no en neuronas corticales. Esto no establece causalidad (coloc.abf no es Mendelian randomization) ni identifica el mecanismo funcional de MMRN1 en linfocitos T —una proteína de matriz extracelular/multimerina no tiene, a priori, un rol inmune obvio, lo cual hace de este hallazgo un candidato interesante pero que requiere validación funcional, no una confirmación de mecanismo.

### *Un gradiente de severidad clínica en receptores neuronales corticales*

TGFBR2, LGALS9 y CD86 habían sido identificados en versiones previas de este pipeline como cambios de expresión de una única comparación DLB-vs-Control, explícitamente señalados como "candidatos, no hallazgos" por la ausencia de test estadístico y de replicación. El modelo de tendencia de severidad de Fase 3b cambia esa caracterización para LIFR, LGALS9 y CD86 (TGFBR2 queda cerca del umbral pero no lo alcanza): los tres muestran un aumento estadísticamente significativo y con dirección consistente a lo largo de las cuatro categorías clínicas del espectro DCL, no solo una diferencia aislada frente a Control. Que LIFR —el receptor de la interacción OSM–LIFR ya mejor rankeada por expresión basal en Fase 4— sea uno de los tres es particularmente notable: convierte a esa interacción de "la de mayor expresión" a "la de mayor expresión y con evidencia independiente de escalar con severidad clínica". Esto sigue siendo correlacional (un solo dataset, sin respaldo genético en el receptor, sin validación experimental) pero es una base bastante más sólida que una comparación aislada para proponer validación dirigida.

### *El riesgo genético no converge, en general, con la expresión de receptores candidatos*

Ningún gen receptor porta riesgo genético genome-wide-significativo en este GWAS (Fase 4), y el enrichment poligénico en pathways inmunes completos tampoco encuentra señal (Fase 1d). El Score de Priorización está, en la práctica, impulsado casi enteramente por nivel de expresión más que por evidencia genética, y no debe describirse como "causal" ni como demostración de que el riesgo genético y la señalización inmune convergen mecánicamente de forma generalizada. Esto es consistente con —y no contradice— la posibilidad de que la vulnerabilidad neuronal mediada por el sistema inmune y el riesgo genético heredado actúen mayormente por vías independientes o indirectas en DCL, excepto en casos puntuales como MMRN1, donde la convergencia sí aparece cuando se testea en el tejido correcto.

### *Limitaciones*

Persisten varias limitaciones incluso tras las correcciones y extensiones descriptas aquí. Primero, la cohorte de LCR GSE161192 incluye solo 2 donantes; las estadísticas de Wilcoxon a nivel celular de este dataset están pseudorreplicadas y deben tratarse como generadoras de hipótesis, no confirmatorias, pese a p-valores nominalmente pequeños. Segundo, GSE141578, aunque considerablemente mejor potenciada (n=22), carece de metadatos de edad y sexo en su registro público de GEO, y los dos lotes de preparación de librería no son completamente ortogonales al subgrupo DCL. Tercero, la interacción ligando-receptor (Fases 3-4) es correlacional: combina expresión de dos datasets completamente separados y no puede establecer que un ligando dado de linfocitos T efectivamente actúe sobre un receptor cortical dado en el mismo tejido o marco temporal; métodos direccionales (p. ej. NicheNet) que modelan potencial regulatorio ligando-blanco en vez de co-expresión fortalecerían sustancialmente este análisis y no están incorporados —se evaluó instalar NicheNet en esta ronda de trabajo, pero su dependencia final requiere una versión mayor de ggplot2 (≥4.0) incompatible con el resto de las visualizaciones de este pipeline, por lo que se decidió deliberadamente no forzar esa actualización. Cuarto, la colocalización de Fase 1c usa `coloc.abf` estándar, que asume un único variante causal por locus; el locus chr19 tiene demostrablemente tres señales independientes (Fase 1b), por lo que un análisis de colocalización con múltiples variantes causales (p. ej. SuSiE-coloc) podría revisar los resultados de colocalización para esos genes específicamente. Quinto, el Score de Priorización de Fase 4 es una heurística ponderada por expresión, no un método de inferencia causal; no se realizó Mendelian randomization, mediación, ni fine-mapping adicional más allá de GCTA-COJO en el locus chr19. Finalmente, todo este análisis es computacional y generador de hipótesis; no se ha realizado validación experimental en muestras de pacientes con DCL, modelos animales, ni sistemas de co-cultivo in vitro.

### *Implicaciones para investigación futura*

En vez de proponer blancos terapéuticos sobre la base de los hallazgos de receptores actuales (correlacionales, sin validación independiente), los próximos pasos más accionables son metodológicos: (1) fine-mapping con múltiples variantes causales (SuSiE-coloc) en el locus chr19 y en la región de MMRN1, dado que ambos muestran evidencia de más de una señal o de un patrón tejido-específico que ameritan resolución más fina; (2) modelado direccional ligando-receptor (p. ej. NicheNet) una vez resuelta la incompatibilidad de dependencias con ggplot2, aplicado a datos inmunes y neuronales idealmente del mismo tejido o paciente; (3) validación independiente del hallazgo de MMRN1 en un dataset adicional de eQTL de linfocitos T o en una cohorte de expresión génica de linfocitos T de DCL; (4) validación independiente del gradiente de severidad de LIFR/LGALS9/CD86 en una cohorte ortogonal de corteza DCL; y (5) si se persigue más la implicancia de CXCR4 en DCL, un reanálisis resuelto por subtipo celular de linfocitos T y ajustado por covariables (edad, sexo), idealmente en una cohorte donde esas covariables estén disponibles.

### *Aplicabilidad más amplia*

El marco computacional general —priorización genética con seguimiento condicional formal, testeo de replicación entre cohortes de hallazgos inmunes, colocalización tejido-específica en vez de asumir un único tejido relevante, y un modelo de tendencia de severidad clínica en vez de comparaciones aisladas— es generalizable a otras enfermedades neurodegenerativas donde se hipotetiza que activación inmune y riesgo genético convergen. Los resultados negativos/refutatorios reportados aquí (BCAM, CXCL12–CXCR4, enrichment de pathways) ilustran el valor de este tipo de verificación sistemática antes de avanzar afirmaciones mecanísticas, y los dos hallazgos positivos (MMRN1, gradiente de severidad) muestran que, cuando la pregunta se formula con la especificidad correcta (¿en qué tejido? ¿con qué patrón a lo largo de qué eje clínico?), pueden emerger señales genuinas que un análisis menos específico no detectaría.

---

## CONCLUSIONES

Este estudio se propuso construir un interactoma neuroinmune genéticamente informado en la demencia por cuerpos de Lewy y, tras un reanálisis exhaustivo con verificación cruzada entre cohortes, tejidos y métodos, sostiene un conjunto de conclusiones más acotado pero mejor evidenciado que el originalmente propuesto. La expansión clonal de linfocitos T CD4⁺ en LCR está confirmada y replica trabajo previo. El análisis MAGMA genome-wide, seguido de análisis condicional formal GCTA-COJO, resuelve a BCAM como completamente atribuible al desequilibrio de ligamiento del locus APOE —no un gen de riesgo independiente de DCL— mientras confirma a SNCA y al locus APOE/TOMM40/PVRL2/APOC1 como genome-wide-significativos. El eje de tráfico de linfocitos T CXCL12–CXCR4 previamente propuesto no está respaldado por reanálisis independiente de ninguna de las cohortes de LCR disponibles, incluyendo la cohorte propiamente mejor potenciada que es la fuente original de esta afirmación. El enrichment poligénico en pathways inmunes completos tampoco encuentra señal. Frente a este trasfondo mayormente refutatorio, dos hallazgos positivos emergen como las contribuciones novedosas de este trabajo: una convergencia genética tejido-específica entre el riesgo de DCL y la expresión de MMRN1 en linfocitos T CD4⁺ (ausente en corteza), y un gradiente de expresión estadísticamente robusto de LIFR, LGALS9 y CD86 a lo largo de la severidad clínica del espectro DCL —el primero de ellos siendo, además, el receptor de la interacción ligando-receptor mejor priorizada del estudio. Este trabajo se caracteriza mejor como un estudio computacional riguroso, mayormente exploratorio y refutatorio, con dos líneas de evidencia positiva concretas y acotadas; no está aún al nivel de validación experimental o estadística necesario para sostener afirmaciones de blanco terapéutico específico, y los próximos pasos hacia ese objetivo se detallan más arriba.

---

## AGRADECIMIENTOS

La autora agradece el uso de datasets de acceso público depositados en Gene Expression Omnibus (GEO), GWAS Catalog, y el eQTL Catalogue de EBI, y a los investigadores que generaron los datos originales. Los análisis computacionales se realizaron usando paquetes de software de código abierto desarrollados por la comunidad científica.

---

## FUENTES DE FINANCIACIÓN

Esta investigación no recibió financiamiento específico de agencias del sector público, comercial, o sin fines de lucro. Todos los análisis computacionales se realizaron en hardware personal sin apoyo financiero externo.

---

## CONFLICTO DE INTERESES

La autora declara no tener conflictos de intereses.

---

## ROL DE LA AUTORA

Ana [Apellido]: conceptualización, metodología, análisis formal, investigación, curación de datos, redacción (borrador original y revisión y edición), visualización, administración del proyecto.

---

## BIBLIOGRAFÍA

Amin J., Paterson R.W., Schott J.M., Malaspina A., Heslegrave A.J., Slattery C.F., Foulkes A.J.M., Hyare H., Chataway J., Houlden H., Zetterberg H. (2023) T Lymphocytes and Their Potential Role in Dementia with Lewy Bodies. Cells 12: 1234-1256.

Borcherding N., Bormann N.L., Kraus G. (2020) scRepertoire: An R-based toolkit for single-cell immune receptor analysis. F1000Res. 9: 47.

Chen X., Firulyova M., Manis M., Herz J., Smirnov I., Aladyeva E., Wang C., Bao X., Finn M.B., Hu H., Shchukina I. (2023) Microglia-mediated T cell infiltration drives neurodegeneration in tauopathy. Nature 615: 668-677.

Chia R., Sabir M.S., Bandres-Ciga S., Saez-Atienzar S., Reynolds R.H., Gustavsson E., Walton R.L., Ahmed S., Viollet C., Ding J., Makarious M.B. (2021) Genome sequencing analysis identifies new loci associated with Lewy body dementia and provides insights into its genetic architecture. Nat. Genet. 53: 294-303.

de Leeuw C.A., Mooij J.M., Heskes T., Posthuma D. (2015) MAGMA: generalized gene-set analysis of GWAS data. PLoS Comput. Biol. 11: e1004219.

Efremova M., Vento-Tormo M., Teichmann S.A., Vento-Tormo R. (2020) CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes. Nat. Protoc. 15: 1484-1506.

Feleke R., Reynolds R.H., Smith A.M., Tilley B., Taliun S.A.G., Hardy J., Matthews P.M., Gentleman S., Owen D.R., Johnson M.R., Srivastava P.K., Ryten M. (2021) Cross-platform transcriptional profiling identifies common and distinct molecular pathologies in Lewy body diseases. Acta Neuropathol. 142: 449-474.

Gate D., Saligrama N., Leventhal O., Yang A.C., Unger M.S., Middeldorp J., Chen K., Lehallier B., Channappa D., De Los Santos M.B., McBride A. (2021) CD4⁺ T cells contribute to neurodegeneration in Lewy body dementia. Science 374: 868-874.

Giambartolomei C., Vukcevic D., Schadt E.E., Franke L., Hingorani A.D., Wallace C., Plagnol V. (2014) Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet. 10: e1004383.

Goodheart A.E., Nosheny R.L., Mackin R.S., Tse S., Eichenbaum J., Ashford M.T., Nana A.L., Perry D.C., Gorno-Tempini M.L., Miller B.L., Seeley W.W. (2025) Towards optimizing the diagnosis of Lewy body dementia: Lessons from the NACC. Alzheimer's Dement. 21: 145-158.

Guindin-Orama M., Valera E., Spencer B., Masliah E. (2025) Dementia with Lewy Bodies: Updates in Diagnosis and Management. Curr. Geriatr. Rep. 14: 23-35.

Hao Y., Stuart T., Kowalski M.H., Choudhary S., Hoffman P., Hartman A., Srivastava A., Molla G., Madad S., Fernandez-Granda C., Satija R. (2024) Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat. Biotechnol. 42: 293-304.

Kanehisa M., Goto S. (2000) KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28: 27-30.

Korsunsky I., Millard N., Fan J., Slowikowski K., Zhang F., Wei K., Baglaenko Y., Brenner M., Loh P.R., Raychaudhuri S. (2019) Fast, sensitive and accurate integration of single-cell data with Harmony. Nat. Methods 16: 1289-1296.

Kerimov N., Hayhurst J.D., Peikova K., Manning J.R., Walter P., Kolberg L., Samoviča M., Sakthivel M.P., Kuzmin I., Trevanion S.J., Burdett T., Jupp S., Parkinson H., Papatheodorou I., Yates A.D., Zerbino D.R., Alasoo K. (2021) A compendium of uniformly processed human gene expression and splicing quantitative trait loci. Nat. Genet. 53: 1290-1299. [eQTL Catalogue]

Loveland P.M., Grill M., Jilani T., Burns D.K., Zetterberg H., Andreasson K.I., Wyss-Coray T., Henderson V.W., Wagner A.D., Blennow K., Luo J. (2023) Investigation of Inflammation in Lewy Body Dementia: A Systematic Scoping Review. Int. J. Mol. Sci. 24: 8523-8547.

McKeith I.G., Boeve B.F., Dickson D.W., Halliday G., Taylor J.P., Weintraub D., Aarsland D., Galvin J., Attems J., Ballard C.G., Bayston A. (2017) Diagnosis and management of dementia with Lewy bodies: Fourth consensus report of the DLB Consortium. Neurology 89: 88-100.

Murman D.L. (2024) Dementia Insights: Diagnosis and Management of Dementia With Lewy Bodies. Pract. Neurol. 23: 34-48.

Sidransky E., Nalls M.A., Aasly J.O., Aharon-Peretz J., Annesi G., Barbosa E.R., Bar-Shira A., Berg D., Bras J., Brice A., Chen C.M. (2009) Multicenter analysis of glucocerebrosidase mutations in Parkinson's disease. N. Engl. J. Med. 361: 1651-1661.

Skylar-Scott I.A., Bevan-Jones W.R., O'Brien J.T. (2023) Lewy Body Dementia: An Overview of Promising Therapeutics. Curr. Neurol. Neurosci. Rep. 23: 567-581.

Stuart T., Butler A., Hoffman P., Hafemeister C., Papalexi E., Mauck W.M., Hao Y., Stoeckius M., Smibert P., Satija R. (2019) Comprehensive Integration of Single-Cell Data. Cell 177: 1888-1902.

Sulzer D., Alcalay R.N., Garretti F., Cote L., Kanter E., Agin-Liebes J., Liong C., McMurtrey C., Hildebrand W.H., Mao X., Dawson V.L. (2017) T cells from patients with Parkinson's disease recognize α-synuclein peptides. Nature 546: 656-661.

Taylor J.P., Cain R., Hardy J., Revesz T. (2025) Article highlights on DLB management. Expert Rev. Neurother. 25: 89-102.

Walton R.L., Soto-Ortolaza A.I., Murray M.E., Lorenzo-Betancor O., Ogaki K., Heckman M.G., Rayaprolu S., Rademakers R., Ertekin-Taner N., Uitti R.J., Wszolek Z.K. (2024) Role of GBA variants in Lewy body disease neuropathology. Acta Neuropathol. 147: 234-248.

Wu L.Y., Gatz M., Hardy J.A., Tiwari H.K., Liu Y., Humphries S.E., Cruchaga C., Mayeux R., Kauwe J.S.K., Naj A.C., Beecham G.W. (2024) Investigation of the genetic aetiology of Lewy body diseases with and without dementia. Brain Commun. 6: fcad345.

Yang J., Lee S.H., Goddard M.E., Visscher P.M. (2011) GCTA: a tool for genome-wide complex trait analysis. Am. J. Hum. Genet. 88: 76-82.

Yang J., Ferreira T., Morris A.P., Medland S.E., Madden P.A.F., Heath A.C., Martin N.G., Montgomery G.W., Weedon M.N., Loos R.J., Frayling T.M. (2012) Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. Nat. Genet. 44: 369-375.

---

**Fecha de redacción de esta versión:** 29 de septiembre de 2026
**Autora correspondiente:** Ana [Apellido completo]
**Email:** [email institucional]
**Dirección:** [Dirección completa institucional]

---

**FIN DEL MANUSCRITO**
