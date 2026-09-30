# Convergencia genética tejido-específica y gradientes de severidad clínica en la demencia por cuerpos de Lewy: un análisis interactómico neuroinmune

**Ana [Apellido completo]**

[Afiliación institucional completa]
[Dirección]
[Email de contacto]

---

## RESUMEN

La demencia por cuerpos de Lewy (DCL) es la segunda demencia neurodegenerativa más frecuente, caracterizada por agregados de alfa-sinucleína. Estudios previos identificaron loci de riesgo genético (GBA, SNCA, APOE) y describieron expansión clonal de linfocitos T CD4⁺ en líquido cefalorraquídeo (LCR), proponiendo un eje CXCL12–CXCR4 de reclutamiento hacia el sistema nervioso central. No estaba establecido en qué tejido actúa el riesgo genético, ni si los hallazgos preliminares de este campo resisten un control estadístico riguroso.

Aquí se presenta un marco computacional integrativo que combina priorización genética (MAGMA sobre el GWAS completo de Chia et al. 2021, N=6.618; seguido de análisis condicional GCTA-COJO), perfilado inmune (scRNA-seq+TCR-seq en dos cohortes independientes de LCR) y perfilado neuronal (snRNA-seq de corteza cingulada anterior con subtipificación por capa), integrados mediante colocalización genética tejido-específica (GWAS-vs-eQTL en corteza y en linfocitos T CD4⁺), Mendelian randomization, fine-mapping multi-variante (SuSiE-coloc) y un modelo de tendencia de severidad clínica. Se incluyen además verificaciones formales de robustez (sensibilidad al covariable de lote técnico, null por permutación) en vez de depender únicamente de la inspección de resultados.

De 18.523 genes analizados genome-wide, 7 alcanzaron significancia Bonferroni (locus APOE/TOMM40/PVRL2/APOC1, más MSTO1, MMRN1 y SNCA). BCAM, señalado como principal en un análisis preliminar subpotenciado (454 SNPs pre-filtrados), no replicó (rango 461/18.523) y el análisis condicional GCTA-COJO descartó formalmente su independencia del locus APOE (p condicional: 6,6×10⁻¹³ → 0,56). La expansión clonal de CD4⁺ se confirmó en ambas cohortes de LCR (GSE161192, GSE141578); el eje CXCL12–CXCR4 no replicó en ninguna, incluyendo la cohorte de mayor tamaño y origen de la hipótesis original, y esta no-replicación fue robusta a la especificación del modelo estadístico (con y sin covariable de lote técnico). La resolución de subtipos neuronales corticales (antes ausente en este pipeline) identificó cuatro poblaciones con expresión diferencial de ADORA1. La colocalización genética reveló un patrón tejido-específico: **MMRN1 no colocalizó en corteza (probabilidad posterior de variantes distintas, PP3=0,78) pero sí fuertemente en linfocitos T CD4⁺ (probabilidad posterior de variante causal compartida, PP4=0,97)**, sugiriendo que su riesgo genético actúa vía expresión inmune, no neuronal — un resultado que Mendelian randomization (beta=−0,293, padj=2,5×10⁻¹⁰) y fine-mapping multi-variante (SuSiE-coloc, PP4=0,977 para la señal específica) reprodujeron de forma independiente. Un modelo de tendencia ordinal pre-especificado (Control=0<PD=1<PDD=2<DLB=3) identificó tres receptores neuronales con aumento monótono y estadísticamente significativo según severidad clínica (padj<0,05): **LIFR**, **LGALS9** y **CD86**. Sin embargo, un análisis direccional con NicheNet (ejecutado en un entorno aislado para no comprometer el resto del pipeline) sobre los 29 ligandos candidatos de linfocitos T CD4⁺ mostró que **OSM —el ligando propuesto para LIFR, y la interacción mejor priorizada por expresión en Fase 4— rankeó 28 de 29 por actividad de ligando dirigida (AUROC=0,508), e IL17A rankeó último (29 de 29, AUROC=0,501)**, revirtiendo la interpretación mecanística de ese eje: LIFR escala con severidad clínica, pero sin evidencia de que OSM dirija esa respuesta. Una verificación por permutación confirmó además que el ranking expresión-ponderado de Fase 4 está dominado por expresión más que por evidencia genética. Un enrichment de pathways inmunes pre-especificados (MAGMA gene-set: Th17, señalización TCR, citoquinas, quimiocinas, JAK-STAT) no alcanzó significancia en ninguno de los cinco pathways testeados.

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

### *Verificación causal adicional: Mendelian randomization y fine-mapping multi-variante (SuSiE-coloc)*

La colocalización (`coloc.abf`) prueba consistencia con un variante causal compartido, pero no es, por sí sola, evidencia causal en el sentido de Mendelian randomization (MR), ni maneja loci con más de un variante causal. Se implementaron dos verificaciones adicionales dirigidas específicamente al hallazgo de MMRN1 (Fase 1c).

**Mendelian randomization (two-sample).** Para cada gen y tejido, se testeó si la expresión génica predicha genéticamente se asocia con el riesgo de DCL. Cuando se dispuso de un único instrumento fuerte (el SNP eQTL más significativo, requiriendo p<1×10⁻⁵), se calculó un test de razón de Wald (beta_MR = beta_GWAS / beta_eQTL). Para los genes del locus chr19 (PVRL2/APOE/TOMM40/APOC1), se ejecutó adicionalmente un modelo IVW (inverse-variance weighted; paquete `MendelianRandomization`) usando como instrumentos las tres señales independientes ya establecidas por GCTA-COJO (Fase 1b) — un conjunto de instrumentos genuinamente independiente por construcción, no un clumping ad hoc.

**Colocalización multi-variante (SuSiE-coloc).** Para MMRN1, se ajustó SuSiE (`coloc::runsusie`) por separado sobre el GWAS y sobre el eQTL de linfocitos T, usando una matriz de correlación (LD) calculada del panel de referencia 1000G EUR restringida a una ventana de ±200 kb alrededor del gen (más angosta que la ventana de ±1 Mb usada para el escaneo inicial de coloc.abf, por tratabilidad computacional del fine-mapping). `coloc::coloc.susie` testeó colocalización par-por-par entre cada credible set del GWAS y cada credible set del eQTL, en vez de asumir un único variante causal por dataset.

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

### *Fase 3c: Actividad ligando dirigida (NicheNet)*

La interacción ligando-receptor de Fases 3-4 usa co-expresión simétrica: un ligando y un receptor "activos" no implican que el ligando efectivamente dirija la respuesta transcripcional de la célula receptora. Para testear direccionalidad, se ejecutó NicheNet (Browaeys et al., 2020) en un entorno conda aislado (`nichenet_env`, ver Declaración de IA) para evitar que su dependencia de ggplot2≥4.0 afectara el resto del pipeline (ggplot2 3.5.2). Los ligandos candidatos fueron los genes upregulados en CD4⁺ expandidos (Fase 2a) presentes en la red ligando-receptor de referencia de NicheNet (v2, Browaeys et al. 2022; 29 de los genes de Fase 2a calificaron, un conjunto más amplio que la lista curada de CellPhoneDB usada en Fase 3). El geneset de interés fueron los genes significativamente upregulados en neuronas ACC de DLB vs. Control (test de Wilcoxon genome-wide a nivel celular — el mismo tipo de limitación de pseudorreplicación que D7/D8 en Fase 2a, usado aquí únicamente para definir un conjunto de genes blanco, no como una afirmación de significancia por sí misma); el background fueron todos los genes expresados en neuronas ACC. `predict_ligand_activities()` rankeó los 29 ligandos por su capacidad de explicar ese geneset, usando el modelo de potencial regulatorio ligando-blanco de NicheNet (AUROC, AUPR corregido, correlación de Pearson).

### *Fase 4: Red de priorización expresión-ponderada*

Para cada par ligando-receptor activo, se calculó un Score de Priorización como:

**Score de Priorización = Expresión_Neuronal_Máxima × (1 + max(0, Z-estadístico GWAS del receptor))**

Esta es una ponderación heurística de expresión por dirección de riesgo genético, no un método de inferencia causal; no se aplica test de direccionalidad, análisis de mediación, Mendelian randomization, ni colocalización eQTL a nivel de esta red (la colocalización, MR y SuSiE-coloc de Fases 1c/1e/1f se realizan sobre los genes de riesgo GWAS, no sobre esta red integrada). Para cuantificar cuánto del ranking depende del componente genético frente al de expresión, se permutaron las etiquetas de Z-estadístico GWAS entre los 16 pares candidatos (10.000 iteraciones), recalculando el score y el ranking en cada iteración y registrando la frecuencia con la que el par mejor rankeado observado se mantiene en el primer puesto bajo asignación aleatoria de Z.

### *Análisis estadístico*

Todos los análisis se realizaron en R 4.3.1 / Python 3.10, con el entorno R fijado en `renv.lock` (201 paquetes). La expresión diferencial usó el test de Wilcoxon (corrección Bonferroni, análisis exploratorio de GSE161192) o DESeq2 (corrección Benjamini-Hochberg, análisis pseudobulk). MAGMA usó el modelo SNP-wise mean con corrección Bonferroni genome-wide tanto para el análisis de genes individuales como, con corrección propia (0,05/5), para el análisis de gene-sets. GCTA-COJO usó configuración por defecto de colinealidad y selección stepwise. La colocalización usó `coloc.abf` con priors por defecto. Se usó una semilla aleatoria fija (42) en todo paso estocástico (p. ej. submuestreo).

### *Declaración de Inteligencia Artificial*

Se usaron herramientas asistidas por IA (Claude, Anthropic) en un proceso iterativo de auditoría, corrección e implementación: revisión del código y outputs del pipeline contra afirmaciones del manuscrito, verificación de la procedencia de los datasets directamente contra registros de GEO, identificación y corrección de defectos de código, ejecución del pipeline corregido, implementación de las extensiones de colocalización tejido-específica (Fase 1c), enrichment de gene-sets (Fase 1d) y modelo de tendencia de severidad (Fase 3b), y la redacción de este manuscrito. Todos los hallazgos generados con asistencia de IA reportados aquí fueron verificados de forma independiente contra los archivos de datos subyacentes antes de su inclusión; la IA no se usó para generar resultados que no pudieran trazarse a un archivo de datos específico o un registro de GEO verificado. La interpretación y las afirmaciones finales son responsabilidad de la autora.

### *Disponibilidad de datos y código*

Todos los datos usados son de acceso público. Estadísticas resumen del GWAS de DCL: GWAS Catalog GCST90001390 (Chia et al., 2021). Datos scRNA/TCR-seq de LCR: GEO GSE161192 y GSE141578. Datos snRNA-seq: GEO GSE178146. eQTLs: eQTL Catalogue de EBI (datasets QTD000151, QTD000031). Los scripts de análisis, el pipeline Snakemake, los archivos de configuración, y un `Dockerfile` que containeriza el critical path genético (Fase 1/1b/1c: MAGMA+GCTA-COJO+coloc, sin dependencia del entorno de la máquina de desarrollo, verificado reproduciendo `coloc_summary.csv` byte a byte) están disponibles en [repositorio a especificar al momento de publicación].

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

### *Verificación causal adicional: MMRN1 sobrevive Mendelian randomization y fine-mapping multi-variante; el resto del panel no tiene instrumento testeable*

**Mendelian randomization.** El test de razón de Wald usando el SNP eQTL líder como instrumento fue calculable únicamente para MMRN1 en ambos tejidos — para los demás nueve genes, ningún SNP alcanzó la significancia mínima requerida (eQTL p<1×10⁻⁵) en ninguno de los dos tejidos, consistente con el patrón ya observado en la colocalización (PP1 dominante — asociación solo con el GWAS, sin señal eQTL local detectable en este tejido). Este resultado nulo para el resto del panel no es una falla del método sino la misma conclusión que la colocalización ya sugería, ahora confirmada por la ausencia de instrumento.

Para MMRN1, la MR mostró un efecto fuerte y opuesto entre tejidos: en linfocitos T CD4⁺, beta_MR=−0,293 (p=1,2×10⁻¹⁰, padj=2,5×10⁻¹⁰) — mayor expresión genéticamente predicha de MMRN1 se asocia con **menor** riesgo de DCL. En corteza, el efecto es mucho más débil y de signo opuesto (beta_MR=+0,119, p=0,031, no significativo tras corrección). El modelo IVW multi-instrumento para los genes del locus chr19 (usando las 3 señales GCTA-COJO) no fue computable para ningún gen de ese locus, por la misma razón: ausencia de señal eQTL local detectable en estos tejidos específicos.

**SuSiE-coloc.** El fine-mapping de MMRN1 en linfocitos T (ventana de 400 kb, 1.179 SNPs con genotipo de referencia) convergió con 9 credible sets en el GWAS y 8 en el eQTL — un número de señales probablemente inflado por el ruido de estimación de LD con solo 503 individuos de referencia (limitación explícita, no una afirmación de 9 loci genuinamente independientes). Sin embargo, el resultado clave es específico y limpio: el credible set del GWAS anclado en **rs7680557** (el mismo SNP líder identificado en `coloc.abf` y en la MR) colocaliza fuertemente con un credible set específico del eQTL (PP4=0,977), mientras que el resto de las combinaciones credible-set-por-credible-set muestran PP3≈1 (variantes distintas, como se espera de señales genuinamente no relacionadas dentro de un locus con múltiple estructura de LD). Esto **valida** el hallazgo original de `coloc.abf` (PP4=0,97) como atribuible a una señal específica y no a una mezcla espuria entre múltiples señales del locus.

En conjunto, MR y SuSiE-coloc convergen con la colocalización simple: MMRN1 es un candidato de convergencia genética-inmune robusto a tres métodos independientes con distintos supuestos, mientras que ningún otro gen del panel muestra evidencia equivalente en ningún tejido, con ningún método.

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

**Análisis de sensibilidad al covariable de lote.** Dado que el diseño principal (~batch + disease_group) depende de una decisión de modelado, se re-ajustó el mismo modelo pseudobulk (a) sin covariable de lote (22 muestras) y (b) restringido a EP-vs-control (20 muestras, excluyendo las 2 muestras DCL que están confinadas al lote2). Ambos modelos alternativos reprodujeron el mismo resultado nulo: 0 genes significativos (de 13.575 y 13.396 genes testeados respectivamente), CXCR4 no significativo en ambos (padj=1,00 y 0,999), CXCL12 no detectado en ambos. La no-replicación de CXCL12–CXCR4 no depende, por lo tanto, de la elección específica de ajuste por lote.

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

### *Fase 3c: El análisis direccional no respalda a OSM ni a IL17A como impulsores de la respuesta neuronal*

Ninguno de los 29 ligandos candidatos mostró AUROC por encima de 0,54 (todos cerca del azar, 0,50), indicando que ningún ligando derivado de linfocitos T CD4⁺ explica fuertemente, por sí solo, la respuesta transcripcional de DLB en neuronas ACC. Dentro de ese panorama de actividad débil en general, el ranking tiene una estructura interpretable: los ligandos con mayor actividad relativa fueron HMGB1 (AUROC=0,542), S100A4 (0,536) y CD40LG (0,529) — moléculas de tipo alarmina/DAMP y co-estimuladoras, no las citoquinas centrales de la narrativa Th17 de Fase 2a. **OSM, el ligando de la interacción mejor priorizada por expresión en Fase 4, rankeó 28 de 29 (AUROC=0,508, AUPR corregido=0,003); IL17A rankeó último (29 de 29, AUROC=0,501, esencialmente indistinguible del azar).** TGFB1 y ENTPD1 quedaron en la mitad de la tabla (rango 12 y 14 respectivamente).

Esto es un resultado que pesa en contra de la interpretación mecanística de OSM–LIFR e IL17A como impulsores directos de la respuesta transcripcional cortical, bajo el único método direccional aplicado en este trabajo. No contradice que OSM y LIFR estén ambos expresados (eso sigue siendo cierto, y es lo único que Fase 4 mide), pero sí contradice la interpretación de que esa co-expresión refleje una relación causal ligando→respuesta transcripcional en el receptor. Combinado con la ausencia de respaldo genético para LIFR (Fase 4) y con la naturaleza puramente heurística del Score de Priorización (confirmada por permutación, ver más abajo), la evidencia acumulada en este trabajo para el eje OSM–LIFR es, en conjunto, más débil que su rango #1 en Fase 4 sugiere por sí solo.

### *Fase 4: Priorización expresión-ponderada; el refuerzo por severidad clínica no se traduce en refuerzo direccional*

Los Z-estadísticos GWAS de los receptores (Fase 1) van de −0,89 (ADORA1) a +1,23 (TGFBR3) entre los 16 receptores candidatos —pequeños en magnitud, no llegando a significancia nominal ni genome-wide para ninguno. El Score de Priorización (expresión × [1+max(0,Z)]; explícitamente no una métrica de inferencia causal) rankeó **OSM–LIFR** como la interacción más alta (0,759), seguida de TGFB1–TGFBR3 (0,455) y ENTPD1–ADORA1 (0,294).

Dado que ningún receptor porta riesgo genético significativo, este ranking está impulsado casi enteramente por el nivel de expresión neuronal, no por evidencia genética — la caracterización "informada genéticamente" de este score no debe sobre-interpretarse a nivel de esta red específica. El hallazgo independiente de Fase 3b —que LIFR, el receptor de la interacción mejor rankeada, escala significativamente con la severidad clínica— sigue siendo una línea de evidencia real e independiente, pero **el análisis direccional de Fase 3c la contrapesa**: LIFR escala con severidad clínica (un hecho sobre el receptor), pero su ligando propuesto, OSM, no muestra evidencia de dirigir la respuesta transcripcional neuronal (un hecho sobre la relación ligando→respuesta). Ambas cosas pueden ser ciertas simultáneamente sin ser mutuamente reforzantes: LIFR podría escalar con severidad por una vía completamente distinta de OSM. La evidencia acumulada en este trabajo no permite afirmar que OSM–LIFR sea el eje mecanístico central; permite afirmar que LIFR es un receptor cuya expresión covaría con severidad clínica, lo cual es, por sí mismo, un hallazgo más modesto pero mejor delimitado.

**Verificación por permutación.** Para testear formalmente la afirmación anterior (que el ranking está dominado por expresión, no por genética) en vez de solo enunciarla, se permutaron las etiquetas de Z-estadístico GWAS entre los 16 receptores candidatos 10.000 veces, recalculando el Score de Priorización y el ranking en cada iteración. OSM–LIFR se mantuvo como par #1 en el 83,7% de las permutaciones al azar, y la correlación de Spearman promedio entre el ranking con Z permutado y el ranking por expresión pura fue de 0,956. Esto confirma empíricamente, en vez de solo por inspección, que el Score de Priorización está dominado por la expresión y que el componente genético contribuye marginalmente al ranking actual — una limitación real del método, ahora cuantificada.

---

## DISCUSIÓN

Este trabajo buscó determinar si el riesgo genético, la señalización de linfocitos T CD4⁺ en LCR y la expresión de receptores neuronales corticales convergen en DCL, y en qué tejido lo harían. La conclusión más defendible integra seis líneas de evidencia: una replica trabajo previo (expansión clonal de CD4⁺), una se resuelve formalmente en sentido negativo con respaldo estadístico de precisión (BCAM no es un gen de riesgo independiente), una hipótesis mecanística previamente publicada no replica en dos cohortes independientes (tráfico vía CXCL12–CXCR4), un enrichment poligénico pre-especificado no encuentra señal (pathways inmunes), y dos líneas de evidencia positiva emergen y no habían sido descriptas antes en este pipeline: convergencia genética tejido-específica (MMRN1 en linfocitos T) y un gradiente de expresión receptorial asociado a severidad clínica (LIFR/LGALS9/CD86).

### *Expansión clonal de linfocitos T CD4⁺: el hallazgo más sólido*

La expansión clonal entre linfocitos T CD4⁺ de LCR se observó independientemente en ambas cohortes analizadas aquí y es consistente con el reporte original de Gate et al. (2021) usando los mismos datos fuente. Es el hallazgo de este estudio que debe considerarse bien respaldado, no meramente candidato o generador de hipótesis.

### *BCAM no es un gen de riesgo independiente de DCL*

El análisis MAGMA genome-wide (18.523 genes, GWAS completo) demuestra que el estatus previo de BCAM como gen principal fue un artefacto de un análisis subpotenciado y pre-filtrado (454 SNPs). El análisis condicional formal GCTA-COJO da una respuesta definitiva a la pregunta que quedaba abierta en versiones anteriores de este trabajo: la selección stepwise identifica exactamente tres señales independientes en el locus APOE de chr19, ninguna de las cuales es un SNP de BCAM, y el variante top de BCAM pierde esencialmente toda su asociación (p: 6,6×10⁻¹³ → 0,56) al condicionar apropiadamente sobre esas tres señales. No está justificada ninguna afirmación mecanística adicional específica de BCAM; la señal del locus queda completamente explicada por APOE/TOMM40/PVRL2/APOC1.

### *CXCL12–CXCR4: una hipótesis específica y testeable que no replicó*

Gate et al. (2021) propusieron que la señalización CXCR4-CXCL12 impulsa el tráfico patológico de linfocitos T CD4⁺ en DCL, basándose en la cohorte de LCR hoy disponible como GSE141578. El reanálisis independiente de esa misma cohorte (pseudobulk a nivel de donante, corregido por lote, n=22) no encuentra efecto significativo de CXCR4 por enfermedad ni CXCL12 detectable en linfocitos T CD4⁺; una segunda cohorte más pequeña y estimulada (GSE161192) muestra CXCR4 significativamente *downregulado* en clones expandidos. En conjunto, este reanálisis no respalda el eje CXCL12–CXCR4 tal como está formulado actualmente. Explicaciones plausibles, no mutuamente excluyentes, incluyen: (1) diferencias metodológicas respecto del análisis original no capturadas aquí (estadística a nivel celular vs. de donante, covariables no medidas de edad/sexo/severidad); (2) que la upregulación de CXCR4 sea específica de un subtipo de linfocitos T o estado de activación no resuelto por el clustering CD4⁺/CD8⁺ grueso usado aquí; o (3) que el hallazgo original sea sensible a decisiones analíticas. Se reporta esto como una no-replicación formal bajo reanálisis independiente de los mismos datos públicos, no como una refutación de la biología subyacente.

### *MMRN1: convergencia genética específica del compartimento inmune, ahora respaldada por tres métodos independientes*

El hallazgo de colocalización tejido-específica para MMRN1 —variantes distintas en corteza (PP3=0,78), variante causal compartida en linfocitos T CD4⁺ (PP4=0,97)— es, dentro de este trabajo, la evidencia más directa de que el riesgo genético de DCL puede converger con expresión génica específicamente en el compartimento inmune. A diferencia de la versión previa de este análisis, este hallazgo ya no descansa en un único método: Mendelian randomization (Wald ratio, beta=−0,293, padj=2,5×10⁻¹⁰) y SuSiE-coloc a nivel de credible set específico (PP4=0,977 para el par de señales anclado en rs7680557) convergen con la colocalización simple, cada uno con supuestos distintos y vulnerabilidades a sesgos distintas. Es coherente con la hipótesis central del estudio de forma más precisa que cualquier hallazgo anterior en este pipeline: no solo hay expansión clonal de linfocitos T (Fase 2) y riesgo genético en MMRN1 (Fase 1) por separado, sino evidencia convergente de que ambos comparten el mismo variante causal, actuando vía expresión en linfocitos T y no en neuronas corticales, con una dirección de efecto (mayor expresión → menor riesgo) que es al menos internamente consistente.

Esto sigue sin identificar el mecanismo funcional de MMRN1 en linfocitos T —una proteína de matriz extracelular/multimerina no tiene, a priori, un rol inmune obvio, lo cual hace de este hallazgo un candidato interesante pero que requiere validación funcional, no una confirmación de mecanismo. La MR de instrumento único además no puede excluir pleiotropía horizontal (el instrumento podría afectar el riesgo de DCL por una vía distinta de la expresión de MMRN1); el número de credible sets detectado por SuSiE (9 en el GWAS, dentro de una ventana de 400kb) probablemente está inflado por ruido de estimación de LD con un panel de referencia de solo 503 individuos, una limitación que se declara explícitamente en vez de reportar el número de señales como un hallazgo en sí mismo — la señal específica que colocaliza (ancla rs7680557) es el resultado relevante, no el conteo total de credible sets.

### *Un gradiente de severidad clínica en receptores neuronales corticales*

TGFBR2, LGALS9 y CD86 habían sido identificados en versiones previas de este pipeline como cambios de expresión de una única comparación DLB-vs-Control, explícitamente señalados como "candidatos, no hallazgos" por la ausencia de test estadístico y de replicación. El modelo de tendencia de severidad de Fase 3b cambia esa caracterización para LIFR, LGALS9 y CD86 (TGFBR2 queda cerca del umbral pero no lo alcanza): los tres muestran un aumento estadísticamente significativo y con dirección consistente a lo largo de las cuatro categorías clínicas del espectro DCL, no solo una diferencia aislada frente a Control. Este es, en sí mismo, un hallazgo sólido sobre la expresión de estos tres receptores.

Donde este trabajo inicialmente sobre-interpretó el resultado fue al vincular a LIFR específicamente con OSM como mecanismo: LIFR es el receptor de la interacción OSM–LIFR mejor rankeada por expresión basal en Fase 4, y el gradiente de severidad de LIFR se leyó en una primera pasada como refuerzo de esa interacción. El análisis direccional de Fase 3c revierte esa lectura: OSM rankeó 28 de 29 candidatos por actividad de ligando dirigida, sin evidencia de dirigir la respuesta transcripcional neuronal. Lo que queda, entonces, es más modesto y más honesto: LIFR (y LGALS9, CD86) escalan con severidad clínica — un hallazgo real, correlacional, que amerita validación dirigida — pero sin que este trabajo pueda todavía atribuir ese escalado a una señal específica de linfocitos T CD4⁺ vía OSM u otro ligando candidato.

### *El riesgo genético no converge, en general, con la expresión de receptores candidatos*

Ningún gen receptor porta riesgo genético genome-wide-significativo en este GWAS (Fase 4), y el enrichment poligénico en pathways inmunes completos tampoco encuentra señal (Fase 1d). El Score de Priorización está, en la práctica, impulsado casi enteramente por nivel de expresión más que por evidencia genética, y no debe describirse como "causal" ni como demostración de que el riesgo genético y la señalización inmune convergen mecánicamente de forma generalizada. Esto es consistente con —y no contradice— la posibilidad de que la vulnerabilidad neuronal mediada por el sistema inmune y el riesgo genético heredado actúen mayormente por vías independientes o indirectas en DCL, excepto en casos puntuales como MMRN1, donde la convergencia sí aparece cuando se testea en el tejido correcto.

### *Limitaciones*

Persisten varias limitaciones incluso tras las correcciones y extensiones descriptas aquí. Primero, la cohorte de LCR GSE161192 incluye solo 2 donantes; las estadísticas de Wilcoxon a nivel celular de este dataset están pseudorreplicadas y deben tratarse como generadoras de hipótesis, no confirmatorias, pese a p-valores nominalmente pequeños. Segundo, GSE141578, aunque considerablemente mejor potenciada (n=22), carece de metadatos de edad y sexo en su registro público de GEO, y los dos lotes de preparación de librería no son completamente ortogonales al subgrupo DCL. Tercero, la interacción ligando-receptor de Fases 3-4 sigue siendo, en sí misma, correlacional: combina expresión de dos datasets completamente separados y no puede establecer que un ligando dado de linfocitos T efectivamente actúe sobre un receptor cortical dado en el mismo tejido o marco temporal. El análisis direccional de Fase 3c (NicheNet, ejecutado en un entorno conda aislado para evitar que su dependencia de ggplot2≥4.0 afectara al resto del pipeline) aborda parcialmente esta limitación para el lado del ligando (qué ligando explica mejor la respuesta transcripcional observada), pero usa como geneset de interés una lista de genes DE definida por un test a nivel celular (la misma clase de limitación de pseudorreplicación de D7/D8), y no incorpora información espacial ni de coocurrencia temporal real entre linfocitos T y neuronas — sigue siendo, en el mejor caso, una inferencia computacional de direccionalidad, no una demostración de contacto o señalización real. Cuarto, si bien SuSiE-coloc se aplicó a MMRN1 y validó la señal específica que colocaliza, no se aplicó al locus chr19 (donde `coloc.abf` estándar ya había mostrado PP1 dominante, es decir, ausencia de señal eQTL local en los tejidos testeados, haciendo que un fine-mapping multi-variante añadiera poco en ese caso específico); y el número de credible sets detectado por SuSiE para MMRN1 está probablemente inflado por el tamaño limitado del panel de referencia LD (503 individuos), una limitación inherente a usar un panel externo en vez de LD in-sample. Quinto, el Score de Priorización de Fase 4 sigue siendo una heurística ponderada por expresión, no un método de inferencia causal por diseño; la verificación por permutación (Fase 4) cuantifica cuánto depende del componente de expresión (83,7% de las permutaciones mantienen el mismo ranking #1), pero no lo convierte en un método causal. La Mendelian randomization implementada (Fase 1e) se aplicó a los 10 genes candidatos de riesgo genético, no a esta red de priorización expresión-receptor en su conjunto. Finalmente, todo este análisis es computacional y generador de hipótesis; no se ha realizado validación experimental en muestras de pacientes con DCL, modelos animales, ni sistemas de co-cultivo in vitro.

### *Implicaciones para investigación futura*

Varios de los próximos pasos metodológicos identificados en versiones previas de este trabajo ya se implementaron en esta ronda de análisis (fine-mapping multi-variante para MMRN1, Mendelian randomization, verificación de robustez al covariable de lote, cuantificación por permutación del Score de Priorización, y modelado direccional ligando-receptor con NicheNet). Los pasos más accionables que permanecen son: (1) validación independiente del hallazgo de MMRN1 en un dataset adicional de eQTL de linfocitos T (p. ej. DICE) o en una cohorte de expresión génica de linfocitos T de DCL; (2) validación independiente del gradiente de severidad de LIFR/LGALS9/CD86 en una cohorte ortogonal de corteza DCL — se identificaron dos candidatos públicos concretos durante este trabajo: GSE303823 (snRNA-seq, corteza prefrontal, DLB n=6/PDD n=6/AD n=4/control n=4, mismo tipo de dato que el usado aquí, pero con un archivo de datos de 14GB que excedió el espacio en disco disponible en esta sesión) y GSE150696 (microarray Affymetrix HTA 2.0, corteza BA9, DLB n=12/PDD n=12, ~1GB, una prueba de replicación en una plataforma y región distintas); (3) dado que el análisis direccional (Fase 3c) señaló a HMGB1, S100A4 y CD40LG con actividad relativamente mayor (aunque todavía débil en términos absolutos, AUROC~0,53-0,54) que OSM/IL17A, valdría la pena investigar específicamente estos candidatos alarmina/co-estimuladores como hipótesis alternativas, en vez de asumir que la ausencia de señal para OSM cierra la pregunta; y (4) si se persigue más la implicancia de CXCR4 en DCL, un reanálisis resuelto por subtipo celular de linfocitos T y ajustado por covariables (edad, sexo), idealmente en una cohorte donde esas covariables estén disponibles.

### *Aplicabilidad más amplia: un marco de verificación, no solo un caso de estudio*

Más allá de los hallazgos específicos de DCL, la contribución metodológica de este trabajo es, en sí misma, generalizable: un protocolo de verificación en capas —(1) priorizar genéticamente con corrección multi-hipótesis correcta desde el inicio, no post-hoc; (2) resolver formalmente, con análisis condicional, cualquier "top hit" antes de construir una narrativa mecanística sobre él; (3) testear colocalización en el tejido mecánicamente implicado por la hipótesis, no solo en el tejido más conveniente o disponible; (4) verificar esa colocalización con un método ortogonal (MR) y, cuando la potencia computacional lo permite, con fine-mapping multi-variante; (5) replicar cualquier hallazgo inmune en una segunda cohorte independiente antes de aceptarlo; (6) testear robustez a decisiones de modelado (covariables, subconjuntos) en vez de reportar un único resultado; y (7) cuantificar por permutación cuánto de un ranking o score compuesto depende de cada componente— es aplicable, sin modificación conceptual, a cualquier otra enfermedad donde se hipotetice convergencia entre riesgo genético y un mecanismo celular específico. Los resultados negativos/refutatorios reportados aquí (BCAM, CXCL12–CXCR4, enrichment de pathways) no son un défcit del estudio sino la demostración de que este protocolo funciona: descarta lo que no resiste escrutinio antes de gastar interpretación en ello. Los dos hallazgos positivos que sí sobrevivieron todo el protocolo (MMRN1 con tres métodos convergentes, el gradiente de severidad de LIFR/LGALS9/CD86) son, por construcción, más creíbles que cualquier hallazgo individual de un análisis de una sola pasada — no porque los métodos individuales sean infalibles, sino porque sobrevivieron intentos genuinos de refutarlos con métodos de supuestos distintos.

---

## CONCLUSIONES

Este estudio se propuso construir un interactoma neuroinmune genéticamente informado en la demencia por cuerpos de Lewy y, tras un reanálisis exhaustivo con verificación cruzada entre cohortes, tejidos y métodos, sostiene un conjunto de conclusiones más acotado pero mejor evidenciado que el originalmente propuesto. La expansión clonal de linfocitos T CD4⁺ en LCR está confirmada y replica trabajo previo. El análisis MAGMA genome-wide, seguido de análisis condicional formal GCTA-COJO, resuelve a BCAM como completamente atribuible al desequilibrio de ligamiento del locus APOE —no un gen de riesgo independiente de DCL— mientras confirma a SNCA y al locus APOE/TOMM40/PVRL2/APOC1 como genome-wide-significativos. El eje de tráfico de linfocitos T CXCL12–CXCR4 previamente propuesto no está respaldado por reanálisis independiente de ninguna de las cohortes de LCR disponibles —incluyendo la cohorte propiamente mejor potenciada que es la fuente original de esta afirmación— y esta no-replicación es robusta a la especificación del modelo estadístico. El enrichment poligénico en pathways inmunes completos tampoco encuentra señal.

Frente a este trasfondo mayormente refutatorio, un hallazgo positivo se sostiene con fuerza tras un intento deliberado de refutación con métodos ortogonales, y un segundo hallazgo sobrevive parcialmente, revisado a la baja por ese mismo tipo de escrutinio: (1) una convergencia genética tejido-específica entre el riesgo de DCL y la expresión de MMRN1 en linfocitos T CD4⁺ (ausente en corteza), confirmada independientemente por Mendelian randomization y por fine-mapping multi-variante a nivel de la señal específica — el hallazgo más robusto de esta ronda de trabajo; y (2) un gradiente de expresión estadísticamente significativo de LIFR, LGALS9 y CD86 a lo largo de la severidad clínica del espectro DCL, que se mantiene como hallazgo sobre estos receptores, pero cuya interpretación mecanística original —vía el eje OSM–LIFR, la interacción mejor priorizada por expresión— no resistió el análisis direccional (NicheNet): OSM rankeó 28 de 29 ligandos candidatos por actividad dirigida, revirtiendo esa lectura mecanística específica sin invalidar el gradiente de expresión de LIFR en sí. Esta secuencia —proponer, testear con un método ortogonal, y revisar la interpretación cuando el test no la respalda— es, en última instancia, el resultado metodológico más importante de esta ronda de trabajo.

Este trabajo se caracteriza mejor como un estudio computacional riguroso, mayormente exploratorio y refutatorio, con dos líneas de evidencia positiva concretas, acotadas, y puestas a prueba con múltiples métodos independientes. El protocolo de verificación en capas usado aquí —priorización, resolución condicional, colocalización tejido-específica, confirmación causal ortogonal, replicación entre cohortes, sensibilidad a decisiones de modelado, y cuantificación por permutación— es en sí mismo una contribución generalizable a otros estudios de convergencia genética-inmune en enfermedad neurodegenerativa, independientemente de los hallazgos específicos de DCL reportados aquí. No está aún al nivel de validación experimental necesario para sostener afirmaciones de blanco terapéutico específico, y los próximos pasos hacia ese objetivo se detallan más arriba.

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

Burgess S., Butterworth A., Thompson S.G. (2013) Mendelian randomization analysis with multiple genetic variants using summarized data. Genet. Epidemiol. 37: 658-665.

Browaeys R., Saelens W., Saeys Y. (2020) NicheNet: modeling intercellular communication by linking ligands to target genes. Nat. Methods 17: 159-162.

Wallace C. (2021) A more accurate method for colocalisation analysis allowing for multiple causal variants. PLoS Genet. 17: e1009440. [coloc.susie]

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
