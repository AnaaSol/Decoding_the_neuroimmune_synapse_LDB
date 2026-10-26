# ¿Dónde actúa el riesgo genético de la demencia por cuerpos de Lewy? Un atlas de colocalización resuelto por tipo celular en 47 conjuntos de eQTL inmunes y cerebrales

**Ana [Apellido completo]**

[Afiliación institucional completa]
[Email de contacto]

> Versión 2 (actual). Este trabajo reemplaza el enfoque de genes candidatos de la versión 1 (`MANUSCRIPT_Neuroimmune_Interactions_LBD.md`) por un análisis sistemático y no sesgado. Todas las cifras están tomadas directamente de los archivos de salida citados entre paréntesis.

---

## RESUMEN

**Introducción.** La demencia por cuerpos de Lewy (DCL) tiene loci de riesgo genético bien establecidos (GBA, BIN1, TMEM175, SNCA, APOE), pero se desconoce a través de qué genes y en qué tipos celulares actúan. Existe evidencia de participación inmune en DCL (expansión clonal de linfocitos T CD4⁺ en LCR), pero no está establecido si el riesgo genético opera en células inmunes, en el cerebro, o en ambos.

**Métodos.** Construimos un atlas de colocalización (coloc.abf) entre el GWAS de DCL (Chia et al. 2021; 2.591 casos, 4.027 controles) y 47 conjuntos de eQTL del eQTL Catalogue: 28 inmunes (24 tipos celulares de sangre periférica por scRNA-seq de OneK1K, n≈1.000; BLUEPRINT; sangre GTEx), microglía y 18 tejidos cerebrales. Se testearon todos los genes con eQTL en ±500 kb de cada uno de los 5 loci genome-wide (1.528 tests, 286 genes), sin preselección. Una colocalización se consideró robusta si PP4>0,8, si el SNP colocalizante estaba en LD (r²≥0,5) con el SNP líder del GWAS y si sobrevivía a un prior conservador (p12=5×10⁻⁶). Cada colocalización fuerte se sometió además a fine-mapping multi-señal (SuSiE-coloc). Los genes efectores se evaluaron además en snRNA-seq de corteza (GSE178146, 28 donantes), y un gradiente de severidad previamente descripto se sometió a réplica en una cohorte independiente (GSE150696, n=42).

**Resultados.** Seis pares locus–gen colocalizan con la señal líder. En el locus **SNCA**, un único haplotipo de riesgo específico de DCL (r²=0,81–0,99 entre SNPs; distinto de la señal 3′ de Parkinson, r²≤0,20) **disminuye** la expresión de *SNCA-AS1* en tres tejidos cerebrales (PP4 hasta 0,95), de *SNCA* en monocitos CD14⁺ (PP4=0,87) y de *CCSER1* en linfocitos B naïve (PP4=0,93). SuSiE identifica una única señal GWAS en el locus y confirma las colocalizaciones con *SNCA-AS1* y *CCSER1* (PP4 de SuSiE 0,87–0,91). Los resultados en cerebro reproducen hallazgos publicados; los de células inmunes resueltas por célula única son nuevos para DCL. En **TMEM175**, el alelo de riesgo **aumenta** *GAK* en corteza (ROSMAP, PP4=0,98) y **disminuye** *IDUA* (corteza DLPFC, PP4=0,99), aunque la señal GWAS es demasiado débil para que SuSiE la resuelva. APOE y BIN1 no mostraron colocalización (PP4 máx. ≈0,55), y las colocalizaciones en GBA correspondieron a una señal secundaria, no a la variante líder. La colocalización previamente reportada de *MMRN1* en linfocitos T (PP4=0,97) apareció en un único dataset bulk y en ninguno de los 24 de célula única, donde *MMRN1* no alcanza niveles testeables; además, *MMRN1* se detecta en ≤0,65% de las células de los clústeres T del LCR. El gradiente de severidad de LIFR/LGALS9/CD86 mantuvo la dirección en la cohorte de réplica, pero no alcanzó el umbral pre-especificado.

**Conclusión.** El riesgo genético de DCL actúa a través de efectores distintos según el tipo celular: el haplotipo SNCA de DCL regula genes tanto en cerebro como en células inmunes innatas y B, y el locus TMEM175 actúa principalmente en tejido cerebral. El aporte nuevo es la dimensión inmune de célula única: el haplotipo SNCA de DCL también actúa en monocitos y linfocitos B. El atlas además corrige una inferencia previa (MMRN1 en linfocitos T) y ofrece una lista priorizada de efectores y tipos celulares para validación experimental.

**Palabras clave:** demencia por cuerpos de Lewy, colocalización, eQTL de célula única, SNCA, TMEM175, GAK, neuroinmunología

---

## INTRODUCCIÓN

La demencia por cuerpos de Lewy (DCL) es la segunda demencia neurodegenerativa más frecuente (McKeith et al., 2017). El GWAS más grande a la fecha (Chia et al., 2021) identificó asociaciones genome-wide en GBA, BIN1, TMEM175, SNCA y APOE. Un GWAS es, sin embargo, un mapa de *regiones* y no de *genes* ni de *tipos celulares*: la mayoría de las variantes de riesgo son no codificantes y ejercen su efecto regulando la expresión de genes que pueden no ser el más cercano, y solo en ciertos tipos celulares.

En paralelo, crece la evidencia de participación del sistema inmune adaptativo en sinucleinopatías: linfocitos T reactivos contra alfa-sinucleína (Sulzer et al., 2017) y expansión clonal de linfocitos T CD4⁺ en LCR de pacientes con DCL (Gate et al., 2021). Esto plantea una pregunta concreta y testeable: **¿el riesgo genético de DCL actúa en células inmunes, en el cerebro o en ambos, y a través de qué genes?**

**Qué se sabe.** Chia et al. (2021) colocalizaron el locus SNCA con un eQTL de *SNCA-AS1* en cerebro (PsychENCODE, PP4=0,96) y el locus TMEM175 con un eQTL de *TMEM175* en sangre (eQTLGen, PP4=0,99), ambos con datos bulk. Alvarado et al. (2023) mapearon *GAK* y *TMEM175* al mismo SNP (rs6599388), con efectos eQTL opuestos en DCL. Junn (2026) caracterizó el haplotipo 5′ de *SNCA* en corteza, asociado a menor expresión de *SNCA-AS1*. Gao et al. (2026) aplicaron randomización mendeliana con eQTL de célula única, pero solo en 8 tipos celulares cerebrales. Wu et al. (2024) y Zhang et al. (2026) usaron eQTL de sangre bulk y de cerebro. **Hasta donde pudimos verificar, ningún estudio evaluó los loci de DCL con eQTL de células inmunes resueltas por célula única.**

En una versión previa de este trabajo abordamos esa pregunta con 10 genes candidatos en dos tejidos y encontramos una colocalización fuerte de *MMRN1* en linfocitos T CD4⁺ (BLUEPRINT). Un enfoque de genes candidatos, sin embargo, no permite saber si ese hallazgo es específico o simplemente el único que se miró, y un único dataset bulk no permite atribuirlo a un tipo celular. Aquí reemplazamos ese enfoque por un **atlas de colocalización sistemático**. Testeamos todos los genes de todos los loci genome-wide en 47 conjuntos de eQTL, incluyendo 24 tipos celulares inmunes resueltos por scRNA-seq (OneK1K; Yazar et al., 2022), con criterios de robustez definidos antes de ver los resultados.

---

## MATERIALES Y MÉTODOS

### *GWAS y definición de loci*

Se usaron las estadísticas resumen armonizadas (GRCh38) del GWAS de DCL de Chia et al. (2021; GWAS Catalog GCST90001390; 2.591 casos, 4.027 controles). Se definieron 5 loci genome-wide (p<5×10⁻⁸) agrupando SNPs significativos a <1 Mb, y se tomó como líder el SNP de menor p de cada uno: GBA (rs2230288, E326K), BIN1 (rs6733839), TMEM175 (rs6599388), SNCA (rs7680557) y APOE (rs769449). Para cada locus se extrajo la ventana ±500 kb alrededor del SNP líder (`scripts/phase5_coloc_atlas/02_extract_atlas_gwas.py`).

### *Selección de conjuntos de eQTL*

Del índice oficial del eQTL Catalogue (Kerimov et al., 2021; `tabix_ftp_paths.tsv`) se seleccionaron, por estudio y no por resultado, todos los datasets de expresión génica ("ge") en condición basal de: OneK1K (24 tipos celulares de PBMC por scRNA-seq, n=534–981), BLUEPRINT (monocitos, neutrófilos, linfocitos T CD4⁺; n=167–196), GTEx (sangre y los 13 tejidos cerebrales), ROSMAP (corteza, n=560), CommonMind (DLPFC, n=586), BrainSeq (DLPFC, n=479), Braineac2 (putamen, sustancia nigra) y Young et al. 2019 (microglía, n=104). En total fueron 47 datasets: 28 inmunes, 1 de microglía y 18 cerebrales (`config/pipeline_config.yaml`, sección `phase5_coloc_atlas`).

Las estadísticas se obtuvieron mediante consultas remotas por región con `tabix`, sin descargar los archivos completos (que ocupan decenas de GB), en forma secuencial y espaciada según la política de acceso de EBI (`01_fetch_atlas_eqtl.py`; 235 consultas, 0 fallas).

### *Colocalización*

Para cada combinación dataset × locus × gen con al menos un eQTL nominal (p<10⁻⁴) en la ventana, se ejecutó `coloc.abf` (coloc v5.2.3; Giambartolomei et al., 2014) con priors por defecto (p1=p2=10⁻⁴, p12=10⁻⁵), idénticos a los de la versión previa (`03_run_coloc_atlas.R`). Los alelos se alinearon por posición y alelos ref/alt, invirtiendo el signo cuando los alelos estaban intercambiados y descartando los discordantes. Para cada test se registró la dirección del efecto: el signo del efecto eQTL del alelo que *aumenta* el riesgo de DCL, evaluado en el SNP con mayor probabilidad posterior H4.

### *Criterios de robustez (pre-especificados)*

Se consideró colocalización fuerte PP4>0,8. Cada una se clasificó en tres niveles (`04_summarize_coloc_atlas.R`):
- **Nivel 1:** el SNP colocalizante está en LD r²≥0,5 con el SNP líder del locus (1000 Genomes EUR, PLINK) **y** PP4 sigue >0,8 con el prior conservador p12=5×10⁻⁶ (Wallace, 2020).
- **Nivel 2:** está en LD con el líder, pero no sobrevive al prior conservador.
- **Señal secundaria:** r²<0,5 con el líder. En ese caso, coloc.abf colocalizó una asociación más débil dentro de la ventana, no la señal que define el locus.

Este último criterio es necesario porque coloc.abf asume una única variante causal por región y puede devolver un PP4 alto para una señal secundaria.

### *Fine-mapping multi-señal (SuSiE-coloc)*

coloc.abf supone una única variante causal por rasgo y región. Para cada colocalización con PP4>0,8 se ajustó SuSiE (susieR v0.14.2, vía `coloc::runsusie`; hasta L=5 señales por rasgo) sobre la ventana de ±250 kb alrededor del SNP líder, en el GWAS y en el eQTL por separado. Luego se evaluó la colocalización credible set por credible set con `coloc.susie` (`07_susie_coloc_atlas.R`). La matriz de LD con signo se calculó con PLINK en 1000 Genomes EUR (n=503). Ambos efectos se alinearon al alelo contado del panel de referencia, porque la matriz de LD con signo tiene que referirse al mismo alelo que los z-scores. La consistencia entre z-scores y LD se evaluó con `susieR::estimate_s_rss`: s≈0 indica LD consistente, y s>0,3 indicaría que los credible sets no son confiables. Si SuSiE no encontraba ningún credible set para alguno de los dos rasgos, la colocalización se reportó como sostenida solo por coloc.abf.

### *Genes efectores en corteza con DCL*

Los genes de nivel 1–2 se evaluaron en snRNA-seq de corteza cingulada anterior (GSE178146; Feleke et al., 2021; 28 donantes: 7 Control, 7 EP, 7 EPD, 7 DCL; 2.000 núcleos por donante). Los tipos celulares se asignaron con la misma regla de marcadores canónicos que en la versión previa. Los genes que el control de calidad había eliminado por detección mínima (*SNCA-AS1*) se recuperaron desde las matrices 10x crudas para exactamente los mismos núcleos (`05a_recover_qc_filtered_genes.R`). Por tipo celular se ajustó un modelo pseudobulk por donante (DESeq2, ~severidad con codificación pre-especificada Control=0, EP=1, EPD=2, DCL=3), con corrección BH sobre todos los genes expresados del tipo celular (`05_effector_genes_in_lbd_cortex.R`).

### *Réplica independiente del gradiente de severidad*

El gradiente de expresión de LIFR, LGALS9 y CD86 con la severidad, descubierto en neuronas de GSE178146, se sometió a réplica en GSE150696 (Affymetrix HTA 2.0, corteza prefrontal BA9; Control n=9, EPD n=12, DCL n=12, Alzheimer n=9). Esta cohorte difiere de la de descubrimiento en donantes, región, plataforma y resolución (tejido bulk). La expresión por gen se calculó como la mediana de los probesets exónicos (PSR). El modelo fue limma ~severidad (Control=0, EPD=2, DCL=3) + edad + sexo + intervalo post-mortem; un donante sin dato de PMI recibió la mediana. El criterio de réplica, fijado de antemano, fue: misma dirección que en el descubrimiento y p unilateral <0,05/3. Se agregaron dos análisis complementarios: un ajuste de sensibilidad por contenido neuronal (score de marcadores pan-neuronales) y el contraste Alzheimer vs Control como control de especificidad de enfermedad (`scripts/phase6_replication/01_GSE150696_receptor_replication.R`).

### *Reproducibilidad*

Todos los pasos están integrados en el pipeline Snakemake (reglas `fetch_atlas_eqtl`, `extract_atlas_gwas`, `run_coloc_atlas`, `summarize_coloc_atlas`, `effector_genes_in_lbd_cortex`, `fetch_gse150696` y `replicate_receptor_trend_gse150696`). Las versiones de los paquetes de R están fijadas con renv. Los símbolos génicos provienen de GENCODE v39, la anotación que usa el eQTL Catalogue.

---

## RESULTADOS

### *1. Panorama del atlas: seis pares locus–gen colocalizan con la señal líder de DCL*

Se ejecutaron 1.528 tests de colocalización (286 genes, 5 loci, 47 datasets; `data/coloc_atlas/results/coloc_atlas_annotated.tsv`). Once tests superaron PP4>0,8, correspondientes a 8 pares locus–gen (Tabla 1, Figura 1). Seis de esos pares colocalizan con la señal líder del GWAS: *SNCA-AS1*, *SNCA*, *CCSER1* y *MMRN1* en el locus SNCA, y *GAK* e *IDUA* en el locus TMEM175. Las 11 colocalizaciones se distribuyen en 6 tests en tejido cerebral y 5 en tejidos inmunes; ninguna se observó en microglía.

**Tabla 1.** Colocalizaciones con PP4>0,8 (`coloc_atlas_strong_PP4gt0.8.tsv`). Dirección: efecto del alelo de riesgo de DCL sobre la expresión.

| Locus | Gen | Dataset (n) | PP4 | PP4 (prior conservador) | SNP | r² con líder | Dirección | Nivel |
|---|---|---|---|---|---|---|---|---|
| SNCA | *SNCA-AS1* | CommonMind DLPFC (586) | 0,954 | 0,911 | rs6532198 | 0,92 | ↓ | 1 |
| SNCA | *SNCA-AS1* | GTEx cerebelo (209) | 0,932 | 0,872 | rs3756063 | 0,97 | ↓ | 1 |
| SNCA | *SNCA-AS1* | GTEx hemisferio cerebeloso (175) | 0,846 | 0,732 | rs10005233 | 0,89 | ↓ | 2 |
| SNCA | *SNCA* | OneK1K monocitos CD14⁺ (959) | 0,871 | 0,772 | rs7668531 | 0,92 | ↓ | 2 |
| SNCA | *CCSER1* | OneK1K linfocitos B naïve (980) | 0,926 | 0,863 | rs7680557 | 1,00 | ↓ | 1 |
| SNCA | *MMRN1* | BLUEPRINT linfocitos T CD4⁺, bulk (167) | 0,967 | 0,937 | rs7680557 | 1,00 | ↓ | 1 |
| TMEM175 | *GAK* | ROSMAP corteza (560) | 0,983 | 0,966 | rs6599388 | 1,00 | ↑ | 1 |
| TMEM175 | *IDUA* | CommonMind DLPFC (586) | 0,991 | 0,982 | rs11248057 | 0,79 | ↓ | 1 |
| TMEM175 | *IDUA* | BLUEPRINT linfocitos T CD4⁺, bulk (167) | 0,853 | 0,744 | rs6599388 | 1,00 | ↓ | 2 |
| GBA | *ZBTB7B* | BLUEPRINT neutrófilos (196) | 0,895 | 0,811 | rs2306124 | 0,03 | ↓ | secundaria |
| GBA | *IL6R* | BrainSeq DLPFC (479) | 0,815 | 0,688 | rs12047986 | 0,05 | ↑ | secundaria |

### *2. El haplotipo SNCA de DCL disminuye la expresión de SNCA-AS1 en cerebro y de SNCA en monocitos*

Los cinco SNPs colocalizantes del locus SNCA (rs7680557, rs6532198, rs3756063, rs10005233, rs7668531) están en LD fuerte entre sí (r²=0,81–0,99, 1000G EUR). Es decir, se trata de **un único haplotipo de riesgo**. Ese haplotipo es distinto de la señal 3′ de *SNCA* asociada a la enfermedad de Parkinson (rs356182; r²=0,17–0,20 con todos ellos). En todos los genes y tejidos donde colocaliza, el alelo de riesgo **disminuye** la expresión:

- ***SNCA-AS1*** (ARN largo no codificante antisentido de *SNCA*) en tres datasets cerebrales independientes: DLPFC de CommonMind (PP4=0,95) y dos de cerebelo de GTEx (0,93 y 0,85). Esto **reproduce** lo reportado por Chia et al. (2021) con PsychENCODE y por Junn (2026), ahora en tres cohortes adicionales.
- ***SNCA*** en monocitos clásicos CD14⁺ (OneK1K, n=959, PP4=0,87). El mismo gen muestra una colocalización parcial en linfocitos T CD4⁺ de memoria central (OneK1K CD4_TCM, PP4=0,42). Un posible artefacto sería la contaminación por plaquetas, que expresan *SNCA* en abundancia. Sin embargo, el eQTL de *SNCA* en plaquetas de OneK1K **no** colocaliza (PP4=0,06), lo que va en contra de esa explicación.
- ***CCSER1*** en linfocitos B naïve (OneK1K, n=980, PP4=0,93).

No encontramos antecedentes publicados de las colocalizaciones en monocitos y linfocitos B para DCL.

En BLUEPRINT, el eQTL de *SNCA* en monocitos corresponde a una variante distinta (PP3=1,00; SNP líder rs356199, próximo a la señal 3′ de Parkinson). Esto concuerda con que DCL y Parkinson comparten el locus pero no el haplotipo causal.

### *3. MMRN1 en linfocitos T: una inferencia previa que el atlas no sostiene*

En la versión previa, *MMRN1* en linfocitos T CD4⁺ (BLUEPRINT, PP4=0,97) era el hallazgo principal. El atlas lo reproduce numéricamente (PP4=0,967), pero lo pone en contexto:

1. Es la **única** colocalización de *MMRN1* entre 19 datasets testeados. En sangre total GTEx el PP4 es 0,37, y en los 17 tejidos cerebrales es ≤0,11 (PP3 dominante en la mayoría).
2. *MMRN1* **no fue testeado en ninguno de los 24 datasets de célula única de OneK1K**, incluidas plaquetas y todos los subtipos de linfocitos T, porque no alcanza el umbral de expresión en ningún tipo celular circulante resuelto individualmente.
3. En scRNA-seq de LCR (GSE141578, 24.478 células) se detecta *MMRN1* en 0–0,65% de las células de los clústeres T (CD3E⁺ >88%). Los dos clústeres con mayor detección no son linfocitos T CD4⁺: el clúster 14 (10,1%) es CD3E-negativo, y el clúster 10 (5,6%) es CD4-negativo.
4. En corteza, *MMRN1* es principalmente endotelial (2,0% de núcleos endoteliales frente a 0,12% de neuronas), lo que concuerda con su función conocida en plaquetas y endotelio.

La colocalización *estadística* es sólida: SuSiE la confirma (PP4=0,91, sección 5b). Lo que no se sostiene es su *atribución celular*. Con los alelos alineados al panel de referencia, SuSiE encuentra una sola señal GWAS en la ventana, y no las 9 del análisis SuSiE previo (Fase 1f), que carecía de esa alineación; por eso el valor anterior (PP4=0,977) queda reemplazado. Nuestra interpretación es que la señal de BLUEPRINT refleja una fuente de *MMRN1* no linfocitaria en la muestra bulk, o una expresión marginal. El hallazgo de fondo sí se sostiene: el haplotipo SNCA de DCL tiene efectos regulatorios en células inmunes. Pero los efectores mejor respaldados son *SNCA* en monocitos y *CCSER1* en linfocitos B, no *MMRN1* en linfocitos T.

### *4. En TMEM175, el riesgo aumenta GAK en cerebro y disminuye IDUA*

En el locus TMEM175, el alelo de riesgo **aumenta** la expresión de ***GAK*** (quinasa asociada a ciclina G, reguladora del tráfico de clatrina). La colocalización es fuerte en corteza de ROSMAP (n=560, PP4=0,98; prior conservador 0,97) y parcial en sustancia nigra (GTEx, 0,74) y cerebelo (0,66). En tejidos inmunes el PP4 de *GAK* es ≤0,43 en todos los casos. El mismo alelo **disminuye** ***IDUA*** (alfa-L-iduronidasa, enzima lisosomal) en DLPFC (PP4=0,99) y en linfocitos T de BLUEPRINT (0,85). Ninguno de los datasets de célula única de OneK1K supera 0,53 para *IDUA*. Ambos genes se relacionan con la vía endolisosomal, igual que *TMEM175*. El hallazgo de *GAK* concuerda con Alvarado et al. (2023), que mapearon *GAK* al mismo SNP (rs6599388) en DCL. En cambio, no reproducimos la colocalización de *TMEM175* en sangre de Chia et al. (2021): en nuestro atlas, el PP4 máximo de *TMEM175* es 0,46 (cerebelo) y en sangre GTEx es 0,09. La diferencia se explica probablemente por potencia, ya que eQTLGen (n≈31.700) es unas 47 veces más grande que GTEx sangre (n=670).

### *5. Loci sin colocalización y el límite de GBA*

- **APOE** (94 genes testeados) y **BIN1** (28 genes): el PP4 máximo fue 0,54 y 0,55, respectivamente. En APOE esto es esperable, porque el riesgo depende del alelo codificante ε4 y no de un mecanismo regulatorio.
- **GBA**: las dos colocalizaciones (*ZBTB7B* en neutrófilos, *IL6R* en DLPFC) están en LD despreciable con el líder (r²=0,03 y 0,05), y sus SNPs tienen p GWAS≈2×10⁻⁵. Corresponden a una señal secundaria. La variante líder rs2230288 (E326K, frecuencia ≈1%) casi no está presente en los paneles de eQTL, así que este diseño no puede evaluar el mecanismo regulatorio de GBA, que en cualquier caso es mayormente codificante.

### *5b. Fine-mapping multi-señal: el locus SNCA tiene una única señal y sus colocalizaciones se confirman*

El diagnóstico de consistencia LD–z fue bueno en todos los casos (s≤0,08; `data/coloc_atlas/susie/susie_coloc_summary.tsv`).

| Locus | Gen | Dataset | coloc.abf PP4 | Credible sets (GWAS/eQTL) | SuSiE-coloc PP4 | Veredicto |
|---|---|---|---|---|---|---|
| SNCA | *SNCA-AS1* | CommonMind DLPFC | 0,954 | 1/2 | 0,907 | confirmado |
| SNCA | *SNCA-AS1* | GTEx cerebelo | 0,932 | 1/1 | 0,875 | confirmado |
| SNCA | *SNCA-AS1* | GTEx hemisferio cerebeloso | 0,845 | 1/1 | 0,788 | sugestivo |
| SNCA | *CCSER1* | OneK1K B naïve | 0,926 | 1/1 | 0,866 | confirmado |
| SNCA | *MMRN1* | BLUEPRINT T (bulk) | 0,967 | 1/1 | 0,913 | confirmado (estadísticamente) |
| SNCA | *SNCA* | OneK1K monocitos CD14⁺ | 0,871 | 1/0 | — | solo coloc.abf (sin credible set eQTL) |
| TMEM175 | *GAK* | ROSMAP | 0,983 | 1/0 | — | solo coloc.abf (sin credible set eQTL) |
| TMEM175 | *IDUA* | CommonMind DLPFC | 0,991 | 0/0 | — | solo coloc.abf (sin credible set GWAS) |
| TMEM175 | *IDUA* | BLUEPRINT T (bulk) | 0,853 | 0/0 | — | solo coloc.abf |
| GBA | *ZBTB7B* | BLUEPRINT neutrófilos | 0,895 | 1/1 | 0,002 | no sostenido |
| GBA | *IL6R* | BrainSeq | 0,815 | 0/0 | — | solo coloc.abf |

En el locus SNCA, SuSiE encuentra **un único credible set GWAS, que contiene rs7680557**, en todas las ventanas analizadas. Es decir, la señal de DCL en SNCA es una sola, y las colocalizaciones con *SNCA-AS1* y *CCSER1* se refieren a esa señal. En GBA, SuSiE separa la señal líder del eQTL de *ZBTB7B* (PP4=0,002), lo que confirma la clasificación como señal secundaria. En TMEM175 (p líder=3,5×10⁻⁸), el GWAS no tiene potencia suficiente para formar credible sets en dos de las tres ventanas, así que los hallazgos de *GAK* e *IDUA* descansan solo en coloc.abf. La colocalización de *SNCA* en monocitos también descansa en coloc.abf: el eQTL es claro, pero SuSiE no forma un credible set (con pureza ≥0,5) para él.

### *6. Los genes efectores no cambian con la severidad en corteza*

En la corteza con DCL (GSE178146), ningún gen efector mostró un gradiente de severidad significativo tras la corrección BH en ningún tipo celular (`effector_severity_trend_by_celltype.csv`). La señal más fuerte fue *GAK* en astrocitos y neuronas (log2FC por nivel de severidad =+0,14 y +0,10; p nominal 0,048 y 0,059; padj ≥0,19). Su dirección coincide con el efecto genético, pero no es significativa. *SNCA-AS1* resultó prácticamente indetectable por snRNA-seq (29 UMIs en 56.000 núcleos), por lo que su tipo celular cortical no puede resolverse con estos datos. Esto es compatible con que los eQTL germinales fijen un riesgo basal y no cambien con el estadio de la enfermedad.

### *7. Réplica del gradiente receptorial: dirección concordante, sin significancia*

En GSE150696, los tres receptores pre-especificados mantuvieron la dirección del descubrimiento (aumento con la severidad), pero ninguno alcanzó el umbral de réplica (`receptor_replication_GSE150696.csv`):

| Receptor | Descubrimiento (log2FC/nivel, padj) | Réplica (logFC/nivel) | p unilateral | Ajustado por contenido neuronal |
|---|---|---|---|---|
| CD86 | +0,62 (0,016) | +0,13 | 0,048 | +0,16 (p=0,13) |
| LGALS9 | +0,69 (0,011) | +0,04 | 0,15 | −0,01 (p=0,90) |
| LIFR | +0,18 (0,007) | +0,08 | 0,16 | 0,00 (p=1,00) |

Con umbral 0,05/3=0,017, el resultado es **no replicado**. El efecto en tejido bulk desaparece al ajustar por contenido neuronal, lo que indica que en esta cohorte refleja en gran medida cambios de composición celular. Ninguno de los tres mostró cambio significativo en Alzheimer vs Control (p≥0,27).

---

## DISCUSIÓN

Este trabajo pasa de preguntar "¿colocaliza este gen candidato?" a preguntar "¿dónde y a través de qué genes actúa cada locus de riesgo de DCL?". El cambio de diseño produce tres resultados que un enfoque de candidatos no habría permitido obtener. Conviene separar qué es réplica y qué es nuevo. *SNCA-AS1* en cerebro y *GAK* en el locus TMEM175 reproducen trabajos previos (Chia et al., 2021; Alvarado et al., 2023; Junn, 2026), y aquí se confirman en cohortes adicionales y, en el caso de SNCA, con fine-mapping multi-señal. **Lo nuevo es la dimensión inmune resuelta por célula única**: el mismo haplotipo de DCL colocaliza con *CCSER1* en linfocitos B (confirmado por SuSiE) y con *SNCA* en monocitos CD14⁺ (solo coloc.abf).

**El locus SNCA de DCL es regulatorio y actúa en varios compartimentos.** El haplotipo de riesgo de DCL es distinto de la señal 3′ de Parkinson, lo que concuerda con la idea de arquitecturas genéticas parcialmente diferentes entre ambas enfermedades (Chia et al., 2021). Ese haplotipo disminuye *SNCA-AS1* en cerebro de forma reproducible en tres datasets, y disminuye *SNCA* en monocitos CD14⁺ con resolución de célula única. Cómo se relaciona la reducción de *SNCA-AS1* con los niveles de alfa-sinucleína requiere estudio experimental, y no hacemos afirmaciones sobre esa dirección. La colocalización en monocitos y linfocitos B conecta el riesgo genético con la inmunidad innata y humoral, más que con los linfocitos T que motivaron la versión anterior.

**El locus TMEM175 señala a GAK en cerebro.** *GAK* ya fue propuesto como gen efector en este locus, tanto en Parkinson como en DCL (Alvarado et al., 2023). Aquí el patrón de colocalización favorece un efecto en tejido cerebral (corteza, sustancia nigra) y no en células inmunes.

**El atlas corrige una inferencia previa.** El hallazgo de *MMRN1* en linfocitos T era numéricamente real, pero dependía de un único dataset bulk y no se sostiene en célula única. Presentarlo como el hallazgo principal habría sido un error. Mostramos que un atlas sistemático sirve también para evaluar resultados obtenidos con genes candidatos.

### *Limitaciones*

1. **Potencia del GWAS.** Con 2.591 casos, solo 5 loci son genome-wide significativos; hay loci de DCL reales que este diseño no puede ver.
2. **Fine-mapping limitado por potencia y por el panel de LD.** SuSiE confirmó el locus SNCA, pero no pudo resolver TMEM175 (sin credible sets GWAS) ni el eQTL de *SNCA* en monocitos. Además, el LD proviene de un panel externo de 503 personas; el diagnóstico s≤0,08 indica consistencia, pero no reemplaza el LD propio de la muestra.
3. **Composición y tamaño de los paneles de eQTL.** El atlas no incluye eQTLGen (sangre, n≈31.700), por eso no reproduce la colocalización de *TMEM175* en sangre de Chia et al. (2021). La microglía está representada por un único dataset pequeño (n=104), y no hay eQTL de células inmunes en LCR ni de linfocitos T residentes en cerebro. La falta de colocalización no demuestra ausencia de efecto.
4. **Ancestría.** Tanto el GWAS como los eQTL y el panel de LD son de ancestría europea.
5. **Colocalización no es causalidad.** PP4 alto indica una variante compartida, no que el gen sea el mediador causal.
6. **snRNA-seq.** Las etiquetas de tipo celular se basan en marcadores (no en un atlas de referencia), no se resolvió un clúster microglial, y los genes de baja expresión (*SNCA-AS1*) no son detectables.

### *Trabajo futuro*

- Replicar las colocalizaciones inmunes (*SNCA* en monocitos, *CCSER1* en B) en otros eQTL de célula única de sangre, y con métodos de colocalización específicos de célula única.
- Réplica del atlas con eQTL de célula única de cerebro (por ejemplo, eQTL de núcleos por tipo celular).
- Validación funcional de *SNCA-AS1*/*SNCA* en monocitos y de *GAK* en astrocitos y neuronas.
- Incorporación de GWAS de DCL más grandes cuando estén disponibles.

---

## CONCLUSIONES

1. El riesgo genético de DCL actúa a través de genes efectores específicos y en tipos celulares identificables: *SNCA-AS1* (cerebro), *SNCA* (monocitos CD14⁺), *CCSER1* (linfocitos B) y *GAK*/*IDUA* (principalmente cerebro).
2. El haplotipo SNCA de DCL, distinto del de Parkinson, reduce la expresión de sus efectores tanto en cerebro como en células inmunes.
3. La hipótesis previa de *MMRN1* en linfocitos T no se sostiene con resolución de célula única.
4. El gradiente de severidad de LIFR/LGALS9/CD86 no se replicó en una cohorte independiente de tejido bulk.

---

## DECLARACIÓN DE INTELIGENCIA ARTIFICIAL

Se usaron herramientas asistidas por IA (Claude, Anthropic) para diseñar e implementar el atlas de colocalización y la réplica, ejecutar el pipeline y redactar este manuscrito. Cada cifra reportada fue verificada contra el archivo de salida correspondiente, que se cita en el texto. La interpretación y las afirmaciones finales son responsabilidad de la autora.

## DISPONIBILIDAD DE DATOS Y CÓDIGO

Todos los datos son públicos: GWAS GCST90001390; eQTL Catalogue (47 datasets listados en `config/pipeline_config.yaml`); GEO GSE178146, GSE141578 y GSE150696. El código, el pipeline Snakemake, el `renv.lock` y el `Dockerfile` están en [repositorio a especificar].

---

## REFERENCIAS

- Chia R, et al. (2021). Genome sequencing analysis identifies new loci associated with Lewy body dementia and provides insights into its genetic architecture. *Nature Genetics* 53:294–303.
- Alvarado CX, et al. (2023). Human brain single nucleus cell type enrichments in neurodegenerative diseases. *Research Square* (preprint), PMC10680930.
- Feleke R, et al. (2021). Cross-platform transcriptional profiling identifies common and distinct molecular pathologies in Lewy body diseases. *Acta Neuropathologica* 142:449–474.
- Gao YL, Xu QH, Ye H, Hu N (2026). Cell-type-specific genetic associations in Lewy body dementia identified using single-cell eQTL-based Mendelian randomization. *Archives of Gerontology and Geriatrics*. PMID 41855783.
- Gate D, et al. (2021). CD4⁺ T cells contribute to neurodegeneration in Lewy body dementia. *Science* 374:868–874.
- Giambartolomei C, et al. (2014). Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. *PLoS Genetics* 10:e1004383.
- Junn E (2026). Noncoding variation and regulatory control of SNCA in Lewy body dementia. *npj Dementia*. doi:10.1038/s44400-026-00079-7.
- Kerimov N, et al. (2021). A compendium of uniformly processed human gene expression and splicing quantitative trait loci. *Nature Genetics* 53:1290–1299.
- McKeith IG, et al. (2017). Diagnosis and management of dementia with Lewy bodies: fourth consensus report of the DLB Consortium. *Neurology* 89:88–100.
- Sulzer D, et al. (2017). T cells from patients with Parkinson's disease recognize α-synuclein peptides. *Nature* 546:656–661.
- Wallace C (2020). Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. *PLoS Genetics* 16:e1008720.
- Wu LY, et al. (2024). Investigation of the genetic aetiology of Lewy body diseases with and without dementia. *Brain Communications* 6:fcae190.
- Yazar S, et al. (2022). Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease. *Science* 376:eabf3041.
- Zhang Y, Zhang Z, Qiu S, Hu Y (2026). Disease-predominant loci across Alzheimer's disease, Parkinson's disease and Lewy body dementia: evidence from the UK Biobank prospective cohort, conditional GWAS and colocalization. *Frontiers in Genetics*. doi:10.3389/fgene.2026.1865257.

---

**Figura 1.** Atlas de colocalización (`data/coloc_atlas/results/coloc_atlas_heatmap.pdf`). Filas: pares locus–gen con PP4>0,5 en al menos un dataset. Columnas: 47 datasets agrupados en inmunes, microglía y cerebro. ▲/▼: el alelo de riesgo aumenta/disminuye la expresión (PP4>0,8). †: colocalización con una señal secundaria (r²<0,5 con el SNP líder). Gris: gen no testeado en ese dataset.

**Figura 2.** Expresión de los genes efectores por tipo celular en corteza con DCL (`data/coloc_atlas/effector_expression/effector_expression_by_celltype.pdf`).

**Figura 3.** Réplica de los receptores pre-especificados en GSE150696 (`data/replication_GSE150696/replication_prespecified_receptors.pdf`).
