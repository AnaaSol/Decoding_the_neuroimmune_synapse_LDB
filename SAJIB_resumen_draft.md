TÍTULO: Convergencia genética tejido-específica y gradientes de severidad clínica en la demencia por cuerpos de Lewy: un análisis interactómico neuroinmune

Autor(es)¹: Ana [Apellido]

¹ [Afiliación institucional completa]

Categorías temáticas sugeridas: Ómicas (principal); Bioinformática clínica y traslacional (secundaria)

PALABRAS CLAVE: demencia por cuerpos de Lewy; neuroinmunología; colocalización genética; linfocitos T CD4+; transcriptómica de célula única

---

**Introducción:** La demencia por cuerpos de Lewy (DCL) es la segunda demencia neurodegenerativa más frecuente. Estudios GWAS identificaron loci de riesgo (GBA, SNCA, APOE) y trabajos recientes describen expansión clonal de linfocitos T CD4+ en líquido cefalorraquídeo (LCR), proponiendo un eje CXCL12–CXCR4 de reclutamiento al sistema nervioso central. No está establecido en qué tejido actúa el riesgo genético, ni si los hallazgos preliminares resisten un control estadístico riguroso.

**Objetivos:** (1) priorizar genes de riesgo mediante MAGMA y resolver por análisis condicional (GCTA-COJO) si las señales "top" son independientes o reflejan desequilibrio de ligamiento; (2) caracterizar la expansión clonal de CD4+ en dos cohortes de LCR; (3) poner a prueba el eje CXCL12–CXCR4; (4) evaluar convergencia genética tejido-específica (corteza cerebral vs. linfocitos T) mediante colocalización GWAS-eQTL; (5) testear si receptores corticales candidatos muestran un gradiente de expresión asociado a la severidad clínica del espectro DCL.

**Materiales y métodos:** GWAS de DCL (Chia et al. 2021; N=6.618) analizado con MAGMA (18.523 genes) y GCTA-COJO sobre el locus chr19. Se procesaron dos cohortes de scRNA-seq+TCR-seq de LCR (GSE161192, GSE141578) y snRNA-seq de corteza cingulada anterior (GSE178146; 28 muestras, Control/DLB/PDD/PD), con subtipificación neuronal por marcadores de capa cortical. Se ejecutó colocalización (coloc.abf) contra eQTLs de corteza (GTEx) y de linfocitos T CD4+ (BLUEPRINT) para los genes Bonferroni-significativos y receptores candidatos. El gradiente de severidad se evaluó con DESeq2 pseudobulk a nivel de muestra (28 muestras), modelando la severidad clínica como variable ordinal pre-especificada (Control=0<PD=1<PDD=2<DLB=3).

**Resultados:** MAGMA identificó 7 genes Bonferroni-significativos (locus APOE y SNCA); BCAM, señalado como principal en un análisis preliminar subpotenciado, fue descartado formalmente por desequilibrio de ligamiento (p condicional=0,56). La expansión clonal de CD4+ se confirmó en ambas cohortes de LCR; el eje CXCL12–CXCR4 no replicó en ninguna. La resolución de subtipos neuronales corticales (inhibitoria, L2/3, L4, L6) reveló expresión diferencial de ADORA1 entre capas. La colocalización mostró un patrón tejido-específico marcado: MMRN1 no colocalizó en corteza (variantes distintas) pero sí fuertemente en linfocitos T CD4+ (probabilidad posterior de variante causal compartida=0,97), sugiriendo que su riesgo genético actúa vía expresión inmune y no neuronal. El modelo de tendencia ordinal identificó tres receptores con aumento monótono y significativo según severidad clínica (p ajustado<0,05): LIFR (receptor de OSM, la interacción mejor priorizada por expresión), LGALS9 y CD86.

**Conclusión:** Además de confirmar la expansión clonal de CD4+ y descartar formalmente dos hipótesis previas (BCAM, CXCL12–CXCR4), este trabajo aporta dos líneas de evidencia positiva no descriptas antes: una convergencia genética tejido-específica entre el riesgo de DCL y la expresión de MMRN1 en linfocitos T, y un gradiente de expresión de LIFR/LGALS9/CD86 que escala con la severidad clínica del espectro DCL. Generados bajo un esquema de verificación cruzada entre cohortes, tejidos y métodos, estos hallazgos constituyen candidatos concretos y priorizados para validación experimental futura.
