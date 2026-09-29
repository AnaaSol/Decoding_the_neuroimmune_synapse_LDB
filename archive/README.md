# Archivo de auditoría (histórico, no vigente)

Estos documentos registran el proceso de auditoría y corrección que llevó al estado actual del proyecto (julio–septiembre 2026): inventario inicial del pipeline, verificación de cada afirmación del manuscrito contra los archivos de salida reales, y el reporte de correcciones aplicadas.

**No reflejan el estado actual del pipeline ni de los resultados** — para eso, ver `README.md`, `MANUSCRIPT_Neuroimmune_Interactions_LBD.md` y `PIPELINE_DIAGRAM.md` en la raíz del proyecto. Se conservan acá como evidencia del proceso de verificación (útil para responder preguntas sobre metodología o mostrar el rigor del proceso), no como documentación de referencia.

| Archivo | Contenido |
|---|---|
| `INVESTIGATION_REVIEW_SUMMARY.md` | Documento de traspaso inicial: qué se sabía del proyecto antes de auditar cada afirmación contra el código y los datos. |
| `AUDIT_INVENTORY.md` | Mapa completo del repositorio, grafo de dependencias del pipeline, y los ~20 defectos encontrados en el código original (incluida la etiqueta de condición invertida en Fase 3, el pre-filtrado que invalidaba MAGMA, etc.). |
| `CLAIMS_CHECK.md` | Verificación afirmación por afirmación del manuscrito original contra los datos reales: 53 afirmaciones revisadas, 41 incorrectas. |
| `REPORT.md` | Reporte de las correcciones aplicadas, con el detalle de cada defecto (antes/después) y el estado de reproducibilidad en distintos puntos del proceso. |
