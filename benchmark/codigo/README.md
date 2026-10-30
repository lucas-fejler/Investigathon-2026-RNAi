# Código de análisis

La entrada al análisis es `python benchmark/reproducir.py`, desde la raíz del
proyecto. Los módulos separan la agregación de datos, la definición de criterios
y la comparación de regiones.

| Módulo | Función |
| --- | --- |
| `parsear_dsrnaengineer.py` | Agrega los alineamientos de dsRNAEngineer a ventanas de 200 nt y agrupa isoformas Trinity. |
| `criterio.py` | Define las prioridades de ordenamiento de si-Fi21 y dsRNAEngineer. |
| `comparar_herramientas.py` | Une si-Fi21 y SIREN por gen y coordenadas; calcula correlaciones y coincidencia de rankings. |
| `comparar_tres.py` | Compara las tres herramientas y distingue ausencia de datos y empate total. |
| `cruce_recomendaciones.py` | Calcula percentiles de las regiones priorizadas en los otros análisis. |
| `sensibilidad_criterio.py` | Compara siete criterios sobre los datos de si-Fi21. |

## Agregación y análisis complementarios

- `ventanas_sifi.py`: construye tablas regionales desde las salidas de si-Fi21.
- `recomputar_siren.py`: aplica la fórmula de SIREN a la rejilla completa de
  regiones, excluyendo el transcripto identificado como blanco.
- `comparar_no_objetivo.py`: resume coincidencias con maíz y abeja, conservando
  la distinción entre las unidades de conteo de cada método.
- `prefiltro_pares.py`: implementación de un prefiltro por pares de consulta y
  referencia para reducir el trabajo de RNAhybrid. Su equivalencia necesita
  comprobarse para el alcance y los parámetros de cada evaluación.

Estos módulos complementarios requieren sus archivos de entrada y herramientas
externas. La [metodología](../METODOLOGIA.md) define cómo interpretar sus
indicadores y la [guía de reproducción](../REPRODUCIBILIDAD.md) especifica el
alcance del análisis ejecutable con los datos disponibles.
