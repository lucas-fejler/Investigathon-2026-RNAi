# Datos y resultados del benchmark

Las coordenadas de las regiones se expresan en base 1, incluyendo ambos
extremos, sobre `genes_9.fasta`. Cada ventana tiene 200 nt y pertenece a la
rejilla con paso de 4 nt definida en la [metodología](../METODOLOGIA.md).

## Comparaciones

| Archivo | Contenido |
| --- | --- |
| [comparacion_tres.tsv](comparacion_tres.tsv) | Resumen de las tres herramientas por secuencia: regiones priorizadas, correlaciones y Jaccard. |
| [comparacion_tres_ventanas.tsv](comparacion_tres_ventanas.tsv) | Indicadores de las tres herramientas para las 2.752 ventanas comunes. |
| [comparacion_resumen.tsv](comparacion_resumen.tsv) | Resumen de si-Fi21 y SIREN por secuencia. |
| [comparacion_ventanas.tsv](comparacion_ventanas.tsv) | Métricas de si-Fi21 y SIREN unidas por gen y coordenadas. |
| [cruce_recomendaciones.tsv](cruce_recomendaciones.tsv) | Regiones priorizadas y percentiles de conteo en los otros análisis. |
| [sensibilidad_criterio.tsv](sensibilidad_criterio.tsv) | Regiones priorizadas con siete criterios de ordenamiento. |

## Indicadores por herramienta

| Archivo | Contenido |
| --- | --- |
| [sifi_ventanas200.tsv](sifi_ventanas200.tsv) | Posiciones de siRNA por ventana: filtros de actividad y coincidencias intraespecie, maíz y abeja. |
| [siren_regiones200.tsv](siren_regiones200.tsv) | Score con y sin el blanco, registros de transcriptos y entradas de siRNA por ventana. |
| [dsrnaengineer_ventanas200.tsv](dsrnaengineer_ventanas200.tsv) | Alineamientos de Bowtie agregados por ventana y categoría. |
| [dsrnaengineer_resumen.tsv](dsrnaengineer_resumen.tsv) | Correspondencia de entradas, grupo on-target y conteos por secuencia. |
| [sifi_mejor_ventana.tsv](sifi_mejor_ventana.tsv), [siren_mejor_region.tsv](siren_mejor_region.tsv) | Regiones priorizadas de cada análisis. |

Las evaluaciones de SIREN contra especies no objetivo se encuentran en
`no_objetivo_ZMAYS/` y `no_objetivo_AMELLIFERA/`. La tabla
[no_objetivo_tres.tsv](no_objetivo_tres.tsv) reúne conteos con unidades distintas;
sus cocientes no representan diferencias de riesgo biológico.

## Alineamientos de dsRNAEngineer

`dsRNAEngineer_crudo/` contiene los archivos comprimidos de las consultas y
sus alineamientos extraídos, junto con `input.fa`. Los nombres de contigs
Trinity permiten seguir las asignaciones operativas de blanco e isoformas.

## Valores ausentes y empates

`NA` indica una correlación sin variación suficiente. `sin datos` y
`empate total` distinguen falta de información de falta de discriminación entre
regiones. Los campos vacíos no deben transformarse automáticamente en cero.
Los conteos son indicadores computacionales, no números de genes afectados
ni mediciones de eficacia o daño.
