# Resultados del algoritmo

## Regiones en primer puesto

El ranking se limita a las tres regiones admitidas de cada gen. Las coordenadas
son de base 1, con ambos extremos incluidos, sobre el FASTA principal.
**Todas las regiones del panel tienen P = 180/180 = 1**. Los tres pesos
considerados conservan el orden de A y los mismos grupos de empate.

| Secuencia | Regiones en primer puesto dentro del panel | A exacta |
| --- | --- | ---: |
| ATP sintasa β | 2319–2518; 2320–2519; 2321–2520 | 103/180 |
| Burs-α | 36–235; 37–236 | 98/180 |
| ETH | 371–570; 372–571; 373–572 | 104/180 |
| Proteasoma (`DM_RPN6`) | 1341–1540; 1342–1541 | 103/180 |
| V-ATPasa B | 1300–1499; 1301–1500; 1302–1501 | 102/180 |
| V-ATPasa D | 262–461; 263–462; 264–463 | 99/180 |
| BicC | 1883–2082; 1884–2083; 1885–2084 | 100/180 |
| CP (MRFV) | 192–391; 193–392; 194–393 | 95/180 |
| MP putativa (MRFV) | 626–825 | 101/180 |

Todas las regiones enumeradas en una fila comparten el primer puesto; el orden
de presentación no elige una entre ellas. Las regiones empatadas excluidas
durante la admisión se registran en `panel_manifest.json`. Los resultados no
identifican el óptimo combinado entre todas las regiones posibles.

## Archivos para explorar los resultados

| Archivo | Contenido |
| --- | --- |
| [release_region_coverage.tsv](panel/release_region_coverage.tsv) | A, P, estado, puntajes y puestos de las 27 regiones. |
| [release_rankings.tsv](panel/release_rankings.tsv) | Los 81 registros de ranking: 27 regiones por tres pesos. |
| [release_gene_coverage.tsv](panel/release_gene_coverage.tsv) | Completitud y disponibilidad del ranking por gen. |
| [release_exact_ties.json](panel/release_exact_ties.json) | Grupos completos de empate del panel. |
| [species_search_summary.tsv](panel/species_search_summary.tsv) | Búsqueda efectiva por especie y trabajo pendiente o no requerido para P global. |
| [species_gene_coverage.tsv](panel/species_gene_coverage.tsv) | Cobertura por gen y especie. |
| [species_region_coverage.tsv](panel/species_region_coverage.tsv) | Cobertura por región y especie. |
| [release_position_status.tsv](panel/release_position_status.tsv) | Estado y bandera de las 1.821 posiciones distintas. |
| [panel_manifest.json](panel/panel_manifest.json) | Admisión fijada, empates de A, identidad de entradas, políticas y referencias. |

Las versiones JSON/JSONL conservan los campos estructurados utilizados por el
recálculo. `selected_windows.jsonl` y `selected_queries.jsonl` fijan las
posiciones y sus dos orientaciones de consulta. `pair_status.jsonl` contiene
las decisiones registradas por par consulta–referencia;
`native_evidence_inventory.jsonl` identifica la evidencia nativa externa.
El inventario no contiene los alineamientos y no sustituye su revisión.

`release_verification.json` conserva las comprobaciones registradas por la
ejecución científica de origen. El alcance del recálculo que se puede realizar
con esta distribución se detalla en [Reproducibilidad](../REPRODUCIBILIDAD.md).
La correspondencia con los archivos incluidos está fijada en
[integridad.json](integridad.json).

## Cómo interpretar P completo

Una posición positiva ya tiene al menos una interacción no deseada calificante
en las referencias. No es necesario terminar los demás pares para resolver ese
indicador global. Esto deja búsquedas parciales o sin iniciar en algunas
especies. Una especie sin testigo positivo no debe interpretarse como libre
de interacciones si su búsqueda no terminó.

La penalización saturada limita la capacidad de discriminar entre las regiones
de este panel. No demuestra daño biológico, igualdad de riesgos entre especies,
ni superioridad o inferioridad frente a las herramientas del benchmark.
