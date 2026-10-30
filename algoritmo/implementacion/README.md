# Implementación científica

Los módulos separan la medición de asimetría, la selección del panel, la
evaluación de interacciones y el ranking. La entrada y las fórmulas se describen
en [Metodología](../METODOLOGIA.md).

| Archivo | Responsabilidad |
| --- | --- |
| `asimetria.py` | Construcción de la guía y cálculo energético de los extremos del dúplex. |
| `puntuacion.py` | Agregación de A y P y cálculo de puntajes con fracciones exactas. |
| `ranking.py` | Ordenamiento dentro de cada gen, puestos compartidos y solapamientos. |
| `seleccion_panel.py` | Admisión por A y manifiesto inmutable de regiones y referencias. |
| `rnahybrid.py` | Lectura estricta de alineamientos, coordenadas y clasificación de interacciones. |
| `evaluar_interacciones.py` | Evaluación por pares, progreso, recuperación y estados incompletos. |
| `priorizar_referencias.cpp` | Priorización heurística de pares; no certifica resultados negativos. |
| `cobertura.py` | Separación entre P global y cobertura de búsqueda por especie. |
| `exportar_resultados.py` | Agregación de las decisiones en resultados regionales. |
| `consolidar_resultados.py` | Comprobación de evidencia y exportación de puntajes, empates y cobertura. |
| `ejecutar.py` | Interfaz de selección, preparación, ejecución y exportación nativa. |
| `archivos.py` | Lectura, escritura y hashes de los archivos científicos. |
| `configuracion.json` | Identidad de candidatos, blancos previstos y parámetros de evaluación. |
| `parametros_energia/` | Tabla energética de dsRIP, procedencia y licencia. |
| `pruebas/` | Pruebas portables de cálculo, coordenadas y estados de cobertura. |

Para el recálculo de los resultados incluidos se utiliza
[`reproducir.py`](../reproducir.py). La ejecución nativa necesita los activos
externos descritos en [Reproducibilidad](../REPRODUCIBILIDAD.md).
