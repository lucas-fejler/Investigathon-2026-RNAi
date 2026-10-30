# Metodología del benchmark

## Secuencias y referencias

El conjunto comprende nueve secuencias, con un total de 12.784 nt. Las fuentes,
accesiones, longitudes, límites codificantes conocidos y hashes se encuentran en
los [metadatos de secuencias](../seleccion_genes/secuencias/metadatos_genes.tsv).
Las [fuentes](../seleccion_genes/secuencias/fuentes/) permiten consultar los
registros empleados para fijar las entradas.

En la evaluación de interacciones, si-Fi21 y SIREN utilizan cDNA de Ensembl para
*D. maidis*, mientras que dsRNAEngineer utiliza un ensamblado Trinity. Las
referencias no objetivo corresponden a maíz (*Zea mays*) y abeja
(*Apis mellifera*). Estas bases representan un panel definido de comparación;
no cubren toda la diversidad de organismos potencialmente expuestos.

## Ventanas y coordenadas

Las regiones tienen 200 nt y comienzan en las posiciones 1, 5, 9 y siguientes,
con un paso de 4 nt. No se incorpora una ventana adicional al extremo 3′.
Se obtienen 2.752 regiones. Un siRNA de 21 nt contribuye a una región cuando
queda contenido por completo entre sus extremos.

Las tablas usan coordenadas de **base 1 con fin inclusivo** sobre
`genes_9.fasta`. Las ventanas de las herramientas se unen por
`(gen, inicio, fin)`, sin interpolación.

## Unidades de conteo

| Herramienta | Unidad utilizada en las tablas |
| --- | --- |
| si-Fi21 | Posiciones distintas de siRNA que cumplen un filtro o presentan al menos una coincidencia. La agregación evita contar varias veces la misma posición por sus múltiples alineamientos. |
| SIREN | `n_off_transcriptos`: registros de transcriptos con interacciones; `n_off_sirnas`: suma de entradas de siRNA asociadas a esos transcriptos. Una misma guía puede contribuir en más de un transcripto. |
| dsRNAEngineer | Filas de alineamiento de Bowtie contenidas en la región. Una misma guía puede aportar varios alineamientos. |

Estas unidades no equivalen a números de genes afectados ni a probabilidades de
efecto biológico. En particular, un conteo de SIREN puede superar las 180
posiciones posibles de 21 nt de una región por sus múltiples asociaciones.

Bowtie con coincidencia exacta y RNAhybrid con apareamientos parciales aplican
condiciones distintas. La comparación de perfiles y rankings conserva esa
diferencia; no convierte sus conteos en una escala biológica común.

## Identificación operativa del blanco

En SIREN se excluye del score el transcripto identificado como blanco. Se utiliza
la identificación de si-Fi cuando está disponible; la alternativa operativa es
el registro con más entradas de siRNA, con desempate por identificador.

En dsRNAEngineer, los nombres Trinity se agrupan por gen retirando el sufijo de
isoforma `_iN`. El grupo con más posiciones distintas de siRNA se asigna como
on-target. La tabla `dsrnaengineer_resumen.tsv` permite examinar esa asignación y
la diferencia respecto del segundo grupo. Es una heurística operativa; su
interpretación requiere evaluar la correspondencia de secuencias.

## Criterios de ordenamiento

El módulo [criterio.py](codigo/criterio.py) define dos modos. El principal es
`control_plaga`; el modo alternativo es `funcional`.

| Análisis | Modo `control_plaga` | Modo `funcional` |
| --- | --- | --- |
| si-Fi21 | Menos coincidencias externas; después más coincidencias intraespecie; después más posiciones que pasan los filtros de actividad. | Menos coincidencias externas; después menos coincidencias intraespecie; después más posiciones que pasan los filtros. |
| dsRNAEngineer | Menos coincidencias externas; después mayor suma de alineamientos on-target e intraespecie. | Menos coincidencias externas; después menos coincidencias intraespecie; después más alineamientos on-target. |

El orden es lexicográfico: cada nivel se utiliza únicamente para desempatar el
anterior. Los rankings de las dos herramientas comparten prioridades generales,
pero utilizan variables diferentes.

Para SIREN se utiliza su fórmula de score, aplicada a la rejilla común y con el
blanco excluido:

```text
Score(R) = −0,1 × U(R) − 30 × Σt max(n_t(R) − 1, 0)
```

`U(R)` es el número de nombres de siRNA distintos penalizados en la región, y
`n_t(R)` el número de entradas contenidas en ella para cada transcripto no
objetivo `t`. Un score mayor recibe prioridad. La fórmula representa una
penalización de interacciones y no contiene una medición de eficacia.

El tratamiento favorable de las coincidencias intraespecie en `control_plaga`
es un supuesto de priorización. Su valor biológico depende de si corresponden
al propio blanco, a isoformas o a otros genes y de la función de esos genes.

## Comparación y empates

Se calcula la correlación de rangos de Spearman entre perfiles de conteo y el
índice de Jaccard entre los conjuntos del 10 % superior. En la comparación de
tres herramientas se toman `max(1, n // 10)` ventanas por ranking. Los empates
pueden cortar ese conjunto según el orden estable de las filas; el Jaccard
debe leerse teniendo en cuenta esa condición.

Una correlación sin variación suficiente se expresa como `NA`. La ausencia de
datos y el empate total se registran explícitamente. El
[cruce de regiones](resultados/cruce_recomendaciones.tsv) contiene 24 filas:
los dos candidatos virales y ATPasa D no tienen una región discriminada por
dsRNAEngineer. Sus valores ausentes no se convierten en recomendaciones.

El [análisis de sensibilidad](resultados/sensibilidad_criterio.tsv) compara siete
criterios sobre si-Fi21. Cambiar el tratamiento de las coincidencias
intraespecie modifica la región priorizada en Burs-α, Rpn6 y BicC. Esta
dependencia no determina por sí sola qué criterio tiene mayor validez biológica.

## Alcance

La evaluación es descriptiva y computacional. No utiliza un conjunto de eficacia
experimental como patrón de verdad. Los resultados dependen de las bases,
parámetros, unidades de conteo y criterios declarados. Las coincidencias
predichas no prueban daño; su ausencia no establece seguridad biológica.
