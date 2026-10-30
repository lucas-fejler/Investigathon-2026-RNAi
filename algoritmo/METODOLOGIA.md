# Metodología del algoritmo

## Entradas, coordenadas y blanco previsto

La entrada es `genes_9.fasta`: siete registros de *Dalbulus maidis* y dos de
MRFV, con 12.784 nt en total. Se conservan las bases y la orientación de las
[referencias seleccionadas](../seleccion_genes/README.md), incluida la UTR
disponible. Se enumeran todas las regiones contiguas de 200 nt, con paso de
1 nt. Las tablas exportan coordenadas desde 1 con ambos extremos incluidos;
los identificadores internos incorporan la secuencia, su hash y coordenadas
desde 0 con extremo final excluido.

El blanco previsto se define mediante registros e intervalos explícitos en
[configuracion.json](implementacion/configuracion.json). Para los genes del
insecto se consideran los transcriptos anotados del locus declarado; para los
blancos virales, la región correspondiente dentro del genoma de MRFV. Una
interacción en otro gen del insecto o fuera de la región viral prevista puede
contribuir a P.

El identificador `DM_RPN6` fija la secuencia `ENSCHMT00000015190.1`. En la
configuración de esta evaluación, su locus es `ENSCHMG00000008527.1`, anotado
como PSMD6/Rpn7, con tres transcriptos previstos: `ENSCHMT00000015190.1`,
`ENSCHMT00000015206.1` y `ENSCHMT00000015227.1`. La
[selección biológica](../seleccion_genes/README.md#correspondencia-de-nombres-y-secuencias)
explica el uso de la denominación Rpn6 en el antecedente. Los nombres no
reemplazan la identidad de secuencia ni resuelven por sí solos la ortología.

## Asimetría termodinámica: A

Para cada posición de 21 nt se modela una guía antisentido, complementaria
inversa de la ventana suministrada. El dúplex incluye un núcleo de 19 pares
y extremos salientes de 2 nt. Se comparan los cuatro pares terminales del
extremo 5′ de la guía con los cuatro del extremo opuesto. La diferencia
energética positiva se clasifica como favorable; cero es neutral y una
diferencia negativa es desfavorable.

Se utiliza la tabla de 256 tetranucleótidos de la implementación de
[dsRIP](https://doi.org/10.1186/s12915-025-02219-6), con energías en centésimas
de kcal/mol, basada en parámetros de
[Zuber et al.](https://doi.org/10.1093/nar/gkac261). La tabla, su licencia y
procedencia se incluyen en
[parametros_energia](implementacion/parametros_energia/provenance.json).

```text
A(región) = número de posiciones clasificadas como favorables / 180
```

A es un descriptor del modelo, no una eficacia medida. Se conserva la tabla
original, incluido su tratamiento de tetranucleótidos autocomplementarios.
El modelo tiene supuestos de temperatura y contexto distintos del organismo
objetivo; la transferencia de sus asociaciones de eficacia a *D. maidis* o
MRFV no está validada por este trabajo.

## Selección del panel

Se ordenan las regiones de cada gen por A exacta descendente y se admiten tres.
Entre regiones empatadas, las coordenadas y el identificador determinan cuáles
entran en el presupuesto de evaluación. No se eliminan solapamientos. El
[manifiesto del panel](resultados/panel/panel_manifest.json) conserva también
los integrantes excluidos de los grupos empatados.

El panel se fija antes de observar P: 27 regiones, 4.860 pertenencias regionales
de posiciones y 1.821 posiciones únicas, con dos orientaciones de consulta.
Las otras 10.966 regiones carecen de un puntaje combinado en esta evaluación.
Una región con menor A podría superar a una admitida si tuviera menor P;
por eso el ranking del panel no establece el óptimo global de cada secuencia.

## Interacciones no deseadas: P

Se usa RNAhybrid 2.1.2 sobre registros completos de referencia, en su
orientación almacenada, con ambas orientaciones de cada consulta de 21 nt.
Una interacción califica cuando satisface:

```text
MFE ≤ −10 kcal/mol
y (al menos 18 pares Watson–Crick o GU, o 16 pares Watson–Crick consecutivos)
```

La configuración permite bucles internos y salientes de hasta 15 nt, no fuerza
una hélice y utiliza `-p 1 -d 0,1`. Los núcleos apareados contenidos en los
intervalos previstos del candidato se excluyen de P; otras interacciones con
el mismo registro siguen siendo evaluables.

Un barrido de semillas exactas de 8 nt y similitud sin huecos prioriza pares
consulta–referencia. **No descarta definitivamente los pares restantes**:
las posiciones aún no decididas continúan con las referencias pendientes.
Un prefijo inicial de cinco resultados se amplía cuando sigue sin decidirse
la posición, comprobando continuidad y progreso del enmascaramiento.

Una posición cuenta una sola vez aunque tenga varios sitios, orientaciones o
especies con interacciones calificantes. Un positivo permite detener el trabajo
pendiente para esa posición. Un resultado negativo exige completar la política
de búsqueda declarada. Fallos, interrupciones y búsquedas incompletas no se
convierten en ceros. La ausencia de interacciones bajo esta política de reportes
de RNAhybrid no equivale a enumerar todos los dúplex físicamente posibles.

```text
P(región) = posiciones con alguna interacción no deseada calificante / 180
```

## Referencias y cobertura

El catálogo fijado contiene 340.503 registros:

| Referencia | Registros | Pares evaluados | Posiciones con testigo positivo en esa especie |
| --- | ---: | ---: | ---: |
| *D. maidis* | 49.881 | 1.671 | 329 |
| *H. sapiens* | 186.185 | 879 | 876 |
| *Z. mays* | 76.549 | 468 | 466 |
| *A. mellifera* | 27.887 | 99 | 97 |
| MRFV | 1 | 53 | 53 |

Los registros preservan las variantes de las fuentes; no se deduplican como si
fueran genes. El manifiesto describe ensamblados, versiones y políticas de
inclusión. Las posiciones positivas pueden coincidir entre especies, por lo
que la última columna no debe sumarse como posiciones únicas.

Las 1.821 posiciones tienen al menos un testigo positivo en la unión de
referencias. Se resolvió ese indicador con 3.170 pares evaluados. Los
1.240.108.756 pares restantes quedaron sin evaluar tras alcanzarse un positivo;
no son negativos. La cobertura por especie continúa siendo parcial y está
documentada separadamente en las tablas.

## Puntaje, ranking e interpretación

```text
Sλ = A − λP, con λ = 0,5; 1; 2
```

Las cuentas y comparaciones utilizan fracciones exactas. Dentro de cada gen,
las regiones con el mismo puntaje comparten puesto; el siguiente puesto
considera el tamaño del empate, por ejemplo 1, 1, 3. Las coordenadas ordenan la
presentación de los empates, sin resolverlos biológicamente.

En este panel, todas las regiones tienen P = 1. Para cualquier λ considerado,
el puntaje resta una constante y conserva el orden de A. El resultado muestra
una **saturación de la penalización bajo esta política y este panel**, sin
demostrar que todas las regiones tengan el mismo riesgo real. Evaluar la
discriminación del método requiere ampliar la cobertura y estudiar criterios
con datos independientes; modificar umbrales después de ver estos resultados
constituiría un análisis nuevo.

## Fuentes del método

- Cedden D et al. (2025). *Optimizing dsRNA sequences for RNAi in pest control and research with the dsRIP web platform*. [BMC Biology](https://doi.org/10.1186/s12915-025-02219-6).
- Cedden D. *dsRIP web platform v1.0*, depósito de código y parámetros, versión 3. [Figshare](https://doi.org/10.6084/m9.figshare.28156787).
- Zuber J et al. (2022). *Nearest neighbor rules for RNA helix folding thermodynamics: improved end effects*. [Nucleic Acids Research](https://doi.org/10.1093/nar/gkac261).
- Rehmsmeier M et al. (2004). *Fast and effective prediction of microRNA/target duplexes*. [RNA](https://doi.org/10.1261/rna.5248604).
