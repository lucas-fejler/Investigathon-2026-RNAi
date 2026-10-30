# Selección computacional de candidatos para RNAi frente a *Dalbulus maidis* y MRFV

**Investigathon 2026 · Equipo B04**

Mateo Riesco · Tomás Gorodetzky · Lucas Fejler

## Resumen

La interferencia por ARN (RNAi) permite reducir la expresión de genes mediante
moléculas de ARN complementarias a sus transcriptos. La selección de blancos
requiere considerar tanto su función biológica como las características de las
secuencias que se analizarán.

Este proyecto aborda la selección computacional de candidatos frente a la
chicharrita del maíz, *Dalbulus maidis*, y el virus del rayado fino del maíz
(MRFV), en el marco del desafío de EPILAB. El trabajo combina selección de genes,
comparación de herramientas existentes y desarrollo de un algoritmo propio para
priorizar regiones de 200 nucleótidos.

## Contenido del proyecto

| Componente | Contenido |
| --- | --- |
| [Selección de genes](seleccion_genes/README.md) | Criterios biológicos, antecedentes de RNAi y justificación de cinco candidatos adicionales. |
| [Secuencias de referencia](seleccion_genes/secuencias/genes_9.fasta) | Nueve registros FASTA: siete de *D. maidis* y dos de MRFV, con metadatos y fuentes. |
| [Benchmark](benchmark/README.md) | Comparación de si-Fi21, SIREN y dsRNAEngineer sobre 2.752 ventanas de 200 nt. |
| [Algoritmo propio](algoritmo/README.md) | Asimetría termodinámica e interacciones no deseadas; evaluación conjunta sobre un panel de 27 regiones. |

## Selección y análisis

A las cuatro referencias iniciales de EPILAB se incorporaron ETH, Burs-α,
ATP sintasa β, Rpn6 y V-ATPasa B. La selección consideró antecedentes de
silenciamiento, función biológica y disponibilidad de secuencias con procedencia
e identidad revisadas.

El benchmark evalúa la concordancia entre regiones priorizadas bajo criterios
explícitos. Los análisis construidos con si-Fi21 y dsRNAEngineer coinciden en
ETH 217–416, mientras que otros candidatos presentan prioridades diferentes.
La interpretación depende de las métricas, las bases de referencia y el
tratamiento de las coincidencias dentro y fuera de la especie objetivo.

El algoritmo propio calcula la asimetría para las 10.993 regiones posibles de
200 nt, con paso de 1 nt. Sobre un panel de tres regiones por gen, fijado según
esa asimetría, combina el descriptor con una penalización por interacciones no
deseadas. En las 27 regiones evaluadas la penalización alcanza el mismo valor,
por lo que no modifica el orden de asimetría ni sus empates.

El benchmark y el algoritmo propio utilizan distintas poblaciones de regiones,
referencias y reglas de evaluación. Sus resultados caracterizan criterios
computacionales complementarios; no establecen qué método predice mejor la
eficacia biológica.

Los resultados corresponden a una evaluación computacional. No constituyen
validación experimental de eficacia o seguridad biológica.

## Datos y reproducibilidad

El [informe del benchmark](benchmark/informe/benchmark_rnai.html) presenta
metodología, resultados, discusión y referencias. El resumen puede consultarse
directamente en [Benchmark](benchmark/README.md), y las tablas están disponibles
en [resultados](benchmark/resultados/README.md).

Para recalcular las comparaciones con Python 3.10 o posterior:

```bash
python benchmark/reproducir.py
```

Para recalcular la asimetría y los rankings del algoritmo propio con Python
3.11 o posterior:

```bash
python algoritmo/reproducir.py
```

El alcance, las entradas y las dependencias se describen en las guías del
[benchmark](benchmark/REPRODUCIBILIDAD.md) y del
[algoritmo](algoritmo/REPRODUCIBILIDAD.md). El recálculo de los puntajes del
algoritmo utiliza las decisiones de interacción registradas; la revisión de
todos los alineamientos originales requiere activos externos no incluidos.

## Licencia

El código del proyecto se distribuye bajo [licencia MIT](LICENSE). Los datos,
publicaciones y herramientas de terceros mantienen sus respectivas condiciones
de uso. Las fuentes científicas se indican en la selección de genes y en la
[bibliografía](benchmark/bibliografia/README.md).
