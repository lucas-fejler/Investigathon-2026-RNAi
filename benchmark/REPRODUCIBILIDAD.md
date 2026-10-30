# Reproducción de los análisis

## Requisitos

Python 3.10 o posterior. El análisis de las tablas utiliza únicamente la
biblioteca estándar y puede ejecutarse en Windows, Linux o macOS.

Desde la raíz del proyecto:

```bash
python benchmark/reproducir.py
```

Para elegir otro directorio de resultados:

```bash
python benchmark/reproducir.py --salida salidas/otra_ejecucion
```

La carpeta debe ser nueva o estar vacía. Los archivos de referencia permanecen
en sus directorios y las nuevas salidas se escriben dentro de `salidas/`.

## Entradas y etapas

| Etapa | Archivos utilizados | Resultado |
| --- | --- | --- |
| Reconstrucción de dsRNAEngineer | Archivos `.bwt`, `input.fa` y FASTA de referencia. | Tablas de alineamientos agregados por ventana y resumen por secuencia. |
| Comparación de dos herramientas | `sifi_ventanas200.tsv`, `siren_regiones200.tsv`. | Tablas unidas por ventana y resumen por secuencia. |
| Comparación de tres herramientas | Las anteriores y `dsrnaengineer_ventanas200.tsv`. | Rankings, correlaciones y solapamiento entre conjuntos superiores. |
| Cruce de regiones priorizadas | Rankings y métricas de las tres herramientas. | Percentiles de cada región en los otros análisis. |
| Sensibilidad | Tabla regional de si-Fi21. | Comparación de siete criterios de ordenamiento. |

Las tablas de referencia están en [resultados](resultados/README.md). Las
secuencias se encuentran en
[seleccion_genes/secuencias](../seleccion_genes/secuencias/README.md).

## Salidas y comprobaciones

El programa guarda las tablas recalculadas, un log por etapa y un archivo
`reproduccion.json` con los comandos, la versión de Python, hashes de las
entradas y comparación del contenido de los TSV. Se comparan ocho tablas.
La comparación de contenido admite diferencias de saltos de línea entre
plataformas, sin alterar campos o valores.

Los códigos de salida son: `0`, resultados coincidentes; `1`, error de ejecución;
`2`, diferencias respecto de las tablas de referencia.

La integridad de las secuencias puede comprobarse por separado:

```bash
python seleccion_genes/secuencias/verificar.py --json
```

Ese control verifica formato, fuentes y hashes del conjunto FASTA. No evalúa
eficacia biológica.

## Repetición de las ejecuciones de las herramientas

El comando principal reconstruye el posprocesamiento y la comparación. La
repetición de si-Fi21 y SIREN desde sus entradas requiere sus entornos, bases de
referencia y salidas crudas completas; esos recursos no están todos incluidos.
dsRNAEngineer dispone de sus archivos de alineamiento en este repositorio.

Los módulos para agregación desde salidas externas y el controlador de si-Fi21
se describen en [código](codigo/README.md) y
[herramientas](herramientas/README.md). Una repetición completa debe fijar las
versiones y referencias de las ejecuciones antes de comparar sus resultados.
