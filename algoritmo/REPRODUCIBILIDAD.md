# Reproducibilidad del algoritmo

## Recálculo con los archivos incluidos

Requiere **Python 3.11 o posterior**, con biblioteca estándar. Funciona sin
RNAhybrid, compilador ni descarga de referencias. Desde la raíz del repositorio:

```bash
python algoritmo/reproducir.py
python -m unittest discover -s algoritmo/implementacion/pruebas -v
```

El recálculo utiliza el FASTA principal, la tabla energética y los registros
del panel para:

1. Comprobar hashes de entradas, código y resultados.
2. Recalcular la asimetría de todas las posiciones de 21 nt y las 10.993 regiones
   de 200 nt posibles.
3. Reconstruir la admisión de las 27 regiones y los empates de A, incluidos
   integrantes excluidos del panel.
4. Comprobar coordenadas, secuencias de consulta y correspondencia entre
   posiciones positivas y decisiones de interacción registradas.
5. Reagrupar las banderas de P por región y recalcular los puntajes y los
   **81 registros de ranking**, conservando sus empates exactos.

Los resultados derivados se escriben en `salidas/algoritmo/`, fuera de los
archivos publicados. El comando rechaza una carpeta de salida que ya tenga
contenido. Para otra ubicación:

```bash
python algoritmo/reproducir.py --salida salidas/otra_ejecucion
```

| Archivo derivado | Contenido |
| --- | --- |
| `asimetria_regiones.tsv` | A para las 10.993 regiones con paso de 1 nt. |
| `puntajes_panel.tsv` | Puntajes, puestos y tamaños de empate para las 27 regiones y tres pesos. |
| `resumen.json` | Cantidades comprobadas y alcance del recálculo. |

Las **43 pruebas portables** comprueban parámetros energéticos, ejemplos
calculados a mano, correspondencia con el cálculo de referencia de dsRIP,
aritmética exacta, empates, cobertura, coordenadas, fallos y búsqueda de respaldo.
Utilizan datos de prueba identificados como sintéticos. No ejecutan RNAhybrid
ni sustituyen pruebas experimentales.

## Alcance de la evidencia de interacción

El repositorio incluye los resultados por posición, las decisiones por par
consulta–referencia y un inventario de evidencia con sus hashes. **No incluye
los alineamientos nativos completos ni los catálogos de referencia**. Por ello,
el recálculo portable verifica consistencia y agregación de las decisiones
registradas, sin volver a determinar si cada alineamiento original cumple
los umbrales. Un hash identifica un archivo; no reemplaza ese archivo.

Los resultados originales del panel se conservan sin modificar sus bytes.
La identidad del código y las adaptaciones de nombres y rutas están registradas
en [integridad.json](resultados/integridad.json). La versión científica de origen
se identifica por el commit `f28a2788f89e3ea1da885113d42976571c6aee05` del
repositorio de desarrollo con acceso restringido. El código científico y las
tablas necesarios para el recálculo descrito están incluidos aquí.
Las rutinas de cálculo conservan sus fórmulas y umbrales. Una nueva ejecución
nativa de esta implementación tiene una nueva identidad de código; no se
reetiquetan los resultados almacenados como si provinieran de esa ejecución.

## Nueva evaluación con RNAhybrid

La ejecución completa requiere activos externos que **no se obtienen al clonar
este repositorio**:

| Activo | Función |
| --- | --- |
| Referencias y sus catálogos del 28/10/2026 | Secuencias completas, anotaciones e identidad de los 340.503 registros. |
| RNAhybrid 2.1.2 y su versión instrumentada | Evaluación nativa, opciones y registro de evidencia. |
| Enumeración preparada | Padres, regiones, posiciones y consultas con sus identificadores. |
| Asimetría preparada y evidencia de orientación | Medidas por posición y región, con las identidades requeridas por el evaluador. |
| Evidencia nativa de la evaluación | Necesaria para revisar los alineamientos o revalidar las decisiones originales. |

La implementación nativa requiere un entorno Linux compatible, Python 3.11+
y `g++` para el priorizador. La preparación de referencias y enumeración no está
incluida como un comando autónomo de esta distribución. Los activos preparados
deben ajustarse al contrato e identidades del manifiesto; no basta con descargar
transcriptomas de otra versión.

Los módulos permiten fijar un panel nuevo, preparar las referencias, ejecutar
la evaluación y exportar resultados mediante `implementacion/ejecutar.py`.
Esos pasos se realizan en una carpeta de ejecución nueva una vez disponibles
y verificados los activos requeridos. No deben ejecutarse sobre las tablas
históricas publicadas. El código conserva los controles que impiden reutilizar
evidencia con entradas o política diferentes.

La evaluación exhaustiva de P para las 10.993 regiones, una comparación homogénea
con las otras herramientas y la validación biológica quedan fuera del alcance
de estos resultados.
