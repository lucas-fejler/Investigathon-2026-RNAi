# Secuencias de genes para diseño computacional y benchmark

**Archivo principal: `genes_9.fasta`.** Reúne los cinco candidatos de *Dalbulus maidis* y las cuatro entradas proporcionadas por EPILAB. Contiene una secuencia de referencia por entrada, con toda la longitud disponible de la fuente elegida: **9 registros y 12.784 nt**. Todas las secuencias usan A/C/G/T, tienen al menos 200 nt y conservan sus bases y orientación originales.

El objetivo es que el equipo pueda implementar la búsqueda y evaluación de regiones de 200 nt y compararla con herramientas existentes sobre un conjunto reproducible. No se eligió una región ganadora ni se impuso un recorte por coincidencia con otras referencias.

## Archivos para trabajar

| Archivo o carpeta | Contenido y uso |
|---|---|
| **`genes_9.fasta`** | Entrada principal. Cinco referencias públicas y cuatro secuencias de EPILAB, con la longitud disponible intacta. |
| `metadatos_genes.tsv` | ID, gen, organismo, fuente, accesión con versión cuando existe, longitud, tipo de secuencia, límites CDS conocidos y SHA256. |
| `comparacion_EPILAB/` | Los cuatro fragmentos originales de 200/300 nt de EPILAB y sus coordenadas. |
| `concordancia/` | Análisis opcional de tramos codificantes compartidos entre fuentes; no sustituye la entrada principal. |
| `auditoria/` | Diferencias, mapas de alineamiento y método usado para identificarlas. |
| `fuentes/` | Registros originales pequeños, extractos de Ensembl, planilla EPILAB, suplemento y procedencia. |
| `verificar.py` y `SHA256SUMS.txt` | Verificación de integridad, formato, correspondencia con fuentes y análisis de concordancia. |

Los encabezados FASTA contienen únicamente un ID ASCII estable, sin espacios. Los nombres y accesiones se consultan en los metadatos.

## Referencias elegidas y alcance

| ID | Gen / entrada | nt | Fuente fijada | Qué representa |
|---|---|---:|---|---|
| DM_ETH | ETH | 615 | PZ371507.1 | CDS completa anotada de D. maidis, sin UTR en el depósito |
| DM_BURS_A | Burs-α | 468 | PZ371508.1 | CDS completa anotada de D. maidis, sin UTR en el depósito |
| DM_ATP5F1B | ATP sintasa β | 2714 | ENSCHMT00000015711.1 | Transcripto predicho de D. maidis, con UTR |
| DM_RPN6 | Rpn6 | 2513 | ENSCHMT00000015190.1 | Transcripto predicho de D. maidis, con UTR; una isoforma fijada |
| DM_VATP_B | V-ATPasa B | 1672 | ENSCHMT00000003421.1 | Transcripto predicho de D. maidis, con UTR |
| EPILAB_BicC_Dm | BicC | 2520 | Planilla EPILAB, Hoja1, D4 | Secuencia original proporcionada; integridad del transcripto completo no establecida |
| EPILAB_ATPaseD_Dm | ATPasa D | 491 | Planilla EPILAB, Hoja1, D5 | Segmento original proporcionado |
| EPILAB_MP_MRFV | MP putativa de MRFV | 1263 | Planilla EPILAB, Hoja1, D2 | Secuencia original bajo la etiqueta MP |
| EPILAB_CP_MRFV | CP de MRFV | 528 | Planilla EPILAB, Hoja1, D3 | Segmento original bajo la etiqueta CP; no se afirma CDS independiente completa |

Los modelos de Ensembl corresponden al ensamblado `GCA_975739785.1`, anotación `2026_02`. Las nueve entradas incluyen siete del vector y dos virales. Las secuencias públicas de los cinco candidatos son de D. maidis; las referencias experimentales de especies emparentadas usadas para justificar algunos genes no se incorporan como blancos.

## Por qué se elige cada referencia y cómo se tratan las discrepancias

- **ETH:** se elige el depósito asociado al estudio, con accesión y versión. Sus 615 nt coinciden con la CDS del modelo genómico y con una subsecuencia del suplemento al considerar la orientación. La diferencia de longitud con el suplemento incluye secuencia adicional fuera de la CDS depositada.
- **Burs-α:** se elige el depósito del estudio `PZ371508.1`. Su CDS de 468 nt coincide exactamente con la CDS y proteína del modelo `ENSCHMT00000035038.1`. El suplemento difiere en 14 sustituciones y 3 columnas con huecos; la causa no está establecida. Se conserva el depósito completo como referencia y se documenta el conflicto, sin mezclar bases. La etiqueta automática `SRSF3` del modelo genómico es inconsistente con la identidad del depósito; no se usa para reasignar el blanco.
- **ATP sintasa β:** se fija `ENSCHMT00000015711.1`, cuyo tramo CDS 115–1794 se traduce exactamente a `ENSCHMP00000011487.1`. La correspondencia del gen tiene respaldo proteico previo frente a dos referencias de Cicadellidae. La secuencia sigue siendo un modelo predicho; no se presume validada experimentalmente en toda una población.
- **Rpn6:** se fija el modelo compacto `ENSCHMT00000015190.1`, apoyado por comparación proteica. Su CDS 300–1514 aparece idéntica dentro de las otras dos isoformas largas anotadas revisadas, `ENSCHMT00000015206.1` y `ENSCHMT00000015227.1`. Esto respalda ese núcleo codificante; no confirma la arquitectura completa ni la expresión de todas las isoformas.
- **V-ATPasa B:** se fija `ENSCHMT00000003421.1` / `ATP6V1B2`, con CDS 163–1665 y traducción exacta a `ENSCHMP00000001802.1`. La consistencia de la anotación y la proteína justifica el uso computacional de ese modelo; no equivale a validación experimental de toda su secuencia ni de su eficacia.
- **Las cuatro entradas EPILAB:** se eligen sus originales porque son las secuencias proporcionadas para el desafío y permiten reproducir sus comparaciones. La planilla no declara accesiones ni aislados exactos. Un modelo o aislado público parecido no basta para reemplazarlas. Se preservan todas sus bases y longitudes.

Las comparaciones adicionales de EPILAB muestran: BicC frente al modelo `BICC1` presenta 3 sustituciones y una diferencia interna de 36 nt; ATPasa D frente a `ATP6V1D`, 2 sustituciones; MP frente a `AF265566.1`, 10 sustituciones y 3 nt de diferencia interna; CP frente al mismo aislado, 6 sustituciones dentro del segmento comparado. No se atribuye su causa a variación natural, errores o cambios de isoforma sin evidencia adicional.

**La conducta ante estas diferencias es fijar una fuente por entrada, conservar sus bases y registrar el cotejo.** Coincidir con otra referencia no es un requisito de validez de un FASTA ni demuestra mayor eficacia. Las regiones discrepantes permanecen disponibles en la entrada principal; pueden examinarse mediante análisis de sensibilidad entre referencias si resulta relevante para los resultados.

## CDS, UTR y coordenadas

Una CDS es la región que codifica la proteína principal. Las UTR forman parte del ARN mensajero fuera de esa región y pueden participar en su regulación. Los tres transcriptos predichos conservan sus UTR; los depósitos de ETH y Burs-α aportan solo la CDS. La falta de UTR en un registro no demuestra que el ARN natural carezca de ellas.

La entrada principal permite explorar toda la longitud disponible. Si el equipo decide restringir un análisis a CDS, sus límites conocidos están en `metadatos_genes.tsv`, en **base 1 e incluyendo ambos extremos y el codón de stop**. Los campos vacíos de EPILAB significan que no se ha establecido una CDS propia completa para esa entrada. La restricción a CDS debe aplicarse y declararse por igual en las condiciones comparadas.

Los cDNA ya están orientados 5′→3′ como transcriptos; no se invierten de nuevo por la hebra genómica del locus. Se usa T en los FASTA nucleotídicos. Las herramientas pueden generar el complemento inverso o convertir T→U según sus requisitos.

Todas las coordenadas de los archivos derivados se refieren a `genes_9.fasta`. Los tramos opcionales incluyen tanto base 1 inclusiva como base 0 con fin exclusivo. Para convertir un resultado local base 0: sumar el inicio del tramo a su inicio y fin locales.

## Fragmentos originales de EPILAB

Las cuatro secuencias de la columna F se conservan intactas en `comparacion_EPILAB/fragmentos_originales.fasta`:

| Entrada | nt | Coordenadas base 1 inclusivas | Base 0, fin exclusivo |
|---|---:|---|---|
| MP_MRFV | 200 | 1–200 | [0, 200) |
| CP_MRFV | 200 | 201–400 | [200, 400) |
| BicC | 300 | 264–563 | [263, 563) |
| ATPasa D | 300 | 151–450 | [150, 450) |

En CP, la planilla declara 200–400. La ubicación exacta del fragmento es 201–400 al contar desde 1 e incluir ambos extremos. Se registran los valores declarados y los calculados, conservando la secuencia. No se infiere una única convención para toda la planilla.

Estos fragmentos sirven como referencias proporcionadas, no como demostración independiente de eficacia. Una comparación de regiones de 200 nt debe considerar que dos fragmentos originales tienen 300 nt.

## Análisis opcional de concordancia

`concordancia/tramos_concordantes.fasta` contiene **12 tramos de al menos 200 nt, procedentes de los mismos nueve genes**, que coinciden exactamente en los cotejos declarados y alinean dentro de sus límites codificantes. Suman 8.541 nt y 6.153 ventanas posibles de 200 nt por hebra. Su mapa y criterios están en esa misma carpeta.

Este conjunto permite explorar un escenario conservador o comprobar cuánto cambian los resultados al restringir la búsqueda. **No es la entrada predeterminada ni una condición para aceptar las secuencias de EPILAB.** Restringir el análisis a esos tramos puede excluir regiones útiles. En ATP sintasa β y V-ATPasa B, el respaldo completo es el modelo anotado y su proteína, no dos secuencias experimentales independientes.

Los fragmentos originales de EPILAB incluyen diferencias respecto de los comparadores utilizados y ninguno contiene una ventana de 200 nt que pase ese filtro estricto. Eso refleja el criterio elegido; no demuestra error o ineficacia de esos fragmentos. No se puede presentar una diferencia de resultados entre ambos espacios de búsqueda como superioridad biológica del filtrado.

Cada tramo opcional se procesa por separado: concatenarlos crearía uniones artificiales. Si se calcula accesibilidad o estructura dependiente de contexto, usar la secuencia fuente de `genes_9.fasta` y mapear las coordenadas, conservando el mismo contexto en todas las condiciones. Las UTR o extremos ausentes no se completan artificialmente.

## Benchmark y verificación

Para comparar las formas de búsqueda o las herramientas, fijar las mismas consultas, referencias, alcance, contexto y parámetros para todas las condiciones. La entrada disponible contiene 10.993 posiciones posibles de inicio de regiones de 200 nt por hebra. Una región de 200 nt contiene 180 ventanas solapadas de 21 nt por hebra; ese número fijo no indica cuántas serán eficaces.

Desde la carpeta completa, con Python 3.10 o posterior:

```sh
python verificar.py
python verificar.py --json
```

El verificador no requiere internet ni bibliotecas externas. Comprueba formato, IDs, longitudes, hashes, igualdad con las fuentes, CDS/proteínas públicas, originales EPILAB, coordenadas y tramos opcionales. Ante una diferencia termina con `FAIL` y código de salida 1.

El análisis de concordancia usó Biopython 1.88: alineamiento local con coincidencia +2, sustitución −1, apertura de hueco −5 y extensión −1; cobertura completa de la consulta y ubicación única. Dentro de esa ubicación se enumeraron todos los alineamientos globales óptimos; solo se conservaron correspondencias y adyacencias invariantes. Burs tuvo 3 alineamientos óptimos y MP tuvo 2. Los mapas y posiciones observadas están en `auditoria/`. No se calcularon E-values ni probabilidades de eficacia.

Este conjunto permite continuar la implementación y comparación computacional con integridad y procedencia verificables. Las referencias predichas, los límites de cobertura y las discrepancias documentadas condicionan la interpretación biológica; no se certifica la secuencia de toda una colonia, todas las isoformas o aislados, ni eficacia o especificidad experimental.

Fuentes públicas: [Andrada et al. y suplementos](https://www.frontiersin.org/journals/insect-science/articles/10.3389/finsc.2026.1811933/full), [ETH PZ371507.1](https://www.ebi.ac.uk/ena/browser/api/embl/PZ371507.1), [Burs-α PZ371508.1](https://www.ebi.ac.uk/ena/browser/api/embl/PZ371508.1), [Ensembl GCA_975739785.1 / 2026_02](https://ftp.ebi.ac.uk/pub/ensemblorganisms/GCA/975/739/785/1/ensembl/2026_02/geneset/), [MRFV AF265566.1](https://www.ebi.ac.uk/ena/browser/api/embl/AF265566.1). Las consultas, checksums y originales están documentados en `fuentes/`.
