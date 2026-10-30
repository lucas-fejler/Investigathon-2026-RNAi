# Selección de genes y secuencias de referencia

La selección de blancos para interferencia por ARN (RNAi) parte de cuatro
referencias de EPILAB: BicC y V-ATPasa D de *Dalbulus maidis*, y MP putativa y CP
del virus del rayado fino del maíz (MRFV). Se incorporaron cinco candidatos de
*D. maidis* para explorar blancos vinculados con la muda, la expansión alar,
la producción de ATP y el funcionamiento celular.

Se consideraron tres criterios:

1. **Antecedentes funcionales de RNAi:** qué se observó al silenciar el gen,
   distinguiendo especie, estadio, vía de administración y variable medida.
2. **Función biológica:** qué proceso podría verse afectado al reducir su
   expresión.
3. **Secuencias trazables:** qué referencia de *D. maidis* se utiliza, cómo se
   revisó su correspondencia con el gen y dónde se registra su origen.

## Cinco candidatos adicionales

| Candidato | Función de interés | Respaldo y límite principal |
| --- | --- | --- |
| **ETH** | Coordinación de la ecdisis. | La microinyección en ninfas de quinto estadio de *D. maidis* produjo alteraciones de la muda y menor supervivencia. Es evidencia de laboratorio por inyección. [1] |
| **Burs-α** | Expansión alar y procesos posteriores a la muda. | Su silenciamiento por microinyección en *D. maidis* alteró la expansión de las alas y redujo la supervivencia. El depósito y el suplemento del estudio presentan diferencias documentadas. [1] |
| **ATP sintasa β** | Producción mitocondrial de ATP. | La inyección produjo esterilidad femenina en *Euscelidius variegatus*. Otro estudio obtuvo reducción del transcripto mediante alimentación sobre hojas tratadas en esa especie y en *Scaphoideus titanus*. La evidencia por alimentación no demuestra el mismo efecto reproductivo. [2, 3] |
| **Rpn6 (`DM_RPN6`, denominación del antecedente)** | Componente de la maquinaria de eliminación controlada de proteínas. | Canuto y colaboradores observan efectos sobre supervivencia y reproducción en *E. variegatus* y *S. titanus* tras microinyección. La referencia de *D. maidis* se eligió mediante comparación proteica con las secuencias citadas en el estudio. [4] |
| **V-ATPasa B** | Subunidad de una bomba que usa ATP para transportar protones y regular la acidez celular. | Antecedente de silenciamiento oral en *D. maidis*. La versión consultada es un preprint con alta mortalidad en los controles. B y D son genes diferentes del mismo complejo; ambos fueron evaluados en ese trabajo. [5] |

Los antecedentes de ATP sintasa β y del candidato del proteasoma proceden de
otros cicadélidos. Las comparaciones de proteínas apoyan la correspondencia de
secuencias; no trasladan automáticamente el efecto experimental a *D. maidis*.

## Secuencias para trabajar

**Entrada principal:** [genes_9.fasta](secuencias/genes_9.fasta).
Contiene nueve registros y 12.784 nt, con A/C/G/T y la longitud disponible de
cada fuente elegida. Los cuatro originales de EPILAB conservan sus bases.

| ID del FASTA | Entrada | nt |
| --- | --- | ---: |
| `DM_ETH` | ETH | 615 |
| `DM_BURS_A` | Burs-α | 468 |
| `DM_ATP5F1B` | ATP sintasa β | 2.714 |
| `DM_RPN6` | Rpn6, denominación del antecedente | 2.513 |
| `DM_VATP_B` | V-ATPasa B | 1.672 |
| `EPILAB_BicC_Dm` | BicC | 2.520 |
| `EPILAB_ATPaseD_Dm` | V-ATPasa D | 491 |
| `EPILAB_MP_MRFV` | MP putativa de MRFV | 1.263 |
| `EPILAB_CP_MRFV` | CP de MRFV | 528 |

La [documentación de secuencias](secuencias/README.md) describe las
referencias, discrepancias, CDS, UTR, orientación y coordenadas. Los
[metadatos](secuencias/metadatos_genes.tsv) fijan el
identificador y el origen de cada registro. Los recortes de `concordancia/`
son un análisis opcional y no reemplazan al FASTA principal.

Los registros de procedencia y los hashes permiten identificar la referencia
utilizada para cada secuencia y comprobar su integridad.

## Correspondencia de nombres y secuencias

Los nombres de la tabla remiten a los antecedentes utilizados para seleccionar
los candidatos. En el caso del proteasoma, Rpn6 es la denominación de Canuto et
al.; sus referencias son `GJQL01135962.1` y `GFTU01006765.1`. La secuencia
seleccionada de *D. maidis* queda fijada por `ENSCHMT00000015190.1` y el ID
`DM_RPN6`. El identificador, la secuencia y la comparación con el antecedente
permiten seguir el blanco aunque las bases empleen nomenclaturas distintas.
La comparación apoya la correspondencia de secuencias; no constituye por sí
sola una reconstrucción formal de ortología ni una validación funcional local.

## Referencias

1. Andrada NL et al. (2026). *ETH and Burs-α are necessary for the normal molting in Dalbulus maidis and Delphacodes kuscheli*. [Frontiers in Insect Science](https://doi.org/10.3389/finsc.2026.1811933).
2. Galetto L et al. (2021). *Silencing of ATP synthase β induces female sterility in a leafhopper phytoplasma vector*. [Entomologia Generalis](https://doi.org/10.1127/entomologia/2021/1267).
3. Rossi M et al. (2024). *Lab-scale method for plant-mediated delivery of dsRNAs to phloem-feeding leafhoppers*. [Journal of Pest Science](https://doi.org/10.1007/s10340-023-01670-0).
4. Canuto F et al. (2025). *RNAi-mediated silencing of a regulatory subunit of the 26S proteasome induces mortality and female sterility in leafhopper vectors of grapevine Flavescence dorée phytoplasmas*. [Journal of Pest Science](https://doi.org/10.1007/s10340-025-01961-8).
5. Jones TKL et al. (2021). *Assessing the Functionality of RNA Interference (RNAi) in the Phloem-feeding Maize pest Dalbulus maidis*. [bioRxiv, versión consultada](https://doi.org/10.1101/2021.09.29.462424).
