# Informe de resultados

[benchmark_rnai.html](benchmark_rnai.html) presenta el objetivo, la metodología,
los resultados y la interpretación del benchmark. Los estilos y la figura están
incrustados en el HTML, que puede descargarse y abrirse por separado en un
navegador. Los enlaces a los datos apuntan al repositorio del proyecto.

El documento se apoya en las [tablas de resultados](../resultados/README.md).
La figura de concordancia utiliza el Jaccard del 10 % superior con una escala
fija de 0 a 1, y diferencia los valores no evaluables de los valores iguales a cero.

## Regenerar la figura

La generación de la figura es independiente del análisis de tablas. Requiere
las dependencias de `requirements_figura.txt`:

```bash
python -m pip install -r benchmark/informe/requirements_figura.txt
python benchmark/informe/generar_figura.py
```

El programa utiliza `comparacion_tres.tsv`, escribe
`figuras/concordancia_top10.svg` y actualiza la figura incrustada en el HTML.
El SVG también está disponible para usarlo como archivo independiente.
El análisis principal del benchmark se ejecuta según la
[guía de reproducibilidad](../REPRODUCIBILIDAD.md).
