# Herramientas y acceso al software

| Herramienta | Fuente | Uso en el benchmark |
| --- | --- | --- |
| si-Fi21 | [Publicación](https://doi.org/10.3389/fpls.2019.01023) · [Código](https://github.com/snowformatics/siFi21-) | Indicadores por siRNA, filtros de actividad y coincidencias con referencias. |
| SIREN | [Código y documentación](https://github.com/pablovargasmejia/SIREN) | Evaluación con RNAhybrid y score de regiones. |
| dsRNAEngineer | [Publicación](https://doi.org/10.1016/j.tibtech.2025.01.002) · [Plataforma](https://www.dsrna-engineer.cn/) | Alineamientos sobre referencias del portal. |

## Controlador de si-Fi21

`sifi_headless/sifi_headless.py` permite ejecutar una instalación de si-Fi21 sin
su interfaz gráfica. Recibe la ruta de esa instalación, las secuencias de
consulta, la base de Bowtie y las rutas de las herramientas externas.

```bash
python benchmark/herramientas/sifi_headless/sifi_headless.py --help
```

El controlador depende del entorno compatible con si-Fi21, Bowtie y RNAplfold;
es independiente del análisis de tablas, que utiliza la biblioteca estándar.
Los módulos `show_plot.py` y `popup.py` proporcionan adaptaciones de interfaz
para la ejecución sin entorno gráfico.

Las versiones, parámetros y referencias de las herramientas deben fijarse al
repetir una ejecución. El alcance reproducible con los archivos de este
repositorio se describe en [REPRODUCIBILIDAD.md](../REPRODUCIBILIDAD.md).
