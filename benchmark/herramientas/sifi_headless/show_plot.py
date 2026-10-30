# Stub sin GUI de show_plot.py (si-Fi).
# El original dibuja con matplotlib/Qt. El driver headless corta antes de graficar:
# toma el json_lst que data_to_json() produce, que es la misma estructura que
# alimenta a los graficos en la version original.


class DrawPlot(object):
    def __init__(self, *args, **kwargs):
        raise RuntimeError(
            "show_plot.DrawPlot requiere Qt: el modo headless consume json_lst "
            "directamente y no grafica."
        )

    def show(self):
        raise RuntimeError("sin GUI en modo headless")
