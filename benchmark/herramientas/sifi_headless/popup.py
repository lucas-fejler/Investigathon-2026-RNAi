# Stub sin GUI de popup.py (si-Fi).
# El original abre dialogos Qt. El driver headless nunca los invoca:
# get_main_target() se reemplaza por una seleccion explicita.


class ListSelection(object):
    def __init__(self, *args, **kwargs):
        raise RuntimeError(
            "popup.ListSelection es un dialogo Qt: el modo headless debe "
            "pasar main_targets explicitamente en lugar de llamar get_main_target()."
        )

    def exec_(self):
        raise RuntimeError("sin GUI en modo headless")


def show_info_message(*args, **kwargs):
    pass
