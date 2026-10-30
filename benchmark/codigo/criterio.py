"""Define lexicographic ranking criteria for the RNAi benchmark.

The control_plaga mode prioritizes coverage within the target species after
minimizing external matches. The funcional mode penalizes matches to other
transcripts. These are operational priorities, not validated biological effects.
Activity-filter counts in si-Fi and alignment coverage in dsRNAEngineer remain
distinct measurements; see METODOLOGIA.md for their units and interpretation.
"""

MODO = "control_plaga"


def clave_sifi(fuera, intra, eficientes, modo=None):
    """Return the ascending sort key for a si-Fi window.

    Args:
        fuera: Sum of positions matching maize and bee references.
        intra: Positions matching other transcripts within the target reference.
        eficientes: Positions passing the si-Fi activity filters.
        modo: Ranking mode; the configured MODO is used when omitted.

    Returns:
        Tuple defining the lexicographic order, with the preferred key first.
    """
    m = modo or MODO
    if m == "control_plaga":
        return (fuera, -intra, -eficientes)
    return (fuera, intra, -eficientes)


def clave_dsrna(fuera, intra, on_target, modo=None):
    """Return the ascending sort key for a dsRNAEngineer window.

    Args:
        fuera: Sum of alignments against maize and bee references.
        intra: Alignments assigned outside the main target group.
        on_target: Alignments assigned to the main target group.
        modo: Ranking mode; the configured MODO is used when omitted.

    Returns:
        Tuple defining the lexicographic order. Alignment counts do not measure
        silencing efficacy.
    """
    m = modo or MODO
    if m == "control_plaga":
        return (fuera, -(on_target + intra))
    return (fuera, intra, -on_target)


_DESC = {
    "control_plaga": (
        "menos coincidencias externas; despues mayor cobertura intraespecie "
        "y, en si-Fi, mas posiciones que pasan los filtros de actividad"
    ),
    "funcional": (
        "menos coincidencias externas; despues menos coincidencias intraespecie; "
        "despues actividad predicha de si-Fi o cobertura on-target de dsRNAEngineer"
    ),
}
DESCRIPCION = "[modo %s] %s" % (MODO, _DESC[MODO])
