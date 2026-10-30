"""Execute panel admission, preparation, assessment and export."""

import argparse
from pathlib import Path

from evaluar_interacciones import prepare, repair, run
from seleccion_panel import freeze
from exportar_resultados import export


def main():
    """Dispatch explicit freeze, preparation, assessment and export phases.

    Args:
        None. Arguments come from the command line.

    Returns:
        None.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=(
        "freeze", "prepare", "repair", "run", "export"))
    repository = Path(__file__).resolve().parents[2]
    parser.add_argument("--output-dir", type=Path,
                        default=repository / "salidas/algoritmo_nativo")
    parser.add_argument("--per-gene", type=int, default=3)
    parser.add_argument("--wall-seconds", type=float, default=3600)
    parser.add_argument("--max-pairs", type=int, default=20000)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    published = repository / "algoritmo/resultados"
    if output.is_relative_to(published):
        parser.error("Las tablas publicadas se preservan; use una salida nueva.")
    implementation = Path(__file__).resolve().parent
    required = (
        "references/2026-10-28/manifest.json",
        "datos_preparados/enumeracion/parents.json",
        "datos_preparados/asimetria/asymmetry.json",
        "sequences/genes_9.fasta",
    )
    missing = [name for name in required if not (implementation / name).is_file()]
    if missing:
        parser.error(
            "La evaluacion nativa requiere activos externos: "
            + ", ".join(missing)
            + ". Consulte algoritmo/REPRODUCIBILIDAD.md. "
            "Para recalcular las tablas incluidas: python algoritmo/reproducir.py.")
    if args.action == "freeze":
        result = freeze(output, args.per_gene)["expected_coverage"]
    elif args.action == "prepare":
        result = prepare(output)
    elif args.action == "repair":
        result = repair(output)
    elif args.action == "run":
        run(output, args.wall_seconds, args.max_pairs)
        result = export(output)
    else:
        result = export(output)
    print(result)


if __name__ == "__main__":
    main()
