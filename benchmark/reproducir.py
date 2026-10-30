"""Recalculate benchmark analyses from the reference data tables."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


BENCHMARK = Path(__file__).resolve().parent
ROOT = BENCHMARK.parent
INPUT_NAMES = (
    "sifi_ventanas200.tsv", "siren_regiones200.tsv",
    "dsrnaengineer_ventanas200.tsv",
)
SUMMARY_NAMES = (
    "comparacion_ventanas.tsv", "comparacion_resumen.tsv",
    "comparacion_tres_ventanas.tsv", "comparacion_tres.tsv",
    "cruce_recomendaciones.tsv",
    "sensibilidad_criterio.tsv",
)


def digest(path):
    """Return a file's SHA-256.

    Args:
        path: File to read.

    Returns:
        Hexadecimal digest.
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare_tables(original, rebuilt):
    """Compare TSV content independently of platform line endings.

    Args:
        original: Reference table.
        rebuilt: Recalculated table.

    Returns:
        Content comparison with row counts and representative differences.
    """
    with original.open(encoding="utf-8", newline="") as handle:
        first = list(csv.DictReader(handle, delimiter="\t"))
    with rebuilt.open(encoding="utf-8", newline="") as handle:
        second = list(csv.DictReader(handle, delimiter="\t"))
    differences = []
    for number, (left, right) in enumerate(zip(first, second), start=1):
        if left != right and len(differences) < 3:
            differences.append({"fila": number, "referencia": left, "recalculado": right})
    return {
        "coincide": first == second,
        "filas_referencia": len(first),
        "filas_recalculadas": len(second),
        "diferencias_ejemplo": differences,
        "sha256_referencia": digest(original),
        "sha256_recalculado": digest(rebuilt),
    }


def run_stage(script, arguments, output, commands):
    """Execute an existing analysis and retain its stdout and stderr.

    Args:
        script: Filename under benchmark/codigo.
        arguments: Explicit CLI arguments.
        output: Destination for logs.
        commands: Mutable list of commands and outcomes.

    Returns:
        None.

    Raises:
        RuntimeError: The analysis returns an unsuccessful status.
    """
    command = [sys.executable, "-B", str(BENCHMARK / "codigo" / script), *arguments]
    environment = os.environ.copy()
    environment["PYTHONUTF8"] = "1"
    result = subprocess.run(command, cwd=ROOT, env=environment, capture_output=True)
    log = output / (Path(script).stem + ".log")
    log.write_bytes(result.stdout + result.stderr)
    commands.append({"script": script, "argumentos": arguments,
                     "codigo_salida": result.returncode,
                     "log": log.relative_to(output).as_posix()})
    if result.returncode:
        raise RuntimeError(f"{script} failed; inspect {log}")


def preflight(output):
    """Validate inputs and reserve a fresh output directory under salidas.

    Args:
        output: Requested output path.

    Returns:
        None.

    Raises:
        ValueError: The destination or an input is unsuitable.
    """
    permitted = (ROOT / "salidas").resolve()
    if output == permitted or not output.is_relative_to(permitted):
        raise ValueError("Use a new subdirectory of salidas/ for output")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("The output directory must be new or empty")
    for name in (*INPUT_NAMES, *SUMMARY_NAMES, "dsrnaengineer_resumen.tsv"):
        if not (BENCHMARK / "resultados" / name).is_file():
            raise ValueError(f"Missing reference table: {name}")
    fasta = ROOT / "seleccion_genes/secuencias/genes_9.fasta"
    names = [line[1:].strip() for line in fasta.read_text("ascii").splitlines()
             if line.startswith(">")]
    if len(names) != 9 or len(set(names)) != 9:
        raise ValueError("Expected nine unique FASTA records")
    for name in names:
        directory = BENCHMARK / "resultados/dsRNAEngineer_crudo/extraidos" / name / name
        for filename in ("input.fa", "dalbulus_maidis.21.bwt",
                         "zea_mays_l.21.bwt", "apis_mellifera.21.bwt"):
            if not (directory / filename).is_file():
                raise ValueError(f"Missing dsRNAEngineer evidence: {name}/{filename}")


def main():
    """Recalculate the saved analyses and report deviations explicitly.

    Args:
        None.

    Returns:
        Zero for matching outputs; one for execution failure; two for differences.
    """
    sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--salida", default="salidas/benchmark")
    args = parser.parse_args()
    requested = Path(args.salida)
    output = (requested if requested.is_absolute() else ROOT / requested).resolve()
    report = {"estado": "FAIL", "python": sys.version.split()[0],
              "alcance": "Reanálisis offline de tablas y salida guardada de dsRNAEngineer.",
              "comandos": [], "entradas": {}, "comparaciones": {}}
    try:
        preflight(output)
        output.mkdir(parents=True, exist_ok=True)
        for name in INPUT_NAMES:
            original = BENCHMARK / "resultados" / name
            shutil.copyfile(original, output / name)
            report["entradas"][name] = digest(original)
        reconstructed = output / "reconstruccion_dsrna"
        reconstructed.mkdir()
        run_stage("parsear_dsrnaengineer.py", [
            "--extraidos", str(BENCHMARK / "resultados/dsRNAEngineer_crudo/extraidos"),
            "--genes", str(ROOT / "seleccion_genes/secuencias/genes_9.fasta"),
            "--out-dir", str(reconstructed),
        ], output, report["comandos"])
        for name in ("dsrnaengineer_ventanas200.tsv", "dsrnaengineer_resumen.tsv"):
            report["comparaciones"]["reconstruccion/" + name] = compare_tables(
                BENCHMARK / "resultados" / name, reconstructed / name)
        for script in ("comparar_herramientas.py", "comparar_tres.py",
                       "cruce_recomendaciones.py", "sensibilidad_criterio.py"):
            run_stage(script, ["--dir", str(output)], output, report["comandos"])
        for name in SUMMARY_NAMES:
            report["comparaciones"][name] = compare_tables(
                BENCHMARK / "resultados" / name, output / name)
        report["estado"] = ("PASS" if all(row["coincide"] for row in
                                           report["comparaciones"].values())
                            else "DIFERENCIAS")
    except (OSError, ValueError, RuntimeError) as exc:
        report["error"] = str(exc)
        print(f"FAIL: {exc}", file=sys.stderr)
    if output.is_dir() and report["comandos"]:
        (output / "reproduccion.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(report["estado"])
    for name, row in report["comparaciones"].items():
        print(f"{'OK' if row['coincide'] else 'DIFF'}: {name}")
    if report["estado"] == "PASS":
        return 0
    return 2 if report["estado"] == "DIFERENCIAS" else 1


if __name__ == "__main__":
    sys.exit(main())
