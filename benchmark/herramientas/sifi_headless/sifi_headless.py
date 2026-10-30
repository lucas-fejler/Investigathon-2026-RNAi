# -*- coding: utf-8 -*-
"""Driver headless de si-Fi21 (Lueck et al. 2019, Front. Plant Sci. 10:1023).

Corre el codigo original de si-Fi (Python 2) sin su interfaz Qt.

Replica paso por paso SifiPipeline.run_pipeline, con dos unicas diferencias,
ambas de interfaz y no de calculo:

  1. get_main_target() abre un dialogo Qt para que el usuario elija cuales de
     los hits de bowtie son el blanco principal. Aca esa eleccion se pasa por
     parametro (--main-target-regex, o 'auto' = el hit con mas siRNAs, que es
     lo que elige un usuario en la GUI).
  2. show_plot.DrawPlot (matplotlib/Qt) no se invoca. El driver consume
     json_lst, la misma estructura que en el original alimenta a los graficos,
     y la escribe a disco.

Todo el calculo -- creacion de siRNAs, bowtie, RNAplfold, seleccion de hebra,
estabilidad de extremos, accesibilidad del sitio y el veredicto de eficiencia --
usa los metodos originales sin modificar.
"""

import argparse
import json
import os
import re
import sys
import tempfile
from collections import Counter

from Bio import SeqIO


def build_parser():
    p = argparse.ArgumentParser(description="si-Fi21 headless")
    p.add_argument("--sifi-dir", required=True, help="carpeta con sifi_pipeline.py")
    p.add_argument("--query", required=True, help="FASTA de consulta (multi-FASTA)")
    p.add_argument("--db-location", required=True, help="carpeta del indice bowtie (con / final)")
    p.add_argument("--bowtie-db", required=True, help="nombre base del indice bowtie")
    p.add_argument("--bowtie-location", required=True, help="carpeta del ejecutable bowtie")
    p.add_argument("--rnaplfold-location", required=True, help="carpeta del ejecutable RNAplfold")
    p.add_argument("--out-json", required=True)
    p.add_argument("--out-tsv", required=True)
    # Parametros de si-Fi. Los valores por defecto son los de la GUI en modo
    # 'RNAi design' (main.py set_design_settings + Resources/sifi2015.ui).
    p.add_argument("--mode", type=int, default=0, help="0=RNAi design, 1=off-target")
    p.add_argument("--sirna-size", type=int, default=21)
    p.add_argument("--mismatches", type=int, default=0)
    p.add_argument("--accessibility-check", type=int, default=2)
    p.add_argument("--accessibility-window", type=int, default=8)
    p.add_argument("--strand-check", type=int, default=2)
    p.add_argument("--end-check", type=int, default=2)
    p.add_argument("--end-stability-threshold", type=float, default=1.0)
    p.add_argument("--accessibility-threshold", type=float, default=0.1)
    p.add_argument("--terminal-check", type=int, default=2)
    p.add_argument("--no-efficience", action="store_true")
    p.add_argument("--main-target-regex", default="auto",
                   help="regex sobre nombres de hits que son blanco principal; "
                        "'auto' = el hit con mas siRNAs")
    return p


def pick_main_targets(bowtie_data_l, regex):
    """Reemplaza el dialogo Qt get_main_target() por una seleccion explicita.

    Devuelve la lista de nombres de secuencias de la base que cuentan como
    blanco principal. Todo hit restante queda como off-target, igual que en
    el original."""
    counts = Counter(hit[2] for hit in bowtie_data_l)
    if not counts:
        return None
    ordered = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    if regex == "auto":
        return [ordered[0][0]]
    rx = re.compile(regex)
    chosen = [name for name, _ in ordered if rx.search(name)]
    return chosen or None


def main():
    args = build_parser().parse_args()

    sys.path.insert(0, args.sifi_dir)
    import sifi_pipeline
    import general_helpers  # noqa: F401  (lo usa el original para tablas)

    temp_location = tempfile.mkdtemp()

    # Asegurar que bowtie y RNAplfold resuelvan por PATH: el original los
    # invoca por nombre desnudo y en Windows confia en el chdir.
    os.environ["PATH"] = (args.bowtie_location + os.pathsep +
                          args.rnaplfold_location + os.pathsep +
                          os.environ.get("PATH", ""))

    sp = sifi_pipeline.SifiPipeline(
        args.bowtie_db,
        args.db_location,
        args.query,
        args.sirna_size,
        args.mismatches,
        args.accessibility_check,
        args.accessibility_window,
        args.rnaplfold_location,
        args.bowtie_location,
        args.mode,
        args.strand_check,
        args.end_check,
        args.end_stability_threshold,
        args.accessibility_threshold,
        temp_location,
        args.terminal_check,
        args.no_efficience,
    )

    all_records = []
    summary = []

    for seq_record in SeqIO.parse(args.query, "fasta"):
        query_name = seq_record.id
        query_sequence = str(seq_record.seq)
        sp.len_seq = len(query_sequence)

        sys.stderr.write("[si-Fi] %s (%d nt)\n" % (query_name, sp.len_seq))

        query_seq_temp_file = sp.create_single_fasta_file(query_name, query_sequence)
        sirna_file, tab_file_name = sp.create_sirnas(query_sequence, sp.sirna_size)

        bowtie_data = sp.run_bowtie(sirna_file, sp.bowtie_db, sp.mismatches)
        sp.bowtie_data_l = sp.bowtie_to_lst(bowtie_data)

        lunp_data = sp.run_rnaplfold(query_name, query_seq_temp_file)

        if sp.mode == 0:
            main_targets = pick_main_targets(sp.bowtie_data_l, args.main_target_regex)
            no_target = main_targets is None
        else:
            main_targets = None
            no_target = False

        sp.sirna_l = [x.split("\t") for x in open(tab_file_name, "r").readlines()]

        if sp.bowtie_data_l:
            input_data = sp.bowtie_data_l
            no_target = False
        else:
            if sp.mode == 0:
                input_data = [x.split("\t") for x in open(tab_file_name, "r").readlines()]
                no_target = True
            else:
                input_data = None

        if not input_data:
            sys.stderr.write("  sin hits en la base; se omite %s\n" % query_name)
            summary.append({"query": query_name, "length": sp.len_seq,
                            "n_json": 0, "main_targets": main_targets,
                            "note": "sin hits en la base"})
            continue

        json_lst = sp.data_to_json(query_name, input_data, no_target, lunp_data, main_targets)

        for row in json_lst:
            row["_query"] = query_name
        all_records.extend(json_lst)
        summary.append({"query": query_name, "length": sp.len_seq,
                        "n_json": len(json_lst), "main_targets": main_targets,
                        "n_bowtie_hits": len(sp.bowtie_data_l)})
        sys.stderr.write("  hits bowtie=%d  filas json=%d  blanco_principal=%s\n"
                         % (len(sp.bowtie_data_l), len(json_lst), main_targets))

    with open(args.out_json, "w") as fh:
        json.dump({"parameters": vars(args), "summary": summary,
                   "records": all_records}, fh, indent=1)

    # TSV plano con todas las claves observadas
    keys = []
    for row in all_records:
        for k in row.keys():
            if k not in keys:
                keys.append(k)
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for row in all_records:
            fh.write("\t".join([unicode(row.get(k, "")).encode("utf-8")
                                for k in keys]) + "\n")

    sys.stderr.write("\n[si-Fi] listo: %d filas -> %s\n" % (len(all_records), args.out_tsv))


if __name__ == "__main__":
    main()
