#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Aggregate dsRNAEngineer Bowtie alignments on the benchmark coordinate grid.

Read input.fa and the alignment files for Dalbulus maidis, Zea mays and
Apis mellifera. Compare each query with the reference FASTA. Group Trinity
isoforms by their gene prefix and assign the dominant group operationally as
the intended target. Alignment rows and unique siRNA positions are different
units and are retained as such in the regional and sequence summaries.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

LARGO_REGION = 200
PASO = 4
ESPECIES = {
    "dalbulus_maidis": "objetivo",
    "zea_mays_l": "zmays",
    "apis_mellifera": "amellifera",
}


def leer_fasta_simple(ruta):
    seqs = {}
    nombre = None
    partes = []
    with open(ruta) as fh:
        for linea in fh:
            linea = linea.strip()
            if linea.startswith(">"):
                if nombre:
                    seqs[nombre] = "".join(partes).upper()
                nombre = linea[1:].split()[0]
                partes = []
            elif linea:
                partes.append(linea)
    if nombre:
        seqs[nombre] = "".join(partes).upper()
    return seqs


def leer_bwt(ruta):
    """[(ini, fin, contig)] con coordenadas base 1 sobre nuestra secuencia."""
    out = []
    if not os.path.exists(ruta) or os.path.getsize(ruta) == 0:
        return out
    with open(ruta) as fh:
        for linea in fh:
            c = linea.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            try:
                out.append((int(c[1]), int(c[2]), c[4]))
            except ValueError:
                continue
    return out


def gen_trinity(contig):
    """
    'TRINITY_DN64958_c0_g1_i1' -> 'TRINITY_DN64958_c0_g1'.
    Trinity nombra <gen>_i<N> para cada isoforma. Si el nombre no sigue ese
    patron se devuelve tal cual, para no romper con otras referencias.
    """
    i = contig.rfind("_i")
    if i > 0 and contig[i + 2:].isdigit():
        return contig[:i]
    return contig


def rejilla(largo):
    return [(i + 1, i + LARGO_REGION)
            for i in range(0, max(0, largo - LARGO_REGION + 1), PASO)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--extraidos", required=True,
                    help="carpeta con <gen>/<gen>/*.bwt")
    ap.add_argument("--genes", required=True, help="genes_9.fasta de referencia")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    ref = leer_fasta_simple(args.genes)
    os.makedirs(args.out_dir, exist_ok=True)

    filas = []
    resumen = []

    for nombre in sorted(os.listdir(args.extraidos)):
        base = os.path.join(args.extraidos, nombre, nombre)
        if not os.path.isdir(base):
            base = os.path.join(args.extraidos, nombre)
        f_in = os.path.join(base, "input.fa")
        if not os.path.exists(f_in):
            continue
        gen = nombre

        # integridad: lo que se subio tiene que ser nuestra secuencia
        subida = next(iter(leer_fasta_simple(f_in).values()), "")
        nuestra = ref.get(gen, "")
        if not nuestra:
            estado_seq = "gen ausente de genes_9.fasta"
        elif subida == nuestra:
            estado_seq = "OK"
        else:
            estado_seq = ("DIFIERE (subida %d nt vs nuestra %d nt)"
                          % (len(subida), len(nuestra)))
        largo = len(subida) or len(nuestra)

        hits = {}
        for esp, etiqueta in ESPECIES.items():
            hits[etiqueta] = leer_bwt(os.path.join(base, "%s.21.bwt" % esp))

        # separar on-target dentro de la especie objetivo, AGRUPANDO POR GEN
        # Trinity y no por contig: los nombres son TRINITY_DNxxxx_c0_g1_i1,
        # donde _i1 es la ISOFORMA. Dos isoformas del mismo gen son el mismo
        # gen, y contarlas por separado hacia que una apareciera como
        # off-target de la otra. En DM_BURS_A eso daba un empate exacto de
        # 419 contra 419, que era justamente el mismo gen partido en dos.
        por_contig = defaultdict(set)
        for (a, b, contig) in hits["objetivo"]:
            por_contig[gen_trinity(contig)].add((a, b))
        orden = sorted(por_contig.items(), key=lambda x: (-len(x[1]), x[0]))
        on_contig = orden[0][0] if orden else None
        n1 = len(orden[0][1]) if orden else 0
        n2 = len(orden[1][1]) if len(orden) > 1 else 0

        on_hits = [(a, b) for (a, b, c) in hits["objetivo"]
                   if gen_trinity(c) == on_contig]
        off_intra = [(a, b) for (a, b, c) in hits["objetivo"]
                     if gen_trinity(c) != on_contig]

        resumen.append({
            "gen": gen, "largo": largo, "secuencia_subida": estado_seq,
            "hits_objetivo": len(hits["objetivo"]),
            "contigs_objetivo": len(por_contig),
            "on_target_contig": on_contig or "",
            "on_target_sirnas": n1, "segundo_contig_sirnas": n2,
            "off_intra_hits": len(off_intra),
            "hits_zmays": len(hits["zmays"]),
            "hits_amellifera": len(hits["amellifera"]),
        })

        for (ini, fin) in rejilla(largo):
            def cuenta(lista):
                return sum(1 for (a, b) in lista if a >= ini and b <= fin)
            filas.append({
                "gen": gen, "inicio": ini, "fin": fin,
                "ds_on_target": cuenta(on_hits),
                "ds_off_intra": cuenta(off_intra),
                "ds_off_zmays": cuenta([(a, b) for (a, b, _) in hits["zmays"]]),
                "ds_off_amellifera": cuenta([(a, b) for (a, b, _) in hits["amellifera"]]),
            })

    if not filas:
        print("no se encontro ninguna salida utilizable", file=sys.stderr)
        return 1

    campos = ["gen", "inicio", "fin", "ds_on_target", "ds_off_intra",
              "ds_off_zmays", "ds_off_amellifera"]
    ruta = os.path.join(args.out_dir, "dsrnaengineer_ventanas200.tsv")
    with open(ruta, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader()
        w.writerows(filas)
    print("escrito %s (%d ventanas, %d genes)"
          % (ruta, len(filas), len(resumen)))

    campos_r = list(resumen[0].keys())
    ruta_r = os.path.join(args.out_dir, "dsrnaengineer_resumen.tsv")
    with open(ruta_r, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos_r, delimiter="\t")
        w.writeheader()
        w.writerows(resumen)
    print("escrito %s" % ruta_r)

    print("\n%-19s %-6s %-10s %-9s %-8s %-7s %-6s %s"
          % ("gen", "largo", "seq", "hits_obj", "on_tgt", "2do", "maiz", "abeja"))
    print("-" * 92)
    for r in resumen:
        marca = ""
        if r["on_target_sirnas"] and r["segundo_contig_sirnas"]:
            if r["on_target_sirnas"] < 3 * r["segundo_contig_sirnas"]:
                marca = "  <-- sin pico dominante: separacion on/off dudosa"
        if r["hits_objetivo"] == 0:
            marca = "  <-- SIN HITS en ninguna especie"
        print("%-19s %-6d %-10s %-9d %-8d %-7d %-6d %-6d%s"
              % (r["gen"], r["largo"], r["secuencia_subida"][:10],
                 r["hits_objetivo"], r["on_target_sirnas"],
                 r["segundo_contig_sirnas"], r["hits_zmays"],
                 r["hits_amellifera"], marca))
    return 0


if __name__ == "__main__":
    sys.exit(main())
