#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Aggregate si-Fi21 indicators on 200-nt windows with a 4-nt step.

Count distinct siRNA positions only when their full 21-nt interval lies inside
the region. Export activity-filter and reference-match indicators, using
one-based inclusive coordinates over the unmodified query sequences.
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from criterio import clave_sifi  # noqa: E402

# si-Fi escribe el TSV desde un dict, asi que el ORDEN de columnas no es
# estable entre corridas. Siempre parsear por nombre de columna.
COL_POS = "sirna_position"
COL_QUERY = "query_name"
COL_EFF = "is_efficient"
COL_OFF = "is_off_target"

SIRNA_LEN = 21
WINDOW = 200
# Paso de la ventana deslizante. El default es 4 para que la rejilla coincida
# EXACTAMENTE con la que usa SIREN en generate_rnai_sequences() (paso 4, base 1).
# Compartir rejilla permite unir las dos herramientas por (gen, inicio, fin) sin
# interpolar ni asignar la ventana "mas cercana", que meteria error propio.
STEP = 4

# Las dos corridas de "design" (modo 0). Cada gen se analizo contra la base
# donde vive su blanco real: los genes del insecto contra el cDNA de
# D. maidis, los dos genes virales contra el genoma de MRFV.
DESIGN_RUNS = ["sifi_design_dmaidis", "sifi_design_mrfv"]

# Las corridas de off-target (modo 1) contra especies que NO queremos
# afectar. Aqui toda fila es, por construccion, un hit no deseado: el campo
# is_off_target no sirve para filtrar porque si-Fi no tiene un blanco
# principal que elegir en esta base.
OFFTARGET_RUNS = {
    "zmays": "sifi_offtarget_zmays",
    "amellifera": "sifi_offtarget_amellifera",
}


def leer_fasta(ruta):
    """Devuelve {nombre: longitud} del FASTA."""
    largos = {}
    nombre = None
    n = 0
    with open(ruta) as fh:
        for linea in fh:
            linea = linea.strip()
            if not linea:
                continue
            if linea.startswith(">"):
                if nombre is not None:
                    largos[nombre] = n
                nombre = linea[1:].split()[0]
                n = 0
            else:
                n += len(linea)
    if nombre is not None:
        largos[nombre] = n
    return largos


def cargar_design(out_dir):
    """
    De las corridas de design saca, por gen:
      eficientes       posiciones (base 1) con is_efficient = True
      off_intra        posiciones con is_off_target = True (parálogos de la
                       propia especie / otro transcripto de la misma base)
    Se usan conjuntos de POSICIONES, no conteo de filas: el TSV tiene una
    fila por hit, asi que un mismo siRNA aparece repetido y contar filas
    sobreponderaria los siRNAs promiscuos.
    """
    eficientes = defaultdict(set)
    off_intra = defaultdict(set)
    for run in DESIGN_RUNS:
        ruta = os.path.join(out_dir, run + ".tsv")
        if not os.path.exists(ruta):
            print("  aviso: falta %s, se omite" % ruta, file=sys.stderr)
            continue
        with open(ruta) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gen = row[COL_QUERY]
                pos = int(row[COL_POS])
                if row[COL_EFF] == "True":
                    eficientes[gen].add(pos)
                if row[COL_OFF] == "True":
                    off_intra[gen].add(pos)
    return eficientes, off_intra


def cargar_offtarget(out_dir):
    """
    Por cada especie no objetivo, {gen: {posiciones con al menos un hit}}.
    Igual que arriba: posiciones, no filas.
    """
    por_especie = {}
    for especie, run in OFFTARGET_RUNS.items():
        ruta = os.path.join(out_dir, run + ".tsv")
        pos_por_gen = defaultdict(set)
        if not os.path.exists(ruta):
            print("  aviso: falta %s, se omite" % ruta, file=sys.stderr)
        else:
            with open(ruta) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    pos_por_gen[row[COL_QUERY]].add(int(row[COL_POS]))
        por_especie[especie] = pos_por_gen
    return por_especie


def genes_sin_base(out_dir):
    """
    Genes que en una corrida de design quedaron SIN hits (main_targets null).
    Importa distinguirlos: cero off-targets ahi no es especificidad, es que
    la herramienta no tenia contra que comparar. Se reporta como tal.
    """
    sin_base = {}
    for run in DESIGN_RUNS:
        ruta = os.path.join(out_dir, run + ".json")
        if not os.path.exists(ruta):
            continue
        with open(ruta) as fh:
            d = json.load(fh)
        for s in d.get("summary", []):
            if s.get("n_json", 0) == 0:
                # puede tener hits en la OTRA corrida de design
                sin_base.setdefault(s["query"], True)
            else:
                sin_base[s["query"]] = False
    return {g: v for g, v in sin_base.items() if v}


def ventanas(largo):
    """
    Ventanas de 200 nt, paso STEP, base 1 inclusive.

    Se replica la rejilla de SIREN tal cual: range(0, largo-200+1, 4), sin
    anclar una ventana extra al extremo 3'. SIREN no la agrega, y agregarla
    aqui dejaria filas que no se pueden unir con las suyas.
    """
    if largo < WINDOW:
        return [(1, largo)]
    return [(i + 1, i + WINDOW) for i in range(0, largo - WINDOW + 1, STEP)]


def contar(posiciones, start, end):
    """Cuantas posiciones dan un siRNA de 21 nt contenido en [start, end]."""
    return sum(1 for p in posiciones if p >= start and p + SIRNA_LEN - 1 <= end)


def rankear(filas):
    """
    Ordena las ventanas con EL criterio del benchmark, que vive en criterio.py.
    No se reimplementa aqui a proposito: tenerlo escrito en varios lugares ya
    produjo dos versiones distintas conviviendo.
    """
    return sorted(
        filas,
        key=lambda f: clave_sifi(f["off_no_objetivo"], f["off_intra"],
                                 f["eficientes"]),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sifi-out", default=os.path.expanduser("~/rnai_bench/out"),
                    help="carpeta con los .tsv/.json de si-Fi")
    ap.add_argument("--genes", default=os.path.expanduser("~/rnai_bench/work/genes_9.fasta"))
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    largos = leer_fasta(args.genes)
    print("genes leidos: %d" % len(largos))

    eficientes, off_intra = cargar_design(args.sifi_out)
    off_por_especie = cargar_offtarget(args.sifi_out)
    sin_base = genes_sin_base(args.sifi_out)

    filas = []
    for gen in sorted(largos):
        largo = largos[gen]
        for (start, end) in ventanas(largo):
            o_zm = contar(off_por_especie["zmays"].get(gen, ()), start, end)
            o_am = contar(off_por_especie["amellifera"].get(gen, ()), start, end)
            filas.append({
                "gen": gen,
                "largo_gen": largo,
                "inicio": start,
                "fin": end,
                "eficientes": contar(eficientes.get(gen, ()), start, end),
                "off_intra": contar(off_intra.get(gen, ()), start, end),
                "off_zmays": o_zm,
                "off_amellifera": o_am,
                "off_no_objetivo": o_zm + o_am,
            })

    campos = ["gen", "largo_gen", "inicio", "fin", "eficientes",
              "off_intra", "off_zmays", "off_amellifera", "off_no_objetivo"]

    ruta_v = os.path.join(args.out_dir, "sifi_ventanas200.tsv")
    with open(ruta_v, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader()
        for f in filas:
            w.writerow(f)
    print("escrito %s (%d ventanas)" % (ruta_v, len(filas)))

    # Mejor ventana por gen
    ruta_m = os.path.join(args.out_dir, "sifi_mejor_ventana.tsv")
    por_gen = defaultdict(list)
    for f in filas:
        por_gen[f["gen"]].append(f)
    with open(ruta_m, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos + ["nota"], delimiter="\t")
        w.writeheader()
        for gen in sorted(por_gen):
            mejor = rankear(por_gen[gen])[0]
            mejor = dict(mejor)
            mejor["nota"] = ("SIN BASE DE COMPARACION: si-Fi no encontro "
                             "hits, los ceros no son especificidad"
                             if sin_base.get(gen) else "")
            w.writerow(mejor)
    print("escrito %s" % ruta_m)

    if sin_base:
        print("\nAVISO - genes sin hits en su base de design: %s"
              % ", ".join(sorted(sin_base)))
        print("  Para esos genes un cero de off-target NO significa que sean")
        print("  especificos: significa que no habia contra que comparar.")


if __name__ == "__main__":
    main()
