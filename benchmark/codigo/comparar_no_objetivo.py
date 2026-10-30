#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarize matches against maize and bee reference data.

Keep the methods' distinct units visible. Ratios between raw counts are
descriptive comparisons of outputs, not ratios of biological risk.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def sifi_por_gen(ruta_tsv):
    """{gen: nro de transcriptos distintos alcanzados}"""
    if not os.path.exists(ruta_tsv):
        return {}
    d = defaultdict(set)
    with open(ruta_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            d[row["query_name"]].add(row["hit_name"])
    return {g: len(v) for g, v in d.items()}


def dsrna_por_gen(dir_extraidos, archivo):
    d = {}
    if not os.path.isdir(dir_extraidos):
        return d
    for nombre in sorted(os.listdir(dir_extraidos)):
        base = os.path.join(dir_extraidos, nombre, nombre)
        if not os.path.isdir(base):
            base = os.path.join(dir_extraidos, nombre)
        f = os.path.join(base, archivo)
        if not os.path.exists(f):
            continue
        contigs = set()
        with open(f) as fh:
            for linea in fh:
                c = linea.rstrip("\n").split("\t")
                if len(c) >= 5:
                    contigs.add(c[4])
        d[nombre] = len(contigs)
    return d


def siren_por_gen(dir_siren, sufijo):
    d = {}
    if not os.path.isdir(dir_siren):
        return d
    for nombre in sorted(os.listdir(dir_siren)):
        if not nombre.endswith(sufijo):
            continue
        f = os.path.join(dir_siren, nombre, "off_targets_summary.tsv")
        if not os.path.exists(f):
            continue
        with open(f) as fh:
            n = sum(1 for _ in csv.DictReader(fh, delimiter="\t"))
        d[nombre[:-len(sufijo)]] = n
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sifi-out", default=os.path.expanduser("~/rnai_bench/out"))
    ap.add_argument("--siren-out", default=os.path.expanduser("~/rnai_bench/siren_out"))
    ap.add_argument("--dsrna", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    datos = {
        "zmays": {
            "si-Fi": sifi_por_gen(os.path.join(args.sifi_out, "sifi_offtarget_zmays.tsv")),
            "dsRNAeng": dsrna_por_gen(args.dsrna, "zea_mays_l.21.bwt"),
            "SIREN": siren_por_gen(args.siren_out, "_ZMAYS"),
        },
        "amellifera": {
            "si-Fi": sifi_por_gen(os.path.join(args.sifi_out, "sifi_offtarget_amellifera.tsv")),
            "dsRNAeng": dsrna_por_gen(args.dsrna, "apis_mellifera.21.bwt"),
            "SIREN": siren_por_gen(args.siren_out, "_AMELLIFERA"),
        },
    }

    genes = sorted(set().union(*[set(h) for e in datos.values() for h in e.values()]))
    if not genes:
        print("sin datos", file=sys.stderr)
        return 1

    filas = []
    for esp in ("zmays", "amellifera"):
        print()
        print("=" * 88)
        print("TRANSCRIPTOS DISTINTOS ALCANZADOS EN %s"
              % ("Zea mays (maiz)" if esp == "zmays" else "Apis mellifera (abeja)"))
        print("=" * 88)
        print("%-20s %-12s %-12s %-12s %s"
              % ("gen", "si-Fi", "dsRNAeng", "SIREN", "cociente SIREN/bowtie"))
        print("-" * 88)
        for g in genes:
            a = datos[esp]["si-Fi"].get(g)
            b = datos[esp]["dsRNAeng"].get(g)
            c = datos[esp]["SIREN"].get(g)
            bowtie = max([x for x in (a, b) if x is not None], default=None)
            if c and bowtie:
                coc = "%.0f x" % (c / bowtie)
            elif c and bowtie == 0:
                coc = "inf (bowtie = 0)"
            else:
                coc = "-"
            print("%-20s %-12s %-12s %-12s %s"
                  % (g,
                     "-" if a is None else a,
                     "sin datos" if b is None else b,
                     "-" if c is None else c,
                     coc))
            filas.append({"especie": esp, "gen": g,
                          "sifi": "" if a is None else a,
                          "dsrnaengineer": "" if b is None else b,
                          "siren": "" if c is None else c,
                          "cociente_siren_bowtie": coc})

    os.makedirs(args.out_dir, exist_ok=True)
    ruta = os.path.join(args.out_dir, "no_objetivo_tres.tsv")
    with open(ruta, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["especie", "gen", "sifi",
                                           "dsrnaengineer", "siren",
                                           "cociente_siren_bowtie"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(filas)
    print("\nescrito %s" % ruta)
    print("\nLas dos de bowtie exigen match exacto de 21 nt; SIREN acepta >=17 de 21")
    print("apareadas con dG <= -25, admitiendo mismatches y wobble G:U. El cociente")
    print("mide cuanto mas permisivo es el criterio termodinamico entre especies")
    print("lejanas, NO cuanto mas riesgo real hay.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
