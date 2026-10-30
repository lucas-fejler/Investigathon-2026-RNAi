#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculate count percentiles for each prioritized region across analyses.

A percentile is the proportion of windows within the same sequence with a
strictly lower count. Percentiles compare relative positions within a profile;
they do not convert counts into efficacy or biological risk probabilities.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def percentil(valores, v):
    """Porcentaje de ventanas con valor ESTRICTAMENTE menor, 0-100."""
    n = len(valores)
    if n <= 1:
        return 0.0
    menores = sum(1 for x in valores if x < v)
    return 100.0 * menores / (n - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    args = ap.parse_args()

    ruta = os.path.join(args.dir, "comparacion_tres_ventanas.tsv")
    if not os.path.exists(ruta):
        print("falta %s" % ruta, file=sys.stderr)
        return 1

    por_gen = defaultdict(list)
    with open(ruta) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            por_gen[r["gen"]].append({
                "ini": int(r["inicio"]), "fin": int(r["fin"]),
                "sifi": int(r["sifi_off"]),
                "siren": int(r["siren_off"]),
                "ds": int(r["ds_off"]),
                "ds_total": int(r["ds_off"]) + int(r["ds_on"]),
            })

    # mejor ventana por herramienta, tal como la reporta comparacion_tres.tsv
    mejores = {}
    rc = os.path.join(args.dir, "comparacion_tres.tsv")
    with open(rc) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mejores[r["gen"]] = {
                "si-Fi": r["mejor_sifi"],
                "SIREN": r["mejor_siren"],
                "dsRNAeng": r["mejor_dsrna"],
            }

    filas = []
    print("QUE DICEN LAS OTRAS SOBRE LA REGION QUE CADA UNA RECOMIENDA")
    print("percentil de off-target dentro del gen: 0 = la mas limpia, 100 = la peor")
    print("=" * 96)
    print("%-19s %-11s %-13s %-11s %-11s %s"
          % ("gen", "recomienda", "region", "si-Fi", "SIREN", "dsRNAeng"))
    print("-" * 96)

    conflictos = []
    for g in sorted(por_gen):
        f = por_gen[g]
        v_sifi = [x["sifi"] for x in f]
        v_siren = [x["siren"] for x in f]
        v_ds = [x["ds"] for x in f]
        sin_ds = len(set(v_ds)) == 1 and sum(v_ds) == 0 and sum(x["ds_total"] for x in f) == 0

        for herr, reg in mejores.get(g, {}).items():
            try:
                a, b = (int(x) for x in reg.split("-"))
            except ValueError:
                continue
            w = next((x for x in f if x["ini"] == a and x["fin"] == b), None)
            if w is None:
                continue
            p_sifi = percentil(v_sifi, w["sifi"])
            p_siren = percentil(v_siren, w["siren"])
            p_ds = None if sin_ds else percentil(v_ds, w["ds"])

            def fmt(p, propia):
                if p is None:
                    return "sin datos"
                marca = " *" if propia else ""
                return "p%3.0f%s" % (p, marca)

            print("%-19s %-11s %-13s %-11s %-11s %s"
                  % (g, herr, reg,
                     fmt(p_sifi, herr == "si-Fi"),
                     fmt(p_siren, herr == "SIREN"),
                     fmt(p_ds, herr == "dsRNAeng")))

            for otra, p in (("si-Fi", p_sifi), ("SIREN", p_siren), ("dsRNAeng", p_ds)):
                if otra == herr or p is None:
                    continue
                if p >= 75:
                    conflictos.append((g, herr, reg, otra, p))
            filas.append({
                "gen": g, "recomienda": herr, "region": reg,
                "percentil_sifi": "" if p_sifi is None else round(p_sifi),
                "percentil_siren": "" if p_siren is None else round(p_siren),
                "percentil_dsrna": "" if p_ds is None else round(p_ds),
            })
        print("-" * 96)

    print("\n* = percentil de la propia herramienta en su propia recomendacion.")
    print("  Un valor ALTO ahi NO es un error: en modo control_plaga el criterio")
    print("  PREMIA los hits dentro de D. maidis, asi que elige ventanas con")
    print("  muchos, y en una escala cruda de 'off-target' eso se lee como pesimo.")
    print("  La columna mezcla hits intra (deseados) con los de maiz y abeja")
    print("  (indeseados): por eso conviene mirar los conflictos de abajo, que son")
    print("  entre herramientas distintas, y no la celda propia.")

    print("\nCONFLICTOS: region recomendada que OTRA marca en el cuartil peor (p>=75)")
    print("-" * 96)
    if not conflictos:
        print("  ninguno")
    for g, herr, reg, otra, p in conflictos:
        print("  %-19s %s recomienda %-12s pero %-9s la pone en p%.0f"
              % (g, herr, reg, otra, p))

    salida = os.path.join(args.dir, "cruce_recomendaciones.tsv")
    with open(salida, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gen", "recomienda", "region",
                                           "percentil_sifi", "percentil_siren",
                                           "percentil_dsrna"], delimiter="\t")
        w.writeheader()
        w.writerows(filas)
    print("\nescrito %s" % salida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
