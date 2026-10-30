#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare seven ranking criteria on the same si-Fi regional data.

A penalizes external and intra-species matches before activity-filter counts.
B uses activity-filter counts only. C uses total match counts only. D and E
combine those counts with weights 10 and 1. F treats intra-species matches as
neutral; G prioritizes them after minimizing external matches. The main ranking
uses G (control_plaga); the table permits inspection of these assumptions.
"""

import argparse
import csv
import os
import sys


def leer(ruta):
    with open(ruta) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    args = ap.parse_args()

    ruta = os.path.join(args.dir, "sifi_ventanas200.tsv")
    if not os.path.exists(ruta):
        print("falta %s" % ruta, file=sys.stderr)
        return 1

    por_gen = {}
    for r in leer(ruta):
        f = {
            "ini": int(r["inicio"]), "fin": int(r["fin"]),
            "efic": int(r["eficientes"]),
            "intra": int(r["off_intra"]),
            "fuera": int(r["off_zmays"]) + int(r["off_amellifera"]),
        }
        f["total"] = f["intra"] + f["fuera"]
        por_gen.setdefault(r["gen"], []).append(f)

    # A-E asumen que el off-target DENTRO de D. maidis es un costo. Eso vale si
    # el objetivo es caracterizacion funcional (poder atribuir el fenotipo al gen
    # apuntado). Si el objetivo es el control de la plaga, pegarle a otros genes de la
    # plaga no es obviamente malo y hasta puede ayudar a demorar resistencia.
    # F y G prueban ese supuesto: F lo trata como neutro, G como bonus.
    criterios = {
        "A lexicografico": lambda x: (x["fuera"], x["intra"], -x["efic"]),
        "B solo eficacia": lambda x: (-x["efic"],),
        "C solo riesgo": lambda x: (x["total"], x["ini"]),
        "D ponderado x10": lambda x: (-(x["efic"] - 10 * x["total"]),),
        "E ponderado x1": lambda x: (-(x["efic"] - 1 * x["total"]),),
        "F intra neutro": lambda x: (x["fuera"], -x["efic"]),
        "G intra bonus": lambda x: (x["fuera"], -x["intra"], -x["efic"]),
    }

    nombres = list(criterios)
    print("MEJOR VENTANA DE si-Fi SEGUN SIETE CRITERIOS")
    print("=" * 104)
    print("%-19s %s" % ("gen", "  ".join("%-15s" % n for n in nombres)))
    print("-" * 104)

    filas = []
    n_coinciden_mixtos = 0
    n_degenerados = 0
    # A, D y E combinan riesgo y eficacia con pesos distintos. B y C miran una
    # sola dimension, asi que se evaluan aparte.
    mixtos = ["A lexicografico", "D ponderado x10", "E ponderado x1"]

    for g in sorted(por_gen):
        f = por_gen[g]
        mejores = {}
        for n in nombres:
            m = sorted(f, key=criterios[n])[0]
            mejores[n] = "%d-%d" % (m["ini"], m["fin"])

        # Un criterio es DEGENERADO en este gen si todas las ventanas empatan:
        # entonces su "mejor region" no la decide el criterio sino el desempate,
        # y contarla como desacuerdo seria enganoso.
        sin_var_riesgo = len({x["total"] for x in f}) == 1
        sin_var_efic = len({x["efic"] for x in f}) == 1
        if sin_var_riesgo:
            n_degenerados += 1

        coinciden_mixtos = len({mejores[n] for n in mixtos}) == 1
        if coinciden_mixtos:
            n_coinciden_mixtos += 1

        notas = []
        if not coinciden_mixtos:
            notas.append("A/D/E DIFIEREN")
        if sin_var_riesgo:
            notas.append("C degenerado (0 off-target en todas las ventanas)")
        if sin_var_efic:
            notas.append("B degenerado")
        marca = ("   <-- " + "; ".join(notas)) if notas else ""
        print("%-19s %s%s"
              % (g, "  ".join("%-15s" % mejores[n] for n in nombres), marca))
        filas.append(dict(zip(["gen"] + nombres, [g] + [mejores[n] for n in nombres])))

    print()
    print(">>> Criterios que combinan riesgo y eficacia (A, D, E), con pesos muy")
    print(">>> distintos entre si: coinciden en %d de %d genes."
          % (n_coinciden_mixtos, len(por_gen)))
    print(">>> En %d genes el criterio C queda degenerado: si-Fi no encuentra"
          % n_degenerados)
    print(">>> ningun off-target en ninguna ventana, asi que todas empatan y gana")
    print(">>> la primera por desempate. Eso NO es un criterio alternativo.")
    print()
    print(">>> LECTURA: la recomendacion es robusta a COMO se pesan riesgo y")
    print(">>> eficacia, pero no a IGNORAR una de las dos. Donde B (solo eficacia)")
    print(">>> difiere, la region cambia de verdad y hay que declarar el criterio.")

    # El supuesto de fondo: intra-especie como costo (A) vs neutro (F) vs bonus (G)
    print()
    print("EL SUPUESTO SOBRE EL OFF-TARGET DENTRO DE D. maidis")
    print("-" * 104)
    print("A lo trata como COSTO (objetivo: caracterizacion funcional).")
    print("F lo trata como NEUTRO y G como BONUS (objetivo: control de la plaga,")
    print("donde pegarle a otros genes de la plaga puede sumar).")
    print()
    n_igual = 0
    for g in sorted(por_gen):
        f = por_gen[g]
        a = sorted(f, key=criterios["A lexicografico"])[0]
        fr = sorted(f, key=criterios["F intra neutro"])[0]
        gg = sorted(f, key=criterios["G intra bonus"])[0]
        ra = "%d-%d" % (a["ini"], a["fin"])
        rf = "%d-%d" % (fr["ini"], fr["fin"])
        rg = "%d-%d" % (gg["ini"], gg["fin"])
        sin_intra = len({x["intra"] for x in f}) == 1
        if ra == rf == rg:
            n_igual += 1
            nota = "  (el supuesto no cambia nada)"
        else:
            nota = "  <-- EL SUPUESTO CAMBIA LA RESPUESTA"
        if sin_intra:
            nota += "  [sin off-target intra en ninguna ventana]"
        print("  %-19s A=%-13s F=%-13s G=%-13s%s" % (g, ra, rf, rg, nota))
    print()
    print(">>> El supuesto sobre el intra-especie NO cambia la region en %d de %d genes."
          % (n_igual, len(por_gen)))

    salida = os.path.join(args.dir, "sensibilidad_criterio.tsv")
    with open(salida, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gen"] + nombres, delimiter="\t")
        w.writeheader()
        w.writerows(filas)
    print("\nescrito %s" % salida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
