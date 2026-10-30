#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare si-Fi21, SIREN and dsRNAEngineer on a shared regional grid.

Use the criteria defined in criterio.py for the si-Fi and dsRNAEngineer data,
and SIREN's regional score for its ordering. Report missing or non-discriminating
dsRNAEngineer results explicitly. References and count units differ by method.
"""

import argparse
import csv
import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from comparar_herramientas import spearman, jaccard  # noqa: E402
from criterio import clave_sifi, clave_dsrna, DESCRIPCION  # noqa: E402


def leer(ruta):
    if not os.path.exists(ruta):
        return []
    with open(ruta) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def clave(r):
    return (r["gen"], int(r["inicio"]), int(r["fin"]))


def top(filas, f_val, mayor_mejor, frac=0.10):
    if not filas:
        return set()
    k = max(1, int(round(frac * len(filas))))
    orden = sorted(filas, key=lambda r: f_val(r), reverse=mayor_mejor)
    return {(r["inicio"], r["fin"]) for r in orden[:k]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    args = ap.parse_args()
    d = args.dir

    sifi = {clave(r): r for r in leer(os.path.join(d, "sifi_ventanas200.tsv"))}
    siren = {clave(r): r for r in leer(os.path.join(d, "siren_regiones200.tsv"))}
    dsr = {clave(r): r for r in leer(os.path.join(d, "dsrnaengineer_ventanas200.tsv"))}

    faltan = [n for n, m in (("si-Fi", sifi), ("SIREN", siren),
                             ("dsRNAEngineer", dsr)) if not m]
    if faltan:
        print("faltan resultados de: %s" % ", ".join(faltan), file=sys.stderr)
        return 1

    comunes = sorted(set(sifi) & set(siren) & set(dsr))
    print("ventanas con las tres herramientas: %d" % len(comunes))

    por_gen = {}
    for k in comunes:
        g = k[0]
        a, b, c = sifi[k], siren[k], dsr[k]
        por_gen.setdefault(g, []).append({
            "inicio": k[1], "fin": k[2],
            # Se guardan 'fuera' e 'intra' POR SEPARADO, no sumados: el criterio
            # unificado (criterio.py) los trata como niveles distintos, porque
            # pegarle al maiz o a la abeja pesa mas que pegarle a un paralogo
            # del propio insecto.
            "sifi_fuera": int(a["off_zmays"]) + int(a["off_amellifera"]),
            "sifi_intra": int(a["off_intra"]),
            "sifi_off": int(a["off_intra"]) + int(a["off_zmays"]) + int(a["off_amellifera"]),
            "sifi_efic": int(a["eficientes"]),
            "siren_off": int(b["n_off_sirnas"]),
            # SIREN es la UNICA de las tres con ranking propio. Para ordenar sus
            # ventanas se usa SU score (formula exacta, validada con dif_max=0)
            # y no un proxy nuestro: su Score pesa la concentracion por
            # transcripto (-30 por siRNA excedente) y contar siRNAs sueltos da
            # otra cosa. Difieren en 2 de 9 genes, y en DM_RPN6 bastante
            # (109-308 por su score contra 1513-1712 por el proxy).
            "siren_score": float(b["score_corregido"]),
            "ds_fuera": int(c["ds_off_zmays"]) + int(c["ds_off_amellifera"]),
            "ds_intra": int(c["ds_off_intra"]),
            "ds_off": int(c["ds_off_intra"]) + int(c["ds_off_zmays"]) + int(c["ds_off_amellifera"]),
            "ds_on": int(c["ds_on_target"]),
        })

    filas_out = []
    print()
    print("CORRELACION DE RANGOS DEL RIESGO OFF-TARGET ENTRE PARES DE HERRAMIENTAS")
    print("(los dos primeros son bowtie; el tercero compara bowtie contra termodinamica)")
    print("-" * 96)
    print("%-19s %-6s %-14s %-14s %-14s %s"
          % ("gen", "n", "siFi~dsRNA", "siFi~SIREN", "dsRNA~SIREN", "veredicto"))
    print("-" * 96)
    for g in sorted(por_gen):
        f = por_gen[g]
        r_bb = spearman([x["sifi_off"] for x in f], [x["ds_off"] for x in f])
        r_bs = spearman([x["sifi_off"] for x in f], [x["siren_off"] for x in f])
        r_ds = spearman([x["ds_off"] for x in f], [x["siren_off"] for x in f])

        def fmt(v):
            return "%.3f" % v if v is not None else "NA"

        if r_bb is not None and r_bs is not None and r_ds is not None:
            otros = max(r_bs, r_ds)
            veredicto = ("bowtie mas parecido entre si" if r_bb > otros
                         else "NO se cumple")
        else:
            veredicto = "sin varianza suficiente"
        print("%-19s %-6d %-14s %-14s %-14s %s"
              % (g, len(f), fmt(r_bb), fmt(r_bs), fmt(r_ds), veredicto))
        filas_out.append({
            "gen": g, "n_ventanas": len(f),
            "rho_sifi_dsrna": fmt(r_bb),
            "rho_sifi_siren": fmt(r_bs),
            "rho_dsrna_siren": fmt(r_ds),
            "veredicto": veredicto,
        })

    print()
    print("REGION RECOMENDADA POR CADA HERRAMIENTA Y SOLAPAMIENTO DEL TOP 10%")
    print("criterio para si-Fi y dsRNAEngineer (nuestro, no de ellas): %s" % DESCRIPCION)
    print("-" * 96)
    print("%-19s %-12s %-12s %-12s %-7s %-7s %s"
          % ("gen", "si-Fi", "SIREN", "dsRNAeng", "J:sF-dS", "J:sF-SI", "J:dS-SI"))
    print("-" * 96)
    for g in sorted(por_gen):
        f = por_gen[g]
        # criterio unificado, definido en criterio.py
        o_sifi = sorted(f, key=lambda x: clave_sifi(x["sifi_fuera"], x["sifi_intra"],
                                                    x["sifi_efic"]))
        o_ds = sorted(f, key=lambda x: clave_dsrna(x["ds_fuera"], x["ds_intra"],
                                                   x["ds_on"]))

        # dsRNAEngineer no devolvio NINGUN hit para este gen (sus tres .bwt
        # salieron vacios, como pasa con los dos genes de MRFV: secuencia viral
        # sin homologia con transcriptomas de insecto o planta). Entonces todas
        # las ventanas empatan en cero y la "mejor" seria solo la primera por
        # desempate. Reportar esa coordenada como si fuera una recomendacion
        # engaña a quien lea el TSV sin el informe al lado.
        ds_sin_datos = all(x["ds_on"] == 0 and x["ds_off"] == 0 for x in f)
        # Caso distinto: hay datos pero la clave de orden es constante, asi que
        # el primer puesto tampoco lo decide el criterio.
        claves_ds = {clave_dsrna(x["ds_fuera"], x["ds_intra"], x["ds_on"]) for x in f}
        ds_empate = (not ds_sin_datos) and len(claves_ds) == 1
        # SIREN con SU propio ranking: score mas alto es mejor.
        o_siren = sorted(f, key=lambda x: -x["siren_score"])
        k = max(1, len(f) // 10)
        t_sifi = {(x["inicio"], x["fin"]) for x in o_sifi[:k]}
        t_siren = {(x["inicio"], x["fin"]) for x in o_siren[:k]}  # por su score
        t_ds = {(x["inicio"], x["fin"]) for x in o_ds[:k]}

        def j(a, b, valido=True):
            if not valido:
                return "sin datos"
            v = jaccard(a, b)
            return "%.2f" % v if v is not None else "NA"

        if ds_sin_datos:
            mejor_ds = "sin datos"
        elif ds_empate:
            mejor_ds = "empate total"
        else:
            mejor_ds = "%d-%d" % (o_ds[0]["inicio"], o_ds[0]["fin"])
        ds_ok = not (ds_sin_datos or ds_empate)

        print("%-19s %-12s %-12s %-12s %-9s %-9s %s"
              % (g,
                 "%d-%d" % (o_sifi[0]["inicio"], o_sifi[0]["fin"]),
                 "%d-%d" % (o_siren[0]["inicio"], o_siren[0]["fin"]),
                 mejor_ds,
                 j(t_sifi, t_ds, ds_ok), j(t_sifi, t_siren),
                 j(t_ds, t_siren, ds_ok)))
        filas_out[[x["gen"] for x in filas_out].index(g)].update({
            "mejor_sifi": "%d-%d" % (o_sifi[0]["inicio"], o_sifi[0]["fin"]),
            "mejor_siren": "%d-%d" % (o_siren[0]["inicio"], o_siren[0]["fin"]),
            "mejor_dsrna": mejor_ds,
            "jaccard_sifi_dsrna": j(t_sifi, t_ds, ds_ok),
            "jaccard_sifi_siren": j(t_sifi, t_siren),
            "jaccard_dsrna_siren": j(t_ds, t_siren, ds_ok),
        })

    campos = ["gen", "n_ventanas", "rho_sifi_dsrna", "rho_sifi_siren",
              "rho_dsrna_siren", "veredicto", "mejor_sifi", "mejor_siren",
              "mejor_dsrna", "jaccard_sifi_dsrna", "jaccard_sifi_siren",
              "jaccard_dsrna_siren"]
    ruta = os.path.join(d, "comparacion_tres.tsv")
    with open(ruta, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader()
        for r in filas_out:
            w.writerow({c: r.get(c, "") for c in campos})
    print("\nescrito %s" % ruta)

    # tabla completa por ventana
    ruta_v = os.path.join(d, "comparacion_tres_ventanas.tsv")
    campos_v = ["gen", "inicio", "fin", "sifi_off", "sifi_efic",
                "siren_off", "ds_off", "ds_on"]
    with open(ruta_v, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos_v, delimiter="\t")
        w.writeheader()
        for g in sorted(por_gen):
            for x in por_gen[g]:
                w.writerow({"gen": g, **{c: x[c] for c in campos_v if c != "gen"}})
    print("escrito %s" % ruta_v)
    return 0


if __name__ == "__main__":
    sys.exit(main())
