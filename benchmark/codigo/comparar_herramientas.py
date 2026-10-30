#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare si-Fi21 and SIREN on the same 200-nt coordinate grid.

Join by sequence, start and end. Export regional metrics, Spearman correlations
and top-decile overlap. si-Fi position counts and SIREN interaction counts have
different units; ranking agreement is a descriptive result, not efficacy.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from criterio import clave_sifi  # noqa: E402


def leer_tsv(ruta):
    if not os.path.exists(ruta):
        return []
    with open(ruta) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def rangos(valores):
    """
    Rangos con empates promediados (midrank), como pide Spearman.
    Devuelve una lista de rangos en el orden original.
    """
    n = len(valores)
    orden = sorted(range(n), key=lambda i: valores[i])
    rg = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and valores[orden[j + 1]] == valores[orden[i]]:
            j += 1
        promedio = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            rg[orden[k]] = promedio
        i = j + 1
    return rg


def spearman(xs, ys):
    """Correlacion de rangos. None si no hay varianza o faltan datos."""
    if len(xs) != len(ys) or len(xs) < 3:
        return None
    rx, ry = rangos(xs), rangos(ys)
    n = len(rx)
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = sum((rx[i] - mx) ** 2 for i in range(n))
    dy = sum((ry[i] - my) ** 2 for i in range(n))
    if dx <= 0 or dy <= 0:
        return None
    return num / (dx * dy) ** 0.5


def top_conjunto(filas, clave, mayor_es_mejor, frac=0.10):
    """Las mejores ceil(frac*n) ventanas segun 'clave', como set de (ini,fin)."""
    if not filas:
        return set()
    k = max(1, int(round(frac * len(filas))))
    ordenadas = sorted(filas, key=lambda f: f[clave], reverse=mayor_es_mejor)
    return {(f["inicio"], f["fin"]) for f in ordenadas[:k]}


def jaccard(a, b):
    if not a and not b:
        return None
    u = len(a | b)
    return len(a & b) / u if u else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="carpeta con los TSV de las dos herramientas")
    args = ap.parse_args()

    sifi = leer_tsv(os.path.join(args.dir, "sifi_ventanas200.tsv"))
    siren = leer_tsv(os.path.join(args.dir, "siren_regiones200.tsv"))

    if not sifi:
        print("falta sifi_ventanas200.tsv", file=sys.stderr)
        return 1
    if not siren:
        print("falta siren_regiones200.tsv (SIREN todavia no termino ningun gen)",
              file=sys.stderr)
        return 1

    idx_siren = {}
    for r in siren:
        idx_siren[(r["gen"], int(r["inicio"]), int(r["fin"]))] = r

    unidas = defaultdict(list)
    solo_sifi = defaultdict(int)
    for r in sifi:
        clave = (r["gen"], int(r["inicio"]), int(r["fin"]))
        s = idx_siren.get(clave)
        if s is None:
            solo_sifi[r["gen"]] += 1
            continue
        unidas[r["gen"]].append({
            "gen": r["gen"],
            "inicio": int(r["inicio"]),
            "fin": int(r["fin"]),
            "sifi_eficientes": int(r["eficientes"]),
            "sifi_off_intra": int(r["off_intra"]),
            "sifi_off_zmays": int(r["off_zmays"]),
            "sifi_off_amellifera": int(r["off_amellifera"]),
            "siren_score": float(s["score_corregido"]),
            "siren_off_transcriptos": int(s["n_off_transcriptos"]),
            "siren_off_sirnas": int(s["n_off_sirnas"]),
        })

    campos = ["gen", "inicio", "fin", "sifi_eficientes", "sifi_off_intra",
              "sifi_off_zmays", "sifi_off_amellifera",
              "siren_score", "siren_off_transcriptos", "siren_off_sirnas"]

    ruta_v = os.path.join(args.dir, "comparacion_ventanas.tsv")
    with open(ruta_v, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader()
        for g in sorted(unidas):
            w.writerows(unidas[g])
    total = sum(len(v) for v in unidas.values())
    print("escrito %s (%d ventanas en comun, %d genes)"
          % (ruta_v, total, len(unidas)))

    campos_r = ["gen", "n_ventanas",
                "rho_offtarget", "rho_eficacia_vs_score",
                "jaccard_top10", "sifi_mejor", "siren_mejor", "consenso",
                "sifi_off_intra_total", "siren_off_sirnas_mediana"]
    resumen = []
    for g in sorted(unidas):
        filas = unidas[g]
        # off-target: si-Fi exacto vs SIREN termodinamico (ambos "menos es mejor")
        rho_off = spearman([f["sifi_off_intra"] for f in filas],
                           [f["siren_off_sirnas"] for f in filas])
        # eficacia de si-Fi contra el score de SIREN: se espera que NO correlacione
        rho_ef = spearman([f["sifi_eficientes"] for f in filas],
                          [f["siren_score"] for f in filas])

        # criterio unificado del benchmark (criterio.py)
        sifi_ord = sorted(filas, key=lambda f: clave_sifi(
            f["sifi_off_zmays"] + f["sifi_off_amellifera"],
            f["sifi_off_intra"],
            f["sifi_eficientes"]))
        siren_ord = sorted(filas, key=lambda f: -f["siren_score"])

        top_sifi = {(f["inicio"], f["fin"]) for f in sifi_ord[:max(1, len(filas) // 10)]}
        top_siren = top_conjunto(filas, "siren_score", True, 0.10)
        inter = top_sifi & top_siren
        consenso = ""
        if inter:
            # dentro del consenso, la que mas eficientes tenga segun si-Fi
            mejor = max(inter, key=lambda c: next(
                f["sifi_eficientes"] for f in filas
                if (f["inicio"], f["fin"]) == c))
            consenso = "%d-%d" % mejor

        resumen.append({
            "gen": g,
            "n_ventanas": len(filas),
            "rho_offtarget": "%.3f" % rho_off if rho_off is not None else "NA",
            "rho_eficacia_vs_score": "%.3f" % rho_ef if rho_ef is not None else "NA",
            "jaccard_top10": "%.3f" % jaccard(top_sifi, top_siren)
                             if jaccard(top_sifi, top_siren) is not None else "NA",
            "sifi_mejor": "%d-%d" % (sifi_ord[0]["inicio"], sifi_ord[0]["fin"]),
            "siren_mejor": "%d-%d" % (siren_ord[0]["inicio"], siren_ord[0]["fin"]),
            "consenso": consenso or "(sin solapamiento en el top 10%)",
            "sifi_off_intra_total": sum(f["sifi_off_intra"] for f in filas),
            "siren_off_sirnas_mediana": sorted(
                f["siren_off_sirnas"] for f in filas)[len(filas) // 2],
        })

    ruta_r = os.path.join(args.dir, "comparacion_resumen.tsv")
    with open(ruta_r, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos_r, delimiter="\t")
        w.writeheader()
        w.writerows(resumen)
    print("escrito %s" % ruta_r)

    faltan = sorted(set(r["gen"] for r in sifi) - set(unidas))
    if faltan:
        print("\nGenes que si-Fi tiene y SIREN todavia no: %s" % ", ".join(faltan))

    # Solo es desalineacion de rejilla si el gen SI esta en SIREN y aun asi
    # quedaron ventanas sin par. Para un gen que todavia no corrio, tener todas
    # las ventanas sin par es lo esperado y no indica ningun problema.
    desalineadas = {g: n for g, n in solo_sifi.items() if g in unidas}
    if desalineadas:
        print("AVISO - ventanas sin par en genes que SIREN si corrio "
              "(rejillas desalineadas): %s" % desalineadas)

    print()
    for r in resumen:
        print("%-20s n=%-4d rho_off=%-7s rho_efic=%-7s J10=%-6s consenso=%s"
              % (r["gen"], r["n_ventanas"], r["rho_offtarget"],
                 r["rho_eficacia_vs_score"], r["jaccard_top10"], r["consenso"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
