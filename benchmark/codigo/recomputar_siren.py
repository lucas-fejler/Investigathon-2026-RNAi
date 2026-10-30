#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculate SIREN scores on the complete common 200-nt grid.

Exclude the operationally identified target transcript from the interaction
penalty. For each region, apply -0.1 per distinct siRNA name and -30 for each
entry beyond the first within a transcript. Keep direct and reverse names
separate. With --validar, compare the score including the target against the
native regional output to check the implementation of the formula.

Requires per-transcript interaction files and the target-identification data.
The score is a computational interaction penalty, not an efficacy measurement.
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict

PASO = 4
LARGO_REGION = 200
PENAL_UNICO = 0.1
PENAL_EXTRA = 30.0


def leer_largo_fasta(ruta):
    """Largo del primer registro del FASTA (target.fa tiene uno solo)."""
    n = 0
    visto = False
    with open(ruta) as fh:
        for linea in fh:
            linea = linea.strip()
            if linea.startswith(">"):
                if visto:
                    break
                visto = True
            elif linea:
                n += len(linea)
    return n


def parse_nombre_sirna(nombre):
    """'sirna_357-377_r' -> (357, 377). Igual que parse_sirna_name() de SIREN."""
    limpio = nombre.replace("sirna_", "").replace("_r", "")
    partes = limpio.split("-")
    return int(partes[0]), int(partes[1])


def leer_off_targets(ruta):
    """
    Devuelve [(transcripto, [(nombre, ini, fin), ...]), ...].

    Se guarda el NOMBRE y no solo las coordenadas, y esto es esencial:
    parse_sirna_name() borra el sufijo "_r", asi que 'sirna_1-21' y
    'sirna_1-21_r' (las hebras directa y reversa) dan las MISMAS coordenadas
    pero son dos entradas distintas para SIREN. Deduplicar por coordenada
    descuenta la mitad de los siRNAs y el score no cierra.
    """
    filas = []
    with open(ruta) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            nombres = [s.strip() for s in row["siRNA names"].split(",") if s.strip()]
            items = []
            for n in nombres:
                try:
                    ini, fin = parse_nombre_sirna(n)
                except (ValueError, IndexError):
                    continue
                items.append((n, ini, fin))
            filas.append((row["Off target"], items))
    return filas


def regiones(largo):
    """Como generate_rnai_sequences(): paso 4, base 1 inclusive."""
    out = []
    for inicio0 in range(0, largo - LARGO_REGION + 1, PASO):
        out.append((inicio0 + 1, inicio0 + LARGO_REGION))
    return out


def puntuar(region, off_targets, excluir=None):
    """
    Score de SIREN para una region, replicando score_chunk() de
    siren_designVIII.py. 'excluir' es el transcripto on-target a sacar del
    calculo (None = incluir todo, que es el modo de validacion).

    Dos detalles del original que hay que respetar o el numero no cierra:

    - penalized_sirnas se inicializa UNA vez, fuera del bucle de transcriptos,
      asi que el -0.1 se cobra una sola vez por nombre de siRNA aunque ese
      siRNA pegue en varios transcriptos distintos.
    - off_target_counts, en cambio, se reinicia en CADA transcripto, asi que el
      castigo de -30 por siRNA excedente se aplica por transcripto.
    """
    ini, fin = region
    score = 0.0
    nombres_penalizados = set()
    for transcripto, items in off_targets:
        if excluir is not None and transcripto == excluir:
            continue
        n = 0
        for (nombre, s_ini, s_fin) in items:
            if s_ini >= ini and s_fin <= fin:
                if nombre not in nombres_penalizados:
                    score -= PENAL_UNICO
                    nombres_penalizados.add(nombre)
                n += 1
        if n > 1:
            score += (n - 1) * -PENAL_EXTRA
    return score


def on_target_de_sifi(ruta_json):
    """{gen: transcripto_blanco} desde el summary de si-Fi."""
    mapa = {}
    if not os.path.exists(ruta_json):
        return mapa
    with open(ruta_json) as fh:
        d = json.load(fh)
    for s in d.get("summary", []):
        mt = s.get("main_targets")
        if mt:
            mapa[s["query"]] = mt[0]
    return mapa


def validar(dir_gen, largo, off_targets):
    """
    Compara el score recalculado (INCLUYENDO on-target) contra el publicado por
    SIREN, sobre las regiones de 200 nt que SIREN reporto.
    """
    ruta = os.path.join(dir_gen, "rna_sequences_with_scores.tsv")
    if not os.path.exists(ruta):
        return None
    publicados = {}
    with open(ruta) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if int(row["Length"]) != LARGO_REGION:
                continue
            coords = row["Name"].split("_")[-2]
            ini, fin = (int(x) for x in coords.split("-"))
            publicados[(ini, fin)] = float(row["Score"])
    difs = []
    for region, esperado in publicados.items():
        obtenido = puntuar(region, off_targets, excluir=None)
        difs.append(abs(obtenido - esperado))
    if not difs:
        return None
    return len(difs), max(difs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--siren-out", default=os.path.expanduser("~/rnai_bench/siren_out"))
    ap.add_argument("--sifi-json", default=os.path.expanduser("~/rnai_bench/out/sifi_design_dmaidis.json"))
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--validar", action="store_true")
    ap.add_argument("--sin-on-target", action="store_true",
                    help="No excluir ningun transcripto. Es lo correcto para "
                         "las corridas contra especies NO objetivo (maiz, "
                         "abeja): ahi no existe ortologo del gen, asi que TODO "
                         "hit es off-target. Excluir el de mas hits, como se "
                         "hace contra D. maidis, descartaria justamente el "
                         "riesgo mas alto.")
    ap.add_argument("--sufijo", default="_PATCH",
                    help="solo los directorios que terminan asi. Hace falta "
                         "porque siren_out guarda tambien las corridas de "
                         "prueba (_H4/_H6/_H8) y mezclarlas duplicaria genes.")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    mapa_on = on_target_de_sifi(args.sifi_json)

    campos = ["gen", "inicio", "fin", "score_corregido", "score_con_on_target",
              "n_off_transcriptos", "n_off_sirnas", "on_target_excluido"]
    filas_out = []
    resumen = []

    for nombre in sorted(os.listdir(args.siren_out)):
        if args.sufijo and not nombre.endswith(args.sufijo):
            continue
        dir_gen = os.path.join(args.siren_out, nombre)
        ruta_sum = os.path.join(dir_gen, "off_targets_summary.tsv")
        ruta_tgt = os.path.join(dir_gen, "other_files", "target.fa")
        if not (os.path.exists(ruta_sum) and os.path.exists(ruta_tgt)):
            continue
        gen = nombre[:-len(args.sufijo)] if args.sufijo else nombre

        largo = leer_largo_fasta(ruta_tgt)
        off_targets = leer_off_targets(ruta_sum)

        # on-target: el que dijo si-Fi; si no esta, el de mas siRNAs
        # On-target: primero el que identifico si-Fi. Si no aparece en la base
        # de SIREN, se cae al de mas siRNAs, desempatando por NOMBRE para que
        # la eleccion sea deterministica: el orden de filas de
        # off_targets_summary.tsv depende del orden en que terminaron los
        # shards de RNAhybrid, que no es reproducible, y un max() a secas
        # elegiria distinto en cada corrida.
        on = None if args.sin_on_target else mapa_on.get(gen)
        if args.sin_on_target:
            pass
        elif on is None or not any(t == on for t, _ in off_targets):
            if off_targets:
                on = sorted(off_targets, key=lambda x: (-len(x[1]), x[0]))[0][0]
            else:
                on = None

        if args.validar:
            v = validar(dir_gen, largo, off_targets)
            if v:
                n, maxdif = v
                estado = "OK" if maxdif < 1e-6 else "DISCREPA"
                print("  validacion %-20s %s  regiones=%d  dif_max=%.9f"
                      % (gen, estado, n, maxdif))
            else:
                print("  validacion %-20s sin datos publicados" % gen)

        regs = regiones(largo)
        mejor = None
        for r in regs:
            sc_corr = puntuar(r, off_targets, excluir=on)
            sc_con = puntuar(r, off_targets, excluir=None)
            n_tr = 0
            n_si = 0
            for transcripto, items in off_targets:
                if transcripto == on:
                    continue
                dentro = [x for x in items if x[1] >= r[0] and x[2] <= r[1]]
                if dentro:
                    n_tr += 1
                    n_si += len(dentro)
            fila = {
                "gen": gen, "inicio": r[0], "fin": r[1],
                "score_corregido": round(sc_corr, 3),
                "score_con_on_target": round(sc_con, 3),
                "n_off_transcriptos": n_tr, "n_off_sirnas": n_si,
                "on_target_excluido": on or "",
            }
            filas_out.append(fila)
            if mejor is None or sc_corr > mejor["score_corregido"]:
                mejor = fila
        if mejor:
            resumen.append(mejor)

    if not filas_out:
        print("No habia corridas de SIREN completas todavia.", file=sys.stderr)
        return

    ruta_t = os.path.join(args.out_dir, "siren_regiones200.tsv")
    with open(ruta_t, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader()
        w.writerows(filas_out)
    print("escrito %s (%d regiones, %d genes)"
          % (ruta_t, len(filas_out), len(resumen)))

    ruta_m = os.path.join(args.out_dir, "siren_mejor_region.tsv")
    with open(ruta_m, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t")
        w.writeheader()
        w.writerows(resumen)
    print("escrito %s" % ruta_m)


if __name__ == "__main__":
    main()
