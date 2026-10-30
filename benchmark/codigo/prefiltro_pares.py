#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Select query-reference pairs before RNAhybrid using an exact seed filter.

This experimental optimization reduces candidate pairs under the implemented
seed rule. Equivalence to an unfiltered search must be validated for the
queries and evaluator settings in use; observed retention on a subset does
not establish universal completeness.
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

LARGO_SIRNA = 21
SEED_INI = 2      # 1-based, inclusive
SEED_FIN = 17     # 1-based, inclusive
SEED_LARGO = SEED_FIN - SEED_INI + 1      # 16

# Opciones de RNAhybrid: IDENTICAS a los defaults de run_rnahybrid_capture()
# en sirenXII.py. Si se cambian, la salida deja de ser comparable.
OPTS_RNAHYBRID = ["-e", "-25", "-v", "0", "-u", "0", "-f", "2,7",
                  "-p", "0.01", "-d", "0.5,0.1", "-m", "60000"]

COD = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}
# base del target que puede aparear con cada base del siRNA
SOCIOS = {0: (3,), 1: (2,), 2: (1, 3), 3: (0, 2)}

MAX_VARIANTES = 200000   # tope por siRNA; arriba de eso no se filtra ese siRNA


def leer_fasta(ruta):
    """[(id, secuencia_mayuscula)] preservando el orden del archivo."""
    out = []
    nombre = None
    partes = []
    with open(ruta) as fh:
        for linea in fh:
            linea = linea.strip()
            if linea.startswith(">"):
                if nombre is not None:
                    out.append((nombre, "".join(partes).upper()))
                nombre = linea[1:].split()[0]
                partes = []
            elif linea:
                partes.append(linea)
    if nombre is not None:
        out.append((nombre, "".join(partes).upper()))
    return out


def generar_sirnas(seq):
    """
    Replica generate_sirnas() de sirenXII.py con sensibilidad 'high' (paso 1):
    por cada ventana de 21 nt emite el siRNA directo y su reverso complementario
    transcripto, con los nombres 'sirna_{ini}-{fin}' y '..._r'.
    El ORDEN importa: define el orden de los nombres en off_targets_summary.tsv.
    """
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "N": "N"}
    out = []
    for i in range(0, len(seq) - LARGO_SIRNA + 1):
        s = seq[i:i + LARGO_SIRNA]
        rc = "".join(comp.get(b, "N") for b in reversed(s)).replace("T", "U")
        nombre = "sirna_%d-%d" % (i + 1, i + LARGO_SIRNA)
        out.append((nombre, s))
        out.append((nombre + "_r", rc))
    return out


def codigos_aceptables(sirna):
    """
    Todos los 16-meros del target (como entero en base 4) que pueden aparear con
    las posiciones SEED_INI..SEED_FIN del siRNA.

    El duplex es antiparalelo: el siRNA en posicion i (1-based) aparea con la
    base del target que, leida 5'->3', cae en el offset m = SEED_FIN - i.
    Asi, el offset m del seed del target corresponde al siRNA en SEED_FIN - m.
    Devuelve None si hay demasiadas variantes o si aparece una N.
    """
    opciones = []
    for m in range(SEED_LARGO):
        pos_sirna = SEED_FIN - m          # 1-based
        base = sirna[pos_sirna - 1]
        c = COD.get(base)
        if c is None:
            return None
        opciones.append(SOCIOS[c])

    total = 1
    for o in opciones:
        total *= len(o)
    if total > MAX_VARIANTES:
        return None

    codigos = np.zeros(1, dtype=np.int64)
    for o in opciones:
        codigos = np.concatenate([(codigos << 2) | b for b in o])
    return codigos


def codificar(seq):
    """Secuencia -> array int8 con A/C/G/T=0..3 y 4 para cualquier otra cosa."""
    arr = np.full(len(seq), 4, dtype=np.int8)
    b = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    for ch, v in (("A", 0), ("C", 1), ("G", 2), ("T", 3), ("U", 3)):
        arr[b == ord(ch)] = v
    return arr


def codigos_de_target(arr):
    """
    Codigos de los 16-meros validos de un target.
    Devuelve (codigos, posiciones) donde posiciones es el inicio 0-based.
    Las ventanas que contienen una base no ACGT quedan excluidas.
    """
    n = arr.size
    if n < SEED_LARGO:
        return None, None
    nv = n - SEED_LARGO + 1
    codigos = np.zeros(nv, dtype=np.int64)
    malo = np.zeros(nv, dtype=bool)
    for m in range(SEED_LARGO):
        tramo = arr[m:m + nv]
        codigos = (codigos << 2) | np.where(tramo < 4, tramo, 0).astype(np.int64)
        malo |= (tramo == 4)
    ok = ~malo
    return codigos[ok], np.nonzero(ok)[0]


ANCHO_FASTA = 60


def escribir_fasta(ruta, registros):
    """
    FASTA con las lineas CORTADAS a 60 caracteres.

    El ancho no es cosmetico: RNAhybrid trunca las lineas de secuencia largas,
    asi que un registro escrito en una sola linea de 2403 nt le llega cortado y
    pierde los hits de la parte final. SIREN no tropieza con esto porque usa
    SeqIO.write(), que envuelve a 60. Escribir "1 linea por secuencia" daba 7
    bloques donde el original da 12.
    """
    with open(ruta, "w") as fh:
        for nombre, seq in registros:
            fh.write(">%s\n" % nombre)
            for i in range(0, len(seq), ANCHO_FASTA):
                fh.write(seq[i:i + ANCHO_FASTA] + "\n")


def correr_rnahybrid(tgt_id, tgt_seq, sirnas, tmpdir, contador):
    """RNAhybrid de UN target contra solo los siRNAs que pasaron el filtro."""
    base = os.path.join(tmpdir, "t%d" % contador)
    f_t = base + ".t.fa"
    f_q = base + ".q.fa"
    escribir_fasta(f_t, [(tgt_id, tgt_seq)])
    escribir_fasta(f_q, sirnas)
    cmd = ["RNAhybrid", "-t", f_t, "-q", f_q] + OPTS_RNAHYBRID
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True, check=False)
        if p.returncode != 0:
            # No silenciar fallos: un error mudo aca se vuelve un off-target
            # ausente mas adelante, imposible de distinguir de un hit real.
            sys.stderr.write("RNAhybrid fallo en %s (codigo %d): %s\n"
                             % (tgt_id, p.returncode, (p.stderr or "")[:200]))
        return p.stdout or ""
    finally:
        for f in (f_t, f_q):
            try:
                os.remove(f)
            except OSError:
                pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--targets", required=True,
                    help="targets_prefiltered.fa de SIREN (o el FASTA completo)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    t0 = time.time()
    registros = leer_fasta(args.targets)

    # El gen blanco se separa igual que en generate_sirnas(): el registro cuyo
    # id CONTIENE el nombre del gen. Si matchea mas de uno, SIREN aborta; aca
    # tambien, porque el resto del analisis seria ambiguo.
    blanco = [r for r in registros if args.gene in r[0]]
    if len(blanco) != 1:
        print("error: %d registros matchean '%s'" % (len(blanco), args.gene),
              file=sys.stderr)
        return 1
    blanco_id, blanco_seq = blanco[0]
    off_targets = [r for r in registros if r[0] != blanco_id]

    sirnas = generar_sirnas(blanco_seq)
    print("gen %s: %d nt, %d siRNAs, %d targets off"
          % (args.gene, len(blanco_seq), len(sirnas), len(off_targets)))

    # Indice: codigo de 16-mero -> indices de siRNA
    sin_filtro = []
    codigos_list = []
    idx_list = []
    for i, (nombre, s) in enumerate(sirnas):
        cods = codigos_aceptables(s)
        if cods is None:
            sin_filtro.append(i)
            continue
        codigos_list.append(cods)
        idx_list.append(np.full(cods.size, i, dtype=np.int32))

    if codigos_list:
        todos_cod = np.concatenate(codigos_list)
        todos_idx = np.concatenate(idx_list)
        orden = np.argsort(todos_cod, kind="stable")
        todos_cod = todos_cod[orden]
        todos_idx = todos_idx[orden]
    else:
        todos_cod = np.zeros(0, dtype=np.int64)
        todos_idx = np.zeros(0, dtype=np.int32)

    print("indice: %d codigos de seed (%.1f por siRNA); %d siRNAs sin filtrar"
          % (todos_cod.size,
             todos_cod.size / max(1, len(sirnas) - len(sin_filtro)),
             len(sin_filtro)))
    if sin_filtro:
        print("  (esos %d se mandan contra TODOS los targets, para no perder hits)"
              % len(sin_filtro))

    # Escaneo
    pares = {}       # tgt_id -> set(indices de siRNA)
    nt_total = 0
    for tgt_id, tgt_seq in off_targets:
        nt_total += len(tgt_seq)
        arr = codificar(tgt_seq)
        cods, _ = codigos_de_target(arr)
        if cods is None or cods.size == 0:
            continue
        if todos_cod.size:
            izq = np.searchsorted(todos_cod, cods, side="left")
            der = np.searchsorted(todos_cod, cods, side="right")
            hits = der > izq
            if np.any(hits):
                sel = set()
                for a, b in zip(izq[hits], der[hits]):
                    sel.update(todos_idx[a:b].tolist())
                if sel:
                    pares.setdefault(tgt_id, set()).update(sel)

    # los siRNAs sin filtro van contra todos los targets
    if sin_filtro:
        for tgt_id, _ in off_targets:
            pares.setdefault(tgt_id, set()).update(sin_filtro)

    n_pares = sum(len(v) for v in pares.values())
    n_pares_totales = len(sirnas) * len(off_targets)
    print("escaneo: %.1fs, %d nt de targets" % (time.time() - t0, nt_total))
    print("pares sobrevivientes: %d de %d (%.4f%%) en %d targets"
          % (n_pares, n_pares_totales,
             100.0 * n_pares / max(1, n_pares_totales), len(pares)))

    # RNAhybrid solo sobre los pares sobrevivientes
    other = os.path.join(args.outdir, "other_files")
    os.makedirs(other, exist_ok=True)
    salida = os.path.join(other, "all_targets.txt")

    mapa_seq = dict(off_targets)
    t1 = time.time()
    tmpdir = tempfile.mkdtemp(prefix="pares_")
    try:
        with open(salida, "w") as out_all, \
             ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
            futs = {}
            for k, (tgt_id, idxs) in enumerate(pares.items()):
                # orden de sirnas.fa, que es el que define el orden de salida
                subset = [sirnas[i] for i in sorted(idxs)]
                futs[ex.submit(correr_rnahybrid, tgt_id, mapa_seq[tgt_id],
                               subset, tmpdir, k)] = tgt_id
            hechos = 0
            for fut in as_completed(futs):
                txt = fut.result()
                if txt:
                    out_all.write(txt)
                hechos += 1
                if hechos % 500 == 0:
                    print("  RNAhybrid %d/%d" % (hechos, len(futs)), flush=True)
    finally:
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass

    print("RNAhybrid: %.1fs; escrito %s (%d bytes)"
          % (time.time() - t1, salida, os.path.getsize(salida)))
    print("TOTAL %.1fs" % (time.time() - t0))
    print("\nAhora correr SIREN con --outdir %s : va a saltear RNAhybrid."
          % args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
