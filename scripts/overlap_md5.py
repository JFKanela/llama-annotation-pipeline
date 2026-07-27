#!/usr/bin/env python3
"""
Solapamiento por identidad exacta de secuencia (MD5) entre los cuatro proteomas.

QUE PRODUCE
  - Numero de secuencias unicas frente al total (dimensiona InterProScan).
  - Reparto por combinacion de brazos: las 14 intersecciones no vacias.
  - Redundancia interna de cada proteoma (secuencias repetidas dentro del mismo brazo).

QUE SOSTIENE EN EL ARTICULO
  - La figura de solapamiento entre proteomas (diagrama UpSet, no Venn de cuatro).
  - El nucleo de confianza maxima: secuencias identicas en los cuatro brazos.
  - La redundancia interna de miniprot, muy superior a la de los demas.

ADVERTENCIA METODOLOGICA IMPORTANTE
  Esta comparacion es por identidad EXACTA de secuencia. Las secuencias exclusivas
  de Helixer NO son genes especificos de llama: entre prediccion ab initio y
  transferencia por homologia la coincidencia byte a byte es rara por construccion.
  Para novedad real hay que usar solapamiento de COORDENADAS a nivel de locus,
  que es lo que hace overlap_coords.py.

USO
  python3 overlap_md5.py
"""

import gzip
import hashlib
from collections import defaultdict

# Ajustar rutas si cambian de sitio. Helixer apunta al proteoma definitivo
# derivado de la ejecucion en la web tool con overlap.
FILES = {
    "liftoff":  "llama_liftoff.faa.gz",
    "miniprot": "llama_miniprot.faa.gz",
    "lifton":   "llama_lifton.faa.gz",
    "helixer":  "Lgla_hx036_helixer.faa",
}
ARMS = list(FILES)


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_fasta_md5(path):
    """Devuelve la lista de MD5 de cada secuencia del fichero, en orden."""
    out = []
    seq = []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    out.append(hashlib.md5("".join(seq).encode()).hexdigest())
                seq = []
            else:
                seq.append(line.strip())
        if seq:
            out.append(hashlib.md5("".join(seq).encode()).hexdigest())
    return out


def main():
    md5_to_arms = defaultdict(set)
    counts = {}

    for arm, path in FILES.items():
        hashes = read_fasta_md5(path)
        counts[arm] = len(hashes)
        for h in hashes:
            md5_to_arms[h].add(arm)

    total = sum(counts.values())
    unique = len(md5_to_arms)

    print("--- recuento por brazo ---")
    for arm in ARMS:
        print(f"  {arm:10s} {counts[arm]:>7d}")
    print(f"  {'TOTAL':10s} {total:>7d}")
    print(f"  {'UNICAS':10s} {unique:>7d}")
    print(f"  AHORRO     {100 * (1 - unique / total):>6.1f} %")

    print("\n--- en cuantos brazos aparece cada secuencia unica ---")
    dist = defaultdict(int)
    for arms in md5_to_arms.values():
        dist[len(arms)] += 1
    for k in sorted(dist):
        print(f"  en {k} brazo(s): {dist[k]:>7d}")

    print("\n--- reparto por combinacion (las intersecciones de la figura UpSet) ---")
    combo = defaultdict(int)
    for arms in md5_to_arms.values():
        key = " + ".join(a for a in ARMS if a in arms)
        combo[key] += 1
    for key, n in sorted(combo.items(), key=lambda x: -x[1]):
        print(f"  {n:>7d}  {key}")

    print("\n--- redundancia interna (secuencias repetidas dentro del mismo brazo) ---")
    print(f"  {'brazo':10s} {'en FASTA':>9s} {'unicas':>8s} {'duplicadas':>11s} {'%':>6s}")
    for arm in ARMS:
        u = sum(1 for arms in md5_to_arms.values() if arm in arms)
        dup = counts[arm] - u
        print(f"  {arm:10s} {counts[arm]:>9d} {u:>8d} {dup:>11d} {100 * dup / counts[arm]:>5.1f}%")


if __name__ == "__main__":
    main()
