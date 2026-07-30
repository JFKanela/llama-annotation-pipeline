#!/usr/bin/env python3
"""
Analisis de cobertura de alineamiento por brazo de anotacion.

QUE PRODUCE
  Por cada proteoma y contra cada referencia:
    - porcentaje de proteinas con hit
    - cobertura de SUJETO: cuanto del ortologo de referencia cubre la proteina.
      Baja indica TRUNCAMIENTO del modelo.
    - cobertura de CONSULTA: cuanto de la proteina alinea con el ortologo.
      Baja, con sujeto alta, indica SOBRE-EXTENSION: el modelo lleva secuencia
      que el ortologo no tiene.
    - identidad media

QUE RESUELVE
  La pregunta abierta sobre Helixer. Sus proteinas son mas largas que las de los
  brazos de homologia (550,6 aa frente a 524-535), y hasta ahora no era posible
  saber si eso indica mejor modelado o sobre-extension. Las dos coberturas juntas
  lo distinguen:

    sujeto ALTA  + consulta ALTA  -> modelos correctos
    sujeto BAJA  + consulta ALTA  -> modelos truncados
    sujeto ALTA  + consulta BAJA  -> SOBRE-EXTENSION
    ambas bajas                   -> alineamiento parcial o modelo dudoso

INTERPRETACION DE LAS DOS REFERENCIAS
  Contra Camelus dromedarius: es la comparacion valida. Ninguno de los cuatro
  brazos deriva de ese proteoma, de modo que los cuatro son comparables.

  Contra Vicugna pacos: NO es evidencia de calidad. Tres de los cuatro brazos son
  proyecciones de esa anotacion, asi que mide fidelidad de transferencia y saldra
  excelente por construccion. Sirve como control interno entre los brazos de
  homologia, y para ver cuanto se aparta Helixer, que es el unico independiente.

LIMITACION DECLARADA
  Se usa un solo HSP por emparejamiento (--max-hsps 1). En proteinas multidominio
  cuyo alineamiento se fragmenta, la cobertura sale artificialmente baja. Es una
  aproximacion conservadora: subestima la cobertura, no la infla.

USO
  python3 analyze_blastp.py
"""

import os
from collections import defaultdict

BASE = os.path.expanduser("~/llama_annotation_pipeline")
WORK = os.path.join(BASE, "blastp")
MAP = os.path.join(BASE, "md5_to_ids.tsv")
ARMS = ["liftoff", "miniprot", "lifton", "helixer"]
REFS = [("dromedario", "REFERENCIA EXTERNA, la comparacion valida"),
        ("alpaca", "CONTROL INTERNO, circular para los brazos de homologia")]


def load_map():
    md5_to_arms = defaultdict(list)
    with open(MAP) as fh:
        next(fh)
        for line in fh:
            h, arm, _ = line.rstrip("\n").split("\t")
            md5_to_arms[h].append(arm)
    return md5_to_arms


def load_hits(path):
    """md5 -> (cobertura_consulta, cobertura_sujeto, identidad)"""
    hits = {}
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            q, pid = c[0], float(c[2])
            qs, qe, ss, se = int(c[4]), int(c[5]), int(c[6]), int(c[7])
            qlen, slen = int(c[10]), int(c[11])
            qcov = 100.0 * (qe - qs + 1) / qlen if qlen else 0
            scov = 100.0 * (se - ss + 1) / slen if slen else 0
            hits[q] = (qcov, scov, pid)
    return hits


def mean(v):
    return sum(v) / len(v) if v else 0.0


def median(v):
    if not v:
        return 0.0
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def main():
    md5_to_arms = load_map()
    total = defaultdict(int)
    for arms in md5_to_arms.values():
        for a in arms:
            total[a] += 1

    for ref, nota in REFS:
        path = os.path.join(WORK, f"hits_{ref}.tsv")
        if not os.path.exists(path):
            print(f"\n[{ref}] no hay resultados todavia ({path})")
            continue

        hits = load_hits(path)
        print(f"\n{'=' * 78}")
        print(f"CONTRA {ref.upper()}   ({nota})")
        print(f"{'=' * 78}")

        per_arm = {a: {"qcov": [], "scov": [], "pid": [], "n": 0} for a in ARMS}
        for h, arms in md5_to_arms.items():
            if h not in hits:
                continue
            qcov, scov, pid = hits[h]
            for a in arms:
                per_arm[a]["qcov"].append(qcov)
                per_arm[a]["scov"].append(scov)
                per_arm[a]["pid"].append(pid)
                per_arm[a]["n"] += 1

        print(f"{'brazo':10s} {'prot':>6s} {'con hit':>8s} {'%':>6s} "
              f"{'ident':>6s} {'cobSUJ':>7s} {'cobCON':>7s} {'suj>=80%':>9s} {'con>=80%':>9s}")
        for a in ARMS:
            d = per_arm[a]
            t = total[a]
            if not d["n"]:
                print(f"{a:10s} {t:>6d}  sin resultados")
                continue
            s80 = 100.0 * sum(1 for x in d["scov"] if x >= 80) / d["n"]
            q80 = 100.0 * sum(1 for x in d["qcov"] if x >= 80) / d["n"]
            print(f"{a:10s} {t:>6d} {d['n']:>8d} {100*d['n']/t:>5.1f}% "
                  f"{mean(d['pid']):>5.1f}% {median(d['scov']):>6.1f}% "
                  f"{median(d['qcov']):>6.1f}% {s80:>8.1f}% {q80:>8.1f}%")

        print("\n  cobSUJ y cobCON son MEDIANAS. Lectura:")
        print("    sujeto alta + consulta alta -> modelos correctos")
        print("    sujeto BAJA + consulta alta -> modelos truncados")
        print("    sujeto alta + consulta BAJA -> SOBRE-EXTENSION")

        # senal de sobre-extension: consulta claramente por debajo de sujeto
        print("\n  --- indicio de sobre-extension (consulta < sujeto - 20 puntos) ---")
        for a in ARMS:
            d = per_arm[a]
            if not d["n"]:
                continue
            over = sum(1 for qc, sc in zip(d["qcov"], d["scov"]) if sc - qc > 20)
            print(f"    {a:10s} {over:>6d} de {d['n']} ({100*over/d['n']:.1f}%)")


if __name__ == "__main__":
    main()
