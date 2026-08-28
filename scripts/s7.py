#!/usr/bin/env python3
"""Tabla suplementaria S7: sensibilidad al umbral de sobre-extension.

QUE PRODUCE
  S7_sensibilidad_umbral.tsv. La fraccion de secuencias que parecen sobre-extendidas
  en cada brazo, contra cada una de las dos referencias, con el umbral a 10, 15, 20
  y 30 puntos porcentuales de diferencia entre cobertura de sujeto y de consulta.

QUE SOSTIENE EN EL ARTICULO
  Que el hallazgo de circularidad NO es un artefacto del umbral elegido. Contra
  alpaca, el contraste entre Helixer y LiftOn se mantiene entre 27 y 51 veces en los
  cuatro umbrales; contra dromedario, entre 1,1 y 1,4 veces en los cuatro. La
  conclusion es la misma se ponga el umbral donde se ponga.

  El umbral de 20 pp que usa analyze_blastp.py y que se reporta en el manuscrito no
  es, por tanto, una eleccion afortunada: es un punto cualquiera de una meseta.

USO
  python3 scripts/s7.py
"""

import os
from collections import defaultdict
BASE = os.path.expanduser("~/llama_annotation_pipeline")
NOM = {"liftoff":"Liftoff","miniprot":"miniprot","lifton":"LiftOn","helixer":"Helixer"}
ARMS = ["Liftoff","miniprot","LiftOn","Helixer"]

md5 = defaultdict(list)
with open(os.path.join(BASE,"md5_to_ids.tsv")) as fh:
    next(fh)
    for line in fh:
        h, arm, _ = line.rstrip("\n").split("\t")
        md5[h].append(NOM[arm])
print("secuencias unicas:", len(md5))

def cov(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12: continue
            qs,qe,ss,se = int(c[4]),int(c[5]),int(c[6]),int(c[7])
            ql,sl = int(c[10]),int(c[11])
            d[c[0]] = (100.0*(qe-qs+1)/ql if ql else 0,
                       100.0*(se-ss+1)/sl if sl else 0)
    return d

REF = {"alpaca": os.path.join(BASE, "blastp", "hits_alpaca.tsv"),
       "dromedario": os.path.join(BASE, "blastp", "hits_dromedario.tsv")}
COV = {k: cov(v) for k, v in REF.items()}
for k, v in COV.items():
    print(k, "alineamientos:", len(v))

SAL = os.path.join(BASE, "tablas", "S7_sensibilidad_umbral.tsv")
with open(SAL, "w") as out:
    out.write("umbral_pp\treferencia\tbrazo\tn\tsobre_extension_pct\n")
    print()
    print("umbral   ref            " + "  ".join("%9s" % a for a in ARMS))
    for U in (10, 15, 20, 30):
        for ref in ("alpaca", "dromedario"):
            c = COV[ref]
            tot = defaultdict(int); over = defaultdict(int)
            for h, arms in md5.items():
                if h not in c: continue
                q, s = c[h]
                for a in arms:
                    tot[a] += 1
                    if s - q > U: over[a] += 1
            fila = []
            for a in ARMS:
                pct = 100.0 * over[a] / tot[a] if tot[a] else 0.0
                fila.append("%9.2f" % pct)
                out.write("%d\t%s\t%s\t%d\t%.2f\n" % (U, ref, a, tot[a], pct))
            print("%4dpp  %-12s  %s" % (U, ref, "  ".join(fila)))
print()
print("Escrito:", SAL)
