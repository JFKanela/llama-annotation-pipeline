#!/usr/bin/env python3
"""Capas de confianza del proteoma candidato combinado.

QUE PRODUCE
  Lgla_combined_proteome_confidence.tsv, una fila por cada una de las 20.708
  secuencias del proteoma combinado, con su origen, si arrastra un stop interno,
  si tiene hit contra Camelus dromedarius y en que capa cae.

LAS TRES CAPAS SON ANIDADAS, NO EXCLUYENTES EN CALIDAD
  high_confidence      LiftOn, pauta de lectura integra Y homologia externa   18.920
  extended_reference   el resto del nucleo LiftOn                              1.313
  extended_candidate   los loci de Helixer incorporados                          475
                                                                        --------------
                                                                              20.708

QUE SOSTIENE EN EL ARTICULO
  Es lo que respalda el "confidence-aware" del titulo. Sin este fichero, el
  proteoma combinado es una lista plana de 20.708 entradas de las que el usuario no
  puede saber cuales son transferencias limpias y cuales predicciones sin verificar.

POR QUE DROMEDARIO Y NO ALPACA
  Alpaca es la fuente de la que derivan los tres brazos de homologia: medir contra
  ella es circular. Dromedario es el patron externo. Ver README.md, seccion de
  exactitud contra referencias externas.

CAUTELA
  extended_candidate lleva external_hit = NA, no "yes". Esos loci entraron
  precisamente por tener ortologo en algun camelido, pero el cruce se hizo sobre la
  tabla de loci noveles y no sobre el mismo alineamiento que el nucleo, de modo que
  la columna no es comparable entre capas.

USO
  python3 scripts/capas.py
"""

import gzip, os, csv
BASE = os.path.expanduser("~/llama_annotation_pipeline")
ZEN  = os.path.expanduser("~/zenodo_data")
SAL  = os.path.join(BASE, "tablas", "Lgla_combined_proteome_confidence.tsv")

def op(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

lifton, stop_int = [], set()
seq, ident = [], None
def cierra():
    if ident is None: return
    s = "".join(seq)
    lifton.append(ident)
    if "." in s and not (s.count(".") == 1 and s.endswith(".")):
        stop_int.add(ident)
with op(os.path.join(ZEN, "llama_lifton.faa.gz")) as fh:
    for line in fh:
        if line.startswith(">"):
            cierra(); ident = line[1:].split()[0]; seq = []
        else:
            seq.append(line.strip())
    cierra()
print("LiftOn:", len(lifton), "proteinas,", len(stop_int), "con stop interno")

md5_hit = set()
with open(os.path.join(BASE, "blastp", "hits_dromedario.tsv")) as fh:
    for l in fh:
        md5_hit.add(l.split("\t")[0])
id_md5 = {}
with open(os.path.join(BASE, "md5_to_ids.tsv")) as fh:
    next(fh)
    for line in fh:
        h, brazo, ident2 = line.rstrip("\n").split("\t")
        if brazo == "lifton":
            id_md5[ident2] = h
print("mapa lifton:", len(id_md5), "identificadores")
con_hit = set()
for i in lifton:
    if id_md5.get(i) in md5_hit:
        con_hit.add(i)
print("LiftOn con hit contra dromedario:", len(con_hit))

red = set()
with open(os.path.join(BASE, "blastp", "rescued_vs_lifton.tsv")) as fh:
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) >= 6 and float(c[2]) >= 95 and int(c[3]) / int(c[4]) >= 0.8:
            red.add(c[0])
rescatados = []
with open(os.path.join(BASE, "blastp", "helixer_novel_loci.tsv")) as fh:
    r = csv.reader(fh, delimiter="\t"); next(r)
    for f in r:
        if (f[6] == "si" or f[7] == "si") and f[0] not in red:
            rescatados.append(f[0])
print("redundantes:", len(red), " loci incorporados:", len(rescatados))

n = {"high_confidence": 0, "extended_reference": 0, "extended_candidate": 0}
with open(SAL, "w") as out:
    out.write("protein_id\tsource\tinternal_stop\texternal_hit\tconfidence_layer\n")
    for i in lifton:
        st = "yes" if i in stop_int else "no"
        hi = "yes" if i in con_hit else "no"
        capa = "high_confidence" if (st == "no" and hi == "yes") else "extended_reference"
        n[capa] += 1
        out.write(i + "\tLiftOn\t" + st + "\t" + hi + "\t" + capa + "\n")
    for i in rescatados:
        n["extended_candidate"] += 1
        out.write(i + "\tHelixer\tno\tNA\textended_candidate\n")

print()
print("=== RESULTADO ===")
for k in ("high_confidence", "extended_reference", "extended_candidate"):
    print("  ", k, n[k])
print("   TOTAL", sum(n.values()), "(debe ser 20708)")
print("Escrito:", SAL)
