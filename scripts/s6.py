#!/usr/bin/env python3
"""Tabla suplementaria S6: correspondencia con la anotacion de referencia.

QUE PRODUCE
  El recuento de genes, mRNA e isoformas del GFF3 NATIVO de Vicugna pacos
  (GCF_048564905.1), que es el que consume LiftOn, frente a la reduccion por AGAT a
  un transcrito por gen que consumen Liftoff y miniprot.

QUE SOSTIENE EN EL ARTICULO
  Cuantifica la asimetria de referencia entre los tres brazos de homologia, que es
  una limitacion declarada y no un detalle menor: LiftOn dispuso de 2,67 mRNA por
  gen codificante y los otros dos de uno solo. Ver README.md, seccion de asimetria
  de referencia.

  De los 37.446 genes del GFF3 nativo, 21.233 tienen al menos un mRNA y 16.213 no
  tienen ninguno: son los no codificantes. Ese 21.233 es exactamente el numero de
  mRNA que quedan tras la reduccion por AGAT, que es lo que cabe esperar de una
  reduccion a isoforma primaria.

USO
  python3 scripts/s6.py
"""

import gzip, os, collections
BASE = os.path.expanduser("~/llama_annotation_pipeline")
NAT = os.path.join(BASE, "ref", "GCF_048564905.1_VicPac4_genomic.gff.gz")

genes, mrna_por_gen = set(), collections.defaultdict(list)
with gzip.open(NAT, "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9: continue
        d = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
        if c[2] == "gene" and d.get("ID"):
            genes.add(d["ID"])
        elif c[2] == "mRNA" and d.get("Parent"):
            mrna_por_gen[d["Parent"]].append(d.get("ID", "?"))

n = sum(len(v) for v in mrna_por_gen.values())
print("genes totales      :", len(genes))
print("genes con mRNA     :", len(mrna_por_gen))
print("mRNA totales       :", n)
print("genes con 1 mRNA   :", sum(1 for v in mrna_por_gen.values() if len(v) == 1))
print("genes con >1 mRNA  :", sum(1 for v in mrna_por_gen.values() if len(v) > 1))
print("maximo isoformas   :", max(len(v) for v in mrna_por_gen.values()))
