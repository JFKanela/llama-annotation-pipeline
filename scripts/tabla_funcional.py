#!/usr/bin/env python3
"""Tabla 2: anotacion funcional por brazo (InterProScan 5.78-109.0).

Sustituye a la figura 5 del plan, que discriminaba poco: los cuatro brazos
quedan entre 86.6 % y 90.2 % de cobertura Pfam.

"Con dominio" excluye MobiDBLite y Coils, que son prediccion de desorden y de
estructura enrollada y anotan casi cualquier proteina.

CAUTELA que debe constar al presentarla: los denominadores NO son equivalentes.
Helixer predice 18765 proteinas frente a 20000-21000 de los demas, y solo emite
modelos donde consigue construir un ORF completo, de modo que su conjunto esta
preseleccionado por plausibilidad estructural y su porcentaje parte con ventaja.
"""
import glob, os
from collections import defaultdict

BASE = os.path.expanduser("~/llama_annotation_pipeline")
OUT = os.path.expanduser(
    "~/Gdrive/Doctorado/llama_annotation_pipeline/MANUSCRITO/02_tablas")
os.makedirs(OUT, exist_ok=True)

NOMBRE = {"liftoff": "Liftoff", "miniprot": "miniprot",
          "lifton": "LiftOn", "helixer": "Helixer"}
ARMS = ["Liftoff", "miniprot", "LiftOn", "Helixer"]
RUIDO = {"MobiDBLite", "Coils"}

md5_arms = defaultdict(list)
with open(os.path.join(BASE, "md5_to_ids.tsv")) as fh:
    next(fh)
    for line in fh:
        h, arm, _ = line.rstrip("\n").split("\t")
        md5_arms[h].append(NOMBRE[arm])

hits = defaultdict(set)
n_lineas = 0
for f in sorted(glob.glob(os.path.join(BASE, "ips_out", "chunk_*.tsv"))):
    with open(f) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            hits[c[0]].add(c[3])
            n_lineas += 1

tot = defaultdict(int); algo = defaultdict(int)
dom = defaultdict(int); pfam = defaultdict(int); ipr = defaultdict(int)
for h, arms in md5_arms.items():
    dbs = hits.get(h, set())
    reales = dbs - RUIDO
    for a in arms:
        tot[a] += 1
        if dbs: algo[a] += 1
        if reales: dom[a] += 1
        if "Pfam" in dbs: pfam[a] += 1

print(f"anotaciones totales: {n_lineas}")
print(f"proteinas unicas con algun resultado: {len(hits)} de {len(md5_arms)}")
print()
md = ["## Tabla 2 · Anotación funcional (InterProScan 5.78-109.0)\n",
      "| Annotation | Proteins | With domain | % | With Pfam | % |",
      "|---|---|---|---|---|---|"]
for a in ARMS:
    t = tot[a]
    print(f"{a:9s} {t:>6d}  dominio {dom[a]:>6d} ({100*dom[a]/t:5.1f}%)  "
          f"Pfam {pfam[a]:>6d} ({100*pfam[a]/t:5.1f}%)")
    md.append(f"| {a} | {t:,} | {dom[a]:,} | {100*dom[a]/t:.1f} | "
              f"{pfam[a]:,} | {100*pfam[a]/t:.1f} |".replace(",", " "))

md.append(f"\n> {n_lineas:,} anotaciones sobre {len(hits):,} de las {len(md5_arms):,} "
          "secuencias únicas. Los 18 análisis por defecto de InterProScan, con el servicio "
          "de coincidencias precalculadas desactivado. «Con dominio» excluye MobiDBLite y "
          "Coils, que anotan casi cualquier proteína. PANTHER se reporta a nivel de familia: "
          "el formato tabular no emite subfamilia.\n".replace(",", " "))
md.append("> Los denominadores no son equivalentes. Helixer solo emite modelos donde "
          "consigue construir un ORF completo, de modo que su conjunto está preseleccionado "
          "por plausibilidad estructural.\n")

f = os.path.join(OUT, "Tabla2_funcional.md")
open(f, "w").write("\n".join(md))
print("\nescrito:", f)
