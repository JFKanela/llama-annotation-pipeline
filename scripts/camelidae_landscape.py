#!/usr/bin/env python3
"""
Panorama de recursos genomicos anotados en Camelidae, consultado a la API de NCBI.

QUE PRODUCE
  - camelidae.json con la respuesta cruda de NCBI Datasets (evidencia archivable).
  - Tabla de ensamblados CON anotacion, por especie y fecha de release.
  - Recuento de ensamblados SIN anotacion por taxon.

QUE SOSTIENE EN EL ARTICULO
  1. La afirmacion central: la llama es el unico camelido domestico sin anotacion
     estructural publica, pese a disponer de un ensamblado a nivel de cromosoma.
     Deja de ser una impresion y pasa a ser una consulta reproducible que cualquier
     revisor puede repetir. CONSERVAR camelidae.json con su fecha.
  2. La lista de proteomas disponibles para la comparativa con OrthoFinder.

DOS CORRECCIONES QUE ESTE ANALISIS OBLIGO A HACER
  a) El argumento para elegir Vicugna pacos como referencia debe ser EXCLUSIVAMENTE
     filogenetico (llama y alpaca son Lamini; los Camelus son Camelini). NO usar el
     argumento de calidad de anotacion: Camelus bactrianus tiene anotacion cromosomica
     de 2025-06-04, dos meses MAS RECIENTE que la de alpaca (2025-04-02), y el
     dromedario tiene mCamDro1.pat de 2024. El argumento de calidad seria falso y un
     revisor lo detectaria con esta misma consulta.
  b) La vicuna NO puede entrar en la comparativa con OrthoFinder: Vicugna vicugna
     mensalis no tiene anotacion publica. El plan inicial la incluia.

USO
  python3 camelidae_landscape.py            # descarga y analiza
  python3 camelidae_landscape.py --local    # analiza camelidae.json ya descargado
"""

import json
import sys
import urllib.request
from collections import defaultdict

API = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/"
       "Camelidae/dataset_report?page_size=1000")
CACHE = "camelidae.json"


def fetch():
    print(f"Consultando {API}")
    with urllib.request.urlopen(API) as r:
        data = r.read()
    with open(CACHE, "wb") as f:
        f.write(data)
    print(f"Guardado {CACHE} ({len(data)} bytes)")
    return json.loads(data)


def main():
    if "--local" in sys.argv:
        d = json.load(open(CACHE))
    else:
        d = fetch()

    reps = d.get("reports", [])
    print(f"\nregistros devueltos: {len(reps)}\n")

    annotated, unannotated = [], []
    for r in reps:
        acc = r.get("accession", "")
        org = r.get("organism", {}).get("organism_name", "")
        ai = r.get("assembly_info", {}) or {}
        an = r.get("annotation_info", {}) or {}
        row = (org, acc, ai.get("assembly_name", ""), ai.get("assembly_level", ""),
               an.get("release_date", ""), an.get("name", ""))
        (annotated if an.get("release_date") else unannotated).append(row)

    annotated.sort(key=lambda x: (x[0], x[4]))
    print("=== CON ANOTACION ===")
    print(f"{'especie':26s} {'accesion':19s} {'ensamblado':24s} {'nivel':11s} release")
    for org, acc, name, lvl, date, _ in annotated:
        print(f"{org[:25]:26s} {acc:19s} {name[:23]:24s} {lvl:11s} {date}")

    print(f"\n=== SIN ANOTACION: {len(unannotated)} ensamblados ===")
    by_taxon = defaultdict(int)
    for org, *_ in unannotated:
        by_taxon[org] += 1
    for k, v in sorted(by_taxon.items()):
        print(f"  {k:30s} {v}")

    print("\n=== ENSAMBLADOS SIN ANOTAR DE INTERES ===")
    for org, acc, name, lvl, *_ in sorted(unannotated):
        if any(k in org for k in ("glama", "guanicoe", "vicugna mensalis")):
            print(f"  {org:28s} {acc:19s} {name:22s} {lvl}")


if __name__ == "__main__":
    main()
