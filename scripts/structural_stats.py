#!/usr/bin/env python3
"""Estadisticas estructurales de los cuatro proteomas. VERSION CORREGIDA.

CORRECCION del 14 de agosto de 2026, importante. La version anterior contaba
las CDS y los exones agrupando por su atributo Parent SIN COMPROBAR que ese
padre fuera un mRNA. En los GFF3 de homologia el atributo exon cuelga tambien
de lnc_RNA, tRNA, snRNA, snoRNA y rRNA: en Liftoff hay 34545 padres de exon y
solo 20306 son mRNA. Eso inflaba el recuento de monoexonicos de 2516 a 6009.

En CDS la contaminacion era menor (49 padres de 20075) y el resultado apenas
cambiaba, pero el metodo estaba mal igualmente. Y Helixer no emite features no
codificantes, de modo que el sesgo NO era simetrico entre brazos: afectaba solo
a los tres de homologia, que es justo lo que invalida una comparacion.

CRITERIO ADOPTADO: una sola CDS, no un solo exon. Es lo coherente con un
articulo de proteomas, donde importa la estructura de la region codificante y
no la del transcrito completo, y es el unico aplicable a los cuatro brazos,
porque miniprot no emite features exon.

Se reportan ambos criterios para poder contrastarlos con AGAT, que usa exones.
"""
import gzip, os
import statistics as st
from collections import defaultdict

ZEN = os.path.expanduser("~/zenodo_data")
BASE = os.path.expanduser("~/llama_annotation_pipeline")

FAA = {"Liftoff":  os.path.join(ZEN, "llama_liftoff.faa.gz"),
       "miniprot": os.path.join(ZEN, "llama_miniprot.faa.gz"),
       "LiftOn":   os.path.join(ZEN, "llama_lifton.faa.gz"),
       "Helixer":  os.path.join(BASE, "Lgla_hx036_helixer.faa")}
GFF = {"Liftoff":  os.path.join(ZEN, "llama_liftoff.gff3.gz"),
       "miniprot": os.path.join(ZEN, "llama_miniprot.gff3.gz"),
       "LiftOn":   os.path.join(ZEN, "llama_lifton.gff3.gz"),
       "Helixer":  os.path.join(BASE, "Lgla_hx036_helixer_FINAL.gff")}
ARMS = list(FAA)


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")


def parse_gff(path):
    """Devuelve (CDS por mRNA, exones por mRNA), FILTRANDO padres a mRNA."""
    mrna = set()
    cds, exon = defaultdict(int), defaultdict(int)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            d = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            if c[2] == "mRNA" and d.get("ID"):
                mrna.add(d["ID"])
            par = d.get("Parent")
            if not par:
                continue
            if c[2] == "CDS":
                cds[par] += 1
            elif c[2] == "exon":
                exon[par] += 1
    # AQUI esta la correccion: solo se conservan los padres que son mRNA
    return ({k: v for k, v in cds.items() if k in mrna},
            {k: v for k, v in exon.items() if k in mrna},
            len(mrna))


def prot_lengths(path):
    out, cur = [], 0
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    out.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
        if cur:
            out.append(cur)
    return out


print(f"{'brazo':10s} {'#prot':>7s} {'#mRNA':>7s} {'longMed':>8s} "
      f"{'CDS/mRNA':>9s} {'mono CDS':>10s} {'mono exon':>11s}")
for a in ARMS:
    L = prot_lengths(FAA[a])
    cds, exon, n_mrna = parse_gff(GFF[a])
    cv = list(cds.values())
    mono_cds = 100 * sum(1 for x in cv if x == 1) / len(cv) if cv else 0
    ev = list(exon.values())
    mono_ex = 100 * sum(1 for x in ev if x == 1) / len(ev) if ev else float("nan")
    ex_txt = f"{mono_ex:9.1f}%" if ev else "      n/d"
    print(f"{a:10s} {len(L):>7d} {n_mrna:>7d} {st.mean(L):>8.1f} "
          f"{st.mean(cv):>9.2f} {mono_cds:>9.1f}% {ex_txt}")
