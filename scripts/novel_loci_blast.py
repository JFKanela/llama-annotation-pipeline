#!/usr/bin/env python3
"""
Caracterizacion de los loci de Helixer sin solape posicional con la homologia.

QUE HACE
  1. Recalcula los loci novel de Helixer conservando los IDENTIFICADORES, cosa que
     overlap_coords.py no hace porque fusiona intervalos y pierde la correspondencia.
  2. Los cruza con los alineamientos contra dromedario y alpaca que YA EXISTEN, sin
     recalcular nada.
  3. Separa dos poblaciones que no significan lo mismo:
       - CON hit en camelido: gen real que la transferencia por homologia perdio.
         La novedad es un fallo de transferencia, no un gen propio de la especie.
       - SIN hit en ningun camelido: candidato a especifico de linaje, o falso
         positivo de Helixer. Solo estos necesitan busqueda contra bases publicas.
  4. Escribe el FASTA de los que no tienen hit, listo para alinear contra Swiss-Prot.

POR QUE IMPORTA
  Sin esta separacion, la cifra de novedad posicional no puede interpretarse. Y esa
  cifra depende de la unidad que se cuente: overlap_coords.py fusiona intervalos y da
  loci, este script evalua cada mRNA y da transcritos. Ninguna de las dos se escribe
  aqui: ambas se imprimen al ejecutar. Ver README.md, seccion de novedad posicional.
  Con ella, el manuscrito puede afirmar cuantos son fallos de transferencia recuperados
  por la prediccion ab initio, que es un resultado util y defendible, y cuantos quedan
  como candidatos genuinos, que es un resultado mas fuerte pero que exige mas evidencia.

DIFERENCIA CON overlap_coords.py
  Alli se fusionan los intervalos de Helixer antes de comparar, de modo que dos mRNA
  solapados cuentan como un solo locus. Aqui se evalua CADA mRNA por separado, porque
  hace falta la correspondencia con su proteina. El recuento puede diferir ligeramente
  del de aquel analisis, y es correcto que difiera: son unidades distintas.

USO
  python3 novel_loci_blast.py
"""

import gzip
import os
from collections import defaultdict

BASE = os.path.expanduser("~/llama_annotation_pipeline")
ZEN = os.path.expanduser("~/zenodo_data")

HELIXER_GFF = os.path.join(BASE, "Lgla_hx036_helixer_FINAL.gff")
HELIXER_FAA = os.path.join(BASE, "Lgla_hx036_helixer.faa")
MAP = os.path.join(BASE, "md5_to_ids.tsv")
HITS = {
    "dromedario": os.path.join(BASE, "blastp", "hits_dromedario.tsv"),
    "alpaca": os.path.join(BASE, "blastp", "hits_alpaca.tsv"),
}
HOMOLOGY_GFF = {
    "liftoff": os.path.join(ZEN, "llama_liftoff.gff3.gz"),
    "miniprot": os.path.join(ZEN, "llama_miniprot.gff3.gz"),
    "lifton": os.path.join(ZEN, "llama_lifton.gff3.gz"),
}
OUT_FASTA = os.path.join(BASE, "blastp", "helixer_novel_nohit.faa")
OUT_TSV = os.path.join(BASE, "blastp", "helixer_novel_loci.tsv")

MIN_LEN = 1000   # umbral de longitud de locus, el mismo de overlap_coords.py


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def load_mrna_with_ids(path):
    """(scaffold, start, end, ID) por cada mRNA."""
    out = []
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "mRNA":
                continue
            attrs = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            out.append((c[0], int(c[3]) - 1, int(c[4]), attrs.get("ID", "")))
    return out


def load_intervals(path):
    d = defaultdict(list)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 5 or c[2] != "mRNA":
                continue
            d[c[0]].append((int(c[3]) - 1, int(c[4])))
    return d


def merge(intervals):
    if not intervals:
        return []
    ints = sorted(intervals)
    out = [list(ints[0])]
    for s, e in ints[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def read_fasta(path):
    seqs, hdr, buf = {}, None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if hdr:
                    seqs[hdr] = "".join(buf)
                hdr, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if hdr:
            seqs[hdr] = "".join(buf)
    return seqs


def main():
    # --- 1. union de la homologia
    acc = defaultdict(list)
    for path in HOMOLOGY_GFF.values():
        for scf, iv in load_intervals(path).items():
            acc[scf].extend(iv)
    homology = {scf: merge(iv) for scf, iv in acc.items()}
    print(f"homologia: {sum(len(v) for v in homology.values())} loci fusionados")

    # --- 2. mRNA de Helixer sin solape, conservando ID
    hx = load_mrna_with_ids(HELIXER_GFF)
    print(f"helixer:   {len(hx)} mRNA")

    novel = []
    for scf, s, e, mid in hx:
        biv = homology.get(scf, [])
        if not any(not (be <= s or bs >= e) for bs, be in biv):
            novel.append((scf, s, e, mid, e - s))

    long_novel = [x for x in novel if x[4] >= MIN_LEN]
    print(f"\nmRNA de Helixer sin solape posicional: {len(novel)}")
    print(f"  de ellos, >= {MIN_LEN} bp: {len(long_novel)}")

    # --- 3. ID -> md5, desde el mapa ya existente
    id_to_md5 = {}
    with open(MAP) as fh:
        next(fh)
        for line in fh:
            h, arm, oid = line.rstrip("\n").split("\t")
            if arm == "helixer":
                id_to_md5[oid] = h
    print(f"\nmapa de helixer: {len(id_to_md5)} identificadores")

    # --- 4. cruce con los alineamientos ya calculados
    hits = {}
    for ref, path in HITS.items():
        s = set()
        if os.path.exists(path):
            with open(path) as fh:
                for line in fh:
                    s.add(line.split("\t", 1)[0])
        hits[ref] = s
        print(f"  hits contra {ref}: {len(s)}")

    con_dro = con_alp = sin_nada = sin_md5 = 0
    filas, sin_hit_ids = [], []

    for scf, s, e, mid, length in long_novel:
        md5 = id_to_md5.get(mid)
        if md5 is None:
            sin_md5 += 1
            continue
        d = md5 in hits.get("dromedario", set())
        a = md5 in hits.get("alpaca", set())
        if d:
            con_dro += 1
        if a:
            con_alp += 1
        if not d and not a:
            sin_nada += 1
            sin_hit_ids.append(mid)
        filas.append((mid, scf, s, e, length, md5, "si" if d else "no", "si" if a else "no"))

    n = len(filas)
    print(f"\n=== LOCI NOVEL DE HELIXER, >= {MIN_LEN} bp: {n} evaluados ===")
    if sin_md5:
        print(f"  ({sin_md5} sin correspondencia en el mapa; revisar si es alto)")
    if n:
        print(f"  con hit en dromedario         : {con_dro:>5d}  ({100*con_dro/n:.1f}%)")
        print(f"  con hit en alpaca             : {con_alp:>5d}  ({100*con_alp/n:.1f}%)")
        print(f"  SIN hit en ningun camelido    : {sin_nada:>5d}  ({100*sin_nada/n:.1f}%)")
    print()
    print("  LECTURA:")
    print("   con hit  -> gen real que la transferencia por homologia perdio.")
    print("               La novedad es un fallo de transferencia recuperado por ab initio.")
    print("   sin hit  -> candidato a especifico de linaje, o falso positivo de Helixer.")
    print("               Requiere busqueda contra bases publicas para decidir.")

    # --- 5. volcados
    with open(OUT_TSV, "w") as out:
        out.write("mrna_id\tscaffold\tstart\tend\tlength\tmd5\thit_dromedario\thit_alpaca\n")
        for f in filas:
            out.write("\t".join(str(x) for x in f) + "\n")
    print(f"\nEscrito {OUT_TSV}")

    seqs = read_fasta(HELIXER_FAA)
    written = 0
    with open(OUT_FASTA, "w") as out:
        for mid in sin_hit_ids:
            seq = seqs.get(mid)
            if seq is None:
                continue
            out.write(f">{mid}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
            written += 1
    print(f"Escrito {OUT_FASTA} con {written} secuencias sin hit en camelidos")
    if written != len(sin_hit_ids):
        print(f"  AVISO: {len(sin_hit_ids) - written} identificadores no se encontraron"
              f" en el FASTA. Revisar la correspondencia de nombres.")


if __name__ == "__main__":
    main()
