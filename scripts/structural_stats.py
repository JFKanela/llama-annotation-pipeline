#!/usr/bin/env python3
"""
Estadisticas estructurales de los cuatro proteomas y sus anotaciones.

QUE PRODUCE
  - Longitud de proteina: media, mediana y maximo por brazo.
  - Numero de CDS por transcrito (proxy de exones codificantes) y porcentaje
    de transcritos monoexonicos.
  - Numero de mRNA por gen, que detecta asimetrias de isoformas entre brazos.

QUE SOSTIENE EN EL ARTICULO
  Eje de calidad INDEPENDIENTE de BUSCO, siguiendo el criterio de Kourelis et al.
  2019 (BMC Genomics 20:722), que es el modelo declarado del articulo. Los proteomas
  con modelos truncados tienen proteinas mas cortas; un exceso de genes monoexonicos
  suele indicar predicciones incompletas o pseudogenes.

HALLAZGO QUE ESTE SCRIPT DESTAPO
  LiftOn da 1.00 mRNA/gen en el fichero depositado, igual que los demas. Es decir,
  el proteoma publico es la version de ISOFORMA PRIMARIA y no la completa. Los cuatro
  brazos estan por tanto en igualdad y las comparaciones de longitud y estructura
  son directas, sin necesidad de normalizar.

CAUTELA AL INTERPRETAR
  "Proteinas mas largas" no es inequivocamente mejor: puede indicar mejor modelado
  o sobre-extension (union de exones que no van juntos, o inclusion de region no
  codificante). Distinguirlo exige cobertura BLASTP contra la referencia de alpaca.
  Lo defendible sin ese paso es la descripcion, no el juicio de calidad.

  miniprot no emite jerarquia gene -> mRNA, de modo que su fila de isoformas sale
  vacia. Es correcto, no es un fallo.

USO
  python3 structural_stats.py
"""

import gzip
import statistics as st
from collections import defaultdict

GFF = {
    "liftoff":  "llama_liftoff.gff3.gz",
    "miniprot": "llama_miniprot.gff3.gz",
    "lifton":   "llama_lifton.gff3.gz",
    "helixer":  "Lgla_hx036_helixer_FINAL.gff",
}
FAA = {
    "liftoff":  "llama_liftoff.faa.gz",
    "miniprot": "llama_miniprot.faa.gz",
    "lifton":   "llama_lifton.faa.gz",
    "helixer":  "Lgla_hx036_helixer.faa",
}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_attrs(attr):
    return dict(kv.split("=", 1) for kv in attr.split(";") if "=" in kv)


def parse_gff(path):
    """Devuelve (CDS por transcrito, longitud CDS por transcrito, mRNA por gen)."""
    cds_count = defaultdict(int)
    cds_len = defaultdict(int)
    mrna_of_gene = defaultdict(set)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            typ, start, end, attr = c[2], int(c[3]), int(c[4]), c[8]
            if typ == "mRNA":
                d = parse_attrs(attr)
                if d.get("Parent"):
                    mrna_of_gene[d["Parent"]].add(d.get("ID", ""))
            elif typ == "CDS":
                d = parse_attrs(attr)
                par = d.get("Parent")
                if par:
                    cds_count[par] += 1
                    cds_len[par] += (end - start + 1)
    return cds_count, cds_len, mrna_of_gene


def prot_lengths(path):
    lens = []
    cur = 0
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    lens.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
        if cur:
            lens.append(cur)
    return lens


def main():
    print("--- longitud de proteina (aa) ---")
    print(f"{'brazo':10s} {'#prot':>7s} {'media':>9s} {'mediana':>8s} {'maximo':>8s}")
    for arm, path in FAA.items():
        v = prot_lengths(path)
        print(f"{arm:10s} {len(v):>7d} {st.mean(v):>9.1f} {int(st.median(v)):>8d} {max(v):>8d}")

    print("\n--- estructura por transcrito (numero de CDS = exones codificantes) ---")
    print(f"{'brazo':10s} {'#tx':>7s} {'CDS/tx media':>13s} {'CDS/tx med':>11s} {'monoexon%':>10s}")
    for arm, path in GFF.items():
        cc, _, _ = parse_gff(path)
        counts = list(cc.values())
        if not counts:
            print(f"{arm:10s}  (sin CDS con Parent legible)")
            continue
        mono = 100 * sum(1 for x in counts if x == 1) / len(counts)
        print(f"{arm:10s} {len(counts):>7d} {st.mean(counts):>13.2f} "
              f"{int(st.median(counts)):>11d} {mono:>9.1f}%")

    print("\n--- isoformas: mRNA por gen ---")
    print(f"{'brazo':10s} {'#genes':>7s} {'mRNA/gen':>9s}")
    for arm, path in GFF.items():
        _, _, mg = parse_gff(path)
        if not mg:
            print(f"{arm:10s}  (sin jerarquia gene -> mRNA; esperable en miniprot)")
            continue
        total = sum(len(v) for v in mg.values())
        print(f"{arm:10s} {len(mg):>7d} {total / len(mg):>9.2f}")


if __name__ == "__main__":
    main()
