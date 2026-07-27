#!/usr/bin/env python3
"""
Solapamiento por COORDENADAS a nivel de locus entre Helixer y los brazos de homologia.

QUE PRODUCE
  - Loci de Helixer que no solapan ningun locus de la homologia (novedad posicional).
  - Loci de la homologia que no solapan ningun locus de Helixer (deficit de Helixer).
  - Distribucion de longitudes de los loci novel, para separar senal de artefacto.

QUE SOSTIENE EN EL ARTICULO
  Es la unica cifra de novedad defendible. Sustituye a la del analisis por MD5,
  que estaba inflada un orden de magnitud por medir identidad exacta de secuencia
  en lugar de posicion genomica.

UNIDAD DE COMPARACION: EL LOCUS, NO LA CDS
  El locus se define como el intervalo genomico completo del mRNA, de inicio a fin,
  intrones incluidos. Un primer intento que fusionaba CDS conto EXONES (mas de
  200.000 por metodo, cuando un genoma de mamifero tiene del orden de 20.000 genes)
  y dio una cifra falsa. Se usa mRNA y no gene porque miniprot no emite feature gene.

  Comprobado que los cuatro GFF3 comparten la nomenclatura de scaffold, de modo que
  los intervalos son cruzables. Si se cambia alguna anotacion de origen, verificarlo
  de nuevo: si los identificadores de la columna 1 no casan, el resultado seria cero
  solapamiento sin dar ningun error.

LIMITE DE LA AFIRMACION
  "Novel por posicion" no equivale a "gen especifico de llama". Un locus que Helixer
  ve y la homologia no puede ser un gen real ausente en la anotacion de alpaca, un
  gen propio de llama, o un falso positivo. Distinguirlos exige BLAST contra bases
  publicas. Lo defendible sin ese paso es: "N loci codificantes predichos por Helixer
  sin solape posicional con la anotacion por homologia".

USO
  python3 overlap_coords.py
"""

import gzip
import statistics
from collections import defaultdict

FILES = {
    "liftoff":  "llama_liftoff.gff3.gz",
    "miniprot": "llama_miniprot.gff3.gz",
    "lifton":   "llama_lifton.gff3.gz",
    "helixer":  "Lgla_hx036_helixer_FINAL.gff",
}
HOMOLOGY_ARMS = ["liftoff", "miniprot", "lifton"]


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def load_mrna_intervals(path):
    """Un intervalo (start, end) por cada feature mRNA, agrupado por scaffold.

    Coordenadas GFF (1-based, inclusivas) convertidas a media abierta 0-based.
    """
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
    """Fusiona intervalos solapados o contiguos."""
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


def novel_loci(a_merged, b_merged):
    """Loci de a que no solapan NINGUN locus de b."""
    out = []
    for scf, aiv in a_merged.items():
        biv = b_merged.get(scf, [])
        for s, e in aiv:
            if not any(not (be <= s or bs >= e) for bs, be in biv):
                out.append((scf, s, e))
    return out


def main():
    raw = {arm: load_mrna_intervals(path) for arm, path in FILES.items()}

    print("--- numero de loci (mRNA) por metodo, sin fusionar ---")
    for arm in FILES:
        n = sum(len(v) for v in raw[arm].values())
        print(f"  {arm:10s} {n:>7d} loci")

    homology_raw = defaultdict(list)
    for arm in HOMOLOGY_ARMS:
        for scf, iv in raw[arm].items():
            homology_raw[scf].extend(iv)
    homology = {scf: merge(iv) for scf, iv in homology_raw.items()}
    n_hom = sum(len(v) for v in homology.values())

    helixer = {scf: merge(iv) for scf, iv in raw["helixer"].items()}
    n_hx = sum(len(v) for v in helixer.values())

    print(f"\n  union homologia fusionada: {n_hom} loci")
    print(f"  helixer fusionado:         {n_hx} loci")

    hx_novel = novel_loci(helixer, homology)
    ho_novel = novel_loci(homology, helixer)

    print("\n=== RESULTADO ===")
    print(f"Loci de Helixer que NO solapan homologia: {len(hx_novel)}"
          f"  ({100 * len(hx_novel) / n_hx:.1f}%)")
    print(f"Loci de homologia que NO solapan Helixer: {len(ho_novel)}"
          f"  ({100 * len(ho_novel) / n_hom:.1f}%)")

    if hx_novel:
        lens = sorted(e - s for _, s, e in hx_novel)
        print("\nLongitud de los loci novel de Helixer (bp):")
        print(f"  minimo {lens[0]}, mediana {int(statistics.median(lens))}, maximo {lens[-1]}")
        print(f"  menores de 300 bp (candidatos a artefacto): {sum(1 for x in lens if x < 300)}")
        print(f"  mayores de 1000 bp (candidatos a gen real): {sum(1 for x in lens if x >= 1000)}")

    # Volcado de los loci novel, por si se quiere hacer BLAST despues
    with open("helixer_novel_loci.bed", "w") as out:
        for scf, s, e in sorted(hx_novel):
            out.write(f"{scf}\t{s}\t{e}\n")
    print("\nEscrito helixer_novel_loci.bed con los loci novel de Helixer.")


if __name__ == "__main__":
    main()
