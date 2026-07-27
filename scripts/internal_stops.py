#!/usr/bin/env python3
"""
Recuento de codones de stop internos por brazo de anotacion.

QUE PRODUCE
  Por cada proteoma: cuantas secuencias contienen algun caracter de stop, cuantas
  lo tienen solo al final (normal) y cuantas lo tienen INTERNO (pauta de lectura rota).

QUE SOSTIENE EN EL ARTICULO
  Es la evidencia mas solida de que LiftOn es el mejor proteoma de referencia, y la
  unica que NO depende de la asimetria de isoformas que contaminaba las cifras de
  BUSCO. LiftOn reduce los stops internos 4.4 veces respecto a Liftoff, que es
  exactamente su funcion de diseno: rescatar modelos de Liftoff mediante miniprot.

  Concuerda con lo que el propio GFF de Liftoff declara en sus atributos y que no se
  habia aprovechado: valid_ORF=False, missing_start_codon=True, matches_ref_protein=False.

CAUTELA CRITICA SOBRE HELIXER
  Helixer da 0.0% y eso NO es un merito. Es consecuencia estructural del metodo: un
  predictor ab initio construye ORFs por definicion y no puede producir stops internos.
  Es imposible que de otra cifra. Presentarlo como calidad seria un error detectable
  de inmediato en revision. La fila de Helixer aqui es NO INFORMATIVA. Solo los tres
  brazos de homologia son comparables entre si.

ORIGEN
  Este analisis surgio al depurar un fallo de InterProScan, que rechazaba el FASTA
  combinado en el paso stepLoadFromFastaIntoDB por contener el caracter ".".

DECISION DE DISENO PENDIENTE
  Si el entregable es un proteoma de referencia, cabe preguntarse si deben incluirse
  los modelos con pauta rota, que argumentablemente no son secuencias codificantes de
  proteina. Kourelis et al. 2019 filtro por cobertura BLASTP contra referencia (umbral
  del 60%) precisamente para eliminarlos. Opciones: filtrar, o publicar sin filtrar y
  declarar la cifra. Consultar con el director.

USO
  python3 internal_stops.py
"""

import gzip

FILES = {
    "liftoff":  "llama_liftoff.faa.gz",
    "miniprot": "llama_miniprot.faa.gz",
    "lifton":   "llama_lifton.faa.gz",
    "helixer":  "Lgla_hx036_helixer.faa",
}
STOP_CHARS = ".*"


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def count_stops(path):
    total = with_stop = only_final = internal = 0
    seq = []
    started = False

    def flush():
        nonlocal total, with_stop, only_final, internal
        if not started:
            return
        s = "".join(seq)
        total += 1
        if any(ch in s for ch in STOP_CHARS):
            with_stop += 1
            core = s.rstrip(STOP_CHARS)
            if any(ch in core for ch in STOP_CHARS):
                internal += 1
            else:
                only_final += 1

    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                started = True
                seq = []
            else:
                seq.append(line.strip())
        flush()

    return total, with_stop, only_final, internal


def main():
    print(f"{'brazo':10s} {'total':>7s} {'con stop':>9s} {'%':>6s} "
          f"{'solo final':>11s} {'INTERNOS':>9s}")
    for arm, path in FILES.items():
        total, with_stop, only_final, internal = count_stops(path)
        pct = 100 * with_stop / total if total else 0
        print(f"{arm:10s} {total:>7d} {with_stop:>9d} {pct:>5.1f}% "
              f"{only_final:>11d} {internal:>9d}")

    print("\nRecordatorio: la fila de helixer es no informativa. Un predictor ab initio")
    print("construye ORFs por definicion y no puede producir stops internos.")


if __name__ == "__main__":
    main()
