#!/usr/bin/env python3
"""
Prepara la entrada de InterProScan a partir de los cuatro proteomas.

QUE HACE, EN TRES PASOS
  1. DEDUPLICA por MD5 los cuatro proteomas en un unico FASTA, y escribe el mapa
     que permite devolver cada resultado a su proteoma de origen.
  2. LIMPIA los caracteres que InterProScan no admite.
  3. TROCEA en lotes para permitir reanudacion tras interrupciones.

POR QUE DEDUPLICAR
  Los cuatro proteomas suman muchas secuencias identicas entre si, porque tres de
  ellos proyectan la misma referencia. Calcular una sola vez cada secuencia unica
  ahorra en torno al 47% del computo. InterProScan no necesita calcular dos veces
  la misma secuencia.

EL MAPA ES IMPRESCINDIBLE
  md5_to_ids.tsv relaciona cada MD5 con el brazo y el identificador original.
  SIN ESE FICHERO LOS RESULTADOS FUNCIONALES NO SE PUEDEN ATRIBUIR a cada proteoma,
  porque el FASTA deduplicado usa el MD5 como unico identificador.

POR QUE LIMPIAR
  InterProScan aborta en el paso stepLoadFromFastaIntoDB si el FASTA contiene
  caracteres que no son aminoacidos. Los proteomas de homologia contienen el
  caracter "." (codon de stop) porque muchos modelos transferidos tienen la pauta
  de lectura rota; ver internal_stops.py para la cuantificacion.

  Se SUSTITUYE por X (residuo desconocido) en lugar de eliminar, para preservar la
  longitud de la secuencia y por tanto las coordenadas de los dominios que
  InterProScan devuelva. El stop final si se elimina.

USO
  python3 make_interproscan_input.py
"""

import gzip
import hashlib
import re

FILES = {
    "liftoff":  "llama_liftoff.faa.gz",
    "miniprot": "llama_miniprot.faa.gz",
    "lifton":   "llama_lifton.faa.gz",
    "helixer":  "Lgla_hx036_helixer.faa",
}
OUT_RAW = "camelid_unique_proteins.faa"
OUT_CLEAN = "camelid_unique_proteins_clean.faa"
MAP = "md5_to_ids.tsv"
CHUNK_DIR = "ips_chunks"
CHUNK_SIZE = 1000

VALID = "ACDEFGHIKLMNPQRSTVWYXBZJUO"
INVALID = re.compile(f"[^{VALID}]")
LINE_WIDTH = 60


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def write_seq(fh, header, seq):
    fh.write(header + "\n")
    for i in range(0, len(seq), LINE_WIDTH):
        fh.write(seq[i:i + LINE_WIDTH] + "\n")


def dedup():
    """Paso 1: deduplica y escribe el mapa."""
    seen = {}
    mapping = []
    for arm, path in FILES.items():
        header, seq = None, []

        def flush():
            if header is None:
                return
            s = "".join(seq)
            h = hashlib.md5(s.encode()).hexdigest()
            mapping.append((h, arm, header[1:].split()[0]))
            seen.setdefault(h, s)

        with opener(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    flush()
                    header, seq = line.strip(), []
                else:
                    seq.append(line.strip())
            flush()

    with open(OUT_RAW, "w") as out:
        for h, s in seen.items():
            write_seq(out, f">{h}", s)

    with open(MAP, "w") as m:
        m.write("md5\tarm\toriginal_id\n")
        for h, arm, oid in mapping:
            m.write(f"{h}\t{arm}\t{oid}\n")

    print(f"1. deduplicacion: {len(seen)} unicas de {len(mapping)} totales "
          f"(ahorro {100 * (1 - len(seen) / len(mapping)):.1f}%)")
    return len(seen)


def clean():
    """Paso 2: sustituye caracteres no validos por X."""
    n_seq = n_fix = 0
    with open(OUT_CLEAN, "w") as out:
        header, seq = None, []

        def flush():
            nonlocal n_seq, n_fix
            if header is None:
                return
            s = "".join(seq).upper().rstrip(".*")
            if INVALID.search(s):
                s = INVALID.sub("X", s)
                n_fix += 1
            n_seq += 1
            write_seq(out, header, s)

        with open(OUT_RAW) as fh:
            for line in fh:
                if line.startswith(">"):
                    flush()
                    header, seq = line.strip(), []
                else:
                    seq.append(line.strip())
            flush()

    print(f"2. limpieza: {n_seq} secuencias, {n_fix} con caracteres sustituidos por X")


def split():
    """Paso 3: trocea en lotes."""
    import os
    os.makedirs(CHUNK_DIR, exist_ok=True)
    buf, idx, count = [], 1, 0

    def write_chunk(b, i):
        with open(f"{CHUNK_DIR}/chunk_{i:03d}.faa", "w") as o:
            o.writelines(b)

    with open(OUT_CLEAN) as fh:
        for line in fh:
            if line.startswith(">"):
                if count == CHUNK_SIZE:
                    write_chunk(buf, idx)
                    idx += 1
                    buf, count = [], 0
                count += 1
            buf.append(line)
    if buf:
        write_chunk(buf, idx)
    print(f"3. troceado: {idx} lotes de hasta {CHUNK_SIZE} secuencias en {CHUNK_DIR}/")


if __name__ == "__main__":
    dedup()
    clean()
    split()
    print("\nListo. Ejecutar despues run_ips.sh")
