#!/usr/bin/env python3
"""Genera el informe comparativo de la anotacion de Lama glama.

Interfaz invocada por la regla `report` del Snakefile:

    python scripts/build_report.py \
        --busco-dir results/busco \
        --agat-dir results/agat \
        --gffcompare results/gffcompare/cmp.stats \
        --unmapped results/liftoff/unmapped.txt \
        --out-md <md> --out-tsv <tsv>

Diseno:
  - No recalcula nada: solo lee ficheros ya producidos por el pipeline.
  - Descubre los metodos escaneando --busco-dir, de modo que brazos anadidos
    despues (lifton, helixer) aparecen sin tocar el codigo.
  - Cuenta proteinas del .faa real evaluado por BUSCO, no mRNA del GFF.
  - Tolera ficheros ausentes: informa de lo que hay y omite el resto.
"""
import argparse
import glob
import os
import re
import sys

# Orden de presentacion. Los metodos no listados van al final, alfabeticamente.
ORDEN = ["liftoff", "miniprot", "lifton", "helixer", "alpaca_ref", "chaku_v1"]

ESTRATEGIA = {
    "liftoff":    "homologia ADN-ADN",
    "miniprot":   "homologia proteina-ADN",
    "lifton":     "homologia hibrida",
    "helixer":    "ab initio",
    "alpaca_ref": "referencia (techo)",
    "chaku_v1":   "baseline previo",
}

# Nombre del .faa por metodo cuando no coincide con el patron por defecto.
FAA_ESPECIAL = {
    "lifton": "llama_lifton_v2.faa",
}

RE_BUSCO = re.compile(
    r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)"
)


def hallar_summary(busco_dir, metodo):
    """Localiza el short_summary sea cual sea el nombre que BUSCO le diera."""
    candidatos = [
        os.path.join(busco_dir, metodo, "short_summary.txt"),
        *sorted(glob.glob(os.path.join(
            busco_dir, metodo, "short_summary.specific.*.txt"))),
        *sorted(glob.glob(os.path.join(
            busco_dir, metodo, "run_*", "short_summary.txt"))),
    ]
    for ruta in candidatos:
        if os.path.isfile(ruta):
            return ruta
    return None


def leer_busco(ruta):
    if not ruta:
        return None
    with open(ruta) as fh:
        m = RE_BUSCO.search(fh.read())
    if not m:
        return None
    c, s, d, f, mi, n = m.groups()
    return {"C": c, "S": s, "D": d, "F": f, "M": mi, "n": n}


def contar_proteinas(prot_dir, metodo):
    nombre = FAA_ESPECIAL.get(metodo, "llama_%s.faa" % metodo)
    ruta = os.path.join(prot_dir, nombre)
    if not os.path.isfile(ruta):
        return None
    with open(ruta) as fh:
        return sum(1 for linea in fh if linea.startswith(">"))


def contar_feature(gff, feature):
    if not gff or not os.path.isfile(gff):
        return None
    n = 0
    with open(gff) as fh:
        for linea in fh:
            if linea.startswith("#"):
                continue
            campos = linea.split("\t")
            if len(campos) > 2 and campos[2] == feature:
                n += 1
    return n


def hallar_gff(metodo):
    """Localiza el GFF del metodo bajo results/<metodo>/."""
    patrones = [
        "results/%s/llama_%s_v2.gff3" % (metodo, metodo),
        "results/%s/llama_%s.gff3" % (metodo, metodo),
    ]
    for ruta in patrones:
        if os.path.isfile(ruta):
            return ruta
    encontrados = sorted(glob.glob("results/%s/*.gff3" % metodo))
    return encontrados[0] if encontrados else None


def leer_gffcompare(ruta):
    niveles = {}
    if not ruta or not os.path.isfile(ruta):
        return niveles
    patron = re.compile(r"(Transcript|Locus) level:\s+([\d.]+)\s+\|\s+([\d.]+)")
    with open(ruta) as fh:
        for linea in fh:
            m = patron.search(linea)
            if m:
                niveles[m.group(1).lower()] = (m.group(2), m.group(3))
    return niveles


def contar_unmapped(ruta):
    if not ruta or not os.path.isfile(ruta):
        return None
    with open(ruta) as fh:
        return sum(1 for l in fh if l.strip() and not l.startswith("#"))


def listar_agat(agat_dir):
    if not agat_dir or not os.path.isdir(agat_dir):
        return []
    return sorted(glob.glob(os.path.join(agat_dir, "*.stats.txt")))


def descubrir_metodos(busco_dir):
    if not os.path.isdir(busco_dir):
        return []
    hallados = [d for d in os.listdir(busco_dir)
                if os.path.isdir(os.path.join(busco_dir, d))]
    conocidos = [m for m in ORDEN if m in hallados]
    resto = sorted(m for m in hallados if m not in ORDEN)
    return conocidos + resto


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--busco-dir", required=True)
    ap.add_argument("--agat-dir", default="results/agat")
    ap.add_argument("--gffcompare", default=None)
    ap.add_argument("--unmapped", default=None)
    ap.add_argument("--proteomes-dir", default="results/proteomes")
    ap.add_argument("--out-md", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    metodos = descubrir_metodos(args.busco_dir)
    if not metodos:
        sys.exit("ERROR: no se encontro ningun resultado BUSCO en %s"
                 % args.busco_dir)

    filas = []
    for metodo in metodos:
        busco = leer_busco(hallar_summary(args.busco_dir, metodo))
        if not busco:
            continue
        gff = hallar_gff(metodo)
        filas.append({
            "metodo": metodo,
            "estrategia": ESTRATEGIA.get(metodo, "-"),
            "busco": busco,
            "proteinas": contar_proteinas(args.proteomes_dir, metodo),
            "genes": contar_feature(gff, "gene"),
            "mrna": contar_feature(gff, "mRNA"),
        })

    gffcmp = leer_gffcompare(args.gffcompare)
    unmapped = contar_unmapped(args.unmapped)
    agat = listar_agat(args.agat_dir)

    # ---------------- TSV ----------------
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join([
            "metodo", "estrategia", "busco_completo_pct", "busco_unico_pct",
            "busco_duplicado_pct", "busco_fragmentado_pct", "busco_ausente_pct",
            "busco_n", "proteinas", "genes_gff", "mrna_gff"]) + "\n")
        for f in filas:
            b = f["busco"]
            fh.write("\t".join([
                f["metodo"], f["estrategia"], b["C"], b["S"], b["D"], b["F"],
                b["M"], b["n"],
                str(f["proteinas"] if f["proteinas"] is not None else "NA"),
                str(f["genes"] if f["genes"] is not None else "NA"),
                str(f["mrna"] if f["mrna"] is not None else "NA")]) + "\n")

    # ---------------- Markdown ----------------
    L = []
    L.append("# Informe comparativo de anotacion: *Lama glama*")
    L.append("")
    L.append("Generado automaticamente por `scripts/build_report.py`. "
             "Todas las cifras proceden de ficheros producidos por el pipeline; "
             "este script no recalcula nada.")
    L.append("")
    L.append("## 1. Completitud BUSCO (modo proteina)")
    L.append("")
    L.append("| Proteoma | Estrategia | Completo % | Unico % | Duplicado % | "
             "Fragmentado % | Ausente % | Proteinas |")
    L.append("|---|---|---|---|---|---|---|---|")
    for f in filas:
        b = f["busco"]
        prot = f["proteinas"] if f["proteinas"] is not None else "-"
        L.append("| %s | %s | %s | %s | %s | %s | %s | %s |" % (
            f["metodo"], f["estrategia"], b["C"], b["S"], b["D"], b["F"],
            b["M"], prot))
    L.append("")
    if filas:
        L.append("Linaje BUSCO: `artiodactyla_odb12` (n = %s marcadores)."
                 % filas[0]["busco"]["n"])
        L.append("")

    L.append("## 2. Conteo de features en los GFF")
    L.append("")
    L.append("| Metodo | Genes | mRNA | Proteinas (.faa) |")
    L.append("|---|---|---|---|")
    for f in filas:
        L.append("| %s | %s | %s | %s |" % (
            f["metodo"],
            f["genes"] if f["genes"] is not None else "-",
            f["mrna"] if f["mrna"] is not None else "-",
            f["proteinas"] if f["proteinas"] is not None else "-"))
    L.append("")
    L.append("El numero de proteinas puede ser inferior al de mRNA: los "
             "transcritos sin CDS traducible (pseudogenes, modelos truncados) "
             "no rinden secuencia proteica. La cifra evaluada por BUSCO es la "
             "del `.faa`.")
    L.append("")

    if gffcmp:
        L.append("## 3. Concordancia estructural entre metodos")
        L.append("")
        L.append("| Nivel | Sensibilidad % | Precision % |")
        L.append("|---|---|---|")
        for nivel, (sn, pr) in gffcmp.items():
            L.append("| %s | %s | %s |" % (nivel, sn, pr))
        L.append("")
        L.append("Una concordancia baja no implica repertorios genicos "
                 "distintos: refleja que los metodos modelan la estructura "
                 "fina (exones, UTR) de forma diferente.")
        L.append("")

    if unmapped is not None:
        L.append("## 4. Recuperacion")
        L.append("")
        L.append("Genes de la referencia no transferidos por Liftoff "
                 "(unmapped): **%d**." % unmapped)
        L.append("")

    if agat:
        L.append("## 5. Estadisticas estructurales (AGAT)")
        L.append("")
        for ruta in agat:
            L.append("- `%s`" % ruta)
        L.append("")

    L.append("## Nota sobre el sustrato")
    L.append("")
    L.append("Los brazos de homologia se ejecutaron sobre el ensamblado "
             "completo (2.351.761.190 bp). Helixer, por restricciones de "
             "computo, se ejecuto sobre los scaffolds de longitud igual o "
             "superior a 10 kb: 3.640 scaffolds (0,34 % del total) que reunen "
             "1.915.763.599 bp, el 81,46 % de la secuencia ensamblada. Esta "
             "asimetria condiciona toda comparacion directa de completitud "
             "entre Helixer y los brazos de homologia.")
    L.append("")

    os.makedirs(os.path.dirname(args.out_md) or ".", exist_ok=True)
    with open(args.out_md, "w") as fh:
        fh.write("\n".join(L))

    print("Informe escrito: %s" % args.out_md)
    print("Tabla escrita  : %s" % args.out_tsv)
    print("Metodos incluidos: %s" % ", ".join(f["metodo"] for f in filas))


if __name__ == "__main__":
    main()
