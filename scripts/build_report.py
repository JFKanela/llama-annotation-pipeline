#!/usr/bin/env python3
"""Agrega BUSCO + AGAT + gffcompare en una tabla y un informe markdown.

Robusto a ficheros ausentes: lo que falte se marca como 'NA', no rompe.

PRECAUCION -- LEER ANTES DE REGENERAR EL INFORME.
(Repuesto en la ronda 7; se comprobo sobre el codigo en la ronda 5 y sigue
siendo cierto tras la reescritura del 14 de agosto de 2026.)

  1. Las cifras de BUSCO NO estan codificadas en el script: se leen de los
     ficheros de --busco-dir. Precisamente por eso, si results/busco/ sigue
     conteniendo la salida de la primera ejecucion de Helixer (helixerlite,
     scaffolds >= 10 kb, sin solapamiento), la tabla se regenerara con el
     79,3 % antiguo en lugar del 85,5 % vigente. Antes de volver a ejecutar
     este script hay que sustituir esa carpeta por la salida BUSCO de la
     anotacion vigente (web tool de Helixer v0.3.6, scaffolds >= 25 kb).

  2. La nota final sobre el sustrato SI lleva valores codificados en el
     cuerpo del script (funcion main, seccion "Nota sobre el sustrato").
     Si el sustrato vuelve a cambiar, hay que editarla a mano.

Invocado por la regla `report` del Snakefile:

    python scripts/build_report.py \
        --busco-dir results/busco \
        --agat-dir results/agat \
        --gffcompare results/gffcompare/cmp.stats \
        --unmapped results/liftoff/unmapped.txt \
        --out-md <md> --out-tsv <tsv>
"""
import argparse
import glob
import os
import re


def parse_busco(busco_dir):
    """Devuelve {nombre_proteoma: {C,S,D,F,M,n}} leyendo short_summary.txt."""
    out = {}
    # Algunas carpetas de BUSCO no contienen short_summary.txt sin sufijo,
    # solo short_summary.specific.<linaje>.<nombre>.txt. Buscando unicamente el
    # primero quedaban fuera lifton y helixer, es decir dos de los cuatro brazos.
    rutas = glob.glob(os.path.join(busco_dir, "*", "short_summary.txt"))
    vistos = {os.path.dirname(x) for x in rutas}
    for extra in glob.glob(os.path.join(busco_dir, "*", "short_summary.specific.*.txt")):
        if os.path.dirname(extra) not in vistos:
            rutas.append(extra)
            vistos.add(os.path.dirname(extra))
    for path in rutas:
        name = os.path.basename(os.path.dirname(path))
        txt = open(path).read()
        m = re.search(
            r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)",
            txt,
        )
        if m:
            out[name] = dict(
                zip(["C", "S", "D", "F", "M", "n"], m.groups())
            )
    return out


def parse_agat(agat_dir):
    """Extrae recuentos clave del informe de AGAT."""
    out = {}
    for path in glob.glob(os.path.join(agat_dir, "*.stats.txt")):
        method = os.path.basename(path).split(".")[0]
        # AGAT emite un bloque por cada tipo de feature (c_gene_segment,
        # lnc_rna, mrna, trna, v_gene_segment...). Solo interesa el bloque
        # "mrna", que son los genes codificantes de proteina. Sin este filtro
        # el bucle se quedaba con el ULTIMO bloque del fichero, v_gene_segment,
        # y reportaba 31 genes para Liftoff en lugar de 20306.
        d = {}
        bloque = None
        for line in open(path):
            line = line.strip()
            b = re.match(r"^-+\s+(\w+)\s+-+$", line)
            if b:
                bloque = b.group(1).lower()
                continue
            if bloque != "mrna":
                continue
            for key in ("gene", "mrna", "cds", "exon"):
                m = re.match(rf"(?i)number of\s+{key}\b\s+(\d+)", line)
                if m:
                    d[key] = int(m.group(1))
        out[method] = d
    return out


def parse_gffcompare(path):
    """Sensibilidad/precision a nivel de transcrito y locus."""
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        m = re.search(r"(Transcript|Locus) level:\s+([\d.]+)\s+\|\s+([\d.]+)", line)
        if m:
            out[m.group(1).lower()] = {"sensitivity": m.group(2), "precision": m.group(3)}
    return out


def count_lines(path):
    if not os.path.exists(path):
        return "NA"
    return sum(1 for ln in open(path) if ln.strip() and not ln.startswith("#"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--busco-dir", required=True)
    ap.add_argument("--agat-dir", required=True)
    ap.add_argument("--gffcompare", required=True)
    ap.add_argument("--unmapped", required=True)
    ap.add_argument("--out-md", required=True)
    ap.add_argument("--out-tsv", required=True)
    a = ap.parse_args()

    busco = parse_busco(a.busco_dir)
    agat = parse_agat(a.agat_dir)
    gffc = parse_gffcompare(a.gffcompare)
    unmapped_n = count_lines(a.unmapped)

    # ---- TSV ----
    order = ["liftoff", "miniprot", "alpaca_ref", "chaku_v1"]
    names = [n for n in order if n in busco] + [n for n in busco if n not in order]
    with open(a.out_tsv, "w") as f:
        f.write("proteome\tBUSCO_complete_pct\tsingle_pct\tdup_pct\tfragmented_pct\t"
                "missing_pct\tn_markers\tgenes\tmRNAs\n")
        for n in names:
            b = busco.get(n, {})
            g = agat.get(n, {})
            f.write("\t".join(str(x) for x in [
                n, b.get("C", "NA"), b.get("S", "NA"), b.get("D", "NA"),
                b.get("F", "NA"), b.get("M", "NA"), b.get("n", "NA"),
                g.get("gene", "NA"), g.get("mrna", "NA"),
            ]) + "\n")

    # ---- Markdown ----
    L = []
    L.append("# Informe comparativo · Anotacion de proteoma de *Lama glama*\n")
    L.append("Cuatro estrategias de anotacion sobre el mismo genoma diana (DNA Zoo `GCA_028534125.1`), "
             "misma referencia (*Vicugna pacos* `GCF_048564905.1`).\n")

    L.append("\n## 1 · Completitud BUSCO (modo proteina)\n")
    L.append("| Proteoma | Completo % | Unico % | Duplicado % | Fragmentado % | Ausente % | n |")
    L.append("|---|---|---|---|---|---|---|")
    for n in names:
        b = busco.get(n, {})
        L.append(f"| {n} | {b.get('C','NA')} | {b.get('S','NA')} | {b.get('D','NA')} | "
                 f"{b.get('F','NA')} | {b.get('M','NA')} | {b.get('n','NA')} |")

    L.append("\n## 2 · Estructura de la anotacion (AGAT)\n")
    L.append("| Anotacion | Genes | mRNAs |")
    L.append("|---|---|---|")
    for method in ("liftoff", "miniprot"):
        g = agat.get(method, {})
        L.append(f"| {method} | {g.get('gene','NA')} | {g.get('mrna','NA')} |")

    L.append("\n## 3 · Concordancia A (Liftoff) vs B (miniprot) · gffcompare\n")
    L.append("Referencia = anotacion Liftoff; consulta = anotacion miniprot.\n")
    L.append("| Nivel | Sensibilidad % | Precision % |")
    L.append("|---|---|---|")
    for lvl in ("transcript", "locus"):
        d = gffc.get(lvl, {})
        L.append(f"| {lvl} | {d.get('sensitivity','NA')} | {d.get('precision','NA')} |")

    L.append("\n## 4 · Recuperacion Liftoff\n")
    L.append(f"Genes de la referencia NO transferidos por Liftoff (unmapped): **{unmapped_n}**\n")

    L.append("\n## Nota sobre el sustrato\n")
    L.append("Los brazos de homologia se ejecutaron sobre el ensamblado completo "
             "(2.351.761.190 bp). Helixer, por restricciones de computo, se ejecuto "
             "sobre los scaffolds de longitud igual o superior a 25 kb: 244 scaffolds, "
             "el 79,6 % de la secuencia ensamblada. Esta asimetria condiciona toda "
             "comparacion directa de completitud entre Helixer y los brazos de "
             "homologia. (Sustrato de la ejecucion vigente con la web tool de Helixer "
             "v0.3.6; la primera ejecucion, con helixerlite y sin solapamiento, uso "
             "scaffolds >= 10 kb y esta superada.)\n")

    L.append("\n---\n")
    L.append("> Recordatorio de encuadre: los tres brazos de homologia son "
         "PROYECCIONES de la referencia de alpaca y no pueden descubrir genes "
         "ausentes en ella; solo el brazo ab initio (Helixer) es independiente "
         "de toda referencia. El valor del recurso es disponer de un proteoma de "
         "llama anotado y evaluado por cuatro estrategias, no el descubrimiento "
         "de genes nuevos. NO usar la palabra 'primer': existe una anotacion "
         "interna previa (chaku_v1) que se emplea como linea de base.\n")

    open(a.out_md, "w").write("\n".join(L) + "\n")
    print(f"Informe escrito: {a.out_md} y {a.out_tsv}")


if __name__ == "__main__":
    main()
