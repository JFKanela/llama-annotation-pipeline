#!/usr/bin/env python3
"""Seccion 2 de la Tabla 1: estructura de la anotacion, los cuatro brazos.

Sustituye a la seccion que genera build_report.py, que solo cubre Liftoff y
miniprot porque LiftOn y Helixer no pasaron por AGAT.

Todos los recuentos filtran los padres a mRNA. Ver structural_stats.py.
Criterio de monoexonico: una sola CDS.
"""
import gzip, os
from collections import defaultdict

ZEN = os.path.expanduser("~/zenodo_data")
BASE = os.path.expanduser("~/llama_annotation_pipeline")
OUT = os.path.expanduser(
    "~/Gdrive/Doctorado/llama_annotation_pipeline/MANUSCRITO/02_tablas")
os.makedirs(OUT, exist_ok=True)

GFF = {"Liftoff":  os.path.join(ZEN, "llama_liftoff.gff3.gz"),
       "miniprot": os.path.join(ZEN, "llama_miniprot.gff3.gz"),
       "LiftOn":   os.path.join(ZEN, "llama_lifton.gff3.gz"),
       "Helixer":  os.path.join(BASE, "Lgla_hx036_helixer_FINAL.gff")}
FAA = {"Liftoff":  os.path.join(ZEN, "llama_liftoff.faa.gz"),
       "miniprot": os.path.join(ZEN, "llama_miniprot.faa.gz"),
       "LiftOn":   os.path.join(ZEN, "llama_lifton.faa.gz"),
       "Helixer":  os.path.join(BASE, "Lgla_hx036_helixer.faa")}
ARMS = list(GFF)


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")


def stats(path):
    mrna, genes = set(), set()
    cds_n, cds_bp = defaultdict(int), defaultdict(int)
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            d = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
            if c[2] == "mRNA":
                if d.get("ID"):
                    mrna.add(d["ID"])
                if d.get("Parent"):
                    genes.add(d["Parent"])
            elif c[2] == "CDS" and d.get("Parent"):
                cds_n[d["Parent"]] += 1
                cds_bp[d["Parent"]] += int(c[4]) - int(c[3]) + 1
    n = {k: v for k, v in cds_n.items() if k in mrna}
    bp = {k: v for k, v in cds_bp.items() if k in mrna}
    v = list(n.values())
    return dict(genes=len(genes) or len(mrna), mrna=len(mrna), con_cds=len(v),
                cds_tot=sum(v), cds_mean=sum(v) / len(v) if v else 0,
                mono=100 * sum(1 for x in v if x == 1) / len(v) if v else 0,
                bp=sum(bp.values()))


def nprot(path):
    with opener(path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


filas = []
for a in ARMS:
    s = stats(GFF[a])
    s["prot"] = nprot(FAA[a])
    filas.append((a, s))
    print(f"{a:9s} genes={s['genes']:>6d} mRNA={s['mrna']:>6d} prot={s['prot']:>6d} "
          f"CDS={s['cds_tot']:>7d} CDS/mRNA={s['cds_mean']:5.2f} "
          f"mono={s['mono']:4.1f}% CDS={s['bp']/1e6:5.1f} Mb")

md = ["## 2 · Estructura de la anotación\n",
      "| Annotation | Genes | mRNA | Proteins | CDS features | CDS per mRNA | "
      "Single-CDS (%) | Total CDS (Mb) |",
      "|---|---|---|---|---|---|---|---|"]
for a, s in filas:
    md.append(f"| {a} | {s['genes']:,} | {s['mrna']:,} | {s['prot']:,} | "
              f"{s['cds_tot']:,} | {s['cds_mean']:.2f} | {s['mono']:.1f} | "
              f"{s['bp']/1e6:.1f} |".replace(",", " "))
md.append("\n> Todos los recuentos filtran los padres a mRNA, de modo que excluyen "
          "lnc_RNA, tRNA, snRNA, snoRNA y rRNA. Sin ese filtro los brazos de homología "
          "quedarían sobrerrepresentados frente a Helixer, que no emite features no "
          "codificantes. Monoexónico se define aquí como transcrito con una sola CDS.\n")

f = os.path.join(OUT, "tabla1_seccion2.md")
open(f, "w").write("\n".join(md))
print("\nescrito:", f)
print("\n".join(md))
