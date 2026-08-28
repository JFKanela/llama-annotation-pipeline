#!/usr/bin/env python3
"""Figura 3. Propiedades estructurales de los modelos de cada brazo.

Datos: los cuatro proteomas y sus GFF3.
Paneles: A numero de CDS por transcrito, B transcritos con una sola CDS.

QUE DEBE VERSE: Helixer produce modelos mas largos, con mas exones codificantes y
la mitad de transcritos con una sola CDS. Un transcrito con una sola CDS puede ser
una prediccion incompleta, un gen monoexonico real o un pseudogen: la diferencia es
marcada pero no admite una lectura unica.

EL PIE DE FIGURA DEBE REMITIR A LA FIGURA 4: esas proteinas mas largas NO son
sobre-extension. Medido contra C. dromedarius, que es el patron externo, los cuatro
brazos son indistinguibles. La explicacion es la preseleccion: Helixer solo emite
modelos donde consigue construir un ORF completo.

ANCHO. A diferencia de las figuras 1, 2 y 4, esta conserva bbox_inches="tight":
con figsize=(6.8, 3.4) el recorte deja 169,8 mm, que es el ancho de pagina
completa de la revista. Si se cambia el contenido hay que volver a medir el PDF,
porque el recorte depende de lo que se dibuje.
"""
import gzip, os
import statistics as st
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Rutas relativas al repositorio: este script vive en scripts/
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(REPO, "results")
OUT = os.path.join(RESULTS, "figures")
os.makedirs(OUT, exist_ok=True)

# Proteomas y anotaciones producidos por el flujo de trabajo
ZEN = os.path.join(RESULTS, "proteomes")
BASE = RESULTS

ARMS = ["Liftoff", "miniprot", "LiftOn", "Helixer"]
FAA = {"Liftoff":  os.path.join(ZEN, "llama_liftoff.faa.gz"),
       "miniprot": os.path.join(ZEN, "llama_miniprot.faa.gz"),
       "LiftOn":   os.path.join(ZEN, "llama_lifton.faa.gz"),
       "Helixer":  os.path.join(BASE, "Lgla_hx036_helixer.faa")}
GFF = {"Liftoff":  os.path.join(ZEN, "llama_liftoff.gff3.gz"),
       "miniprot": os.path.join(ZEN, "llama_miniprot.gff3.gz"),
       "LiftOn":   os.path.join(ZEN, "llama_lifton.gff3.gz"),
       "Helixer":  os.path.join(BASE, "Lgla_hx036_helixer_FINAL.gff")}
# color fijo por brazo, coherente con el resto de figuras
COL = {"Liftoff": "#0072B2", "miniprot": "#009E73",
       "LiftOn": "#CC79A7", "Helixer": "#E69F00"}

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

def longitudes(p):
    out, cur = [], 0
    with opener(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur: out.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
        if cur: out.append(cur)
    return out

def cds_por_tx(p):
    """CDS por transcrito, FILTRANDO los padres a mRNA.

    Sin ese filtro se cuelan padres que no son mRNA (lnc_RNA, tRNA, snRNA...),
    presentes solo en los brazos de homologia y no en Helixer, lo que sesga la
    comparacion. Ver la cabecera de structural_stats.py.
    """
    mrna, c = set(), defaultdict(int)
    with opener(p) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            d = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "mRNA" and d.get("ID"): mrna.add(d["ID"])
            if f[2] == "CDS" and d.get("Parent"): c[d["Parent"]] += 1
    return [v for k, v in c.items() if k in mrna]

L = {a: longitudes(FAA[a]) for a in ARMS}
C = {a: cds_por_tx(GFF[a]) for a in ARMS}
for a in ARMS:
    mono = 100 * sum(1 for x in C[a] if x == 1) / len(C[a])
    print(f"{a:9s} n={len(L[a]):>6d}  long media={st.mean(L[a]):6.1f}  "
          f"CDS/tx={st.mean(C[a]):5.2f}  monoexon={mono:4.1f}%")

fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.4))

# B. CDS por transcrito
ax = axes[0]
bp = ax.boxplot([C[a] for a in ARMS], patch_artist=True, showfliers=False, whis=(5, 95),
                widths=0.55, medianprops=dict(color="black", lw=1.3))
for patch, a in zip(bp["boxes"], ARMS):
    patch.set_facecolor(COL[a]); patch.set_alpha(0.85); patch.set_edgecolor("black")
ax.set_ylabel("Coding exons per transcript", fontsize=9)
ax.set_title("A", loc="left", fontsize=11, fontweight="bold")

# B. monoexonicos
ax = axes[1]
vals = [100 * sum(1 for x in C[a] if x == 1) / len(C[a]) for a in ARMS]
ax.bar(range(len(ARMS)), vals, color=[COL[a] for a in ARMS], width=0.62)
for i, v in enumerate(vals):
    ax.text(i, v + 0.3, f"{v:.1f}", ha="center", fontsize=8.5)
ax.set_ylim(0, max(vals) * 1.22)
ax.set_ylabel("Single-exon transcripts (%)", fontsize=9)
ax.set_title("B", loc="left", fontsize=11, fontweight="bold")

for ax in axes:
    ax.set_xticks(range(1, len(ARMS) + 1) if ax is axes[0] else range(len(ARMS)))
    ax.set_xticklabels(ARMS, fontsize=8.5, rotation=30, ha="right")
    ax.tick_params(axis="y", labelsize=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

plt.tight_layout()
for ext in ("pdf", "png"):
    f = os.path.join(OUT, f"fig3_estructura.{ext}")
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("escrito:", f)
