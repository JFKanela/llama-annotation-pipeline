#!/usr/bin/env python3
"""Figura 1. Completitud BUSCO de los cuatro proteomas, con techo y linea de base.

Datos: results/busco_summaries/*.txt (linaje artiodactyla_odb12, n=12594).
Salida: PDF vectorial para el manuscrito, mas PNG para revisar.

ANCHO FIJO 170 mm (6.69 pulgadas), que es el ancho de pagina completa de la revista.
NO anadir bbox_inches="tight" al guardar: expande el lienzo y baja el tamano de
letra por debajo del minimo de 7 pt exigido.
Paleta Okabe-Ito, segura para daltonismo.
"""
import re, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Fuentes embebidas como TrueType (Type 42), no Type 3: texto seleccionable y
# sin problemas en los flujos de produccion editoriales.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

# Rutas relativas al repositorio: este script vive en scripts/
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(REPO, "results")
OUT = os.path.join(RESULTS, "figures")
os.makedirs(OUT, exist_ok=True)

# Resumenes BUSCO, un fichero por brazo con el nombre de la clave de ORDEN.
# Se obtienen de results/busco/{brazo}/short_summary.specific.*.txt
BUSCO = os.path.join(RESULTS, "busco_summaries")

# orden de abajo arriba en la figura
ORDEN = [
    ("chaku_v1",    "chaku_v1 (previous in-house)"),
    ("helixer_web", "Helixer"),
    ("liftoff",     "Liftoff"),
    ("miniprot",    "miniprot"),
    ("lifton",      "LiftOn"),
    ("alpaca_ref",  "V. pacos reference (ceiling)"),
]
PAT = re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%")
COL = {"S": "#0072B2", "D": "#56B4E9", "F": "#E69F00", "M": "#D55E00"}

datos, labels = [], []
for clave, etiqueta in ORDEN:
    p = os.path.join(BUSCO, clave + ".txt")
    m = None
    with open(p) as fh:
        for line in fh:
            m = PAT.search(line)
            if m:
                break
    if not m:
        raise SystemExit(f"no encuentro la linea C: en {p}")
    C = float(m.group(1))
    S, Dd, F, M = (float(m.group(i)) for i in (2, 3, 4, 5))
    datos.append((S, Dd, F, M, C))
    labels.append(etiqueta)
    print(f"{etiqueta:32s} S={S:5.1f} D={Dd:4.1f} F={F:4.1f} M={M:5.1f}")

fig, ax = plt.subplots(figsize=(6.69, 3.16))
y = range(len(datos))
izq = [0.0] * len(datos)
for i, (k, nombre) in enumerate([("S", "Single-copy"),
                                 ("D", "Duplicated"),
                                 ("F", "Fragmented"),
                                 ("M", "Missing")]):
    v = [d[i] for d in datos]
    # chaku_v1 (indice 0) se atenua: es comparacion contextual, no un brazo evaluable
    alfas = [0.35 if j == 0 else 1.0 for j in range(len(datos))]
    for j, (yy, vv, ll) in enumerate(zip(y, v, izq)):
        ax.barh(yy, vv, left=ll, color=COL[k], alpha=alfas[j],
                label=nombre if j == 1 else None,
                height=0.62, edgecolor="white", linewidth=0.6)
    izq = [a + b for a, b in zip(izq, v)]

# porcentaje completo al final de cada barra
for i, d in enumerate(datos):
    ax.text(101, i, f"{d[4]:.1f}%", va="center", fontsize=8.5)

ax.axhline(0.5, color="#9A9A9A", linewidth=0.8, linestyle=(0, (4, 3)))
ax.set_yticks(list(y))
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlim(0, 108)
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("% of BUSCO groups (artiodactyla_odb12, n = 12,594)", fontsize=9)
ax.tick_params(axis="x", labelsize=8.5)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=4,
          frameon=False, fontsize=8.5)
plt.tight_layout()

for ext in ("pdf", "png"):
    f = os.path.join(OUT, f"fig1_busco.{ext}")
    fig.savefig(f, dpi=300)
    print("escrito:", f)
