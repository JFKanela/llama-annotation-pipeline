#!/usr/bin/env python3
"""Figura 2. Solapamiento entre los cuatro proteomas, por identidad exacta de secuencia.

Datos: md5_to_ids.tsv, el MISMO mapa que produjo las cifras del manuscrito.
Diagrama UpSet dibujado con matplotlib puro, sin la libreria upsetplot, que da
problemas de compatibilidad con las versiones actuales de pandas y matplotlib.

ADVERTENCIA para el pie de figura: la comparacion es por identidad EXACTA de
secuencia. Las secuencias exclusivas de Helixer NO son genes especificos de llama.
La novedad real por posicion genomica se trata aparte.

ANCHO FIJO 170 mm (6.69 pulgadas). NO anadir bbox_inches="tight" al guardar.
"""
import os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Fuentes embebidas como TrueType (Type 42), no Type 3: texto seleccionable y
# sin problemas en los flujos de produccion editoriales.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
from matplotlib.gridspec import GridSpec

# Rutas relativas al repositorio: este script vive en scripts/
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(REPO, "results")
OUT = os.path.join(RESULTS, "figures")
os.makedirs(OUT, exist_ok=True)

MAP = os.path.join(RESULTS, "md5_to_ids.tsv")

NOMBRE = {"liftoff": "Liftoff", "miniprot": "miniprot",
          "lifton": "LiftOn", "helixer": "Helixer"}
ORDEN = ["Liftoff", "miniprot", "LiftOn", "Helixer"]
AZUL, GRIS = "#0072B2", "#cccccc"

md5_arms = defaultdict(set)
with open(MAP) as fh:
    next(fh)
    for line in fh:
        h, arm, _ = line.rstrip("\n").split("\t")
        md5_arms[h].add(NOMBRE[arm])
print(f"secuencias unicas: {len(md5_arms)}")

combo = defaultdict(int)
for arms in md5_arms.values():
    combo[frozenset(arms)] += 1
items = sorted(combo.items(), key=lambda x: -x[1])
for k, v in items:
    print(f"  {v:>6d}  {' + '.join(a for a in ORDEN if a in k)}")

totales = {a: sum(v for k, v in items if a in k) for a in ORDEN}
n = len(items)

fig = plt.figure(figsize=(6.69, 5.22))
gs = GridSpec(2, 2, width_ratios=[1, 4.2], height_ratios=[3.6, 1.15],
              hspace=0.06, wspace=0.06, figure=fig)

# --- barras de intersecciones, arriba a la derecha
ax = fig.add_subplot(gs[0, 1])
x = range(n)
vals = [v for _, v in items]
ax.bar(x, vals, color=AZUL, width=0.62)
for i, v in enumerate(vals):
    ax.text(i, v + max(vals) * 0.02, f"{v:,}",
            ha="center", va="bottom", fontsize=7.2, rotation=90)
ax.set_ylim(0, max(vals) * 1.24)
ax.set_xlim(-0.7, n + 2.6)
ax.set_ylabel("Sequences in\nintersection", fontsize=8.5)
ax.tick_params(axis="y", labelsize=8)
ax.set_xticks([])
for s in ("top", "right", "bottom"):
    ax.spines[s].set_visible(False)

# --- matriz de puntos, abajo a la derecha
axm = fig.add_subplot(gs[1, 1], sharex=ax)
for j, arm in enumerate(ORDEN):
    y = len(ORDEN) - 1 - j
    if j % 2 == 0:
        axm.axhspan(y - 0.5, y + 0.5, color="#f5f5f5", zorder=0)
    for i, (k, _) in enumerate(items):
        axm.plot(i, y, "o", ms=7,
                 color=AZUL if arm in k else GRIS, zorder=3)
for i, (k, _) in enumerate(items):
    ys = [len(ORDEN) - 1 - j for j, a in enumerate(ORDEN) if a in k]
    if len(ys) > 1:
        axm.plot([i, i], [min(ys), max(ys)], "-", color=AZUL, lw=1.6, zorder=2)
for j, arm in enumerate(ORDEN):
    axm.text(n - 0.15, len(ORDEN) - 1 - j, "  " + arm,
             ha="left", va="center", fontsize=9.5)
axm.set_xlim(-0.7, n + 2.6)
axm.set_yticks(range(len(ORDEN)))
axm.set_yticklabels([])
axm.set_ylim(-0.6, len(ORDEN) - 0.4)
axm.set_xticks([])
for s in ("top", "right", "bottom", "left"):
    axm.spines[s].set_visible(False)
axm.tick_params(length=0)

# --- barras de totales, abajo a la izquierda
axt = fig.add_subplot(gs[1, 0], sharey=axm)
MAXT = max(totales.values())
for j, arm in enumerate(ORDEN):
    y = len(ORDEN) - 1 - j
    axt.barh(y, totales[arm], color=AZUL, height=0.62)
    axt.text(-0.02, y, f"{totales[arm]:,}",
             transform=axt.get_yaxis_transform(),
             ha="right", va="center", fontsize=8.5)
axt.set_xlim(MAXT * 1.02, 0)
axt.set_xlabel("Sequences\nper set", fontsize=8.5)
axt.tick_params(axis="x", labelsize=7.5)
axt.set_yticks([])
for s in ("top", "right", "left"):
    axt.spines[s].set_visible(False)

fig.add_subplot(gs[0, 0]).axis("off")

for ext in ("pdf", "png"):
    f = os.path.join(OUT, f"fig2_upset.{ext}")
    fig.savefig(f, dpi=300)
    print("escrito:", f)
