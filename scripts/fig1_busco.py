#!/usr/bin/env python3
"""Figura 1. Completitud BUSCO de los cuatro proteomas, con techo y linea de base.

Datos: MANUSCRITO/03_datos/busco/*.txt (linaje artiodactyla_odb12, n=12594).
Salida: PDF vectorial para el manuscrito, mas PNG para revisar.
Paleta Okabe-Ito, segura para daltonismo.
"""
import re, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = os.path.expanduser(
    "~/Gdrive/Doctorado/llama_annotation_pipeline/MANUSCRITO")
BUSCO = os.path.join(D, "03_datos", "busco")
OUT = os.path.join(D, "01_figuras")
os.makedirs(OUT, exist_ok=True)

# orden de abajo arriba en la figura
ORDEN = [
    ("chaku_v1",    "chaku_v1 (previous in-house)"),
    ("helixer_web", "Helixer"),
    ("liftoff",     "Liftoff"),
    ("miniprot",    "miniprot"),
    ("lifton",      "LiftOn"),
    ("alpaca_ref",  "V. pacos reference (ceiling)"),
]
PAT = re.compile(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%')
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

fig, ax = plt.subplots(figsize=(7.2, 3.4))
y = range(len(datos))
izq = [0.0] * len(datos)
for i, (k, nombre) in enumerate([("S", "Complete, single-copy"),
                                 ("D", "Complete, duplicated"),
                                 ("F", "Fragmented"),
                                 ("M", "Missing")]):
    v = [d[i] for d in datos]
    ax.barh(list(y), v, left=izq, color=COL[k], label=nombre,
            height=0.62, edgecolor="white", linewidth=0.6)
    izq = [a + b for a, b in zip(izq, v)]

# porcentaje completo al final de cada barra
for i, d in enumerate(datos):
    ax.text(101, i, f"{d[4]:.1f}%", va="center", fontsize=8.5)

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
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("escrito:", f)
