#!/usr/bin/env python3
"""Figura 4. Exactitud de los modelos medida contra una referencia externa.

Panel A: fraccion de proteinas con cobertura >=80 %, de sujeto y de consulta,
         contra Camelus dromedarius, que es el patron de medida valido porque
         ninguno de los cuatro brazos deriva de el.

Panel B: EL PANEL IMPORTANTE. Indicio de sobre-extension medido contra las dos
         referencias. Contra alpaca, Helixer da 2,02 % frente al 0,04 % de LiftOn;
         contra dromedario los cuatro quedan en un margen estrecho, entre 1,67 y
         2,15 %. La causa es la circularidad: los tres brazos de homologia son
         proyecciones de la anotacion de alpaca y alinean contra su propia fuente.
         Usar alpaca como unico patron habria llevado a una conclusion que es un
         artefacto del patron de medida elegido (el orden de los brazos no se
         invierte; cambia la magnitud, de ~57x a 1,3x).
         El contraste se mantiene a umbrales de 10, 15, 20 y 30 puntos.

NO se reportan las medianas de cobertura: dan 100 % en los cuatro brazos y en
ambas referencias, de modo que saturan y no discriminan.

ANCHO FIJO 170 mm (6.69 pulgadas). NO anadir bbox_inches="tight" al guardar.
"""
import os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Fuentes embebidas como TrueType (Type 42), no Type 3: texto seleccionable y
# sin problemas en los flujos de produccion editoriales.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

# Rutas relativas al repositorio: este script vive en scripts/
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(REPO, "results")
OUT = os.path.join(RESULTS, "figures")
os.makedirs(OUT, exist_ok=True)

MAP = os.path.join(RESULTS, "md5_to_ids.tsv")
HITS = {"dromedario": os.path.join(RESULTS, "blastp", "hits_dromedario.tsv"),
        "alpaca": os.path.join(RESULTS, "blastp", "hits_alpaca.tsv")}
NOMBRE = {"liftoff": "Liftoff", "miniprot": "miniprot",
          "lifton": "LiftOn", "helixer": "Helixer"}
ARMS = ["Liftoff", "miniprot", "LiftOn", "Helixer"]
COL = {"Liftoff": "#0072B2", "miniprot": "#009E73",
       "LiftOn": "#CC79A7", "Helixer": "#E69F00"}

md5_arms = defaultdict(list)
with open(MAP) as fh:
    next(fh)
    for line in fh:
        h, arm, _ = line.rstrip("\n").split("\t")
        md5_arms[h].append(NOMBRE[arm])

def cobertura(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12: continue
            qs, qe, ss, se = int(c[4]), int(c[5]), int(c[6]), int(c[7])
            ql, sl = int(c[10]), int(c[11])
            d[c[0]] = (100.0 * (qe - qs + 1) / ql if ql else 0,
                       100.0 * (se - ss + 1) / sl if sl else 0)
    return d

res = {}
for ref, path in HITS.items():
    cov = cobertura(path)
    est = {a: {"q80": 0, "s80": 0, "over": 0, "n": 0} for a in ARMS}
    for h, arms in md5_arms.items():
        if h not in cov: continue
        q, sj = cov[h]
        for a in arms:
            e = est[a]; e["n"] += 1
            if q >= 80: e["q80"] += 1
            if sj >= 80: e["s80"] += 1
            if sj - q > 20: e["over"] += 1
    res[ref] = est
    print(f"\n--- {ref} ---")
    for a in ARMS:
        e = est[a]
        print(f"  {a:9s} n={e['n']:>6d}  suj>=80 {100*e['s80']/e['n']:5.1f}%  "
              f"con>=80 {100*e['q80']/e['n']:5.1f}%  sobre-ext {100*e['over']/e['n']:5.2f}%")

fig, axes = plt.subplots(1, 2, figsize=(6.69, 3.30))
x = np.arange(len(ARMS)); w = 0.36

# --- Panel A: cobertura contra dromedario
ax = axes[0]
e = res["dromedario"]
s80 = [100 * e[a]["s80"] / e[a]["n"] for a in ARMS]
q80 = [100 * e[a]["q80"] / e[a]["n"] for a in ARMS]
ax.bar(x - w/2, s80, w, label="Subject coverage", color="#0072B2")
ax.bar(x + w/2, q80, w, label="Query coverage", color="#56B4E9")
for i in range(len(ARMS)):
    ax.text(x[i] - w/2 - 0.06, s80[i] + 1.5, f"{s80[i]:.1f}", ha="center", fontsize=7.5)
    ax.text(x[i] + w/2, q80[i] + 1.5, f"{q80[i]:.1f}", ha="center", fontsize=7.5)
ax.set_ylim(0, 100)
# El denominador son las secuencias con hit (e["n"]), no el total del brazo (6.7).
ax.set_ylabel("Proteins $\\geq$80% coverage\n(% of proteins with a hit)", fontsize=9)
ax.set_title("A   coverage vs. $\\it{C.\\ dromedarius}$",
             loc="left", fontsize=9.5, fontweight="bold")
ax.legend(fontsize=8, frameon=False, loc="upper center",
          bbox_to_anchor=(0.5, -0.30), ncol=2)

# --- Panel B: el hallazgo metodologico
ax = axes[1]
dro = [100 * res["dromedario"][a]["over"] / res["dromedario"][a]["n"] for a in ARMS]
alp = [100 * res["alpaca"][a]["over"] / res["alpaca"][a]["n"] for a in ARMS]
ax.bar(x - w/2, alp, w, label="vs. $\\it{V.\\ pacos}$", color="#D55E00")
ax.bar(x + w/2, dro, w, label="vs. $\\it{C.\\ dromedarius}$", color="#0072B2")
for i in range(len(ARMS)):
    ax.text(x[i] - w/2, alp[i] + 0.05, f"{alp[i]:.2f}", ha="center", fontsize=7.5)
    ax.text(x[i] + w/2, dro[i] + 0.05, f"{dro[i]:.2f}", ha="center", fontsize=7.5)
ax.set_ylim(0, max(max(dro), max(alp)) * 1.35)
ax.set_ylabel("Apparent over-extension (%)", fontsize=9)
ax.set_title("B   circular vs. external reference",
             loc="left", fontsize=9.5, fontweight="bold")
ax.legend(fontsize=8, frameon=False, loc="upper left")

for ax in axes:
    ax.set_xticks(x); ax.set_xticklabels(ARMS, fontsize=8.5, rotation=20, ha="right")
    ax.tick_params(axis="y", labelsize=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

plt.tight_layout()
# Margen derecho para que el titulo del panel B no se corte en el borde de la pagina.
fig.subplots_adjust(right=0.99)
for ext in ("pdf", "png"):
    f = os.path.join(OUT, f"fig4_cobertura.{ext}")
    fig.savefig(f, dpi=300)
    print("escrito:", f)
