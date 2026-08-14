#!/usr/bin/env python3
"""Figura 4. Exactitud de los modelos medida contra una referencia externa.

Panel A: fraccion de proteinas con cobertura >=80 %, de sujeto y de consulta,
         contra Camelus dromedarius, que es el patron de medida valido porque
         ninguno de los cuatro brazos deriva de el.

Panel B: EL PANEL IMPORTANTE. Indicio de sobre-extension medido contra las dos
         referencias. Contra alpaca, Helixer parece sobre-extender veinte veces
         mas que LiftOn; contra dromedario los cuatro son indistinguibles.
         La causa es la circularidad: los tres brazos de homologia son
         proyecciones de la anotacion de alpaca y alinean contra su propia fuente.
         USAR ALPACA COMO PATRON HABRIA PRODUCIDO UNA CONCLUSION FALSA.

NO se reportan las medianas de cobertura: dan 100 % en los cuatro brazos y en
ambas referencias, de modo que saturan y no discriminan.
"""
import os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.expanduser("~/llama_annotation_pipeline")
OUT = os.path.expanduser(
    "~/Gdrive/Doctorado/llama_annotation_pipeline/MANUSCRITO/01_figuras")
os.makedirs(OUT, exist_ok=True)

MAP = os.path.join(BASE, "md5_to_ids.tsv")
HITS = {"dromedario": os.path.join(BASE, "blastp", "hits_dromedario.tsv"),
        "alpaca": os.path.join(BASE, "blastp", "hits_alpaca.tsv")}
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
              f"con>=80 {100*e['q80']/e['n']:5.1f}%  sobre-ext {100*e['over']/e['n']:4.1f}%")

fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.9))
x = np.arange(len(ARMS)); w = 0.36

# --- Panel A: cobertura contra dromedario
ax = axes[0]
e = res["dromedario"]
s80 = [100 * e[a]["s80"] / e[a]["n"] for a in ARMS]
q80 = [100 * e[a]["q80"] / e[a]["n"] for a in ARMS]
ax.bar(x - w/2, s80, w, label="Subject coverage", color="#0072B2")
ax.bar(x + w/2, q80, w, label="Query coverage", color="#56B4E9")
for i in range(len(ARMS)):
    ax.text(x[i] - w/2, s80[i] + 0.6, f"{s80[i]:.1f}", ha="center", fontsize=7.5)
    ax.text(x[i] + w/2, q80[i] + 0.6, f"{q80[i]:.1f}", ha="center", fontsize=7.5)
ax.set_ylim(78, 96)
ax.set_ylabel("Proteins with coverage $\\geq$ 80 % (%)", fontsize=9)
ax.set_title("A   vs. $\\it{C.\\ dromedarius}$ (external)",
             loc="left", fontsize=10, fontweight="bold")
ax.legend(fontsize=8, frameon=False, loc="upper center",
          bbox_to_anchor=(0.5, -0.30), ncol=2)

# --- Panel B: el hallazgo metodologico
ax = axes[1]
dro = [100 * res["dromedario"][a]["over"] / res["dromedario"][a]["n"] for a in ARMS]
alp = [100 * res["alpaca"][a]["over"] / res["alpaca"][a]["n"] for a in ARMS]
ax.bar(x - w/2, alp, w, label="vs. $\\it{V.\\ pacos}$ (circular)", color="#D55E00")
ax.bar(x + w/2, dro, w, label="vs. $\\it{C.\\ dromedarius}$ (external)", color="#0072B2")
for i in range(len(ARMS)):
    ax.text(x[i] - w/2, alp[i] + 0.05, f"{alp[i]:.1f}", ha="center", fontsize=7.5)
    ax.text(x[i] + w/2, dro[i] + 0.05, f"{dro[i]:.1f}", ha="center", fontsize=7.5)
ax.set_ylim(0, max(max(dro), max(alp)) * 1.35)
ax.set_ylabel("Apparent over-extension (%)", fontsize=9)
ax.set_title("B   the circular reference misleads",
             loc="left", fontsize=10, fontweight="bold")
ax.legend(fontsize=8, frameon=False, loc="upper left")

for ax in axes:
    ax.set_xticks(x); ax.set_xticklabels(ARMS, fontsize=8.5, rotation=20, ha="right")
    ax.tick_params(axis="y", labelsize=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

plt.tight_layout()
for ext in ("pdf", "png"):
    f = os.path.join(OUT, f"fig4_cobertura.{ext}")
    fig.savefig(f, dpi=300, bbox_inches="tight")
    print("escrito:", f)
