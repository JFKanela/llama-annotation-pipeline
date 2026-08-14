#!/usr/bin/env python3
"""Ensambla la Tabla 1 definitiva del manuscrito.

Toma el informe de build_report.py, ordena la seccion 1 por completitud BUSCO
para que coincida con la figura 1, sustituye la seccion 2 por la version de los
cuatro brazos, y retitula las secciones 3 y 4 para dejar claro que solo cubren
los dos brazos de homologia puros.
"""
import os, re

T = os.path.expanduser(
    "~/Gdrive/Doctorado/llama_annotation_pipeline/MANUSCRITO/02_tablas")
src = os.path.join(T, "tabla1.md")
sec2 = os.path.join(T, "tabla1_seccion2.md")
dst = os.path.join(T, "Tabla1_FINAL.md")

s = open(src).read()

# --- 1. ordenar la seccion 1 por completitud, de mayor a menor
filas = re.findall(r"^\| (\w+) \| ([\d.]+) \|.*$", s, flags=re.M)
bloque = re.search(r"(\| Proteoma \|.*?\n\|---.*?\n)((?:\|.*\n)+)", s)
cab, cuerpo = bloque.group(1), bloque.group(2)
lineas = [l for l in cuerpo.strip().split("\n") if l.startswith("|")]
lineas.sort(key=lambda l: float(l.split("|")[2]), reverse=True)
s = s.replace(bloque.group(0), cab + "\n".join(lineas) + "\n")

# --- 2. sustituir la seccion de estructura
nueva = open(sec2).read()
ini = s.index("## 2 · Estructura")
fin = s.index("## 3 ·")
s = s[:ini] + nueva.split("## 2 · Estructura de la anotación\n", 1)[1].lstrip() + "\n" + s[fin:]

# --- 3. retitular las secciones de dos brazos
s = s.replace("## 3 · Concordancia A (Liftoff) vs B (miniprot) · gffcompare",
              "## 3 · Concordancia entre los dos brazos de homología pura (gffcompare)")
s = s.replace("## 4 · Recuperacion Liftoff",
              "## 4 · Genes de la referencia no transferidos por Liftoff")
s = s.replace("Referencia = anotacion Liftoff; consulta = anotacion miniprot.",
              "Referencia = anotación Liftoff; consulta = anotación miniprot. Esta sección "
              "compara únicamente los dos brazos de homología pura; LiftOn es un híbrido de "
              "ambos y Helixer no deriva de ninguna referencia.")

open(dst, "w").write(s)
print("escrito:", dst)
print()
print(s)
