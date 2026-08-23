#!/usr/bin/env python3
"""Ensambla la Tabla 1 definitiva del manuscrito.

Toma el informe de build_report.py, ordena la seccion 1 por completitud BUSCO
para que coincida con la figura 1, sustituye la seccion 2 por la version de los
cuatro brazos, y retitula las secciones 3 y 4 para dejar claro que solo cubren
los dos brazos de homologia puros.
"""
import os, re, sys

T = os.path.expanduser(
    "~/Gdrive/Doctorado/llama_annotation_pipeline/MANUSCRITO/02_tablas")
src = os.path.join(T, "tabla1.md")
sec2 = os.path.join(T, "tabla1_seccion2.md")
dst = os.path.join(T, "Tabla1_FINAL.md")

def fatal(msg):
    """Aborta con un mensaje que diga QUE hay que arreglar, no un traceback.

    Este script localiza las secciones de la salida de build_report.py por texto
    exacto del encabezado. Es fragil por diseno: si alli cambia un titulo, aqui
    deja de funcionar. Lo que no es aceptable es que lo haga con un ValueError
    pelado que no dice donde mirar.
    """
    sys.exit("ERROR en tabla1_final.py: " + msg)


for f in (src, sec2):
    if not os.path.exists(f):
        fatal("no existe %s.\n"
              "  tabla1.md lo produce build_report.py y tabla1_seccion2.md lo\n"
              "  produce tabla_estructura.py. Ejecuta ambos antes que este." % f)

s = open(src, encoding="utf-8").read()

# --- 1. ordenar la seccion 1 por completitud, de mayor a menor
bloque = re.search(r"(\| Proteoma \|.*?\n\|---.*?\n)((?:\|.*\n)+)", s)
if bloque is None:
    fatal("no encuentro la tabla que empieza por '| Proteoma |' en %s.\n"
          "  Si build_report.py ha cambiado ese encabezado, hay que actualizar\n"
          "  la expresion regular de este script." % src)
cab, cuerpo = bloque.group(1), bloque.group(2)
lineas = [l for l in cuerpo.strip().split("\n") if l.startswith("|")]
try:
    lineas.sort(key=lambda l: float(l.split("|")[2]), reverse=True)
except (IndexError, ValueError) as e:
    fatal("la segunda columna de la tabla de proteomas no es un numero (%s).\n"
          "  Se ordena por completitud BUSCO para que la Tabla 1 coincida con la\n"
          "  Figura 1; revisa el formato que emite build_report.py." % e)
s = s.replace(bloque.group(0), cab + "\n".join(lineas) + "\n")

# --- 2. sustituir la seccion de estructura
nueva = open(sec2, encoding="utf-8").read()
MARCA = "## 2 · Estructura de la anotación\n"
try:
    ini = s.index("## 2 · Estructura")
    fin = s.index("## 3 ·")
except ValueError:
    fatal("no encuentro los encabezados '## 2 · Estructura' y/o '## 3 ·' en\n"
          "  %s. Este script los localiza por texto exacto: si cambia un titulo\n"
          "  en build_report.py, hay que actualizarlo aqui tambien." % src)
if MARCA not in nueva:
    fatal("no encuentro '%s' en %s.\n"
          "  Lo escribe tabla_estructura.py; comprueba que se ejecuto y que el\n"
          "  titulo coincide caracter a caracter." % (MARCA.strip(), sec2))
s = s[:ini] + nueva.split(MARCA, 1)[1].lstrip() + "\n" + s[fin:]

# --- 3. retitular las secciones de dos brazos
s = s.replace("## 3 · Concordancia A (Liftoff) vs B (miniprot) · gffcompare",
              "## 3 · Concordancia entre los dos brazos de homología pura (gffcompare)")
s = s.replace("## 4 · Recuperacion Liftoff",
              "## 4 · Genes de la referencia no transferidos por Liftoff")
s = s.replace("Referencia = anotacion Liftoff; consulta = anotacion miniprot.",
              "Referencia = anotación Liftoff; consulta = anotación miniprot. Esta sección "
              "compara únicamente los dos brazos de homología pura; LiftOn es un híbrido de "
              "ambos y Helixer no deriva de ninguna referencia.")

open(dst, "w", encoding="utf-8").write(s)
print("escrito:", dst)
print()
print(s)
