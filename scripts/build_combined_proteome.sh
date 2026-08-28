#!/usr/bin/env bash
# build_combined_proteome.sh - proteoma candidato combinado de referencia.
#
# QUE CONSTRUYE
#   El nucleo LiftOn (20.233 secuencias) extendido con los loci de Helixer sin
#   solape posicional que SI tienen ortologo en otro camelido (555), menos los
#   redundantes con el nucleo (80). Neto: 475 anadidos, 20.708 en total.
#
#     20.233 + 555 - 80 = 20.708
#
#   Es el fichero depositado como llama_combined_reference_proteome.faa.gz en el
#   registro de datos, version 1.3.0, DOI 10.5281/zenodo.22150587. Su nombre de
#   trabajo durante el analisis fue Lgla_combined_reference_proteome_v1.faa.
#
# PROCEDENCIA Y VERIFICACION
#   Escrito el 26 de agosto de 2026, a posteriori, a partir del registro de la
#   ejecucion original del 14 de agosto. NO es la transcripcion de un script
#   ejecutado entonces. Pero, a diferencia de la version que se subio en la ronda
#   8, esta SI esta verificada: reproduce el fichero depositado byte a byte.
#
#     MD5 del fichero de referencia: dc69fa820facc0f087696bdd4885e1c1
#     Tamano: 11.337.876 bytes, 20.708 secuencias
#
#   La comprobacion final aborta con codigo de error si el recuento no da 20.708.
#
# POR QUE EL NUCLEO ES LIFTON
#   No por ser el mas completo, sino por ser el brazo de homologia mejor modelado:
#   la mayor fraccion de proteinas que cubren >= 80 % de la consulta contra la
#   referencia externa (92,6 %) y los menos stops internos de los tres brazos de
#   homologia (2,7 %).
#
# SECUENCIAS, NO PROTEINAS
#   El recuento se expresa en secuencias. No son proteinas validadas: parte del
#   nucleo arrastra pautas de lectura rotas, y esa es justamente la informacion que
#   recoge Lgla_combined_proteome_confidence.tsv (scripts/capas.py).
#
# REQUISITOS
#   diamond y python3 en el PATH. No necesita seqkit.
#   blastp/helixer_novel_loci.tsv lo produce novel_loci_blast.py.
#
# USO
#   bash scripts/build_combined_proteome.sh
#   BASE=/otra/ruta bash scripts/build_combined_proteome.sh
#
# Ver REPRODUCIBILITY.md, seccion 8.4.

set -euo pipefail

# Directorio de trabajo. Por defecto, el padre del que contiene el script,
# que es la disposicion de la ejecucion original (<base>/scripts/).
BASE="${BASE:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$BASE"

echo "== 1. loci de Helixer con ortologo en algun camelido =="
awk -F'\t' 'NR>1 && ($7=="si" || $8=="si") {print $1}' \
    blastp/helixer_novel_loci.tsv | sort -u > /tmp/ids_rescatados.txt
echo -n "   rescatados (esperado 555): "; wc -l < /tmp/ids_rescatados.txt

echo "== 2. extraer sus secuencias del proteoma de Helixer =="
python3 - <<'EOF'
ids = set(open("/tmp/ids_rescatados.txt").read().split())
keep, n = False, 0
with open("blastp/helixer_rescued.faa", "w") as out:
    for line in open("Lgla_hx036_helixer.faa"):
        if line.startswith(">"):
            keep = line[1:].split()[0] in ids
            if keep: n += 1
        if keep: out.write(line)
print("   secuencias extraidas:", n)
EOF

echo "== 3. base de DIAMOND sobre el nucleo SANEADO =="
# DIAMOND aborta ante el caracter '.' que gffread emite para los codones de stop,
# de modo que la BASE se construye sobre una copia con '.' sustituido por 'X'.
# El fichero de SALIDA se construye despues sobre el proteoma ORIGINAL, para no
# perder la marca de los stops internos (decision declarada en el apartado 3.4).
zcat ~/zenodo_data/llama_lifton.faa.gz > /tmp/lifton_orig.faa
python3 - <<'EOF'
import re
VALID = "ACDEFGHIKLMNPQRSTVWYXBZJUO"
pat = re.compile("[^%s]" % VALID)
hdr, seq, n = None, [], 0
with open("/tmp/lifton_clean.faa", "w") as out:
    def vuelca():
        global n
        if hdr is None: return
        s = pat.sub("X", "".join(seq).upper().rstrip(".*"))
        n += 1
        out.write(hdr + "\n")
        for i in range(0, len(s), 60):
            out.write(s[i:i+60] + "\n")
    for line in open("/tmp/lifton_orig.faa"):
        if line.startswith(">"):
            vuelca(); hdr, seq = line.strip(), []
        else:
            seq.append(line.strip())
    vuelca()
print("   secuencias saneadas:", n)
EOF

echo "== 4. detectar redundancia con el nucleo =="
diamond makedb --in /tmp/lifton_clean.faa -d /tmp/lifton_clean --threads 4
diamond blastp -q blastp/helixer_rescued.faa -d /tmp/lifton_clean.dmnd \
  -o blastp/rescued_vs_lifton.tsv --very-sensitive \
  --max-target-seqs 1 --max-hsps 1 \
  --outfmt 6 qseqid sseqid pident length qlen slen --threads 4

echo "== 5. redundantes al 95 % de identidad y 80 % de cobertura =="
LC_ALL=C awk -F'\t' '($3+0)>=95 && ($4+0)/($5+0)>=0.8 {print $1}' \
  blastp/rescued_vs_lifton.tsv | sort -u > /tmp/redundantes.txt
echo -n "   redundantes (esperado 80): "; wc -l < /tmp/redundantes.txt

echo "== 6. ensamblar el proteoma combinado =="
# Sobre el proteoma ORIGINAL, con los codones de stop conservados.
python3 - <<'EOF'
red = set(open("/tmp/redundantes.txt").read().split())
n = 0
with open("Lgla_combined_reference_proteome_v1.faa", "w") as out:
    for line in open("/tmp/lifton_orig.faa"):
        out.write(line)
        if line.startswith(">"): n += 1
    keep = False
    for line in open("blastp/helixer_rescued.faa"):
        if line.startswith(">"):
            keep = line[1:].split()[0] not in red
            if keep: n += 1
        if keep: out.write(line)
print("   secuencias en el proteoma combinado:", n)
EOF

echo "== comprobacion final =="
N=$(grep -c '^>' Lgla_combined_reference_proteome_v1.faa)
echo "   secuencias: $N (esperado 20708)"
[ "$N" -eq 20708 ] && echo "   OK" || { echo "   DISCREPANCIA"; exit 1; }
