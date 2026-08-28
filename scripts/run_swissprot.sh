#!/usr/bin/env bash
# run_swissprot.sh - alineamiento contra Swiss-Prot de los loci de Helixer sin
# ortologo en camelidos.
#
# QUE HACE
#   Descarga Swiss-Prot, construye la base de DIAMOND y alinea contra ella el FASTA
#   que escribe novel_loci_blast.py con los loci huerfanos: los transcritos de
#   Helixer posicionalmente noveles (>= 1 kb) que no dieron hit ni contra
#   Camelus dromedarius ni contra Vicugna pacos.
#
# POR QUE IMPORTA
#   Es la segunda busqueda independiente sobre los huerfanos. De los 235, solo 10
#   tienen homologo en Swiss-Prot. Que dos busquedas independientes no devuelvan
#   nada es lo que sostiene leer los 225 restantes como falsos positivos de la
#   prediccion y no como genes especificos de linaje.
#
# PROCEDENCIA
#   Escrito el 26 de agosto de 2026, a posteriori, a partir del registro de trabajo
#   de la ejecucion original del 14 de agosto. Reproduce los parametros empleados;
#   no es la transcripcion de un script ejecutado entonces.
#
#   Los parametros de abajo NO son los que REPRODUCIBILITY.md documentaba hasta la
#   ronda 8. Aquellos (--max-target-seqs 1 --max-hsps 1 y doce columnas) son
#   incompatibles con el fichero que la ejecucion produjo: novel_vs_sprot.tsv tiene
#   SIETE columnas y 29 filas para 10 consultas, lo que exige --max-target-seqs 5 y
#   el outfmt de abajo. Corregido en la ronda 9.
#
# LA RELEASE IMPORTA
#   La release usada tenia 575.503 secuencias curadas. UniProt sobrescribe
#   current_release/ en cada version, de modo que una reejecucion posterior baja una
#   base DISTINTA y puede devolver otro numero de hits. El script imprime el
#   recuento de lo que ha descargado: comparalo antes de tratar el resultado como
#   una reproduccion y no como una medida nueva.
#
# USO
#   bash scripts/run_swissprot.sh
#
# Ver REPRODUCIBILITY.md, seccion 8.2.

set -euo pipefail

# Directorio de trabajo. Por defecto, el padre del que contiene el script.
BASE="${BASE:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$BASE"
cd blastp

URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

if [ ! -f uniprot_sprot.fasta.gz ]; then
  echo "Descargando Swiss-Prot..."
  wget -q --show-progress "$URL"
fi

echo -n "Secuencias en la base descargada: "
zcat uniprot_sprot.fasta.gz | grep -c '^>'

diamond makedb --in uniprot_sprot.fasta.gz -d sprot

diamond blastp \
  -q helixer_novel_nohit.faa \
  -d sprot.dmnd \
  -o novel_vs_sprot.tsv \
  --very-sensitive \
  --max-target-seqs 5 \
  --evalue 1e-5 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
  --threads 4

echo -n "Consultas con hit: "
cut -f1 novel_vs_sprot.tsv | sort -u | wc -l
