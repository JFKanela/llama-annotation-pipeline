#!/usr/bin/env bash
# run_swissprot.sh - Swiss-Prot search of the Helixer loci with no camelid hit.
#
# WHAT IT DOES
#   Downloads the Swiss-Prot release, builds a DIAMOND database and aligns
#   against it the FASTA that novel_loci_blast.py writes with the orphan loci:
#   the positionally novel Helixer transcripts (>= 1 kb) that returned no hit
#   against either Camelus dromedarius or Vicugna pacos.
#
# WHY IT MATTERS
#   It is the second independent test of the orphan loci. Of the 235 orphans,
#   only 10 have a Swiss-Prot homologue. Two independent searches returning
#   nothing is what supports reading the remaining 225 as false positives of the
#   prediction rather than lineage-specific genes.
#
# PROVENANCE - READ THIS
#   This script was written AFTER the run it documents, in round 8, by
#   reconstructing the commands from the working record. It is NOT a transcript
#   of an executed script and it has NOT been re-run against the original
#   inputs. Its parameters mirror those of REPRODUCIBILITY.md section 7.3, which
#   is what the run used. See REPRODUCIBILITY.md section 8.2.
#
# THE RELEASE MATTERS
#   The Swiss-Prot release used contained 575,503 sequences. UniProt overwrites
#   current_release/ on every release, so re-running this later downloads a
#   DIFFERENT database and may return a different number of hits. The script
#   prints the sequence count it actually downloaded; compare it before treating
#   the result as a reproduction rather than a new measurement.
#
# USAGE
#   bash scripts/run_swissprot.sh
set -euo pipefail

BASE="${BASE:-$HOME/llama_annotation_pipeline}"
WORK="${WORK:-$BASE/blastp}"
QUERY="${QUERY:-$WORK/helixer_novel_nohit.faa}"
DB_FASTA="$WORK/uniprot_sprot.fasta.gz"
DB="$WORK/swissprot"
OUT="$WORK/hits_swissprot.tsv"
THREADS="${THREADS:-4}"

URL="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

mkdir -p "$WORK"

if [[ ! -s "$QUERY" ]]; then
    echo "ERROR: no existe la consulta $QUERY" >&2
    echo "       Ejecuta antes scripts/novel_loci_blast.py, que la produce." >&2
    exit 1
fi

# --- 1. base de datos
if [[ ! -s "$DB_FASTA" ]]; then
    echo "[1/3] descargando Swiss-Prot..."
    wget -q --show-progress -O "$DB_FASTA.partial" "$URL"
    mv "$DB_FASTA.partial" "$DB_FASTA"
else
    echo "[1/3] Swiss-Prot ya descargada, se omite."
fi

N_SEQ=$(gunzip -c "$DB_FASTA" | grep -c '^>')
echo "      secuencias en la release descargada: $N_SEQ"
echo "      la release usada en el manuscrito tenia: 575503"
if [[ "$N_SEQ" != "575503" ]]; then
    echo "      AVISO: release distinta de la usada. El resultado NO es una" >&2
    echo "             reproduccion, es una medida nueva." >&2
fi

# --- 2. indice
if [[ ! -s "$DB.dmnd" ]]; then
    echo "[2/3] diamond makedb..."
    diamond makedb --in "$DB_FASTA" -d "$DB" --threads "$THREADS"
else
    echo "[2/3] indice ya presente, se omite."
fi

# --- 3. alineamiento
if [[ -s "$OUT" ]]; then
    echo "[3/3] $OUT ya existe y no esta vacio, se omite."
else
    echo "[3/3] diamond blastp contra Swiss-Prot..."
    diamond blastp \
      -q "$QUERY" \
      -d "$DB" \
      -o "$OUT.partial" \
      --outfmt 6 qseqid sseqid pident length qstart qend sstart send \
                  evalue bitscore qlen slen \
      --very-sensitive \
      --max-target-seqs 1 \
      --max-hsps 1 \
      --evalue 1e-5 \
      --threads "$THREADS"
    mv "$OUT.partial" "$OUT"
fi

echo
echo "consultas enviadas : $(grep -c '^>' "$QUERY")"
echo "consultas con hit  : $(cut -f1 "$OUT" | sort -u | wc -l)"
echo "escrito: $OUT"
