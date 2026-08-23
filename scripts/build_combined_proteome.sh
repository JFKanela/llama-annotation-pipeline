#!/usr/bin/env bash
# build_combined_proteome.sh - combined reference proteome for Lama glama.
#
# WHAT IT BUILDS
#   The LiftOn core (20,233 proteins) extended with the positionally novel
#   Helixer loci that DO have an orthologue in another camelid (555), minus the
#   ones redundant with the core (80). Net: 475 added, 20,708 total.
#
#     20,233 + 555 - 80 = 20,708
#
#   That is the file deposited as llama_combined_reference_proteome.faa.gz in the
#   data record, version 1.1.0, DOI 10.5281/zenodo.22072343. Its working name
#   during the analysis was Lgla_combined_reference_proteome_v1.faa.
#
# PROVENANCE - READ THIS BEFORE REUSING
#   This script is a RECONSTRUCTION. It was written in round 8, after the fact,
#   from the working record; it is NOT a transcript of the commands that produced
#   the deposited file, and it has NOT been re-executed against the original
#   inputs. The deposited FASTA is the authoritative artefact. The composition of
#   the deposited file WAS verified by direct counting (20,233 LiftOn-origin
#   identifiers + 475 Helixer-origin = 20,708); the route below is the best
#   available account of how it got there, not a validated build.
#   See REPRODUCIBILITY.md section 8.4.
#
# WHY LIFTON IS THE CORE
#   Not because it is the most complete, but because it is the best-modelled
#   homology branch: highest fraction of proteins covering >= 80 % of the query
#   against the external reference (92.6 %), and fewest internal stop codons of
#   the three homology branches (2.7 %).
#
# PREREQUISITES
#   novel_loci_blast.py must have produced blastp/helixer_novel_loci.tsv.
#   Needs diamond and seqkit on PATH.
#
# USAGE
#   bash scripts/build_combined_proteome.sh
set -euo pipefail

BASE="${BASE:-$HOME/llama_annotation_pipeline}"
WORK="${WORK:-$BASE/blastp}"
LIFTON="${LIFTON:-$BASE/results/proteomes/llama_lifton_v2.faa}"
HELIXER="${HELIXER:-$BASE/Lgla_hx036_helixer.faa}"
NOVEL_TSV="${NOVEL_TSV:-$WORK/helixer_novel_loci.tsv}"
OUT="${OUT:-$BASE/Lgla_combined_reference_proteome_v1.faa}"
THREADS="${THREADS:-4}"

# umbrales de redundancia contra el nucleo
MIN_PIDENT="${MIN_PIDENT:-95}"
MIN_QCOV="${MIN_QCOV:-80}"

for f in "$LIFTON" "$HELIXER" "$NOVEL_TSV"; do
    [[ -s "$f" ]] || { echo "ERROR: falta $f" >&2; exit 1; }
done
mkdir -p "$WORK"

# --- 1. los loci rescatados: con hit en dromedario o en alpaca
awk -F'\t' 'NR>1 && ($7=="si" || $8=="si") {print $1}' "$NOVEL_TSV" \
    > "$WORK/rescued_ids.txt"
N_RESCUED=$(wc -l < "$WORK/rescued_ids.txt")
echo "[1/5] loci con ortologo camelido: $N_RESCUED   (esperado 555)"

# --- 2. sus secuencias
seqkit grep -f "$WORK/rescued_ids.txt" "$HELIXER" > "$WORK/rescued.faa"
echo "[2/5] secuencias extraidas: $(grep -c '^>' "$WORK/rescued.faa")"

# --- 3. redundancia contra el nucleo LiftOn
if [[ ! -s "$WORK/lifton_core.dmnd" ]]; then
    diamond makedb --in "$LIFTON" -d "$WORK/lifton_core" --threads "$THREADS"
fi
diamond blastp \
  -q "$WORK/rescued.faa" \
  -d "$WORK/lifton_core" \
  -o "$WORK/rescued_vs_core.tsv" \
  --outfmt 6 qseqid sseqid pident length qstart qend sstart send \
              evalue bitscore qlen slen \
  --very-sensitive --max-target-seqs 1 --max-hsps 1 --evalue 1e-5 \
  --threads "$THREADS"
echo "[3/5] alineados contra el nucleo."

# --- 4. descartar los redundantes
awk -F'\t' -v pid="$MIN_PIDENT" -v qc="$MIN_QCOV" \
    '$3 >= pid && (100*($6-$5+1)/$11) >= qc {print $1}' \
    "$WORK/rescued_vs_core.tsv" | sort -u > "$WORK/redundant_ids.txt"
N_RED=$(wc -l < "$WORK/redundant_ids.txt")
echo "[4/5] redundantes con el nucleo (>= ${MIN_PIDENT} % id, >= ${MIN_QCOV} % cob): $N_RED   (esperado 80)"

grep -v -x -F -f "$WORK/redundant_ids.txt" "$WORK/rescued_ids.txt" \
    > "$WORK/added_ids.txt" || true
N_ADD=$(wc -l < "$WORK/added_ids.txt")

# --- 5. concatenar
seqkit grep -f "$WORK/added_ids.txt" "$HELIXER" > "$WORK/added.faa"
cat "$LIFTON" "$WORK/added.faa" > "$OUT"

N_CORE=$(grep -c '^>' "$LIFTON")
N_TOTAL=$(grep -c '^>' "$OUT")
echo "[5/5] nucleo $N_CORE + anadidos $N_ADD = $N_TOTAL"
echo
echo "escrito: $OUT"
if [[ "$N_TOTAL" != "20708" ]]; then
    echo "AVISO: el total no coincide con las 20708 proteinas del fichero" >&2
    echo "       depositado. Revisa entradas y umbrales antes de usarlo." >&2
fi
