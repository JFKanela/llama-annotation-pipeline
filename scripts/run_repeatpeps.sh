#!/usr/bin/env bash
# run_repeatpeps.sh - alineamiento de los proteomas contra la base de proteinas
# de elementos transponibles de RepeatMasker (RepeatPeps.lib).
#
# QUE HACE
#   Descarga RepeatMasker 4.1.7-p1, extrae Libraries/RepeatPeps.lib (18.011
#   proteinas de elementos transponibles, 16,1 millones de residuos), construye la
#   base de DIAMOND y alinea contra ella siete conjuntos: los cuatro proteomas de
#   llama, el proteoma candidato combinado y los dos proteomas de referencia de
#   RefSeq (Camelus dromedarius y Vicugna pacos), que sirven de calibracion.
#
# POR QUE IMPORTA
#   El ensamblado no se enmascaro antes de anotar (Helixer convierte la secuencia
#   a mayusculas al importarla, de modo que un soft-masking no le afecta). Esta es
#   la comprobacion a posteriori de si esa omision introdujo ORF de elementos
#   transponibles en el brazo ab initio, y en particular en los 790 loci
#   posicionalmente noveles y en los 225 sin ortologo. El resumen lo hace
#   repeatpeps_summary.py.
#
# CRITERIO
#   El mismo que aplica funannotate a los modelos de EVidenceModeler
#   (funannotate/library.py, RepeatBlast + RemoveBadModels): DIAMOND blastp contra
#   RepeatPeps, e-value 1e-10, un solo hit, y cualquier alineamiento bajo ese
#   umbral cuenta como coincidencia, sin exigir cobertura minima. Aqui se usa
#   --very-sensitive en lugar de --sensitive, el modo de todo el trabajo, y se
#   guardan qlen y slen para poder estudiar la sensibilidad a la cobertura.
#
# SANEADO
#   DIAMOND rechaza el caracter '.' con que gffread escribe los stops. Se aplica el
#   mismo saneado que en 6.6 del manuscrito: '.' interno -> 'X', '.' terminal
#   eliminado. Los proteomas de referencia no contienen '.'.
#
# LA RELEASE IMPORTA
#   RepeatPeps.lib cambia entre versiones de RepeatMasker. La usada es la de
#   4.1.7-p1 (fichero fechado el 13-09-2024, MD5 9d055c4370ac3a40dfecc2759e2d78d4).
#   El script comprueba el MD5 y aborta si no coincide.
#
# USO
#   conda activate blast          # diamond 2.2.4
#   bash scripts/run_repeatpeps.sh
#
# SALIDA
#   repeatpeps/RepeatPeps.lib, repeatpeps/repeatpeps.dmnd
#   repeatpeps/hits_<conjunto>.tsv   (12 columnas, ver FMT)
#
# Escrito el 3 de septiembre de 2026.

set -euo pipefail

BASE="${BASE:-$HOME/llama_annotation_pipeline}"
ZEN="${ZEN:-$HOME/zenodo_data}"
DRO="${DRO:-/mnt/d/VetAI_Research_Hub/00_COMMON_DB/00_Reference_Genomes/Camelus_dromedarius/proteome.faa}"
ALP="${ALP:-/mnt/d/VetAI_Research_Hub/00_COMMON_DB/00_Reference_Genomes/Vicugna_pacos/proteome.faa}"
THREADS="${THREADS:-4}"
WORK="$BASE/repeatpeps"
mkdir -p "$WORK/q"
cd "$WORK"

RM_URL="https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.7-p1.tar.gz"
RM_MD5="9d055c4370ac3a40dfecc2759e2d78d4"

command -v diamond >/dev/null 2>&1 || { echo "diamond no esta en el PATH"; exit 1; }
echo "diamond: $(diamond version 2>&1 | head -1)"

# 1. Base de proteinas de elementos transponibles
if [ ! -s RepeatPeps.lib ]; then
  echo "Descargando RepeatMasker 4.1.7-p1 (46 MB) para extraer RepeatPeps.lib..."
  curl -sSL -o RepeatMasker-4.1.7-p1.tar.gz "$RM_URL"
  tar -xzf RepeatMasker-4.1.7-p1.tar.gz RepeatMasker/Libraries/RepeatPeps.lib RepeatMasker/Libraries/RepeatPeps.readme
  mv RepeatMasker/Libraries/RepeatPeps.lib RepeatMasker/Libraries/RepeatPeps.readme .
  rm -rf RepeatMasker RepeatMasker-4.1.7-p1.tar.gz
fi
got=$(md5sum RepeatPeps.lib | cut -d' ' -f1)
[ "$got" = "$RM_MD5" ] || { echo "MD5 de RepeatPeps.lib distinto del esperado: $got"; exit 1; }
echo "RepeatPeps.lib: $(grep -c '^>' RepeatPeps.lib) proteinas, MD5 correcto"
[ -s repeatpeps.dmnd ] || diamond makedb --in RepeatPeps.lib -d repeatpeps --quiet

# 2. Consultas saneadas
sanea() {  # '.' interno -> X, '.' terminal fuera
  awk 'BEGIN{seq=""} /^>/{if(seq!=""){s=seq; sub(/\.$/,"",s); gsub(/\./,"X",s); print s}; print; seq=""; next}
       {seq=seq $0} END{if(seq!=""){s=seq; sub(/\.$/,"",s); gsub(/\./,"X",s); print s}}'
}
zcat "$ZEN/llama_liftoff.faa.gz"  | sanea > q/liftoff.faa
zcat "$ZEN/llama_miniprot.faa.gz" | sanea > q/miniprot.faa
zcat "$ZEN/llama_lifton.faa.gz"   | sanea > q/lifton.faa
sanea < "$BASE/Lgla_hx036_helixer.faa"  > q/helixer.faa
zcat "$ZEN/llama_combined_reference_proteome.faa.gz" | sanea > q/combined.faa
cp "$DRO" q/cdro_ref.faa
cp "$ALP" q/vpac_ref.faa
for f in q/*.faa; do printf "%-18s %6d secuencias\n" "$(basename "$f" .faa)" "$(grep -c '^>' "$f")"; done

# 3. Alineamiento
FMT="6 qseqid sseqid pident length evalue bitscore qstart qend sstart send qlen slen"
for s in helixer liftoff miniprot lifton combined vpac_ref cdro_ref; do
  out="hits_${s}.tsv"
  [ -s "$out" ] && { echo "SALTO $s (ya existe)"; continue; }
  diamond blastp -q "q/$s.faa" -d repeatpeps -o "$out" --outfmt $FMT \
    --very-sensitive --evalue 1e-10 --max-target-seqs 1 --max-hsps 1 --threads "$THREADS" --quiet
  printf "%-18s %5d con hit\n" "$s" "$(wc -l < "$out")"
done
echo "Hecho. Resumen: python3 scripts/repeatpeps_summary.py"
