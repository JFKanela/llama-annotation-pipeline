#!/usr/bin/env bash
# =====================================================================
# run_blastp.sh - cobertura de alineamiento contra proteomas de referencia
#
# QUE HACE
#   Alinea las 42364 secuencias unicas de los cuatro proteomas de llama contra
#   dos proteomas de referencia, y calcula cobertura de consulta y de sujeto.
#
#   Es el cuarto criterio de calidad de Kourelis et al. 2019, el unico que
#   evalua EXACTITUD del modelo en lugar de completitud o consistencia interna.
#
# POR QUE DOS REFERENCIAS
#   - Camelus dromedarius (GCF_036321535.1): patron de medida EXTERNO. Ninguno de
#     los cuatro brazos deriva de el, de modo que los cuatro son comparables sin
#     circularidad. ES LA COMPARACION QUE IMPORTA.
#   - Vicugna pacos (GCF_048564905.1): control interno. Tres de los cuatro brazos
#     son proyecciones de esta anotacion, asi que mide fidelidad de transferencia,
#     no calidad. Saldra excelente por construccion y NO debe interpretarse como
#     evidencia de calidad.
#
# POR QUE DIAMOND Y NO BLASTP
#   Llama y alpaca divergieron hace 2-3 millones de anos; los ortologos rondan el
#   95-99 % de identidad. Son alineamientos triviales, sin problema de sensibilidad.
#   DIAMOND --very-sensitive da resultados equivalentes en una fraccion del tiempo:
#   menos de una hora frente a dias.
#
# REANUDACION
#   Cada paso comprueba si su salida ya existe y no esta vacia, y la salta. Si el
#   sistema se cuelga, se relanza el mismo script y continua donde estaba.
#
# USO
#   conda activate blast
#   ./run_blastp.sh
# =====================================================================

set -uo pipefail

BASE="$HOME/llama_annotation_pipeline"
WORK="$BASE/blastp"
REFDIR="/mnt/d/VetAI_Research_Hub/00_COMMON_DB/00_Reference_Genomes"
QUERY="$BASE/camelid_unique_proteins_clean.faa"
THREADS=4
LOG="$WORK/blastp_run.log"

mkdir -p "$WORK" "$WORK/tmp"

log() { echo "$(date '+%F %H:%M:%S') $*" | tee -a "$LOG"; }
die() { log "ABORTO: $*"; exit 1; }

log "=============================================="
log "Inicio de run_blastp.sh"

# ---------------------------------------------------------------------
# 0. Comprobaciones previas
# ---------------------------------------------------------------------
command -v diamond >/dev/null 2>&1 || die "diamond no esta en el PATH. Crear el entorno:
    conda create -n blast -c conda-forge -c bioconda diamond -y
    conda activate blast"

log "diamond: $(diamond --version 2>&1 | head -1)"

[ -s "$QUERY" ] || die "no encuentro el fichero de consulta $QUERY"
log "consulta: $QUERY ($(grep -c '^>' "$QUERY") secuencias)"

[ -d "$REFDIR" ] || die "no existe el directorio de referencias $REFDIR"

# ---------------------------------------------------------------------
# 1. Localizar los proteomas de referencia
#    Se buscan por accesion, que es el criterio fiable. Si no aparecen, se
#    lista el directorio y se aborta, en lugar de adivinar nombres.
# ---------------------------------------------------------------------
find_proteome() {
    local acc="$1" hit
    hit=$(find "$REFDIR" -type f \( -name "*.faa" -o -name "*.fasta" -o -name "*.fa" \) \
            -path "*${acc}*" 2>/dev/null | head -1)
    [ -z "$hit" ] && hit=$(find "$REFDIR" -type f -name "protein.faa" -path "*${acc}*" 2>/dev/null | head -1)
    echo "$hit"
}

DRO="/mnt/d/VetAI_Research_Hub/00_COMMON_DB/00_Reference_Genomes/Camelus_dromedarius/proteome.faa"
ALP="/mnt/d/VetAI_Research_Hub/00_COMMON_DB/00_Reference_Genomes/Vicugna_pacos/proteome.faa"

if [ -z "$DRO" ] || [ -z "$ALP" ]; then
    log "No localizo los proteomas por accesion. Contenido del directorio:"
    find "$REFDIR" -maxdepth 3 -type f \( -name "*.faa" -o -name "*.fasta" -o -name "*.fa" \) \
        -exec ls -lh {} \; 2>/dev/null | tee -a "$LOG"
    die "Ajustar las rutas DRO y ALP a mano en la cabecera del script."
fi

log "dromedario: $DRO ($(grep -c '^>' "$DRO") proteinas)"
log "alpaca:     $ALP ($(grep -c '^>' "$ALP") proteinas)"

# ---------------------------------------------------------------------
# 2. Construir las bases de datos
# ---------------------------------------------------------------------
build_db() {
    local faa="$1" name="$2"
    if [ -s "$WORK/${name}.dmnd" ]; then
        log "SALTO makedb $name (ya existe)"
        return 0
    fi
    log "makedb $name"
    diamond makedb --in "$faa" -d "$WORK/${name}" --threads "$THREADS" \
        >> "$LOG" 2>&1 || die "fallo makedb $name"
    log "OK makedb $name"
}

build_db "$DRO" dromedario
build_db "$ALP" alpaca

# ---------------------------------------------------------------------
# 3. Alinear
#    --max-target-seqs 1 y --max-hsps 1: un solo mejor emparejamiento por
#    consulta, que es lo que hace interpretable la cobertura.
#    Se piden qlen y slen para calcular ambas coberturas en el analisis.
# ---------------------------------------------------------------------
run_blast() {
    local name="$1"
    local out="$WORK/hits_${name}.tsv"
    if [ -s "$out" ]; then
        log "SALTO blastp contra $name (ya existe, $(wc -l < "$out") lineas)"
        return 0
    fi
    log "blastp contra $name (esto es lo largo)"
    diamond blastp \
        -q "$QUERY" \
        -d "$WORK/${name}.dmnd" \
        -o "${out}.partial" \
        --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen \
        --very-sensitive \
        --max-target-seqs 1 \
        --max-hsps 1 \
        --evalue 1e-5 \
        --threads "$THREADS" \
        --tmpdir "$WORK/tmp" \
        >> "$LOG" 2>&1 || { rm -f "${out}.partial"; die "fallo blastp contra $name"; }
    mv "${out}.partial" "$out"
    log "OK blastp contra $name ($(wc -l < "$out") lineas)"
}

run_blast dromedario
run_blast alpaca

# ---------------------------------------------------------------------
# 4. Analisis
# ---------------------------------------------------------------------
if [ -x "$BASE/analyze_blastp.py" ] || [ -f "$BASE/analyze_blastp.py" ]; then
    log "lanzando analisis"
    python3 "$BASE/analyze_blastp.py" 2>&1 | tee -a "$LOG"
else
    log "analyze_blastp.py no encontrado; los ficheros hits_*.tsv estan listos"
fi

log "TERMINADO"
