#!/bin/bash
# Estado de la ejecución de DIAMOND blastp. Uso: ~/estado_blastp.sh  o alias  blp
BASE=$HOME/llama_annotation_pipeline
WORK=$BASE/blastp
LOG=$WORK/blastp_run.log

echo "=== DIAMOND blastp  $(date '+%F %H:%M') ==="

if pgrep -f "diamond blastp" >/dev/null; then
    echo "estado   : CORRIENDO"
elif [ -s "$WORK/hits_dromedario.tsv" ] && [ -s "$WORK/hits_alpaca.tsv" ]; then
    echo "estado   : TERMINADO"
else
    echo "estado   : PARADO (relanza ./run_blastp.sh)"
fi

done_n=0
for r in dromedario alpaca; do
    f="$WORK/hits_${r}.tsv"
    if [ -s "$f" ]; then
        printf "%-11s: LISTO  %s hits\n" "$r" "$(wc -l < "$f")"
        done_n=$((done_n + 1))
    elif [ -s "${f}.partial" ]; then
        printf "%-11s: en curso, %s hits escritos\n" "$r" "$(wc -l < "${f}.partial")"
    else
        printf "%-11s: pendiente\n" "$r"
    fi
done
echo "progreso : $done_n de 2 comparaciones"

cur=$(grep "blastp contra" "$LOG" 2>/dev/null | tail -1)
if [ -n "$cur" ]; then
    t=$(echo "$cur" | awk '{print $1" "$2}')
    mins=$(python3 -c "
import datetime,sys
try:
    d=datetime.datetime.strptime(sys.argv[1],'%Y-%m-%d %H:%M:%S')
    print(f'{(datetime.datetime.now()-d).total_seconds()/60:.0f}')
except Exception: print('?')" "$t" 2>/dev/null)
    echo "ultimo   : ${cur#* } (${mins} min)"
fi

echo "cpu      : $(ps -o %cpu= -C diamond 2>/dev/null | awk '{s+=$1} END {printf "%.0f%%", s}')"
echo "disco    : $(df -h "$HOME" | tail -1 | awk '{print $4" libres ("$5")"}')"
echo "--- log ---"
tail -3 "$LOG" 2>/dev/null
