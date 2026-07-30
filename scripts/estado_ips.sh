#!/bin/bash
BASE=$HOME/llama_annotation_pipeline
DONE=$(ls "$BASE"/ips_out/chunk_*.tsv 2>/dev/null | wc -l)
TOTAL=$(ls "$BASE"/ips_chunks/chunk_*.faa 2>/dev/null | wc -l)
LEFT=$((TOTAL - DONE))

echo "=== InterProScan  $(date '+%F %H:%M') ==="
pgrep -f interproscan >/dev/null && echo "estado   : CORRIENDO" \
                                 || echo "estado   : PARADO (relanza ./run_ips.sh)"
echo "progreso : $DONE de $TOTAL lotes  ($((DONE*100/TOTAL))%)"
awk -v l="$LEFT" 'BEGIN{printf "restante : %.1f horas (%.1f dias)\n", l*80/60, l*80/1440}'
echo "disco    : $(df -h "$HOME" | tail -1 | awk '{print $4" libres ("$5")"}')"

CUR_LOG=$(ls -t "$BASE"/ips_out/chunk_*.log 2>/dev/null | head -1)
if [ -n "$CUR_LOG" ]; then
  NAME=$(basename "$CUR_LOG" .log)
  PCT=$(grep -oE '[0-9]+% completed' "$CUR_LOG" | tail -1)
  ELAPSED=$(python3 -c 'import sys,re,datetime
t=re.findall(r"(\d{2}/\d{2}/\d{4} \d{2}:\d{2}:\d{2})",open(sys.argv[1],errors="ignore").read())
if t:
    d=datetime.datetime.strptime(t[0],"%d/%m/%Y %H:%M:%S")
    print(f"{(datetime.datetime.now()-d).total_seconds()/60:.0f} min")' "$CUR_LOG")
  echo "lote     : $NAME  ${PCT:-arrancando}  (${ELAPSED:-?} en curso)"
fi
