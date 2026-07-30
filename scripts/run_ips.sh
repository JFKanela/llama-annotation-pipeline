#!/bin/bash
IPS=~/interproscan/interproscan-5.78-109.0/interproscan.sh
CPUS=4
if ! ls ips_chunks/chunk_*.faa >/dev/null 2>&1; then
  echo "ERROR: no hay lotes en ips_chunks/. Ejecuta primero el troceado."; exit 1
fi
mkdir -p ips_out ips_tmp
MIN_GB=15

for f in ips_chunks/chunk_*.faa; do
  base=$(basename "$f" .faa)
  out="ips_out/${base}.tsv"

  [ -s "$out" ] && { echo "$(date '+%F %H:%M') SALTO $base"; continue; }

  free=$(df -BG --output=avail "$HOME" | tail -1 | tr -dc '0-9')
  if [ "$free" -lt "$MIN_GB" ]; then
    echo "$(date '+%F %H:%M') ABORTO: solo ${free} GB libres, minimo ${MIN_GB}"
    exit 1
  fi

  echo "$(date '+%F %H:%M') INICIO $base  (${free} GB libres)"
  nice -n 19 ionice -c3 "$IPS" -i "$f" -f TSV -o "$out" \
       -cpu $CPUS -dp -T ips_tmp >> "ips_out/${base}.log" 2>&1

  if [ -s "$out" ]; then
    echo "$(date '+%F %H:%M') OK    $base  ($(wc -l < "$out") lineas)"
    rm -rf ips_tmp/* 2>/dev/null
  else
    echo "$(date '+%F %H:%M') FALLO $base"
  fi
done
echo "$(date '+%F %H:%M') TERMINADO"
