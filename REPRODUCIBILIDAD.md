# Reproducibilidad y ejecución del pipeline

Este documento describe con precisión **qué parte del análisis está automatizada
en Snakemake y qué parte se ejecutó manualmente**, con los comandos exactos, de
modo que cualquier persona pueda reproducir el resultado completo.

La distinción es deliberada y se declara abiertamente: el pipeline de Snakemake
no reproduce los cuatro brazos de anotación de principio a fin. Dos de ellos
(LiftOn y Helixer) se incorporaron al proyecto después de escribir el flujo de
trabajo y se ejecutaron con los comandos que aquí se documentan.

---

## 1. Resumen del alcance

| Componente | Automatizado en Snakemake | Ejecución |
|---|---|---|
| Descarga de genomas y anotación de referencia | Sí | `snakemake` |
| Brazo Liftoff (homología ADN-ADN) | Sí | `snakemake` |
| Brazo miniprot (homología proteína-ADN) | Sí | `snakemake` |
| Brazo LiftOn (homología híbrida) | **No** | Manual, sección 3 |
| Brazo Helixer (*ab initio*) | **No** | Externa (GPU), sección 4 |
| Extracción de proteomas (gffread) | Parcial | Manual para LiftOn y Helixer |
| Control de calidad BUSCO | Sí para los brazos automatizados | Manual para LiftOn y Helixer |
| Informe comparativo | Sí | `scripts/build_report.py` |

`scripts/build_report.py` descubre los métodos escaneando `results/busco/`, de
modo que incorpora al informe cualquier brazo presente, esté o no automatizado.

---

## 2. Entorno base y ejecución automatizada

```bash
conda activate smk        # Snakemake 9.23.1
snakemake --use-conda --cores 6
```

Los entornos por herramienta se declaran en `envs/` y Snakemake los construye de
forma aislada. El control de calidad emplea BUSCO 6.1.0 con el linaje
`artiodactyla_odb12` (n = 12.594). Conviene advertir que el linaje
`cetartiodactyla_odb12` no existe en OrthoDB v12: el clado fue renombrado.

Sobre requisitos de memoria: con 20 GB de RAM, la ejecución concurrente de
varios BUSCO puede agotar la memoria y provocar la intervención del *OOM
killer*. Se recomienda `--cores 6` en lugar de 8 y disponer de espacio de
intercambio (en el sistema de desarrollo se configuró un `swapfile` de 16 GB).

---

## 3. Brazo LiftOn (ejecución manual)

### 3.1. Construcción del entorno

LiftOn 1.0.9 no está disponible en bioconda, solo en PyPI, y su instalación
requiere resolver una cadena de dependencias en un orden concreto:

```bash
conda create -n lifton -c conda-forge -c bioconda python=3.10 liftoff miniprot -y
conda activate lifton

# 1. cigar: empaquetado obsoleto que requiere pkg_resources
pip install --no-cache-dir setuptools-scm pkg_resources
pip install --no-cache-dir --no-build-isolation cigar

# 2. mappy: no compila con gcc reciente, se instala por conda
conda install -n lifton -c bioconda -c conda-forge mappy -y

# 3. LiftOn sin dependencias, ya resueltas arriba
pip install --no-cache-dir --no-build-isolation --no-deps lifton

# 4. dependencias restantes
pip install --no-cache-dir intervaltree "duckdb>=1.0" "pyarrow>=14"
```

### 3.2. Anotación de referencia en formato nativo

LiftOn requiere el GFF **nativo de NCBI**. El uso de una anotación previamente
procesada con AGAT produce un proteoma inválido (véase la nota 3.4). Descarga:

```bash
BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/564/905"
DIR=$(curl -s "$BASE/" | grep -oP 'GCF_048564905\.1_[^/"]+' | head -1)
wget "${BASE}/${DIR}/${DIR}_genomic.gff.gz" -O resources/reference/annotation_ncbi.gff.gz
gunzip resources/reference/annotation_ncbi.gff.gz
```

### 3.3. Ejecución

LiftOn reutiliza las anotaciones de Liftoff y miniprot ya generadas por el
pipeline (opciones `-L` y `-M`), lo que reduce sustancialmente el tiempo de
cálculo. Se emplean los parámetros por defecto: rescate mediante miniprot,
selección del mejor resultado con verificación y transferencia de elementos
génicos. No se utiliza `--legacy-merge`.

```bash
conda activate lifton
lifton \
    -g resources/reference/annotation_ncbi.gff \
    -ad RefSeq \
    -L results/liftoff/llama_liftoff.gff3 \
    -M results/miniprot/llama_miniprot.gff3 \
    -o results/lifton/llama_lifton_v2.gff3 \
    -t 6 \
    resources/target/genome.fasta \
    resources/reference/genome.fasta
```

Obsérvese que el genoma diana precede al de referencia en la línea de órdenes.

### 3.4. Nota sobre una ejecución fallida previa

Una primera ejecución empleó como referencia la anotación reducida a isoforma
primaria con AGAT. LiftOn generó 33.045 modelos génicos con coordenadas
aparentemente correctas, pero solo 2.533 proteínas traducibles, junto a 575.167
avisos de validación del GFF de entrada. El problema se resolvió íntegramente
utilizando el GFF nativo de NCBI. Se documenta aquí porque el fallo es
silencioso: el GFF resultante parece correcto y solo el recuento de proteínas
revela el problema.

---

## 4. Brazo Helixer (ejecución externa en GPU)

Helixer se ejecutó en infraestructura externa (Kaggle, GPU Tesla P100 de 16 GB)
porque el equipo de desarrollo dispone de una GPU con 4 GB de VRAM, insuficiente
para la inferencia. Por esta razón **no se automatiza en Snakemake**: hacerlo
sugeriría una reproducibilidad que en la práctica no existe, ya que requiere
hardware y una cuenta de plataforma específicos. El GFF3 resultante se trata
como entrada externa documentada y se aporta en el repositorio.

El script de ejecución es `scripts/run_helixer.py`.

### 4.1. Entorno

```bash
uv venv --python 3.11 /tmp/hxenv
uv pip install --python /tmp/hxenv/bin/python helixerlite "tensorflow[and-cuda]"
```

Modelo: `vertebrate_v0.3_m_0080.h5` (Zenodo), longitud de subsecuencia 213.840.
Se emplea la API de Python de `helixerlite` (`fasta2hdf5`, `HybridModel`,
`preds2gff3`) junto con `gfftk` para la conversión final.

### 4.2. Sustrato y parámetros

Se retuvieron los scaffolds de longitud igual o superior a 10 kb: **3.640
scaffolds** (0,34 % del total) que reúnen **1.915.763.599 bp**, el **81,46 %** de
la secuencia ensamblada. Los brazos de homología se ejecutaron sobre el
ensamblado completo (2.351.761.190 bp). Esta asimetría de sustrato condiciona
toda comparación directa de completitud y debe tenerse presente al interpretar
los resultados.

Parámetros: sin solapamiento, tamaño de lote 16. El sustrato se procesó en siete
bloques consecutivos de 300 Mb (`CHUNK_BP = 300_000_000`) debido al límite de
tiempo de sesión de la plataforma (unas 12 horas). Los ficheros intermedios se
escribieron en almacenamiento temporal amplio y los GFF3 parciales en un
directorio de trabajo persistente, de modo que la ejecución fuera reanudable.

### 4.3. Reanudación entre sesiones

La ejecución completa (unas 14 horas) excede el límite de sesión. Al reanudar,
los bloques ya calculados deben recuperarse del resultado de la sesión anterior,
incorporado como conjunto de datos de entrada:

```python
import os, shutil, glob
os.makedirs("/kaggle/working/hxchunks", exist_ok=True)
for f in glob.glob("/kaggle/input/**/chunk*.gff3", recursive=True):
    dst = "/kaggle/working/hxchunks/" + os.path.basename(f)
    if not os.path.exists(dst):
        shutil.copy(f, dst)
```

El script detecta los bloques presentes y solo procesa los pendientes.

### 4.4. Verificación de integridad del GFF3

La concatenación de los siete bloques se verificó antes de su uso:

```bash
# los siete bloques representados
awk -F'\t' '!/^#/ && $3=="gene"' llama_helixer.gff3 | grep -oP 'ID=c\d+' | sort | uniq -c
# ausencia de identificadores duplicados
awk -F'\t' '!/^#/ && $3=="gene"' llama_helixer.gff3 | grep -oP 'ID=[^;]+' | sort | uniq -d | wc -l
# features con coordenadas inválidas, por tipo
awk -F'\t' '!/^#/ && $5 < $4 {c[$3]++} END{for(f in c) print f, c[f]}' llama_helixer.gff3
```

Resultado: los siete bloques presentes, cero identificadores duplicados y 42
features de tipo exón con coordenadas inválidas (fin anterior a inicio), un
artefacto conocido de la conversión a GFF3. Ninguna CDS resultó afectada, por lo
que el proteoma no se vio comprometido: `gffread` descarta esas líneas y traduce
a partir de las CDS.

---

## 5. Extracción de proteomas y control de calidad manuales

Para los dos brazos no automatizados, tras obtener el GFF3:

```bash
conda activate smk

# LiftOn
gffread -y results/proteomes/llama_lifton_v2.faa \
        -g resources/target/genome.fasta \
        results/lifton/llama_lifton_v2.gff3

# Helixer
gffread -y results/proteomes/llama_helixer.faa \
        -g resources/target/genome.fasta \
        results/helixer/llama_helixer.gff3
```

El control de calidad debe emplear **el mismo entorno de BUSCO** que construyó
Snakemake, para garantizar identidad de versión y linaje. Su localización:

```bash
for d in .snakemake/conda/*/; do
  if [ -x "${d}bin/busco" ]; then echo "BUSCO en: $d"; fi
done
```

Y la ejecución, sustituyendo `RUTA` por el resultado anterior y `METODO` por
`lifton` o `helixer`:

```bash
conda run -p RUTA busco \
    -i results/proteomes/llama_METODO.faa \
    -m proteins -l artiodactyla_odb12 -c 6 \
    -o METODO --out_path results/busco
```

Cada ejecución requiere en torno a 100 minutos con seis núcleos. Conviene
lanzarla en una sesión persistente (`tmux`). BUSCO escribe el resumen con nombre
extendido (`short_summary.specific.artiodactyla_odb12.METODO.txt`);
`build_report.py` contempla todas las variantes de nombre.

---

## 6. Generación del informe comparativo

```bash
snakemake --use-conda --cores 4 results/report/comparison_report.md
```

O directamente:

```bash
python scripts/build_report.py \
    --busco-dir results/busco \
    --agat-dir results/agat \
    --gffcompare results/gffcompare/cmp.stats \
    --unmapped results/liftoff/unmapped.txt \
    --out-md results/report/comparison_report.md \
    --out-tsv results/report/comparison_table.tsv
```

El script no recalcula nada: únicamente lee los ficheros ya producidos y compone
la tabla comparativa a partir de ellos.

---

## 7. Datos de partida

| Recurso | Identificador |
|---|---|
| Genoma diana (*Lama glama*) | GCA_028534125.1 (DNA Zoo, Hi-C, ejemplar Fiesta) |
| Genoma y anotación de referencia (*Vicugna pacos*) | GCF_048564905.1 (VicPac4, RefSeq RS_2025_04) |
| Linaje BUSCO | `artiodactyla_odb12` |
| Modelo Helixer | `vertebrate_v0.3_m_0080.h5` |

No se dispuso de evidencia de RNA-seq: los registros públicos disponibles para
*Lama glama* corresponden a amplicones de VHH o nanocuerpos y no constituyen
datos de transcriptoma completo. Esta es la razón de que la estrategia principal
sea la anotación por homología.
