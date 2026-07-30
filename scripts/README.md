# Scripts

Contenido de esta carpeta, agrupado por función. Cada script de análisis lleva en su
cabecera qué produce, qué sostiene en el artículo y sus limitaciones metodológicas.

Para el detalle de ejecución, ver `REPRODUCIBILITY.md` en la raíz del repositorio.

## 1. Ejecución de las anotaciones

Los brazos Liftoff y miniprot se ejecutan desde el `Snakefile`. LiftOn y Helixer se
añadieron después del workflow y se ejecutaron manualmente, de ahí estos dos scripts.

| Script | Papel |
|---|---|
| `run_lifton.sh` | Brazo LiftOn (homología híbrida). Encapsula los comandos de `REPRODUCIBILITY.md` sección 3. No está integrado en el `Snakefile`: se añadió después y reutiliza las anotaciones de Liftoff y miniprot ya producidas. Requiere los entornos conda `lifton` (LiftOn 1.0.9) y `smk` (gffread) |
| `run_helixer.py` | Brazo *ab initio* en GPU externa (Kaggle, Tesla P100), en entorno `uv` aislado. `helixerlite 25.5.27`, modelo `vertebrate_v0.3_m_0080.h5` (DOI 10.5281/zenodo.10836346) **(SUPERADO, ver aviso)** |

> **ASIMETRÍA DE REFERENCIA ENTRE BRAZOS DE HOMOLOGÍA.** `run_lifton.sh` documenta en su
> cabecera un punto que hay que leer antes de comparar cifras de completitud: LiftOn
> consume el **GFF3 nativo de NCBI**, mientras que Liftoff y miniprot reciben la
> anotación **reducida por AGAT a un transcrito por gen**, que es la que produce el
> workflow con `keep_longest_isoform: true`. Los tres brazos **no comparten la misma
> anotación de referencia**: difieren en profundidad de isoformas, no en repertorio
> génico. Los recuentos exactos están en la sección correspondiente del `README.md` de
> la raíz.

> **AVISO SOBRE `run_helixer.py`.** Documenta la primera ejecución de Helixer, con
> scaffolds ≥10 kb, **sin solapamiento de subsecuencias** y cargando el modelo por ruta
> en lugar de por linaje. Su propia cabecera advierte de que esas decisiones hacen que
> la completitud reportada sea una **cota inferior**, no una medida del método.
>
> **Ese no es el método del manuscrito final.** La anotación definitiva se obtuvo con la
> **web tool oficial de Helixer v0.3.6** (plabipd.de), con solapamiento activado por
> defecto para el linaje `vertebrate` y sustrato de scaffolds ≥25 kb (244 scaffolds,
> 79,6 % del ensamblado). El cambio elevó la completitud BUSCO de 79,3 % a 85,5 %.
>
> El script se conserva por trazabilidad histórica. La ejecución definitiva no se hizo
> por script sino por interfaz web, y sus parámetros están documentados en la sección
> 4.A de `REPRODUCIBILITY.md`.
>
> Cita de Helixer: Nature Methods (2025), DOI 10.1038/s41592-025-02939-1.

## 2. Análisis que producen cifras del artículo

| Script | Salida | Dónde va en el artículo |
|---|---|---|
| `overlap_md5.py` | 42.364 únicas de 80.331; 14 intersecciones | Figura UpSet de solapamiento. Redundancia interna de miniprot |
| `overlap_coords.py` | 828 loci novel, 3.477 perdidos | La afirmación de novedad posicional de Helixer |
| `structural_stats.py` | Longitud, CDS/transcrito, monoexónicos | Eje de calidad independiente de BUSCO |
| `internal_stops.py` | 12,0 % / 8,8 % / 2,7 % / 0 % | Justifica LiftOn como proteoma de referencia |
| `camelidae_landscape.py` | `camelidae.json` y tabla | Evidencia de la afirmación central. Especies para OrthoFinder |

Son independientes entre sí y pueden correrse en cualquier orden, siempre que existan
los cuatro proteomas y sus GFF3.

## 3. Evaluación contra referencias externas

| Script | Papel |
|---|---|
| `run_blastp.sh` | Alinea las secuencias únicas contra dos proteomas de referencia con DIAMOND, con reanudación por pasos |
| `analyze_blastp.py` | Calcula cobertura de consulta y de sujeto por brazo, y reparte los resultados con `md5_to_ids.tsv` |
| `estado_blastp.sh` | Monitorización del progreso |

Cuarto criterio de calidad de Kourelis et al. 2019, y el único de todo el trabajo que
evalúa **exactitud del modelo** en lugar de completitud o consistencia interna.

> **POR QUÉ DOS REFERENCIAS, Y POR QUÉ IMPORTA.** *Camelus dromedarius*
> (`GCF_036321535.1`) es el patrón de medida válido: ninguno de los cuatro brazos deriva
> de él, de modo que los cuatro son comparables. *Vicugna pacos* (`GCF_048564905.1`) es
> un control **circular**: tres de los cuatro brazos son proyecciones de esa anotación y
> alinean casi perfectamente contra su propia fuente.
>
> La diferencia no es teórica. Contra alpaca, Helixer parece sobre-extender sus modelos
> veinte veces más que LiftOn (2,0 % frente a 0,0 %). Contra dromedario, los cuatro
> brazos están en el mismo 2 % y son indistinguibles. **Usar alpaca como patrón habría
> producido una conclusión falsa en el manuscrito.**

> **NO REPORTAR LAS MEDIANAS DE COBERTURA.** Dan 100 % en los cuatro brazos y en ambas
> referencias: saturan y no discriminan. Usar las fracciones con cobertura ≥80 %.

## 4. Anotación funcional

Cadena secuencial:

```
make_interproscan_input.py  ->  run_ips.sh  ->  (reasignación por md5_to_ids.tsv)
```

| Script | Papel |
|---|---|
| `make_interproscan_input.py` | Deduplica por MD5, limpia caracteres no válidos y trocea en lotes |
| `run_ips.sh` | Ejecuta InterProScan por lotes, con reanudación tras interrupciones |
| `estado_ips.sh` | Monitorización del progreso |

Ejecución realizada: InterProScan 5.78-109.0, los 18 análisis por defecto, sobre 42.364
secuencias únicas en 43 lotes de 1.000, con el servicio de coincidencias precalculadas
**desactivado** (`-dp`), porque las proteínas de llama no están en UniProtKB. Requiere
OpenJDK 11: la versión 17 no sirve.

> **PANTHER se reporta a nivel de familia.** El formato TSV no emite subfamilia
> (`PTHR12345:SF6`) en ningún caso, ni siquiera cuando el emplazamiento filogenético
> funciona. Comprobado con el fichero de prueba oficial, que da cero subfamilias en TSV y
> las cinco esperadas en XML. La resolución de subfamilia exige formatos estructurados.

## 5. Informes

| Script | Papel |
|---|---|
| `build_report.py` | Invocado por la regla `report` del `Snakefile`. Lee las salidas de BUSCO, AGAT, gffcompare y los genes no transferidos por Liftoff, y compone las tablas de resumen |

Interfaz:

```
python scripts/build_report.py \
    --busco-dir results/busco \
    --agat-dir results/agat \
    --gffcompare results/gffcompare/cmp.stats \
    --unmapped results/liftoff/unmapped.txt \
    --out-md <md> --out-tsv <tsv>
```

Su diseño es deliberadamente conservador: no recalcula nada, descubre los métodos
escaneando `--busco-dir` (de modo que brazos añadidos después aparecen sin tocar el
código), cuenta proteínas del `.faa` real evaluado por BUSCO en lugar de mRNA del GFF, y
tolera ficheros ausentes.

> **DOS RIESGOS DISTINTOS AL REGENERAR EL INFORME.** Comprobados sobre el código en la
> ronda 5 de correcciones, no supuestos.
>
> **1. Las cifras de BUSCO NO están codificadas**: se parsean de los ficheros de
> `--busco-dir`. Pero `results/busco/` conserva la salida de la ejecución con
> `helixerlite` (79,3 %). Hay que sustituir esa carpeta por la de la ejecución vigente
> (85,5 %) antes de regenerar, o la tabla saldrá con la cifra antigua.
>
> **2. El script SÍ lleva otros valores codificados.** La nota final sobre el sustrato
> contenía en duro «10 kb», «3.640 scaffolds», «0,34 %», «1.915.763.599 bp» y «81,46 %».
> Se corrigió en la ronda 5 a ≥25 kb, 244 scaffolds y 79,6 %. Antes de regenerar el
> informe conviene revisar si queda algún otro valor fijo que haya quedado obsoleto.

## Rutas

Los scripts de análisis esperan encontrar en el directorio de trabajo:

- `llama_liftoff.{faa.gz,gff3.gz}`, `llama_miniprot.*`, `llama_lifton.*`
  (del depósito de datos, DOI 10.5281/zenodo.21445840)
- `Lgla_hx036_helixer.faa` y `Lgla_hx036_helixer_FINAL.gff`
  (web tool de Helixer v0.3.6, con overlap)

Si cambian de ubicación, ajustar los diccionarios `FILES` de cada script.

## Tres cautelas que conviene no perder

**El solapamiento por MD5 no mide novedad.** Mide identidad exacta de secuencia. Las
secuencias exclusivas de Helixer en ese análisis no son genes específicos de llama. La
novedad real la da `overlap_coords.py`, y sale un orden de magnitud menor.

**La fila de Helixer en `internal_stops.py` es no informativa.** Un predictor
*ab initio* construye ORFs por definición y no puede producir stops internos. Solo los
tres brazos de homología son comparables entre sí.

**Las proteínas más largas de Helixer no son inequívocamente mejores.** Puede ser mejor
modelado o sobre-extensión. Distinguirlo exigía cobertura por alineamiento contra una
referencia, y **ya está hecho**: es la sección 3. La respuesta, sin embargo, invierte el
planteamiento original. Medido contra **alpaca** —que es lo que aquí se proponía— Helixer
parecería sobre-extender, pero esa lectura es un artefacto de circularidad. Medido contra
**dromedario**, que es la referencia válida, los cuatro brazos son indistinguibles. La
referencia frente a la que se mide no es un detalle de implementación: decide el signo de
la conclusión.
