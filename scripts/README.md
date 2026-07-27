# Scripts

Contenido de esta carpeta, agrupado por función. Cada script de análisis lleva en su
cabecera qué produce, qué sostiene en el artículo y sus limitaciones metodológicas.

## 1. Ejecución de las anotaciones

Los brazos de homología Liftoff y miniprot se ejecutan desde el `Snakefile` del
repositorio. LiftOn y Helixer se ejecutaron manualmente, y de ahí estos dos scripts.

| Script | Papel |
|---|---|
| `run_lifton.sh` | Ejecución manual del brazo LiftOn |
| `run_helixer.py` | Ejecución del brazo *ab initio* con `helixerlite` **(SUPERADO, ver aviso)** |

> **AVISO SOBRE `run_helixer.py`.** Este script documenta la primera ejecución de
> Helixer, con `helixerlite 25.5.27`, **sin solapamiento de subsecuencias** y sobre
> sustrato de scaffolds ≥10 kb. **Ese no es el método del manuscrito final.**
>
> La anotación *ab initio* definitiva se obtuvo con la **web tool oficial de Helixer
> v0.3.6** (plabipd.de), con solapamiento activado por defecto para el linaje
> `vertebrate` y sustrato de scaffolds ≥25 kb (244 scaffolds, 79,6 % del ensamblado).
> El cambio elevó la completitud BUSCO de 79,3 % a 85,5 %.
>
> El script se conserva por trazabilidad histórica, no como método vigente. La
> ejecución definitiva no se hizo por script, sino por interfaz web, y sus parámetros
> están documentados en la ficha de Methods del manuscrito.

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

## 3. Anotación funcional

Cadena secuencial:

```
make_interproscan_input.py  ->  run_ips.sh  ->  (reasignación por md5_to_ids.tsv)
```

| Script | Papel |
|---|---|
| `make_interproscan_input.py` | Deduplica por MD5, limpia caracteres no válidos y trocea en lotes |
| `run_ips.sh` | Ejecuta InterProScan por lotes, con reanudación tras interrupciones **(pendiente de subir)** |
| `estado_ips.sh` | Monitorización del progreso **(pendiente de subir)** |

## 4. Informes

| Script | Papel |
|---|---|
| `build_report.py` | Genera las tablas de resumen a partir de las salidas de BUSCO y de las estadísticas estructurales |

> **REVISAR ANTES DE USAR.** Comprobado sobre el propio script (julio 2026):
>
> - Las cifras de BUSCO **no** están codificadas: se leen de los ficheros de
>   `--busco-dir`. Aun así, si `results/busco/` conserva la salida de la ejecución
>   con `helixerlite`, el informe se regenerará con el **79,3 %** antiguo en lugar
>   del **85,5 %** vigente. Hay que sustituir esa carpeta antes de reejecutar.
> - La nota final sobre el sustrato **sí** llevaba valores codificados en el cuerpo
>   del script (`≥10 kb`, 3.640 scaffolds, 81,46 %). Se han actualizado al sustrato
>   vigente (`≥25 kb`, 244 scaffolds, 79,6 %). Si el sustrato vuelve a cambiar, hay
>   que editarla a mano.

## Rutas

Los scripts de análisis esperan encontrar en el directorio de trabajo:

- `llama_liftoff.{faa.gz,gff3.gz}`, `llama_miniprot.*`, `llama_lifton.*`
  (del depósito de datos, DOI 10.5281/zenodo.21445840)
- `Lgla_hx036_helixer.faa` y `Lgla_hx036_helixer_FINAL.gff`
  (web tool de Helixer v0.3.6, con overlap)

Si cambian de ubicación, ajustar los diccionarios `FILES` de cada script.

## Tres cautelas que conviene no perder

**El solapamiento por MD5 no mide novedad.** Mide identidad exacta de secuencia. Las
secuencias exclusivas de Helixer en ese análisis no son genes específicos de llama.
La novedad real la da `overlap_coords.py`, y sale un orden de magnitud menor.

**La fila de Helixer en `internal_stops.py` es no informativa.** Un predictor
*ab initio* construye ORFs por definición y no puede producir stops internos. Solo los
tres brazos de homología son comparables entre sí.

**Las proteínas más largas de Helixer no son inequívocamente mejores.** Puede ser mejor
modelado o sobre-extensión. Distinguirlo exige cobertura BLASTP contra la referencia de
alpaca, que es trabajo pendiente.

---

### Nota sobre este README

Las descripciones de `run_lifton.sh`, `run_helixer.py` y `build_report.py` se
redactaron originalmente a partir del registro de trabajo del proyecto, **no
leyendo los propios scripts**.

Contraste realizado en julio de 2026 sobre `build_report.py` y `run_helixer.py`:
la nota de la sección 4 se ha corregido con lo que el código hace realmente, y
`run_helixer.py` lleva ya el aviso de SUPERSEDED en su cabecera. **`run_lifton.sh`
sigue sin contrastar.**
