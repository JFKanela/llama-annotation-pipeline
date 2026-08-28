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
| `overlap_coords.py` | 828 **loci** novel (intervalos fusionados), 3.477 perdidos | La afirmación de novedad posicional de Helixer |
| `structural_stats.py` | Longitud, CDS/transcrito, monoexónicos | Eje de calidad independiente de BUSCO |
| `novel_loci_blast.py` | 891 **transcritos** sin solape, 790 de >1 kb; de esos, **555** con ortólogo camélido y 235 sin él | La caracterización de la novedad posicional de Helixer, y los 555 que entran en el proteoma combinado |
| `internal_stops.py` | 12,0 % / 8,8 % / 2,7 % / 0 % | Justifica LiftOn como proteoma de referencia |
| `camelidae_landscape.py` | `camelidae.json` y tabla | Evidencia de la afirmación central. Especies para OrthoFinder |

> **828 Y 891 NO SE CONTRADICEN: CUENTAN COSAS DISTINTAS.** `overlap_coords.py`
> fusiona los intervalos solapados antes de comparar, de modo que dos mRNA superpuestos
> son **un locus**. `novel_loci_blast.py` evalúa **cada mRNA por separado**, porque
> necesita la correspondencia con una proteína. Loci frente a transcritos. Ninguna de
> las dos cifras se recalculó en la ronda 8: ambas proceden de la salida en ejecución de
> su propio script. La que sí se verificó contando sobre los ficheros es la de 790 y su
> desglose 555 / 235 / 10 / 225.

> **CORRECCIÓN DE LA RONDA 8.** Hasta el commit `06eb6be` aquí, en el `README.md` de la
> raíz y en `REPRODUCIBILITY.md` §8.3 se publicaba **545 con ortólogo y 245 sin él**.
> Es incorrecto. El recuento sobre `helixer_novel_loci.tsv` da **555 y 235**, y solo con
> 555 cuadra el proteoma combinado: 20.233 + 555 − 80 = 20.708, que son las proteínas
> que efectivamente tiene el fichero depositado. El error pasó tres rondas porque
> 545+245 también suma 790.

Son independientes entre sí y pueden correrse en cualquier orden, siempre que existan
los cuatro proteomas y sus GFF3. `novel_loci_blast.py` es la excepción: necesita que
`run_blastp.sh` haya producido antes las tablas de hits, porque **no realinea nada**,
solo cruza tablas existentes. Escribe además el FASTA de los loci huérfanos, que es
lo que después se alineó manualmente contra Swiss-Prot.

> **BUG CORREGIDO EN `structural_stats.py` (14 de agosto de 2026).** Contaba CDS y
> exones agrupando por el atributo `Parent` **sin comprobar que el padre fuera un
> mRNA**. En los GFF3 de homología el `exon` cuelga también de `lnc_RNA`, `tRNA`,
> `snRNA`, `snoRNA` y `rRNA`: en Liftoff hay **34.545 padres de exon y solo 20.306 son
> mRNA**, y los monoexónicos pasaban de **2.516 a 6.009**.
>
> Lo grave no es la magnitud sino la **asimetría**: Helixer no emite features no
> codificantes, así que la contaminación afectaba solo a los tres brazos de homología,
> que es justo el eje sobre el que se comparan. En CDS la contaminación era menor
> (49 padres de 20.075) y las cifras publicadas apenas cambian. Mismo fallo corregido
> en `fig3_estructura.py`.
>
> **Cualquier cifra estructural anterior al 14 de agosto de 2026 hay que recalcularla.**

## 3. Evaluación contra referencias externas

| Script | Papel |
|---|---|
| `run_blastp.sh` | Alinea las secuencias únicas contra dos proteomas de referencia con DIAMOND, con reanudación por pasos |
| `analyze_blastp.py` | Calcula cobertura de consulta y de sujeto por brazo, y reparte los resultados con `md5_to_ids.tsv` |
| `estado_blastp.sh` | Monitorización del progreso |
| `run_swissprot.sh` | Alinea contra Swiss-Prot los loci huérfanos que escribe `novel_loci_blast.py`. Reconstruido a posteriori; **sus parámetros se corrigieron en la ronda 9**, ver aviso |

Cuarto criterio de calidad de Kourelis et al. 2019, y el único de todo el trabajo que
evalúa **exactitud del modelo** en lugar de completitud o consistencia interna.

> **POR QUÉ DOS REFERENCIAS, Y POR QUÉ IMPORTA.** *Camelus dromedarius*
> (`GCF_036321535.1`) es el patrón de medida válido: ninguno de los cuatro brazos deriva
> de él, de modo que los cuatro son comparables. *Vicugna pacos* (`GCF_048564905.1`) es
> un control **circular**: tres de los cuatro brazos son proyecciones de esa anotación y
> alinean casi perfectamente contra su propia fuente.
>
> La diferencia no es teórica. El criterio de sobre-extensión es el de
> `analyze_blastp.py`: cobertura de sujeto por encima de la de consulta en más de 20
> puntos. Contra alpaca, Helixer lo cumple en **364 de 18.016** proteínas alineadas
> (**2,0 %**) frente a **7 de 19.726** de LiftOn (**0,04 %**, que se imprime como 0,0 %
> a un decimal): **más de cincuenta veces**. Contra dromedario los cuatro brazos quedan
> entre **1,7 % y 2,2 %** y son indistinguibles. **Usar alpaca como patrón habría
> producido una conclusión falsa en el manuscrito.**
>
> *(Hasta la ronda 8 aquí ponía «veinte veces». La proporción real es de unas 57. Las
> dos cifras porcentuales sí eran correctas; el multiplicador no.)*

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

> **DOS BUGS CORREGIDOS EL 14 DE AGOSTO DE 2026, ambos con efecto sobre cifras.**
>
> **1. El `glob` de BUSCO dejaba fuera dos brazos.** Buscaba solo `short_summary.txt`
> sin sufijo, y las carpetas de `lifton` y `helixer` solo contienen
> `short_summary.specific.<linaje>.<nombre>.txt`. **Dos de los cuatro brazos no
> llegaban al informe.** Ahora se aceptan ambos patrones sin duplicar.
>
> **2. El parser de AGAT leía el bloque equivocado.** AGAT emite un bloque por tipo de
> feature; el bucle sobrescribía y se quedaba con el **último**, `v_gene_segment`,
> reportando **31 genes para Liftoff en lugar de 20.306**. Ahora rastrea el bloque
> activo y solo lee `mrna`.
>
> El pie del informe se reescribió además para dejar constancia de que el valor del
> recurso no es el descubrimiento de genes nuevos, y lleva la instrucción explícita de
> **no usar la palabra «primer»**, porque existe `chaku_v1` como línea de base.

> **DOS RIESGOS DISTINTOS AL REGENERAR EL INFORME.** Comprobados sobre el código en la
> ronda 5 y repuestos en la ronda 7 sobre la versión reescrita, que los había perdido.
>
> **1. Las cifras de BUSCO NO están codificadas**: se parsean de los ficheros de
> `--busco-dir`. Pero `results/busco/` conserva la salida de la ejecución con
> `helixerlite` (79,3 %). Hay que sustituir esa carpeta por la de la ejecución vigente
> (85,5 %) antes de regenerar, o la tabla saldrá con la cifra antigua.
>
> **2. El script SÍ lleva otros valores codificados.** La nota final sobre el sustrato
> lleva en duro «≥25 kb», «244 scaffolds» y «79,6 %». Si el sustrato vuelve a cambiar,
> hay que editarla a mano.

## 6. Figuras y tablas del manuscrito

Escriben en `MANUSCRITO/01_figuras/` (PDF y PNG a 300 dpi) y `MANUSCRITO/02_tablas/`.

| Script | Salida |
|---|---|
| `fig1_busco.py` | Figura 1 · completitud BUSCO de los cuatro proteomas, con techo y línea de base |
| `fig2_upset.py` | Figura 2 · solapamiento entre proteomas, 14 intersecciones |
| `fig3_estructura.py` | Figura 3 · CDS por transcrito y monoexónicos |
| `fig4_cobertura.py` | Figura 4 · exactitud contra referencia externa y hallazgo de circularidad |
| `tabla_estructura.py` | Sección 2 de la Tabla 1 |
| `tabla_funcional.py` | Tabla 2 · anotación funcional por brazo |
| `tabla1_final.py` | Ensambla la Tabla 1 definitiva |

Dependencias externas: solo **matplotlib** (las cuatro figuras) y **numpy**
(`fig4_cobertura.py`). Las tres de tablas usan solo la biblioteca estándar.

> **`upsetplot` NO se usa.** No es compatible con las versiones actuales de pandas y
> matplotlib. `fig2_upset.py` dibuja el diagrama con matplotlib puro.

> **Los tres defectos que arrastraba esta carpeta, resueltos en la ronda 8.**
>
> - `fig3_estructura.py` anunciaba tres paneles y dibujaba dos. **Corregido el
>   docstring, no añadido el panel**: la longitud proteica se retiró a propósito porque
>   en escala logarítmica la diferencia entre 524 y 551 aminoácidos no se distingue. Se
>   sigue calculando e imprimiendo por consola. De paso, los dos paneles llevaban ambos
>   el comentario `# B.` en el código.
> - `tabla1_final.py` localizaba las secciones por texto exacto sin `try` y reventaba
>   con un `ValueError` pelado. **Sigue localizándolas por texto exacto** —es frágil por
>   diseño y no hay forma barata de evitarlo— pero ahora comprueba que existan los dos
>   ficheros de entrada y aborta con un mensaje que dice qué encabezado no encuentra y
>   dónde hay que tocarlo.
> - `tabla_funcional.py` llevaba 86,6 %, 90,2 % y 18.765 escritos a mano en el
>   docstring. **Se contrastaron contra `Tabla2_funcional.md` y eran correctos**; el
>   problema era que estuvieran escritos a mano y pudieran divergir. Se han sustituido
>   por una remisión a la tabla generada.

## 7. El proteoma candidato combinado

| Script | Papel |
|---|---|
| `build_combined_proteome.sh` | Ensambla el proteoma candidato combinado: núcleo LiftOn + loci de Helixer rescatados − redundantes |

```
20.233 (núcleo LiftOn)  +  555 (rescatados)  −  80 (redundantes)  =  20.708
```

Es el **producto** del trabajo, no una de las cuatro anotaciones que se comparan. Se
deposita como `llama_combined_reference_proteome.faa.gz` en el registro de datos 1.3.0
(DOI 10.5281/zenodo.22150587). Contado sobre el fichero depositado: 20.233
identificadores de origen LiftOn y 475 de origen Helixer, que es el mismo 20.708 leído
al revés.

**Secuencias, no proteínas.** El recuento se expresa en secuencias en toda la
documentación desde la ronda 9. No son proteínas validadas: parte del núcleo arrastra
pautas de lectura rotas, y eso es justamente lo que clasifica `capas.py`.

**Por qué LiftOn es el núcleo.** No por ser el más completo, sino por ser el brazo de
homología mejor modelado: la mayor fracción de proteínas que cubren ≥80 % de la consulta
contra la referencia externa (92,6 %) y los menos stops internos de los tres brazos de
homología (2,7 %).

**Los 80 excluidos no son 80 duplicaciones.** `S4_modelos_excluidos.tsv` los lista:
**76 corresponden a un solo modelo del núcleo cada uno**, y los cuatro restantes forman
**dos pares** que apuntan al mismo modelo desde posiciones contiguas del mismo scaffold.
Eso es un locus fragmentado en dos modelos de Helixer, no dos genes duplicando uno.

> **UN PASO DEL ENSAMBLADO NO SE PUEDE SIMPLIFICAR.** DIAMOND aborta ante el carácter
> `.` que gffread emite para un stop en pauta, de modo que la **base** se construye
> sobre una copia saneada del núcleo con `.` sustituido por `X`, mientras que el fichero
> de **salida** se arma sobre el proteoma **original**, con las marcas de stop intactas.
> Sanear también la salida borraría en silencio la evidencia de pautas rotas que después
> clasifica `capas.py`. La versión que subió la ronda 8 no hacía el saneado y **habría
> fallado al ejecutarse**.

## 8. Capas de confianza y análisis de sensibilidad

| Script | Salida | Qué sostiene |
|---|---|---|
| `capas.py` | `Lgla_combined_proteome_confidence.tsv` | Las tres capas del proteoma combinado. Es lo que respalda el «confidence-aware» del título |
| `s7.py` | `S7_sensibilidad_umbral.tsv` | Que el hallazgo de circularidad no depende del umbral elegido |
| `s6.py` | `S6_correspondencia_referencia.tsv` | La asimetría de isoformas entre los brazos de homología |

Capas, contadas sobre el fichero: `high_confidence` 18.920, `extended_reference` 1.313,
`extended_candidate` 475. Suman 20.708. Ninguna entrada de la primera capa tiene stop
interno ni le falta hit externo: el criterio se aplica sin excepción.

Sensibilidad al umbral de sobre-extensión, contraste Helixer/LiftOn contra alpaca:
**27× a 10 pp, 35× a 15, 51× a 20 y 49× a 30**. Contra dromedario, los cuatro brazos
quedan dentro de un factor de 1,4 en los cuatro umbrales. El 20 pp del manuscrito es un
punto de una meseta, no una elección afortunada.

> **`stops.py` NO se ha incorporado.** Existe en el directorio de trabajo y hace lo
> mismo que `internal_stops.py`, que ya está aquí y además documenta la cautela de que
> la fila de Helixer es no informativa por construcción. Añadirlo habría creado un
> duplicado funcional con dos nombres distintos.

> **AVISO DE PROCEDENCIA DE LOS DOS SCRIPTS RECONSTRUIDOS.** `run_swissprot.sh` y
> `build_combined_proteome.sh` se escribieron **después** de las ejecuciones que
> documentan, a partir del registro de trabajo. No son transcripciones. La diferencia
> entre los dos, desde la ronda 9, importa:
>
> - **`build_combined_proteome.sh` está verificado por su salida.** Ejecutado sobre las
>   entradas originales reproduce el fichero depositado **byte a byte**, MD5
>   `dc69fa820facc0f087696bdd4885e1c1`, 11.337.876 bytes, 20.708 secuencias. Su
>   comprobación final aborta con código de error si el recuento no cuadra.
> - **`run_swissprot.sh` no tiene esa comprobación.** Lo único que respalda sus
>   parámetros es la forma del fichero que produjeron: `novel_vs_sprot.tsv` tiene siete
>   columnas y 29 filas para 10 consultas, lo que exige `--max-target-seqs 5` y ese
>   `outfmt`. Los parámetros que documentaba `REPRODUCIBILITY.md` §8.2 desde la ronda 7
>   —`--max-target-seqs 1 --max-hsps 1` y doce columnas— habrían dado diez filas y doce
>   columnas. Estaban mal, y la ronda 8 los copió al script. Corregidos en la ronda 9.

## Rutas

Los scripts de análisis esperan encontrar en el directorio de trabajo:

- `llama_liftoff.{faa.gz,gff3.gz}`, `llama_miniprot.*`, `llama_lifton.*`
  (del depósito de datos, versión 1.3.0, DOI 10.5281/zenodo.22150587; DOI de concepto
  10.5281/zenodo.21445839, que resuelve siempre a la última)
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

**Las proteínas más largas de Helixer NO son sobre-extensión.** Fue una duda abierta
durante un tiempo y la resolvió el análisis de la sección 3. Medido contra *V. pacos*,
Helixer parecería sobre-extender más de cincuenta veces lo que LiftOn, pero es un
artefacto de la circularidad: los brazos de homología alinean contra su propia fuente.
Medido contra *C. dromedarius*, que es el patrón externo, los cuatro quedan entre 1,7 %
y 2,2 % y son indistinguibles. La explicación de la mayor longitud media es la
preselección: Helixer solo emite modelos donde consigue construir un ORF completo.
