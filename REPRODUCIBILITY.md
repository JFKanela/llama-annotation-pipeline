# Reproducibility and pipeline execution

This document describes precisely **which part of the analysis is automated in
Snakemake and which part was executed manually**, giving the exact commands so
that anyone can reproduce the complete result.

The distinction is deliberate and stated openly: the Snakemake pipeline does not
reproduce the four annotation branches end to end. Two of them (LiftOn and
Helixer) were incorporated into the project after the workflow was written and
were executed with the commands documented here.

---

## 1. Scope summary

| Component | Automated in Snakemake | Execution |
|---|---|---|
| Genome download and reference annotation | Yes | `snakemake` |
| Liftoff branch (DNA-DNA homology) | Yes | `snakemake` |
| miniprot branch (protein-DNA homology) | Yes | `snakemake` |
| LiftOn branch (hybrid homology) | **No** | Manual, section 3 |
| Helixer branch (*ab initio*) | **No** | External (GPU), section 4 |
| Proteome extraction (gffread) | Partial | Manual for LiftOn and Helixer |
| BUSCO quality control | Yes for the automated branches | Manual for LiftOn and Helixer |
| Comparative report | Yes | `scripts/build_report.py` |

`scripts/build_report.py` discovers the methods by scanning `results/busco/`, so
it adds to the report any branch present, whether or not it is automated.

---

## 2. Base environment and automated execution

```bash
conda activate smk        # Snakemake 9.23.1
snakemake --use-conda --cores 6
```

The per-tool environments are declared in `envs/` and Snakemake builds them in
isolation. Quality control uses BUSCO 6.1.0 with the `artiodactyla_odb12`
lineage (n = 12,594). Note that the `cetartiodactyla_odb12` lineage does not
exist in OrthoDB v12: the clade was renamed.

On memory requirements: with 20 GB of RAM, running several BUSCO jobs
concurrently can exhaust memory and trigger the OOM killer. Using `--cores 6`
instead of 8 is recommended, together with swap space (a 16 GB swapfile was
configured on the development system).

---

## Reference asymmetry between homology branches (important)

The three homology branches did not start from the same reference annotation.
This is a property of the results reported here, not an incident, and it must be
stated whenever their completeness figures are compared.

`keep_longest_isoform: true` (the setting used) reduces the *Vicugna pacos*
reference to one transcript per gene with AGAT. It is that reduced annotation
(`annotation_primary.gff3`) and its reduced proteome (`protein_primary.faa`)
that feed **Liftoff** and **miniprot**.

**LiftOn** cannot use it: it requires the native NCBI GFF3 (`-ad RefSeq`),
because an AGAT-processed annotation yields an invalid proteome (section 3.4).
LiftOn therefore built its reference dictionary from the full multi-isoform
annotation: **86,028 transcripts and 56,808 proteins, of which 349 truncated**.
It consumed the Liftoff and miniprot outputs only as evidence (`-L`, `-M`).

Visible trace of the asymmetry in the log: **31 loci out of roughly 33,300
processed** failed identifier lookup, **29 of them carrying AGAT-generated
identifiers** (`agat-gene-N`) absent from the native NCBI database. The effect on
the resulting proteome is marginal.

Consequences:

- LiftOn had access to a richer reference than the other two homology branches.
- To remove the asymmetry, all three homology branches would have to be run
  against the same reference annotation. **This was not done for the results
  reported here.**

A second, independent asymmetry affects the *ab initio* branch: Helixer was run
on scaffolds of at least 10 kb, whereas the homology branches were run on the
complete assembly (see section 4.2 for the exact figures). Completeness is
therefore not directly comparable between Helixer and the homology branches.

---

## 3. LiftOn branch (manual execution)

The steps below (sections 3.2 and 3.3, plus the proteome extraction of section 5)
are wrapped in `scripts/run_lifton.sh`. The environment setup in 3.1 is a
one-time prerequisite and is not scripted.

### 3.1. Building the environment

LiftOn 1.0.9 is not available on bioconda, only on PyPI, and its installation
requires resolving a dependency chain in a specific order:

```bash
conda create -n lifton -c conda-forge -c bioconda python=3.10 liftoff miniprot -y
conda activate lifton

# 1. cigar: obsolete packaging that requires pkg_resources
pip install --no-cache-dir setuptools-scm pkg_resources
pip install --no-cache-dir --no-build-isolation cigar

# 2. mappy: does not compile with a recent gcc, install via conda
conda install -n lifton -c bioconda -c conda-forge mappy -y

# 3. LiftOn without dependencies, already resolved above
pip install --no-cache-dir --no-build-isolation --no-deps lifton

# 4. remaining dependencies
pip install --no-cache-dir intervaltree "duckdb>=1.0" "pyarrow>=14"
```

### 3.2. Reference annotation in native format

LiftOn requires the **native NCBI** GFF. Using an annotation previously
processed with AGAT produces an invalid proteome (see note 3.4). Download:

```bash
BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/048/564/905"
DIR=$(curl -s "$BASE/" | grep -oP 'GCF_048564905\.1_[^/"]+' | head -1)
wget "${BASE}/${DIR}/${DIR}_genomic.gff.gz" -O resources/reference/annotation_ncbi.gff.gz
gunzip resources/reference/annotation_ncbi.gff.gz
```

### 3.3. Execution

LiftOn reuses the Liftoff and miniprot annotations already produced by the
pipeline (options `-L` and `-M`), which substantially reduces computation time.
Default parameters are used: rescue via miniprot, best-result selection with
verification, and transfer of gene features. `--legacy-merge` is not used.

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

Note that the target genome precedes the reference on the command line.

### 3.4. Note on a previous failed run

A first run used the annotation reduced to the primary isoform with AGAT as the
reference. LiftOn generated 33,045 gene models with apparently correct
coordinates, but only **2,533 translatable proteins**. The problem was fully
resolved by using the native NCBI GFF, which yielded **20,233 proteins**.

**The distinguishing criterion is the result, not the number of validation
warnings.** It is tempting to read the warning count as the symptom, and it is
wrong: the *valid* run logged **1,578,333** input-GFF validation warnings,
almost three times more than the failed one. Those warnings are duplicate CDS
identifiers, a normal artifact of RefSeq GFF3 files, which LiftOn resolves with a
unique-identifier transformation (`create_unique`). A high warning count is
therefore expected and diagnoses nothing.

This is documented because the failure is silent: the resulting GFF looks
correct, and only the translatable-protein count reveals the problem. Check the
protein count of the resulting `.faa`, not the log verbosity.

---

## 4. Helixer branch (external GPU execution)

Helixer was run on external infrastructure (Kaggle, Tesla P100 GPU, 16 GB)
because the development machine has a GPU with 4 GB of VRAM, insufficient for
inference. For this reason it is **not automated in Snakemake**: doing so would
suggest a reproducibility that does not exist in practice, since it requires
specific hardware and a platform account. The resulting GFF3 is treated as a
documented external input and is provided in the repository. The runner script
is `scripts/run_helixer.py`.

### 4.1. Environment

```bash
uv venv --python 3.11 /tmp/hxenv
uv pip install --python /tmp/hxenv/bin/python helixerlite "tensorflow[and-cuda]"
```

Model: `vertebrate_v0.3_m_0080.h5`, subsequence length 213,840. The model file is
**not** tracked in this repository (`.gitignore` excludes `*.h5`); it is
retrievable from its persistent identifier:

- Record: *Helixer–ab initio Prediction of Primary Eukaryotic Gene Models
  Combining Deep Learning and a Hidden Markov Model. Trained models.*, version
  v0.3 (Denton, Holst & Bolger; published 19 March 2024).
- DOI: **10.5281/zenodo.10836346** — <https://doi.org/10.5281/zenodo.10836346>
- File used: `vertebrate_v0.3_m_0080.h5` (36.3 MB, md5
  `acedf94d7c4f811e877da07844bc58f4`), downloaded from
  <https://zenodo.org/records/10836346/files/vertebrate_v0.3_m_0080.h5?download=1>

Note that the DOI above resolves to this specific version; a newer version of
the record exists. The `helixerlite` Python API is used (`fasta2hdf5`,
`HybridModel`, `preds2gff3`) together with `gfftk` for the final conversion.

### 4.2. Substrate and parameters

Scaffolds of length 10 kb or greater were retained: **3,640 scaffolds** (0.34 %
of the total) that account for **1,915,763,599 bp**, i.e. **81.46 %** of the
assembled sequence. The homology branches were run on the complete assembly
(2,351,761,190 bp). This substrate asymmetry conditions any direct comparison of
completeness and must be kept in mind when interpreting the results.

Parameters: no overlap, batch size 16. The substrate was processed in seven
consecutive blocks of 300 Mb (`CHUNK_BP = 300_000_000`) due to the platform's
session time limit (about 12 hours). Intermediate files were written to ample
temporary storage and the partial GFF3s to a persistent working directory, so
that the run was resumable.

### 4.3. Resuming across sessions

The complete run (about 14 hours) exceeds the session limit. On resume, the
already-computed blocks must be recovered from the previous session's output,
attached as an input dataset:

```python
import os, shutil, glob
os.makedirs("/kaggle/working/hxchunks", exist_ok=True)
for f in glob.glob("/kaggle/input/**/chunk*.gff3", recursive=True):
    dst = "/kaggle/working/hxchunks/" + os.path.basename(f)
    if not os.path.exists(dst):
        shutil.copy(f, dst)
```

The script detects the blocks present and processes only the pending ones.

### 4.4. GFF3 integrity check

The concatenation of the seven blocks was verified before use:

```bash
# all seven blocks represented
awk -F'\t' '!/^#/ && $3=="gene"' llama_helixer.gff3 | grep -oP 'ID=c\d+' | sort | uniq -c
# absence of duplicate identifiers
awk -F'\t' '!/^#/ && $3=="gene"' llama_helixer.gff3 | grep -oP 'ID=[^;]+' | sort | uniq -d | wc -l
# features with invalid coordinates, by type
awk -F'\t' '!/^#/ && $5 < $4 {c[$3]++} END{for(f in c) print f, c[f]}' llama_helixer.gff3
```

Result: the seven blocks present, zero duplicate identifiers, and 42 exon-type
features with invalid coordinates (end before start), a known artifact of the
GFF3 conversion. No CDS was affected, so the proteome was not compromised:
`gffread` discards those lines and translates from the CDS.

---

## 5. Manual proteome extraction and quality control

For the two non-automated branches, after obtaining the GFF3:

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

Quality control must use **the same BUSCO environment** that Snakemake built, to
guarantee identical version and lineage. Its location:

```bash
for d in .snakemake/conda/*/; do
  if [ -x "${d}bin/busco" ]; then echo "BUSCO at: $d"; fi
done
```

And the execution, replacing `PATH` with the previous result and `METHOD` with
`lifton` or `helixer`:

```bash
conda run -p PATH busco \
    -i results/proteomes/llama_METHOD.faa \
    -m proteins -l artiodactyla_odb12 -c 6 \
    -o METHOD --out_path results/busco
```

Each run takes about 100 minutes with six cores. It is best launched in a
persistent session (`tmux`). BUSCO writes the summary with an extended name
(`short_summary.specific.artiodactyla_odb12.METHOD.txt`); `build_report.py`
handles all naming variants.

---

## 6. Generating the comparative report

```bash
snakemake --use-conda --cores 4 results/report/comparison_report.md
```

Or directly:

```bash
python scripts/build_report.py \
    --busco-dir results/busco \
    --agat-dir results/agat \
    --gffcompare results/gffcompare/cmp.stats \
    --unmapped results/liftoff/unmapped.txt \
    --out-md results/report/comparison_report.md \
    --out-tsv results/report/comparison_table.tsv
```

The script does not recompute anything: it only reads the files already produced
and assembles the comparative table from them.

---

## 7. Input data

| Resource | Identifier |
|---|---|
| Target genome (*Lama glama*) | GCA_028534125.1 (DNA Zoo, Hi-C, specimen Fiesta) |
| Reference genome and annotation (*Vicugna pacos*) | GCF_048564905.1 (VicPac4, RefSeq RS_2025_04) |
| BUSCO lineage | `artiodactyla_odb12` |
| Helixer model | `vertebrate_v0.3_m_0080.h5` |

No RNA-seq evidence was available: the public records available for *Lama glama*
correspond to VHH or nanobody amplicons and do not constitute whole-transcriptome
data. This is the reason the main strategy is homology-based annotation.
