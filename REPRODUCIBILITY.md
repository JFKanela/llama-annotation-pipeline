# Reproducibility and pipeline execution

This document describes precisely **which part of the analysis is automated in
Snakemake and which part was executed manually**, giving the exact commands so
that anyone can reproduce the complete result.

The distinction is deliberate and stated openly: the Snakemake pipeline does not
reproduce the four annotation branches end to end. Two of them (LiftOn and
Helixer) were incorporated into the project after the workflow was written and
were executed with the commands documented here. The same applies to the two
downstream analyses added later: functional annotation with InterProScan
(section 6) and model accuracy against external references with DIAMOND
(section 7).

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
| Functional annotation (InterProScan) | **No** | Manual, section 6 |
| Model accuracy against external references (DIAMOND) | **No** | Manual, section 7 |

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

### Exact tool versions

Versions were read from the conda environment metadata on the compute machine
after the fact, except where explicitly noted as reconstructed.

| Tool | Version | Build | Environment | Pinned at run time |
|---|---|---|---|---|
| Liftoff | 1.6.3 | pyhdfd78af_2 | `envs/liftoff.yaml`, `lifton` | yes |
| minimap2 | 2.24 | h7132678_1 | `envs/liftoff.yaml`, `lifton` | yes |
| miniprot (standalone branch) | 0.18 | h577a1d6_0 | `envs/miniprot.yaml` | yes |
| miniprot (inside LiftOn) | 0.13-r248 | he4a0461_1 | `lifton` | no |
| miniprot (BUSCO dependency) | 0.18 | h577a1d6_0 | `envs/busco.yaml` | n/a |
| LiftOn | 1.0.9 | pypi | `lifton` | n/a |
| gffread (Liftoff, miniprot proteomes) | 0.12.7 | h077b44d_6 | `envs/gfftools.yaml` | yes |
| gffread (LiftOn, Helixer proteomes) | 0.12.9 | - | `smk` | no |
| gffcompare | 0.12.10 | h9948957_0 | `envs/gfftools.yaml` | no |
| AGAT | 1.7.0 | pl5321hdfd78af_0 | `envs/agat.yaml` | no |
| BUSCO | 6.1.0 | pyhdfd78af_1 | `envs/busco.yaml` | yes |
| Snakemake | 9.23.1 | - | host | n/a |
| **InterProScan** | **5.78-109.0** | official distribution, `~/interproscan/` | standalone, outside conda | n/a |
| **OpenJDK** | **11** | conda environment, **not the system JDK** | InterProScan dependency | yes |
| **DIAMOND** | **2.2.4** | - | conda | n/a |
| **Helixer (web tool)** | **0.3.6** | official service, plabipd.de | external, managed by the service | n/a |
| helixerlite *(superseded run only)* | 25.5.27 | commit 04c8086, cp311 wheel | Kaggle (Tesla P100), isolated `uv` venv | no |
| gfftk *(superseded run only)* | 26.5.22 | - | same venv (helixerlite dependency) | no |
| TensorFlow (`tensorflow[and-cuda]`) *(superseded run only)* | 2.15.1 | - | same venv | no |
| tensorflow-addons *(superseded run only)* | 0.23.0 | - | resolved dependency | no |
| CPython *(superseded run only)* | 3.11.15 | - | `uv venv --python 3.11` | minor only |

The Helixer annotation reported in the manuscript was produced with the
**official Helixer web tool v0.3.6**. The `helixerlite` stack listed below it
belongs to a **first, superseded run** and is kept in the table only so that the
provenance of `scripts/run_helixer.py` remains legible. See section 4.

Two gffread versions were used. The Liftoff and miniprot proteomes were
extracted by the Snakemake workflow with gffread 0.12.7, pinned in
`envs/gfftools.yaml`. The LiftOn and Helixer proteomes were extracted manually
with gffread 0.12.9 from the working environment `smk`, the only gffread
available outside the Snakemake-managed environments. The difference is
patch-level; no equivalence test between the two versions was performed for this
dataset.

The same pattern applies to miniprot: version 0.18 is pinned in
`envs/miniprot.yaml` for the standalone branch, whereas the manually created
`lifton` environment was built without pinning it, so the solver installed
0.13-r248 alongside LiftOn's dependencies. LiftOn's internal miniprot-based
rescue therefore used an older miniprot than the standalone branch. Liftoff
(1.6.3) and minimap2 (2.24) are identical in both environments.

At the time of execution, gffcompare and AGAT were not pinned in their
environment files, and gffread in the working environment was installed without
a version constraint. The versions listed in the table are those actually
installed, read from the conda environment metadata after the fact.

The version constraints for AGAT and gffcompare were added to their environment
files after the analyses had been run, so that the pinned versions match those
actually used. The results reported in the accompanying manuscript were produced
with the versions listed in the table above, which the solver resolved at run
time in the absence of a constraint. Adding these constraints changes the hash of
the affected environments, so Snakemake will rebuild them on the next run; that
is expected and does not invalidate anything already executed.

**Determined versus reconstructed** *(applies to the superseded `helixerlite`
run only; the current annotation was produced with the Helixer web tool, whose
environment is managed by the service and is not under the user's control).*
`helixerlite` was installed without a
version constraint and its version is not recorded in the output GFF3. It is
*determined* as 25.5.27: only three releases exist on PyPI and that one has been
current since 27 May 2025, more than a year before the run. The same
date-versioning argument determines `gfftk` 26.5.22. By contrast, **TensorFlow
2.15.1, tensorflow-addons 0.23.0 and CPython 3.11.15 are a reconstruction**,
obtained by re-running the same unconstrained installation command one day after
the original run; they are not a record of what was installed. The model,
`vertebrate_v0.3_m_0080`, is documented and is the determining factor for
reproducibility: the code version without the model reproduces nothing, whereas
the model largely determines the prediction.

---
## Reference asymmetry between homology branches (important)

The three homology branches did not start from the same reference annotation.
The difference is one of isoform depth, not of gene repertoire. All figures
below were verified by direct counting on the reference files themselves; they
do not come from the LiftOn log.

`keep_longest_isoform: true` (the setting used) reduces the *Vicugna pacos*
reference annotation to one transcript per gene with AGAT, and it is that
reduced annotation (`annotation_primary.gff3`) and the corresponding reduced
proteome (`protein_primary.faa`) that feed **Liftoff** and **miniprot**.

**LiftOn**, by contrast, requires the native NCBI GFF3 (`annotation_ncbi.gff`,
option `-ad RefSeq`), because an AGAT-processed annotation yields an invalid
proteome (section 3.4). LiftOn therefore built its reference dictionary from the
full multi-isoform annotation. It consumed the Liftoff and miniprot outputs only
as evidence (`-L`, `-M`).

| Reference used by | gene | mRNA | proteins |
|---|---|---|---|
| LiftOn (`annotation_ncbi.gff` / `protein.faa`) | 37,446 | 56,758 | 56,808 |
| Liftoff and miniprot (`annotation_primary.gff3` / `protein_primary.faa`) | 37,446 | 21,233 | 21,283 |

Both references contain the **same 37,446 gene models**. They differ in isoform
depth: 2.67 mRNA per protein-coding gene versus one. The AGAT reduction removed
no gene. (In both files the protein count exceeds the mRNA count by 50, because
32 `V_gene_segment` and 18 `C_gene_segment` features carry CDS without being
`mRNA`; both are retained by the reduction.)

Consequences for anyone reusing or extending this work:

- LiftOn had more isoforms per locus available than the other two homology
  branches. Any direct comparison of completeness between the three should say
  so. The effect is bounded but not zero: on the reference proteome itself the
  primary-isoform set already recovers 98.5 % of BUSCO groups, so additional
  isoforms can add at most about 1.5 points there. That bound applies to the
  reference, not to the transfer step, where additional isoforms also provide
  additional opportunities for a clean mapping.
- The reduction generated 28 new gene identifiers (`agat-gene-N`) for gene
  features that AGAT had to recreate. All 28 failed lookup against the native
  NCBI database during the LiftOn run, and no other gene identifier of that kind
  failed: 28 out of 37,446 gene models (0.075 %). The log records 31 failure
  events because one of those entries also appears as a duplicated copy
  (suffix `_1`), and because two further failures involve native identifiers and
  have a different, unrelated cause. The remaining AGAT-generated identifiers
  (29,430 five-prime UTR, 22,526 three-prime UTR, 3,505 pseudogene) are not
  queried against the reference database and caused no failure.
- To remove the asymmetry, run all three homology branches against the same
  reference annotation. This was not done for the results reported here.

A second, independent asymmetry affects the *ab initio* branch: **Helixer** was
run on scaffolds of at least 25 kb (244 scaffolds, holding 79.6 % of the
assembled sequence), whereas the three homology branches were run on the
complete assembly (see section 4.2). Completeness figures are therefore not
directly comparable between Helixer and the homology branches.

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

## 4. Helixer branch (*ab initio*, external execution)

Helixer was run on external infrastructure because the development machine has a
GPU with 4 GB of VRAM, insufficient for inference. For this reason it is **not
automated in Snakemake**: doing so would suggest a reproducibility that does not
exist in practice. The resulting GFF3 is treated as a documented external input.

**The branch was run twice.** Section 4.A describes the **current** run, which
produced the annotation reported in the manuscript. Section 4.B describes the
**first, superseded** run, kept because `scripts/run_helixer.py` documents it and
because the comparison between the two is itself informative.

**Citation.** Helixer and HelixerPost are covered by the same reference:
*Helixer: ab initio prediction of primary eukaryotic gene models combining deep
learning and a hidden Markov model*, Nature Methods (2025), DOI
[10.1038/s41592-025-02939-1](https://doi.org/10.1038/s41592-025-02939-1). Cite
the published article, not the earlier preprint.

---

### 4.A. Current run — official Helixer web tool v0.3.6

The definitive *ab initio* annotation was produced through the official Helixer
web service, <https://www.plabipd.de/helixer_main.html>, running Helixer
**v0.3.6**. There is no run script: the execution is a form submission, so the
parameters below, not a command line, are the record of what was done.

| Parameter | Value |
|---|---|
| Interface | official Helixer web tool, plabipd.de |
| Helixer version | 0.3.6 |
| Input assembly | `GCA_028534125.1_Lama_glama_HiC_genomic.fna` |
| Scaffold length filter | **>= 25,000 bp** — **244 scaffolds**, **79.6 %** of the assembled sequence |
| Lineage | `vertebrate` |
| Subsequence length | 213,840 bp (lineage default) |
| Subsequence overlap | **enabled** (lineage defaults) |
| `overlap-offset` | 106,920 bp |
| `overlap-core-length` | 160,380 bp |
| `window-size` (HelixerPost) | 100 |
| `edge-threshold` (HelixerPost) | 0.1 |
| `peak-threshold` (HelixerPost) | 0.8 |
| `min-coding-length` (HelixerPost) | 60 |
| Output | `Lgla_hx036_helixer_FINAL.gff`, proteome `Lgla_hx036_helixer.faa` |

The four HelixerPost values are the tool's defaults. The web interface exposes no
control over `peak-threshold`, so alternative values could not be explored; no
parameter tuning was performed.

**Result versus the superseded run.**

| | First run (4.B) | Current run (4.A) |
|---|---|---|
| BUSCO complete | 79.3 % | **85.5 %** |
| BUSCO fragmented | 6.6 % | **5.9 %** |
| BUSCO missing | 14.1 % | **8.6 %** |
| Proteins | 18,933 | **18,765** |

The gain is attributable mainly to overlapping prediction, which removes the
accuracy loss at subsequence boundaries that the first run incurred every
213,840 bp.

#### Non-determinism (must be recorded)

The web tool does not activate Helixer's `--deterministic` mode, which is what
forces deterministic cuDNN and cuBLAS kernels. GPU kernel non-determinism is
therefore not suppressed.

Two consecutive submissions with **identical input and identical settings**
returned **18,765** and **18,763** genes, with micro-differences at the
boundaries of individual models. The 18,765-gene run was the one adopted and
deposited.

**This annotation is not reproducible bit for bit.** Anyone re-running it should
expect differences of this order and should not treat a mismatch of a few genes
as evidence of an error. This is a limitation of the execution route, not of
Helixer: the command-line tool run with `--deterministic` on the same GPU
architecture would not have it.

---

### 4.B. First run (SUPERSEDED) — `helixerlite` on Kaggle

**Everything in this subsection describes a run that is no longer the method of
the manuscript.** It is retained for traceability and to document
`scripts/run_helixer.py`, which is marked as superseded in its own header. The
figures below (scaffolds >= 10 kb, no overlap, BUSCO 79.3 %) must not be quoted
as current.

It was run on Kaggle (NVIDIA Tesla P100 GPU, 16 GB).

### 4.1. Environment *(superseded run)*

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

### 4.2. Substrate and parameters *(superseded run)*

Scaffolds of length 10 kb or greater were retained: **3,640 scaffolds** (0.34 %
of the total) that account for **1,915,763,599 bp**, i.e. **81.46 %** of the
assembled sequence. **The current run uses a different substrate: scaffolds
>= 25 kb, 244 scaffolds, 79.6 % of the assembly (section 4.A).** The homology
branches were run on the complete assembly (2,351,761,190 bp). This substrate
asymmetry conditions any direct comparison of completeness in either case.

Parameters: no overlap, batch size 16. The substrate was processed in seven
consecutive blocks of 300 Mb (`CHUNK_BP = 300_000_000`) due to the platform's
session time limit (about 12 hours). Intermediate files were written to ample
temporary storage and the partial GFF3s to a persistent working directory, so
that the run was resumable.

### Helixer run parameters *(superseded run — for the current one see 4.A)*

Helixer was run through `helixerlite` on Kaggle (NVIDIA Tesla P100) inside an
isolated virtual environment created with `uv`.

| Parameter | Value |
|---|---|
| Input assembly | `GCA_028534125.1_Lama_glama_HiC_genomic.fna.gz` (NCBI FTP) |
| Scaffold length filter | >= 10,000 bp |
| Model | `vertebrate_v0.3_m_0080.h5` (DOI 10.5281/zenodo.10836346) |
| Subsequence length | 213,840 bp |
| Subsequence overlap | **none** |
| Batch size | 16 |
| CPUs | 4 |
| Genome split | ~300 Mb blocks, whole scaffolds only |
| `window_size` | 100 |
| `edge_threshold` | 0.1 |
| `peak_threshold` | 0.8 |
| `min_coding_length` | 60 |
| Final GFF3 normalisation | `gfftk` 26.5.22 (`gff2dict` / `dict2gff3`) |

All four HelixerPost parameters were left at the helixerlite defaults
(`window_size` 100, `edge_threshold` 0.1, `peak_threshold` 0.8,
`min_coding_length` 60); no parameter tuning was performed. The two values that
appear explicitly in the run script are identical to the defaults. These
defaults were read from the signature of `preds2gff3` in the installed version
(helixerlite 25.5.27), not from upstream Helixer documentation.

The model was loaded by direct path (`--load-model-path`) rather than by
lineage, so the lineage-based parameter defaults that Helixer applies when
`--lineage` is used were not in effect.

The genome was processed in seven consecutive blocks of approximately 300 Mb
because of the session time limit of the external platform. Blocks were formed
by accumulating **whole scaffolds**, so no scaffold was split across blocks and
block boundaries introduce no edge artefact beyond the one already present at
subsequence boundaries. The per-block GFF3 files were concatenated and the
result verified: all seven blocks represented, no duplicate identifiers, and 42
exon features with invalid coordinates discarded by gffread, none of them
affecting a CDS.

### 4.3. Resuming across sessions *(superseded run)*

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

### 4.4. GFF3 integrity check *(superseded run)*

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

For the two non-automated branches, after obtaining the GFF3.

The `lifton` environment does **not** contain gffread, so `conda run -n lifton
gffread ...` fails. Use the working environment `smk`, which holds the only
gffread available outside the Snakemake-managed environments:

```bash
conda activate smk    # gffread 0.12.9

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

## 6. Functional annotation (InterProScan)

Run outside the `Snakefile`, like the LiftOn and Helixer branches. Scripts:
`scripts/make_interproscan_input.py`, `scripts/run_ips.sh`,
`scripts/estado_ips.sh`.

### 6.1. Preparing the input

The four proteomes hold 80,331 sequences, of which **42,364 are distinct**.
`make_interproscan_input.py` deduplicates by MD5 digest, cleans characters
InterProScan rejects, and splits the result into batches:

```bash
python3 scripts/make_interproscan_input.py
# -> camelid_unique_proteins_clean.faa   42,364 sequences
# -> md5_to_ids.tsv                      md5 -> original identifiers, per branch
# -> ips_chunks/chunk_XXX.faa            43 batches of 1,000 sequences
```

Two properties of this step are not optional.

**`md5_to_ids.tsv` is indispensable.** The deduplicated FASTA uses the MD5 digest
as the sequence identifier, so InterProScan's output carries no branch
information at all. Without this mapping the 754,835 annotations cannot be
attributed to liftoff, miniprot, lifton or helixer, and the per-branch table
cannot be rebuilt.

**In-frame stops are substituted, not deleted.** The three homology proteomes
contain `.` characters, gffread's rendering of an in-frame stop, because many
transferred models have a broken reading frame. InterProScan aborts at
`stepLoadFromFastaIntoDB` on encountering one. Each `.` is replaced by `X`.
Substitution preserves sequence length and therefore the domain coordinates;
deletion would shift every position downstream of the stop.

### 6.2. Running InterProScan

`scripts/run_ips.sh` iterates over the batches:

```bash
nice -n 19 ionice -c3 \
  ~/interproscan/interproscan-5.78-109.0/interproscan.sh \
    -i  ips_chunks/chunk_XXX.faa \
    -f  TSV \
    -o  ips_out/chunk_XXX.tsv \
    -cpu 4 \
    -dp \
    -T  ips_tmp \
  > ips_out/chunk_XXX.log 2>&1
```

**Version: 5.78-109.0**, running the 18 analyses enabled by default.

**`-dp` disables the precalculated match lookup.** That service resolves a query
by exact MD5 against UniProtKB. Llama proteins are not in UniProtKB, which is the
premise of this work, so the lookup can only ever miss. It was switched off
rather than left to fail silently over 42,364 queries.

**Java: OpenJDK 11 is required.** OpenJDK 17 does not work with this release. The
JDK used came from a conda environment, not from the system installation.

**Resumption and housekeeping.** The script skips any batch whose `.tsv` already
exists and is non-empty, so an interrupted run resumes where it stopped. It
aborts if less than 15 GB remain free in `$HOME`, and clears `ips_tmp/` after
each successful batch; InterProScan's temporary files otherwise accumulate to
tens of gigabytes over 43 batches.

`scripts/estado_ips.sh` is a read-only monitor: it counts completed `.tsv` files
against available batches, detects a live process with `pgrep -f interproscan`,
and extrapolates the remaining time from an observed rate of roughly 80 minutes
per batch. That figure is an empirical observation on the development machine,
not a specification.

### 6.3. Result

**754,835 annotations over 40,402 of the 42,364 unique proteins (95.4 %).** The
per-branch breakdown, obtained by expanding the results through
`md5_to_ids.tsv`, is in the `README.md`. Two points that belong with the numbers:

- The "with domain" counts **exclude MobiDBLite and Coils**. Both predict generic
  properties — disorder and coiled-coil structure — and annotate almost any
  protein, so including them makes every branch look near-perfect and removes the
  discrimination the table is meant to provide.
- **PANTHER is reported at family level only.** The TSV output never emits the
  subfamily (`PTHR12345:SF6`), not even when phylogenetic placement succeeds.
  Verified against InterProScan's own official test file, which yields zero
  subfamilies in TSV and the five expected ones in XML, and whose run produced no
  abort. Subfamily resolution requires a structured output format and was not
  used.

---

## 7. Model accuracy against external references (DIAMOND)

Run outside the `Snakefile`. Scripts: `scripts/run_blastp.sh`,
`scripts/analyze_blastp.py`, `scripts/estado_blastp.sh`.

This is the fourth quality criterion of Kourelis et al. (2019) and the only step
in this work that evaluates **model accuracy** rather than completeness or
internal consistency.

### 7.1. Why DIAMOND rather than BLASTP

Llama and alpaca diverged a few million years ago and their orthologues sit at
95-99 % identity. These are trivial alignments: there is no remote-homology
sensitivity problem for BLASTP's exhaustive search to solve. **DIAMOND 2.2.4** in
`--very-sensitive` mode returns the complete comparison in about 18 minutes,
against days for BLASTP, and the accelerated heuristic costs nothing that matters
at this evolutionary distance.

### 7.2. References

| Reference | Accession | Role |
|---|---|---|
| *Camelus dromedarius* | `GCF_036321535.1` (mCamDro1.pat, RS_2024_04, 50,982 proteins) | **External yardstick. The valid comparison** |
| *Vicugna pacos* | `GCF_048564905.1` | Internal control. **Circular** |

### 7.3. Commands

```bash
# one database per reference
diamond makedb --in <reference>/proteome.faa -d blastp/<name> --threads 4

# alignment, identical parameters for both references
diamond blastp \
  -q camelid_unique_proteins_clean.faa \
  -d blastp/<name> \
  -o blastp/hits_<name>.tsv \
  --outfmt 6 qseqid sseqid pident length qstart qend sstart send \
              evalue bitscore qlen slen \
  --very-sensitive \
  --max-target-seqs 1 \
  --max-hsps 1 \
  --evalue 1e-5 \
  --threads 4 \
  --tmpdir blastp/tmp

python3 scripts/analyze_blastp.py
```

`qlen` and `slen` must be in the output format: the coverage figures are derived
from them, as `100 * (send - sstart + 1) / slen` for the subject and
`100 * (qend - qstart + 1) / qlen` for the query.

**Resumption.** `run_blastp.sh` skips `makedb` and `blastp` when the output
exists and is non-empty, writes to a `.partial` file and renames it atomically on
success, deleting it on failure. An interrupted run therefore never leaves a
truncated hit table that a later step would silently treat as complete.

**Paths are machine-specific.** The script points at the reference proteomes
under a local research-data mount. Anyone re-running it must edit those paths;
the accessions in the table above are what identifies the data.

`scripts/estado_blastp.sh` is a read-only monitor of the run.

### 7.4. Result, and the methodological finding

The figures against *C. dromedarius* are in the `README.md`. The result worth
recording here is not any single number but the effect of the choice of
reference:

> Measured against alpaca, Helixer appears to over-extend its models twenty times
> more often than LiftOn (2.0 % versus 0.0 %). Measured against dromedary, all
> four branches sit at the same 2 % and are indistinguishable.
>
> Three of the four branches are projections of the alpaca annotation. They align
> almost perfectly against their own source, so their 0.0 % measures identity
> with themselves rather than model quality. Only the *ab initio* branch is
> independent of that reference, and it is the only one the comparison penalises.
>
> **Had alpaca been used as the yardstick, this work would have reported that the
> *ab initio* branch over-extends its models, and that claim would have been
> false.**

The alpaca comparison is retained as an internal control and must not be read as
a measure of accuracy.

Two limitations, declared:

- **Median coverages are not reported.** They come out at 100 % for all four
  branches against both references: they saturate and discriminate nothing.
  The reported statistic is the fraction of proteins with coverage >= 80 %.
  `analyze_blastp.py` still computes and prints the medians; that output is
  diagnostic, not a result.
- **One HSP per pair.** `--max-hsps 1` means that in multidomain proteins whose
  alignment fragments into several local matches, only one is counted and
  coverage is recorded as artificially low. The approximation errs in the
  conservative direction: it underestimates coverage, it does not inflate it.


---

## 8. Generating the comparative report

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

## 9. Input data

| Resource | Identifier |
|---|---|
| Target genome (*Lama glama*) | GCA_028534125.1 (DNA Zoo, Hi-C, specimen Fiesta) |
| Reference genome and annotation (*Vicugna pacos*) | GCF_048564905.1 (VicPac4, RefSeq RS_2025_04) |
| BUSCO lineage | `artiodactyla_odb12` |
| Helixer model | `vertebrate_v0.3_m_0080.h5` |

No RNA-seq evidence was available: the public records available for *Lama glama*
correspond to VHH or nanobody amplicons and do not constitute whole-transcriptome
data. This is the reason the main strategy is homology-based annotation.
