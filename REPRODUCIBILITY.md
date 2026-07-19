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

### Exact tool versions

Entries marked *to be confirmed* must be completed from the corresponding
environment on the compute machine. They are deliberately left open rather than
filled in by inference.

| Tool | Version | Environment |
|---|---|---|
| Liftoff | 1.6.3 | `envs/liftoff.yaml` (pinned) and `lifton` |
| minimap2 | 2.24 | `envs/liftoff.yaml` (pinned) and `lifton` |
| miniprot (standalone branch) | 0.18-r281 | `envs/miniprot.yaml` (pins `miniprot=0.18`) |
| miniprot (inside LiftOn) | 0.13-r248 | `lifton` |
| LiftOn | 1.0.9 | `lifton` (pip) |
| Snakemake | 9.23.1 | host |
| BUSCO | 6.1.0 | `envs/busco.yaml` (pinned) |
| gffread | 0.12.9 — to be confirmed, see note below | `envs/gfftools.yaml` pins `gffread=0.12.7` |
| AGAT | to be confirmed | `envs/agat.yaml` (unpinned) |
| gffcompare | to be confirmed | `envs/gfftools.yaml` (unpinned) |
| helixerlite | to be confirmed | Kaggle, Tesla P100 |

The `lifton` environment was created without pinning miniprot, so the solver
installed the version compatible with LiftOn's dependencies (0.13-r248), whereas
the standalone miniprot branch pins `miniprot=0.18` in `envs/miniprot.yaml` and
resolves to 0.18-r281. LiftOn's internal miniprot-based rescue therefore used an
older miniprot than the standalone branch. Liftoff (1.6.3) and minimap2 (2.24)
are identical in both environments, so the discrepancy is limited to miniprot.

**Unresolved: gffread.** The version recorded during the audit was 0.12.9, but
`envs/gfftools.yaml` pins `gffread=0.12.7`. These cannot both describe the same
environment, so one of the two must be corrected: either the audited binary came
from an environment other than `llama_gfftools`, or the pin was changed after the
run. This must be settled on the compute machine before the `v1.0.0` tag.

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
run on scaffolds of at least 10 kb (3,640 scaffolds, 0.34 % of the total,
holding 1,915,763,599 bp or 81.46 % of the assembled sequence), whereas the
three homology branches were run on the complete assembly (see section 4.2).
Completeness figures are therefore not directly comparable between Helixer and
the homology branches.

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
