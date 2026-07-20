# llama-annotation-pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21456816.svg)](https://doi.org/10.5281/zenodo.21456816)

Reproducible Snakemake workflow for homology-based and *ab initio* annotation
of the *Lama glama* genome.

This repository contains the workflow, environment definitions and helper
scripts used to annotate the protein-coding gene complement of a *Lama glama*
(llama) genome assembly and to benchmark the resulting proteomes. **No sequence
or annotation data are tracked in this repository** — all `.fasta`, `.gff3`,
`.faa`, BUSCO downloads and pipeline outputs are produced locally and are
excluded via `.gitignore`.

## Data

The annotation sets and proteomes produced with this workflow are deposited
separately:

**DOI: [10.5281/zenodo.21445840](https://doi.org/10.5281/zenodo.21445840)** (CC0-1.0)

The deposit contains the four annotations (GFF3), the four derived proteomes
(FASTA), and the quality-control output for six protein sets: the four produced
here, the *Vicugna pacos* reference proteome as a ceiling, and a previous
unpublished in-house annotation as a baseline.

The files are under embargo until the associated manuscript is published; the
record metadata, including the full description of contents and caveats, is
publicly visible in the meantime.

Note that code and data are deposited as **separate records with different
licences**: MIT for this workflow, CC0-1.0 for the data. This is deliberate.

Source assembly and reference annotation are not redistributed in either deposit.
They are available from INSDC as GCA_028534125.1 and GCF_048564905.1
respectively.

## Overview

The project compares **four** protein-coding annotation strategies on the same
*Lama glama* target assembly, benchmarked against the *Vicugna pacos* reference
proteome (ceiling) and a previous llama proteome (baseline):

- **Liftoff** (`envs/liftoff.yaml`) — homology, DNA-to-DNA lift-over of the
  reference annotation, with optional `-polish` and `-copies`.
- **miniprot** (`envs/miniprot.yaml`) — homology, protein-to-genome alignment of
  the reference proteome.
- **LiftOn** — hybrid homology, reconciling the Liftoff and miniprot models.
- **Helixer** — *ab initio* deep-learning gene prediction.

### Automation boundary (important)

The Snakemake workflow automates the data acquisition, the **Liftoff** and
**miniprot** branches, and their BUSCO / AGAT / gffcompare QC. **LiftOn** and
**Helixer** were added later and are run **manually** (LiftOn) or on **external
GPU infrastructure** (Helixer); they are *not* wired into the `Snakefile`. The
exact commands to reproduce them are documented in
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md). The comparative report
(`scripts/build_report.py`) discovers whichever branches are present under
`results/busco/`, so all four appear in the final table regardless of how they
were produced.

### Reference asymmetry between homology branches (important)

The three homology branches did not start from the same reference annotation.
The difference is one of isoform depth, not of gene repertoire.

With `keep_longest_isoform: true` (the setting used here), the workflow reduces
the *Vicugna pacos* reference annotation to one transcript per gene with AGAT,
and it is that reduced annotation (`annotation_primary.gff3`) and the
corresponding reduced proteome (`protein_primary.faa`) that feed the **Liftoff**
and **miniprot** branches.

**LiftOn**, by contrast, requires the native NCBI GFF3 (`annotation_ncbi.gff`,
option `-ad RefSeq`), because an AGAT-processed annotation yields an invalid
proteome (see `REPRODUCIBILITY.md`). LiftOn therefore built its reference
dictionary from the full multi-isoform annotation. It consumed the Liftoff and
miniprot outputs only as evidence (`-L`, `-M`).

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

A second, independent set of asymmetries affects the *ab initio* branch. They are
described in full in the next section, which should be read before comparing
Helixer with the homology branches.

### Helixer run conditions (important)

Three aspects of how the *ab initio* branch was run must be taken into account
before comparing its completeness figures with those of the homology branches.

**Substrate.** Helixer was run on scaffolds of at least 10 kb: 3,640 scaffolds,
0.34 % of the total by count, holding 1,915,763,599 bp or 81.46 % of the
assembled sequence. The three homology branches were run on the complete
assembly.

**No subsequence overlap.** Prediction was performed on 213,840 bp subsequences
with no overlap between them. Overlapping prediction exists in Helixer to
mitigate the loss of accuracy at subsequence boundaries; running without it
leaves one such boundary every 213,840 bp across 1.92 Gb. The completeness
reported here is therefore a lower bound on what Helixer can achieve on this
assembly, not a representative measure of the tool's performance.

**Model loaded by path.** The model was supplied with `--load-model-path` rather
than selected with `--lineage`, so the lineage-based parameter defaults were not
applied.

These are properties of this particular run, not of the method. Users comparing
the four annotation sets should read the Helixer figures as descriptive of the
data released here rather than as a benchmark of *ab initio* prediction.

### Genomes (configurable in `workflow/config/config.yaml`)

| Role      | Species          | Accession           |
|-----------|------------------|---------------------|
| Target    | *Lama glama*     | `GCA_028534125.1`   |
| Reference | *Vicugna pacos*  | `GCF_048564905.1`   |

## Quality control & reporting

- **BUSCO** in protein mode (default lineage `artiodactyla_odb12`) on every
  proteome, with an optional secondary lineage.
- **AGAT** structural statistics per annotation.
- **gffcompare** concordance between Branch A and Branch B.
- A comparative report (`results/report/comparison_report.md` +
  `comparison_table.tsv`).

## Repository layout

```
.
├── workflow/
│   ├── Snakefile              # pipeline definition
│   ├── config/
│   │   └── config.yaml        # accessions, BUSCO lineage, tool parameters
│   └── envs/                  # per-rule conda environments
│       ├── datasets.yaml      # NCBI datasets (data acquisition)
│       ├── liftoff.yaml       # Liftoff
│       ├── miniprot.yaml      # miniprot
│       ├── agat.yaml          # AGAT
│       ├── busco.yaml         # BUSCO
│       └── gfftools.yaml      # gffread + gffcompare
├── scripts/
│   ├── build_report.py        # builds the comparative report (called by the `report` rule)
│   ├── run_lifton.sh          # LiftOn branch wrapper (manual; see REPRODUCIBILITY.md sec. 3)
│   └── run_helixer.py         # Helixer ab initio runner (external GPU; see REPRODUCIBILITY.md sec. 4)
├── REPRODUCIBILITY.md         # exact commands: automated branches + manual LiftOn / Helixer
├── CITATION.cff
├── .zenodo.json
├── LICENSE
└── .gitignore
```

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) 9.23.1 (the version used to
  produce the results reported in the accompanying manuscript; later versions
  are likely to work but have not been tested)
- [conda](https://docs.conda.io/) / mamba — per-rule environments are resolved
  automatically with `--use-conda`
- Network access for the data-acquisition rules (NCBI `datasets`)

## Usage

```bash
# from the workflow/ directory
cd workflow

# dry run (plan only)
snakemake -n

# full run with conda environments
snakemake --use-conda --cores 8
```

Input genomes are downloaded automatically from NCBI using the accessions in
`config.yaml`; edit that file to point the workflow at different assemblies,
BUSCO lineages or parameters.

## Data availability

Raw and intermediate data (genomes, proteomes, annotations, BUSCO downloads,
reports) are **not** included in this repository by design. They are regenerated
by running the workflow, or obtained from the deposit listed in the [Data](#data)
section above and from INSDC.

## Citation

If you use this workflow, please cite it via the metadata in
[`CITATION.cff`](CITATION.cff), or through its Zenodo DOI:

- Concept DOI, always resolving to the latest version:
  [10.5281/zenodo.21456816](https://doi.org/10.5281/zenodo.21456816)
- Version 1.0.0:
  [10.5281/zenodo.21456817](https://doi.org/10.5281/zenodo.21456817)

Cite the version DOI when referring to the exact state of the code used for a
given analysis, and the concept DOI otherwise.

## License

Released under the [MIT License](LICENSE).
