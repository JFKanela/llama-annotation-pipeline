# llama-annotation-pipeline

Reproducible Snakemake workflow for homology-based and *ab initio* annotation
of the *Lama glama* genome.

This repository contains the workflow, environment definitions and helper
scripts used to annotate the protein-coding gene complement of a *Lama glama*
(llama) genome assembly and to benchmark the resulting proteomes. **No sequence
or annotation data are tracked in this repository** — all `.fasta`, `.gff3`,
`.faa`, BUSCO downloads and pipeline outputs are produced locally and are
excluded via `.gitignore`.

## Overview

The workflow annotates a target llama assembly by transferring / projecting a
well-annotated reference and then compares two independent strategies on the
**same** target genome:

- **Branch A — Liftoff** (`envs/liftoff.yaml`): DNA-to-DNA lift-over of the
  reference annotation, with optional `-polish` and `-copies`.
- **Branch B — miniprot** (`envs/miniprot.yaml`): protein-to-genome alignment
  of the reference proteome.

Both branches are reduced to a single protein per gene (longest isoform) and
evaluated side by side.

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
│   └── fix_report_4brazos.py  # report post-processing helper
├── docs/
├── CITATION.cff
├── .zenodo.json
├── LICENSE
└── .gitignore
```

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) (≥ 7)
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
by running the workflow, or are available from the sources cited in the
accompanying manuscript.

## Citation

If you use this workflow, please cite it via the metadata in
[`CITATION.cff`](CITATION.cff). A Zenodo DOI will be minted from a tagged
release once the repository is made public (see the manuscript for the
version-of-record DOI).

## License

Released under the [MIT License](LICENSE).
