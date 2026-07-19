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

The three homology branches did not start from the same reference annotation,
and this must be taken into account when comparing them.

With `keep_longest_isoform: true` (the setting used here), the workflow reduces
the *Vicugna pacos* reference annotation to one transcript per gene with AGAT,
and it is that reduced annotation (`annotation_primary.gff3`) and the
corresponding reduced proteome (`protein_primary.faa`) that feed the **Liftoff**
and **miniprot** branches.

**LiftOn**, by contrast, requires the native NCBI GFF3 (`annotation_ncbi.gff`,
option `-ad RefSeq`), because an AGAT-processed annotation yields an invalid
proteome (see `REPRODUCIBILITY.md`, section 3.4). LiftOn therefore built its
reference dictionary from the full multi-isoform annotation: 86,028 transcripts
and 56,808 proteins, of which 349 truncated. It consumed the Liftoff and
miniprot outputs only as evidence (`-L`, `-M`).

Consequences for anyone reusing or extending this work:

- LiftOn had access to a richer reference than the other two homology branches.
  Any direct comparison of completeness between the three must state this.
- A small number of identifier lookups failed for this reason: 31 loci out of
  roughly 33,300 processed, 29 of them carrying AGAT-generated identifiers
  (`agat-gene-N`) absent from the native NCBI database. The effect on the
  resulting proteome is marginal, but it is the visible trace of the asymmetry.
- To remove the asymmetry, run all three homology branches against the same
  reference annotation. This was not done for the results reported here.

A second, independent asymmetry affects the *ab initio* branch: **Helixer** was
run on scaffolds of at least 10 kb (3,640 scaffolds, 0.34 % of the total,
holding 1,915,763,599 bp or 81.46 % of the assembled sequence), whereas the
three homology branches were run on the complete assembly. Completeness figures
are therefore not directly comparable between Helixer and the homology
branches.

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
by running the workflow, or are available from the sources cited in the
accompanying manuscript.

## Citation

If you use this workflow, please cite it via the metadata in
[`CITATION.cff`](CITATION.cff). A Zenodo DOI will be minted from a tagged
release once the repository is made public (see the manuscript for the
version-of-record DOI).

## License

Released under the [MIT License](LICENSE).
