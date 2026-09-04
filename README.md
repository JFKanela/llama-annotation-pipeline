# llama-annotation-pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21456816.svg)](https://doi.org/10.5281/zenodo.21456816)

Reproducible Snakemake workflow for homology-based and *ab initio* annotation
of the *Lama glama* genome.

This repository contains the workflow, environment definitions and helper
scripts used to annotate the protein-coding gene complement of a *Lama glama*
(llama) genome assembly, to benchmark the resulting proteomes, and to build a
combined candidate reference proteome from them. **No sequence or annotation
data are tracked in this repository** — all `.fasta`, `.gff3`, `.faa`, BUSCO
downloads and pipeline outputs are produced locally and are excluded via
`.gitignore`.

## Identifiers

| | |
|---|---|
| Code, concept DOI | [10.5281/zenodo.21456816](https://doi.org/10.5281/zenodo.21456816) |
| Data, concept DOI | [10.5281/zenodo.21445839](https://doi.org/10.5281/zenodo.21445839) |
| WorkflowHub | [10.48546/workflowhub.workflow.2250.1](https://doi.org/10.48546/workflowhub.workflow.2250.1) |
| RRID | `SCR_028889` |

## Data

The annotation sets, proteomes and analysis outputs produced with this workflow
are deposited separately, **openly accessible**, under CC0-1.0:

**Data, version 1.3.0: [10.5281/zenodo.22150587](https://doi.org/10.5281/zenodo.22150587)**

The deposit contains the four annotations (GFF3), the four derived proteomes and
the combined candidate reference proteome (FASTA), two derived files of that
product (the high-confidence layer on its own and a copy with the internal stop
characters masked), the `chaku_v1` baseline proteome, the confidence-layer
assignment table, the quality-control output for six protein sets, the functional
annotation and alignment results, the transposable-element protein check, and the
nine additional files of the accompanying manuscript.

Note that code and data are deposited as **separate records with different
licences**: MIT for this workflow, CC0-1.0 for the data. This is deliberate.

Source assembly and reference annotation are not redistributed in either
deposit. They are available from INSDC as `GCA_028534125.1` and
`GCF_048564905.1` respectively.

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

The counting unit throughout is the **sequence**, understood as one entry per
gene model. None of the deposited files contains alternative isoforms. The
entries are predicted models, not experimentally validated proteins.

### Automation boundary (important)

The Snakemake workflow automates the data acquisition, the **Liftoff** and
**miniprot** branches, and their BUSCO / AGAT / gffcompare QC. **LiftOn**,
**Helixer**, **InterProScan** and **DIAMOND** are run **manually** (LiftOn,
InterProScan, DIAMOND) or on **external infrastructure** (Helixer, through the
official web tool); they are *not* wired into the `Snakefile`. The exact
commands to reproduce them are documented in `REPRODUCIBILITY.md`. The
comparative report (`scripts/build_report.py`) discovers whichever branches are
present under `results/busco/`, so all four appear in the final table regardless
of how they were produced.

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

| Reference used by | gene | mRNA | sequences |
|---|---|---|---|
| LiftOn (`annotation_ncbi.gff` / `protein.faa`) | 37,446 | 56,758 | 56,808 |
| Liftoff and miniprot (`annotation_primary.gff3` / `protein_primary.faa`) | 37,446 | 21,233 | 21,283 |

Both references contain the **same 37,446 gene models**, of which 21,233 have at
least one coding transcript. They differ in isoform depth: 2.67 mRNA per
protein-coding gene versus one. The AGAT reduction removed no gene. (In both
files the sequence count exceeds the mRNA count by 50, because 32
`V_gene_segment` and 18 `C_gene_segment` features carry CDS without being
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
  failed: 28 out of 37,446 gene models (0.075 %).
- To remove the asymmetry, run all three homology branches against the same
  reference annotation. This was not done for the results reported here.

### Helixer run conditions (important)

The *ab initio* branch was re-run after the first version of this repository was
released. **The annotation reported in the manuscript is the second run**, and
the two differ in method, substrate and results. `scripts/run_helixer.py`
documents the first run only and is retained for traceability, marked as
superseded.

| | First run (superseded) | Current run (manuscript) |
|---|---|---|
| Tool | `helixerlite` 25.5.27 on Kaggle | **official Helixer web tool v0.3.6** |
| Subsequence overlap | none | **enabled** (`vertebrate` defaults) |
| Substrate | scaffolds >= 10 kb (3,640; 81.46 %) | scaffolds >= 25 kb (**244**; **79.6 %**) |
| BUSCO complete | 79.3 % | **85.5 %** |
| Sequences | 18,933 | **18,765** |

**Substrate.** The three homology branches were run on the complete assembly and
place around **11 %** of their models on sequences shorter than 25 kb, which
Helixer did not analyse. Completeness is therefore not directly comparable
between Helixer and the homology branches.

**Subsequence overlap.** Subsequence length 213,840 bp, overlap offset
106,920 bp, overlap core length 160,380 bp. Overlapping is the main reason for
the 6.2-point gain in BUSCO completeness over the first run.

**Post-processing.** HelixerPost at its defaults: `window-size` 100,
`edge-threshold` 0.1, `peak-threshold` 0.8, `min-coding-length` 60. The web
interface does not allow other values of `peak-threshold`.

**Independence is one of inference, not of training.** The vertebrate model
`vertebrate_v0.3_m_0080.h5` was trained on 936 genomes that include *Vicugna
pacos* (`GCF_000164845.3`) and *Camelus dromedarius* (`GCF_000803125.2`) in the
training and validation split, in earlier assembly versions. No annotation is
consulted while the target genome is annotated, but the model learned
gene-structure patterns from references that included camelids.

**Non-determinism.** Two consecutive runs on identical input yielded **18,765
and 18,763 genes**, a difference of 0.011 % confined to model boundaries,
because the web interface does not enable deterministic mode. The 18,765-gene
run was adopted and is the one deposited; the GFF3 is identified by its MD5.

### Functional annotation

InterProScan 5.78-109.0 with the 18 default analyses under OpenJDK 11. The
80,331 sequences of the four proteomes were reduced to **42,364 unique** ones by
MD5, a 47.3 % reduction in sequences submitted, and `md5_to_ids.tsv` maps each
result back to its proteome of origin.

The precalculated match lookup (`-dp`) was disabled: that service searches by
exact MD5 in UniProtKB and the sequences of this species are not there. Invalid
characters were replaced by `X` rather than removed, to preserve sequence length
and with it the domain coordinates. PANTHER is reported at family level, because
the tabular output never emits subfamily identifiers.

Result: 754,835 annotations over 40,402 of the 42,364 unique sequences (95.4 %).

### Model accuracy against an external reference

DIAMOND 2.2.4 in `blastp` mode, `--very-sensitive --evalue 1e-5`, one best match
and one HSP per pair, against *Camelus dromedarius* `GCF_036321535.1`.

**This is where the methodological finding of the accompanying manuscript comes
from.** Measured against *V. pacos*, the reference from which the homology
branches derive, the *ab initio* branch appears to over-extend its models
(2.02 % against 0.04 % for LiftOn). Measured against the external reference, the
four branches are indistinguishable (1.67 to 2.15 %).

The three homology branches align almost perfectly against their own source:
their near-zero values do not measure quality, they measure identity with
themselves. The contrast holds at thresholds of 10, 15, 20 and 30 percentage
points.

**Evaluating a derived annotation against its own reference does not measure
quality. The measuring standard must be independent of the method evaluated.**

### Positional novelty

Of the 18,765 mRNA predicted by Helixer, 891 do not overlap any locus of the
union of the three homology branches, and 790 exceed 1 kb. Of those 790, **555**
have an orthologue in some camelid and are homology-supported candidates that
transfer failed to place; **235** have none, and of these only **10** match
Swiss-Prot, so the remaining 225 are low-confidence candidates.

Positional novelty is not biological novelty.

### The combined candidate reference proteome

**This is the product of the work.** LiftOn core (20,233) plus the 555 rescued
loci, minus **80** redundant with the core at 95 % identity and 80 % query
coverage: **20,708 sequences**.

Of the 80 excluded, **76 correspond to a single core model each**. The remaining
four form two pairs of contiguous or overlapping loci that align with the same
model, which indicates **fragmentation** of a gene in the *ab initio* prediction
rather than duplication.

The threshold is more permissive than the exact identity used by Kourelis et al.
(2019), because that work compared proteomes of different species whereas here
two methods are compared on the same genome. At exact identity only 42 redundant
models are detected; the product varies by 70 sequences between the strictest and
the most permissive threshold, 0.34 %.

`scripts/build_combined_proteome.sh` reproduces the deposited file **byte for
byte**, MD5 `dc69fa820facc0f087696bdd4885e1c1`, and fails with a non-zero exit
code if the result does not contain 20,708 sequences.

### Confidence layers

The deposit includes `Lgla_combined_proteome_confidence.tsv`, which classifies
each of the 20,708 entries into three nested layers:

| Layer | Criterion | n |
|---|---|---|
| `high_confidence` | LiftOn model with intact reading frame and homology in *C. dromedarius* | 18,920 |
| `extended_reference` | Remaining LiftOn models | 1,313 |
| `extended_candidate` | Loci from the *ab initio* prediction | 475 |

Within the intermediate layer, 545 models have a broken reading frame and 956
lack detectable external homology, 188 meeting both conditions. The external
homology criterion therefore discards more than twice as many models as the
reading-frame criterion. Of the 545 with an internal stop, 357 still retain
homology with *C. dromedarius*.

**Internal stop codons were not filtered.** Liftoff marks them in its own GFF3
with `valid_ORF=False` and `missing_start_codon=True`, so a filtered subset can
be derived in one command.

### Genomes (configurable in `workflow/config/config.yaml`)

| Role | Species | Accession |
|---|---|---|
| Target | *Lama glama* | `GCA_028534125.1` |
| Reference | *Vicugna pacos* | `GCF_048564905.1` |
| External evaluation | *Camelus dromedarius* | `GCF_036321535.1` |

## Quality control & reporting

- **BUSCO** in protein mode (default lineage `artiodactyla_odb12`, n = 12,594) on
  every proteome, with an optional secondary lineage.
- **AGAT** structural statistics per annotation.
- **gffcompare** concordance between the Liftoff and miniprot branches.
- **InterProScan** functional coverage.
- **DIAMOND** coverage against the external and the circular reference.
- **DIAMOND** against RepeatMasker's `RepeatPeps.lib`: transposable-element protein
  check of every proteome and of the novel loci (the assembly was not masked).
- **BUSCO with the substrate equalised**: the homology arms recomputed on the
  scaffolds the ab initio arm analysed, so that the two are comparable
  (`scripts/busco_substrate_restricted.py`).
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
│   ├── README.md                    # what each script produces and its caveats
│   ├── build_report.py              # comparative report (called by the `report` rule)
│   ├── run_lifton.sh                # LiftOn branch wrapper (manual)
│   ├── run_helixer.py               # first Helixer run (SUPERSEDED; kept for traceability)
│   ├── run_ips.sh                   # InterProScan run
│   ├── estado_ips.sh                # InterProScan progress monitor
│   ├── make_interproscan_input.py   # MD5 deduplication and batching for InterProScan
│   ├── run_blastp.sh                # DIAMOND against both reference proteomes
│   ├── estado_blastp.sh             # DIAMOND progress monitor
│   ├── analyze_blastp.py            # coverage and identity statistics
│   ├── run_swissprot.sh             # Swiss-Prot search for loci without camelid orthologue
│   ├── run_repeatpeps.sh            # TE-protein check against RepeatPeps.lib (seven sets)
│   ├── repeatpeps_summary.py        # summary of the TE-protein check; additional file 5
│   ├── busco_substrate_restricted.py # BUSCO with the substrate equalised; additional file 1
│   ├── overlap_md5.py               # exact-sequence overlap between the four proteomes
│   ├── overlap_coords.py            # positional novelty against the homology branches
│   ├── novel_loci_blast.py          # characterisation of the novel loci
│   ├── build_combined_proteome.sh   # builds the combined candidate reference proteome
│   ├── structural_stats.py          # structural statistics per annotation
│   ├── internal_stops.py            # internal stop codons per proteome
│   ├── camelidae_landscape.py       # annotation landscape across Camelidae
│   ├── capas.py                     # confidence layers of the combined proteome
│   ├── s6.py                        # supplementary table S6, reference correspondence
│   ├── s7.py                        # supplementary table S7, over-extension threshold
│   ├── tabla1_final.py              # manuscript Table 1
│   ├── tabla_estructura.py          # manuscript Table 2
│   ├── tabla_funcional.py           # manuscript Table 4
│   ├── fig1_busco.py                # Figure 1, BUSCO completeness
│   ├── fig2_upset.py                # Figure 2, sequence overlap
│   ├── fig3_estructura.py           # Figure 3, structural properties
│   └── fig4_cobertura.py            # Figure 4, coverage and over-extension
├── .gitattributes             # forces LF endings; a CRLF .sh does not run on Linux
├── REPRODUCIBILITY.md         # exact commands, automated and manual branches
├── CITATION.cff
├── .zenodo.json
├── LICENSE
└── .gitignore
```

**Figures.** The four figure scripts write at 170 mm width (6.69 inches) so that
the declared font sizes are the final ones. **Do not reintroduce
`bbox_inches="tight"`**: it expands the canvas on save, which scales the text
below the 7 pt minimum required by the target journal.

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) 9.23.1 (the version used to
  produce the results reported in the accompanying manuscript)
- [conda](https://docs.conda.io/) / mamba — per-rule environments are resolved
  automatically with `--use-conda`
- Network access for the data-acquisition rules (NCBI `datasets`)
- InterProScan requires OpenJDK 11; system version 17 is not compatible

Versions of every tool used are listed in `additional_file_4.csv` within the data
deposit.

## Usage

```
# from the workflow/ directory
cd workflow

# dry run: resolves the DAG and reports every rule that would run,
# without downloading anything or executing a single tool
snakemake -n

# full run with conda environments
snakemake --use-conda --cores 8
```

The dry run is the check to make before committing to a full run: it fails fast on
a malformed `config.yaml`, a missing input or a broken rule graph. It does not
validate the tools themselves, and it does not cover the manual stages listed in
[Automation boundary](#automation-boundary-important), which are run from
`scripts/`.

Input genomes are downloaded automatically from NCBI using the accessions in
`config.yaml`; edit that file to point the workflow at different assemblies,
BUSCO lineages or parameters.

## Data availability

Raw and intermediate data (genomes, proteomes, annotations, BUSCO downloads,
reports) are **not** included in this repository by design. They are regenerated
by running the workflow, or obtained from the deposit listed in the
[Data](#data) section above and from INSDC.

## Citation

If you use this workflow, please cite it via the metadata in `CITATION.cff`, or
through its Zenodo DOI:

- Concept DOI, always resolving to the latest version:
  [10.5281/zenodo.21456816](https://doi.org/10.5281/zenodo.21456816)

Cite the version DOI when referring to the exact state of the code used for a
given analysis, and the concept DOI otherwise. Version DOIs are listed on the
Zenodo record.

## License

Released under the [MIT License](LICENSE).
