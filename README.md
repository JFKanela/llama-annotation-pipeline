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
infrastructure** (Helixer, through the official web tool); they are *not* wired
into the `Snakefile`. The
exact commands to reproduce them are documented in
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md). The comparative report
(`scripts/build_report.py`) discovers whichever branches are present under
`results/busco/`, so all four appear in the final table regardless of how they
were produced.

The two downstream analyses added later — **functional annotation** with
InterProScan and **model accuracy** with DIAMOND — are outside the `Snakefile`
for the same reason, and are documented in their own sections below and in
`REPRODUCIBILITY.md`.

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

The *ab initio* branch was re-run after the first version of this repository was
released. **The annotation reported in the manuscript is the second run**, and
the two differ in method, substrate and results. `scripts/run_helixer.py`
documents the first run only and is retained for traceability, marked as
superseded.

| | First run (superseded) | Current run (manuscript) |
|---|---|---|
| Tool | `helixerlite` 25.5.27 on Kaggle | **official Helixer web tool v0.3.6** (plabipd.de) |
| Subsequence overlap | none | **enabled** (`vertebrate` lineage defaults) |
| Substrate | scaffolds >= 10 kb (3,640; 81.46 % of assembly) | scaffolds >= 25 kb (**244**; **79.6 %** of assembly) |
| BUSCO complete | 79.3 % | **85.5 %** |
| BUSCO fragmented | 6.6 % | **5.9 %** |
| BUSCO missing | 14.1 % | **8.6 %** |
| Proteins | 18,933 | **18,765** |

**Substrate.** The current run used scaffolds of at least 25 kb: 244 scaffolds,
holding 79.6 % of the assembled sequence. The three homology branches were run on
the complete assembly, so completeness is still not directly comparable between
Helixer and the homology branches.

**Subsequence overlap.** Prediction used the overlapping mode that the web tool
applies by default for the `vertebrate` lineage: subsequence length 213,840 bp,
overlap offset 106,920 bp, overlap core length 160,380 bp. Overlapping mitigates
the loss of accuracy at subsequence boundaries and is the main reason for the
6.2-point gain in BUSCO completeness over the first run, which had no overlap.

**Post-processing.** HelixerPost was left at its defaults: `window-size` 100,
`edge-threshold` 0.1, `peak-threshold` 0.8, `min-coding-length` 60. The web
interface does not allow other values of `peak-threshold` to be explored, so no
parameter tuning was performed.

**Non-determinism (important).** The web tool does not enable Helixer's
deterministic mode, so GPU kernel non-determinism is not suppressed. Two
consecutive runs on identical input yielded **18,765 and 18,763 genes**, with
micro-differences at model boundaries. The 18,765-gene run was adopted. **This
annotation is therefore not reproducible bit for bit**, and any attempt to
regenerate it should expect differences of this order.

Full parameters and provenance are in
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md), section 4.

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

## Functional annotation (InterProScan)

The four proteomes were annotated with **InterProScan 5.78-109.0**, running the
18 analyses enabled by default. Like LiftOn and Helixer, this step is **not
wired into the `Snakefile`**; the exact commands are in
[`REPRODUCIBILITY.md`](REPRODUCIBILITY.md).

**Input: deduplication by MD5.** The four proteomes hold 80,331 sequences
between them, but only **42,364 are distinct**. Annotating the unique set
instead of the union saves **47.3 %** of the compute. The deduplicated FASTA
uses the MD5 digest as the sequence identifier, which is what makes the mapping
file **`md5_to_ids.tsv` indispensable**: without it the results cannot be
attributed back to each proteome, because the annotated identifiers no longer
carry any branch information.

**Precalculated match lookup disabled (`-dp`).** The lookup service resolves a
query by exact MD5 against UniProtKB. Llama proteins are not in UniProtKB —
which is the premise of this work — so the service can only ever miss. It was
switched off rather than left to fail silently.

**Java.** InterProScan 5.78-109.0 requires **OpenJDK 11**. OpenJDK 17 does not
work.

**Input cleaning.** The three homology proteomes contain the `.` character,
gffread's rendering of an in-frame stop, because many transferred models have a
broken reading frame. InterProScan aborts at `stepLoadFromFastaIntoDB` when it
finds one. Each `.` is therefore **replaced by `X`, not deleted**: substitution
preserves sequence length and hence the domain coordinates, whereas deletion
would shift every downstream position.

| Branch | Proteins | With domain | % | With Pfam | % |
|---|---|---|---|---|---|
| liftoff | 20,073 | 18,127 | 90.3 | 17,379 | 86.6 |
| miniprot | 21,260 | 19,921 | 93.7 | 19,056 | 89.6 |
| lifton | 20,233 | 18,996 | 93.9 | 18,254 | 90.2 |
| helixer | 18,765 | 17,735 | 94.5 | 16,913 | 90.1 |

**754,835 annotations** over **40,402** of the 42,364 unique proteins
(**95.4 %**). "With domain" **excludes MobiDBLite and Coils**, which predict
disorder and coiled-coil structure respectively and annotate almost any protein;
counting them would make every branch look near-perfect.

> **PANTHER is reported at family level only.** The TSV output never emits the
> subfamily (`PTHR12345:SF6`), not even when phylogenetic placement succeeds.
> Verified against InterProScan's own official test file, which yields zero
> subfamilies in TSV and the five expected ones in XML, and whose run produced no
> abort. Subfamily resolution requires a structured output format and was not
> used here.

**Caution — the denominators are not equivalent.** Helixer predicts 18,765
proteins against the 20,000–21,000 of the other three, and it only emits a model
where it can build a complete ORF. Its protein set is therefore pre-filtered for
structural plausibility, so its percentage starts from an advantage. The
column to compare across branches is the count, not the ratio, unless that
asymmetry is stated alongside it.

## Model accuracy against external references (DIAMOND)

The four proteomes were aligned against two reference proteomes with **DIAMOND
2.2.4** in `blastp --very-sensitive` mode, with `--max-target-seqs 1
--max-hsps 1 --evalue 1e-5`. This is the fourth quality criterion of Kourelis
et al. (2019), and the only one in this work that evaluates **model accuracy**
rather than completeness or internal consistency. Like the previous step, it is
not automated in the `Snakefile`.

**Why DIAMOND and not BLASTP.** Llama and alpaca diverged a few million years
ago and their orthologues sit at 95–99 % identity. These are trivial alignments
with no sensitivity problem to solve, so the accelerated heuristic costs nothing
that matters here and returns the whole comparison in **18 minutes instead of
days**.

**Two references, and the distinction between them is the point of this
section.**

| Reference | Accession | Role |
|---|---|---|
| *Camelus dromedarius* | `GCF_036321535.1` (mCamDro1.pat, RS_2024_04, 50,982 proteins) | **External yardstick. The valid comparison** |
| *Vicugna pacos* | `GCF_048564905.1` | Internal control. **Circular** |

Results against *C. dromedarius*:

| Branch | With hit | % | Identity | Subject >=80 % | Query >=80 % |
|---|---|---|---|---|---|
| liftoff | 18,423 | 91.8 | 95.6 | 85.5 | 89.8 |
| miniprot | 20,519 | 96.5 | 94.5 | 82.7 | 88.7 |
| lifton | 19,277 | 95.3 | 95.5 | 84.2 | 92.6 |
| helixer | 17,945 | 95.6 | 92.6 | 81.7 | 91.4 |

> **The reference chosen decides the conclusion.** Measured against alpaca,
> Helixer appears to over-extend its models twenty times more often than LiftOn
> (2.0 % versus 0.0 %). Measured against dromedary, all four branches sit at the
> same 2 % and are indistinguishable.
>
> The cause is circularity. Three of the four branches are projections of the
> alpaca annotation and align almost perfectly against their own source, so their
> 0.0 % measures identity with themselves, not model quality. Only the *ab
> initio* branch is independent of that reference, and it is the only one the
> comparison penalises.
>
> **Had alpaca been used as the yardstick, this work would have reported that the
> *ab initio* branch over-extends its models, and that claim would have been
> false.** The alpaca figures are retained as an internal control, and must not be
> read as a measure of accuracy.

Two further caveats:

- **Do not report the median coverages.** They come out at 100 % for all four
  branches against both references: they saturate and discriminate nothing. Use
  the fraction of proteins with coverage >=80 %, which is what the tables above
  report. (`scripts/analyze_blastp.py` still computes and prints the medians;
  they are diagnostic output, not a result.)
- **One HSP per pair.** `--max-hsps 1` means that in multidomain proteins whose
  alignment fragments, coverage is recorded as artificially low. The approximation
  is conservative in the right direction: it underestimates, it does not inflate.

## Positional novelty of the *ab initio* branch

The three homology branches are projections of the *Vicugna pacos* annotation and
cannot, by construction, place a gene that is absent from it. Only Helixer can.
The question this section answers is how many of those Helixer-only loci are real
genes and how many are prediction artefacts — because the honest answer to that
is what makes the rest of the resource credible.

Of the Helixer transcripts with **no positional overlap** with any homology
model, **790 loci exceed 1 kb**. Cross-referencing them against the DIAMOND
alignments already computed for section "Model accuracy":

| | Loci | % |
|---|---|---|
| Total, > 1 kb, no positional overlap | 790 | 100 |
| **With an orthologue in another camelid** | **545** | **69.0** |
| Without a camelid orthologue | 245 | 31.0 |
| — of those, with a Swiss-Prot homologue | 10 | 1.3 |
| — of those, with no homologue anywhere | ~235 | 29.7 |

Read plainly:

- **545 loci (69.0 %) are real genes that homology transfer failed to place.**
  They have an orthologue in another camelid, so their existence does not depend
  on trusting the *ab initio* predictor. This is the substantive gain of running
  a reference-independent method alongside the projections.
- **~235 loci have no homologue in any camelid and none in Swiss-Prot.** With
  two independent searches returning nothing, the most probable explanation is
  that they are **false positives of the prediction**, not lineage-specific
  genes. They are reported as such.

Only **10** of the 245 orphans matched Swiss-Prot, which is what makes the
false-positive reading hard to avoid: a genuinely novel protein-coding gene in a
camelid would be unlikely to have no detectable homologue in either a sister
species or a manually curated reference database.

Two caveats belong with these numbers:

- **Positional novelty is not sequence novelty.** A locus counts here because no
  homology model overlaps its coordinates, which can also happen when the
  transfer placed the same gene elsewhere. The MD5 overlap analysis
  (`scripts/overlap_md5.py`) measures a different thing — exact sequence identity
  — and its counts are an order of magnitude larger. They are not
  interchangeable.
- **The 1 kb threshold is a choice, not a property of the data.** Shorter loci
  were excluded because below that length the prediction is dominated by
  fragments; a different threshold gives a different denominator.

The commands are in [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md), section 8.

## Figures and tables

The four manuscript figures and the summary tables are generated by scripts in
`scripts/`, outside the `Snakefile`, from the outputs of the analyses above.
They write to `MANUSCRITO/01_figuras/` (PDF and PNG at 300 dpi) and
`MANUSCRITO/02_tablas/` (markdown).

| Script | Output |
|---|---|
| `fig1_busco.py` | BUSCO completeness of the four proteomes, with the reference ceiling and the in-house baseline |
| `fig2_upset.py` | UpSet diagram of proteome overlap, 14 intersections |
| `fig3_estructura.py` | Coding exons per transcript and single-exon fraction |
| `fig4_cobertura.py` | Accuracy against the external reference, including the circularity finding |
| `tabla_estructura.py` | Section 2 of Table 1 |
| `tabla_funcional.py` | Table 2, functional annotation per branch |
| `tabla1_final.py` | Assembles the final Table 1 |

The only external dependencies are **matplotlib** (the four figures) and
**numpy** (`fig4_cobertura.py` alone). The table scripts use the standard library
only. See `REPRODUCIBILITY.md` section 9 for the environment and for why
`upsetplot` is not used.

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
│   ├── README.md              # what each script produces and its caveats
│   ├── build_report.py        # builds the comparative report (called by the `report` rule)
│   ├── run_lifton.sh          # LiftOn branch wrapper (manual; see REPRODUCIBILITY.md sec. 3)
│   ├── run_helixer.py         # first Helixer run (SUPERSEDED; kept for traceability)
│   ├── overlap_md5.py         # exact-sequence overlap between the four proteomes
│   ├── overlap_coords.py      # positional novelty / loss of Helixer against the homology branches
│   ├── structural_stats.py    # structural statistics per annotation (length, CDS/transcript, monoexonic)
│   ├── internal_stops.py      # internal stop codons per proteome
│   ├── camelidae_landscape.py # annotation landscape across Camelidae (species selection)
│   ├── make_interproscan_input.py  # MD5 deduplication and batching for InterProScan
│   ├── run_ips.sh             # InterProScan in resumable batches
│   ├── estado_ips.sh          # InterProScan progress monitor
│   ├── run_blastp.sh          # DIAMOND alignment against two reference proteomes
│   ├── analyze_blastp.py      # coverage analysis and per-branch assignment
│   ├── estado_blastp.sh       # DIAMOND progress monitor
│   ├── novel_loci_blast.py    # characterises Helixer loci with no positional overlap
│   ├── fig1_busco.py          # Figure 1 — BUSCO completeness
│   ├── fig2_upset.py          # Figure 2 — proteome overlap (UpSet, matplotlib only)
│   ├── fig3_estructura.py     # Figure 3 — coding exons per transcript, single-exon fraction
│   ├── fig4_cobertura.py      # Figure 4 — accuracy against the external reference
│   ├── tabla_estructura.py    # Table 1, section 2
│   ├── tabla_funcional.py     # Table 2 — functional annotation per branch
│   └── tabla1_final.py        # assembles the final Table 1
├── .gitattributes             # forces LF endings; a CRLF .sh does not run on Linux
├── REPRODUCIBILITY.md         # exact commands: automated branches, manual LiftOn / Helixer,
│                              #   InterProScan (sec. 6) and DIAMOND (sec. 7)
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
