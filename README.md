# Human miRNA Analysis Pipeline

A Nextflow-based workflow for preprocessing small RNA sequencing data, identifying known and candidate novel microRNAs (miRNAs), quantifying expression, performing differential expression analysis, and conducting downstream target, network, and functional enrichment analyses.

## Overview

This project implements a modular workflow for human miRNA sequencing analysis using Nextflow DSL2 and containerized bioinformatics tools.

The workflow is organized into three main stages:

1. **Preprocessing:** Raw read quality assessment, adapter trimming, post-trimming QC, and MultiQC reporting.
2. **miRNA identification and quantification:** Known miRNA analysis and novel miRNA prediction using miRDeep2.
3. **Downstream analysis:** Expression analysis, target prediction, network analysis, and functional enrichment.

The pipeline is intended to support reproducible analysis of small RNA sequencing datasets. Its current implementation uses human reference resources, including GRCh38 and miRNA reference FASTA files.

> **Development status:** The workflow structure and modules are under review. Individual stages, their dependencies, and end-to-end reproducibility should be validated before the pipeline is considered production-ready.

## Workflow

```text
Single-end FASTQ files
          |
          v
   Raw read QC
     FastQC
          |
          v
  Adapter trimming
     cutadapt
          |
          v
  Post-trimming QC
     FastQC
          |
          v
     MultiQC
          |
          v
  Reference alignment
      Bowtie
          |
          v
       miRDeep2
     /          \
    v            v
Known miRNAs   Novel miRNA candidates
    |            |
    v            v
Quantification  Quantification
    |            |
    v            v
Differential    Differential
expression      expression
    |            |
    v            v
Target prediction and downstream analyses
    |
    +--> Interaction/network analysis
    +--> GO enrichment
    +--> KEGG enrichment
    +--> Reactome enrichment
```

*The diagram represents the intended analytical flow inferred from the current module structure. Exact dependencies and execution order must be verified against `main.nf`.*

## Features

### 1. Quality control and preprocessing

- Raw read quality assessment using FastQC.
- Adapter trimming using cutadapt.
- Post-trimming quality assessment.
- Aggregated QC reporting using MultiQC.

### 2. Known miRNA analysis

The known miRNA branch includes modules intended for:

- Known miRNA expression quantification and count matrix generation.
- Differential expression analysis using R.
- miRNA target prediction.
- Interaction/network analysis.
- Functional enrichment analysis.

### 3. Novel miRNA discovery

The novel miRNA branch uses miRDeep2-related modules for alignment, prediction, and quantification.

Downstream modules are provided for:

- Novel miRNA count matrix generation.
- Differential expression analysis.
- Target prediction.
- Network analysis.
- Gene Ontology (GO) enrichment.
- KEGG pathway enrichment.
- Reactome pathway enrichment.

Predicted novel miRNAs should be treated as candidates requiring appropriate quality assessment and biological validation, rather than as experimentally confirmed miRNAs.

### 4. Containerized execution

The project includes Singularity container images and definition files for selected tools and analysis environments.

Containers are intended to help standardize software dependencies across compatible Linux systems. For reproducible deployment, tool versions, container sources, and image checksums should be documented.

## Repository structure

```text
.
├── main.nf
├── nextflow.config
├── metadata.csv
├── container/
│   ├── bowtie.sif
│   ├── cutadapt.sif
│   ├── fastqc.sif
│   ├── mirdeep2.sif
│   ├── multiqc.sif
│   ├── python.def
│   ├── python.sif
│   ├── r_mirna_analysis.def
│   └── r_mirna_analysis.sif
├── modules/
│   ├── preprocessing/
│   ├── mirdeep2/
│   ├── known_mirna/
│   └── novel_mirna/
├── reference/
│   ├── GRCh38.fa
│   ├── GRCh38_clean.fa
│   ├── hairpin.fa
│   ├── hairpin_hsa.fa
│   ├── hairpin_hsa_clean.fa
│   ├── mature.fa
│   ├── mature_hsa.fa
│   └── mature_hsa_clean.fa
├── r_scripts/
│   ├── 01_diff_exp_analysis.R
│   ├── 02_target_prediction.R
│   ├── 03_network.R
│   └── 04_enrichment.R
├── database/
│   └── stringdb/
└── test/
```

The repository may use externally stored references, databases, and containers in the final version. Large binary files and raw sequencing data should not be committed directly to Git.

## Requirements

- Linux environment.
- Nextflow.
- SingularityCE or Apptainer, according to the configured execution profile.
- Compatible container images for the required tools.
- Human genome reference compatible with the alignment and annotation workflow.
- Mature and precursor miRNA reference sequences.
- Experimental sample metadata.
- R and the packages required by the analysis scripts, available through the configured container environment.

Tool versions and package dependencies should be pinned and documented before release.

## Input data

### FASTQ files

The pipeline currently uses single-end FASTQ files in its example dataset. Paired-end support has not been established.

FASTQ files should be accessible from the execution environment, with valid paths supplied through the sample metadata or configured input parameters.

### Sample metadata

The project includes `metadata.csv`. Its exact required columns must match the parsing logic in `main.nf`.

A possible sample sheet design is:

```csv
sample,fastq,condition
sample_01,data/sample_01.fastq,control
sample_02,data/sample_02.fastq,control
sample_03,data/sample_03.fastq,treatment
sample_04,data/sample_04.fastq,treatment
```

This is an illustrative format only. Confirm the required column names, delimiters, and condition labels against the implemented pipeline before using it.

For differential expression analysis, biological replicates and clearly defined experimental groups are important. The metadata should identify the biological condition associated with each sample.

## Reference resources

### Reference Data

This directory contains reference files required by the human miRNA analysis pipeline. Large reference files are not tracked in this Git repository.

### Required references

| Expected filename | Description | Source |
|---|---|---|
| `GRCh38.fa` | Human reference genome, GRCh38 assembly | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) |
| `mature_hsa.fa` | Human mature miRNA sequences | [miRBase](https://www.mirbase.org/download/) |
| `hairpin_hsa.fa` | Human precursor miRNA sequences | [miRBase](https://www.mirbase.org/download/) |

These filenames are examples based on the current project structure. Confirm the exact filenames and formats required by the workflow before running it.

### Setup

1. Download the appropriate reference files from the official sources.
2. Select the correct genome assembly and human miRNA sequences.
3. Place the files in this directory using the expected filenames.
4. Check the configured paths in `nextflow.config`.
5. Confirm that any required alignment indexes are generated or available.

### Reproducibility

Record the genome assembly accession, miRNA database release, download date, and checksums for each reference file used in an analysis.

Reference versions should remain consistent across samples and analysis runs.

### Note

Reference files are not included in this repository because of their size and to keep cloning and version control efficient. Users are responsible for obtaining the required files from their respective sources and following applicable data-use terms.

### Installation

Clone the repository:

```bash
git clone <YOUR_GITHUB_REPOSITORY_URL>
cd <YOUR_REPOSITORY_DIRECTORY>
```

Check the execution environment:

```bash
nextflow -version
singularity --version
```

If your system uses Apptainer instead of SingularityCE:

```bash
apptainer --version
```

Review `nextflow.config` and update the relevant paths and execution settings for your system.

Container image locations, reference paths, and output directory settings must be configured according to the actual implementation.

## Running the pipeline

The exact parameter names and execution profiles should be confirmed from `main.nf` and `nextflow.config`.

A typical Nextflow invocation has the following form:

```bash
nextflow run main.nf \
    -profile singularity \
    -resume
```

This example assumes that a `singularity` profile exists and that required inputs and paths are configured elsewhere. Add the actual input, metadata, reference, and output parameters required by this project before running it on a new dataset.

For example, a project may expose parameters such as:

```bash
--input <FASTQ_OR_SAMPLE_SHEET>
--genome <GENOME_FASTA>
--outdir <OUTPUT_DIRECTORY>
```

These are placeholders, not verified parameter names for this repository.

For troubleshooting, inspect the Nextflow execution log and the work directory associated with the failed process. The `-resume` option can reuse cached results when task inputs and execution conditions remain compatible.

## Expected outputs

The following output categories are intended by the current module structure. Exact filenames and directories should be documented after a successful test run.

| Analysis stage | Expected output category |
|---|---|
| Raw QC | FastQC reports |
| Trimming | Filtered FASTQ files and trimming statistics |
| Post-trimming QC | FastQC reports |
| Aggregated QC | MultiQC report |
| Alignment and miRDeep2 | Alignment files and miRNA prediction/quantification results |
| Known miRNA analysis | Count matrix and differential expression results |
| Novel miRNA analysis | Candidate predictions, count matrix and differential expression results |
| Target analysis | Predicted miRNA–target associations |
| Network analysis | Interaction/network tables and visualizations, if implemented |
| Functional enrichment | GO, KEGG and Reactome results, where supported |

### Interpreting results

- **QC reports** help identify read quality issues and assess the effect of trimming.
- **Count matrices** summarize expression measurements across samples.
- **Differential expression results** identify miRNAs with evidence of expression differences between defined conditions.
- **Novel miRNA predictions** require assessment of sequence, precursor structure, read support, and relevant evidence.
- **Target predictions** are computational hypotheses and are not proof of direct regulation.
- **Network and enrichment analyses** depend on appropriate identifiers, database versions, and organism-specific mappings.

## Reproducibility and limitations

For reproducible analyses:

- Record Nextflow and software versions.
- Use consistent genome and miRNA reference releases.
- Document all parameters and thresholds.
- Preserve the sample metadata used for each run.
- Validate each workflow stage using a small test dataset.
- Keep raw data, generated outputs, and large container images out of the source repository unless a deliberate data-distribution strategy is in place.

The results depend on read quality, sequencing depth, sample design, reference completeness, software settings, and the statistical methods used.

A computational prediction should not be interpreted as biological validation. Differential expression, target prediction, network analysis, and pathway enrichment should be evaluated in the context of the study design and available evidence.

## Contributing

Contributions that improve reproducibility, testing, documentation, and workflow portability are welcome. Please include a description of the change and, where possible, a test demonstrating the expected behavior.
