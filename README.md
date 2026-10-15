# ProkASM - Prokaryotic Genome Assembly Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![DSL2](https://img.shields.io/badge/DSL-2-blue.svg)](https://www.nextflow.io/docs/latest/dsl2.html)

A modular Nextflow DSL2 pipeline for de novo bacterial genome assembly from **Illumina**, **PacBio**, **Nanopore**, and **hybrid** sequencing data.

## Pipeline Overview

```
                  ┌─────────────────────────────────────────────────────┐
                  │                    ProkASM                          │
                  │     Prokaryotic Genome Assembly Pipeline            │
                  └─────────────────────────────────────────────────────┘
                                          │
           ┌──────────────────────────────┼──────────────────────────────┐
           │                              │                              │
           ▼                              ▼                              ▼
    ┌─────────────┐              ┌─────────────┐              ┌─────────────────┐
    │  Illumina   │              │  Long Read  │              │     Hybrid      │
    │  Workflow   │              │  Workflows  │              │    Workflow     │
    └─────────────┘              └─────────────┘              └─────────────────┘
           │                     ┌──────┴──────┐                      │
           │                     │             │                      │
           ▼                     ▼             ▼                      ▼
    ┌─────────────┐       ┌──────────┐  ┌──────────┐          ┌─────────────┐
    │   FastQC    │       │  PacBio  │  │ Nanopore │          │ FastQC +    │
    │   Fastp     │       │ [BAM2FQ] │  │ [BAM2FQ] │          │ Fastp +     │
    │   MEGAHIT/  │       │ NanoPlot │  │ NanoPlot │          │ NanoPlot +  │
    │  MetaSPAdes │       │ Hifiasm/ │  │   Flye   │          │  Unicycler  │
    └─────────────┘       │   Flye   │  │          │          └─────────────┘
           │              └──────────┘  └──────────┘                │
           └──────────────────────┴─────────────┴───────────────────┘
                                          │
                                          ▼
                              ┌─────────────────────┐
                              │    Post-Assembly    │
                              │  - Reorientation    │
                              │    (dnaapler/       │
                              │     circlator)      │
                              │  - BUSCO            │
                              └─────────────────────┘
```

## Features

- **Multi-platform support**: Illumina, PacBio, Nanopore, and hybrid assemblies
- **Flexible tool selection**: Multiple options for trimming and assembly
- **Unified samplesheet**: Single CSV format for all data types
- **Multi-lane/file support**: Automatically merges multiple files per sample (both Illumina and long reads)
- **Input flexibility**: Supports FASTQ and BAM formats for long reads
- **Skip options**: Skip individual pipeline steps as needed
- **Modular design**: Easy to extend with additional tools
- **Resource management**: Dynamic memory scaling and error retry
- **Comprehensive reporting**: Timeline, trace, and DAG generation

## Quick Start

```bash
# Clone the repository
git clone <repository-url>
cd prokasm

# Run with Illumina data
nextflow run main.nf \
    --data_type illumina \
    --samplesheet illumina.csv \
    --outdir results \
    -profile cluster

# Run with PacBio HiFi data (BAM format)
nextflow run main.nf \
    --data_type pacbio \
    --samplesheet pacbio.csv \
    --longread_format bam \
    --pacbio_assembler hifiasm \
    --outdir results \
    -profile cluster

# Run with Nanopore data (FASTQ format)
nextflow run main.nf \
    --data_type nanopore \
    --samplesheet nanopore.csv \
    --longread_format fastq \
    --flye_mode nano-hq \
    --outdir results \
    -profile cluster

# Run with Nanopore data (BAM format)
nextflow run main.nf \
    --data_type nanopore \
    --samplesheet nanopore_bam.csv \
    --longread_format bam \
    --flye_mode nano-hq \
    --outdir results \
    -profile cluster

# Run hybrid assembly (Illumina + Nanopore)
nextflow run main.nf \
    --data_type hybrid \
    --samplesheet hybrid.csv \
    --longread_type nanopore \
    --longread_format fastq \
    --outdir results \
    -profile cluster
```

## Samplesheet Format

The pipeline uses a unified CSV samplesheet with 4 columns:

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier (alphanumeric, underscores, hyphens) |
| `R1` | Path to forward Illumina reads (FASTQ/FASTQ.GZ) |
| `R2` | Path to reverse Illumina reads (FASTQ/FASTQ.GZ) |
| `LReads` | Path to long reads (BAM or FASTQ/FASTQ.GZ) |

### Multi-lane / Multi-file Support
To merge multiple files for the same sample, list the `sample_id` multiple times:

```csv
sample_id,R1,R2,LReads
sample_D,data/sample_D_L001_R1.fq.gz,data/sample_D_L001_R2.fq.gz,data/flowcell1.nanopore.bam
sample_D,data/sample_D_L002_R1.fq.gz,data/sample_D_L002_R2.fq.gz,data/flowcell2.nanopore.bam
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--data_type` | Type of data: `illumina`, `pacbio`, `nanopore`, or `hybrid` |
| `--samplesheet` | Path to CSV samplesheet |
| `--outdir` | Output directory (default: `./results`) |

### Platform-Specific Parameters

| Parameter | Description | Required When |
|-----------|-------------|---------------|
| `--longread_type` | Long read platform (`pacbio` or `nanopore`) | `--data_type hybrid` |
| `--longread_format` | Long read format (`bam` or `fastq`) | `pacbio`, `nanopore`, `hybrid` types |
| `--flye_mode` | Flye mode (`nano-raw`, `nano-hq`, `pacbio-hifi`, etc.) | Flye assembler |

### Tool Selection

| Parameter | Options | Default |
|-----------|---------|---------|
| `--illumina_trimmer` | `fastp`, `trimmomatic`, `trimgalore` | `fastp` |
| `--illumina_assembler` | `megahit`, `metaspades` | `megahit` |
| `--pacbio_assembler` | `hifiasm`, `flye` | `hifiasm` |
| `--reorientation_tool` | `dnaapler`, `circlator` | `dnaapler` |

### Skip Options

| Parameter | Description |
|-----------|-------------|
| `--skip_qc` | Skip QC steps (FastQC/NanoPlot) |
| `--skip_trimming` | Skip Illumina read trimming |
| `--skip_reorientation` | Skip contig reorientation |
| `--skip_busco` | Skip BUSCO assessment |

### Tool-Specific Parameters

Each tool has dedicated parameters:

```bash
# Fastp
--fastp_opts "--cut_front --cut_tail -l 50"

# Hifiasm
--hifiasm_k 51
--hifiasm_a 4
--hifiasm_opts "-l 3"

# Flye
--flye_mode nano-hq
--flye_opts "--genome-size 5m"

# Unicycler
--unicycler_mode normal
--unicycler_min_fasta_length 100
--unicycler_opts "--no_correct"

# BUSCO
--busco_mode genome
--busco_lineage bacteria_odb10
--busco_opts "--offline"
```

## Output Directory Structure

```
results/
├── qc/
│   ├── fastqc/           # FastQC reports
│   └── nanoplot/         # NanoPlot reports
├── trimming/
│   ├── fastp/            # Fastp outputs
│   ├── trimmomatic/      # Trimmomatic outputs
│   └── trimgalore/       # TrimGalore outputs
├── assembly/
│   ├── megahit/          # MEGAHIT assemblies
│   ├── metaspades/       # MetaSPAdes assemblies
│   ├── hifiasm/          # Hifiasm assemblies
│   ├── flye/             # Flye assemblies
│   └── unicycler_hybrid/ # Unicycler hybrid assemblies
├── dnaapler/             # Reoriented contigs (dnaapler)
├── circlator/            # Reoriented contigs (circlator)
├── busco/                # BUSCO results
├── alignment/            # Read alignments
├── bam2fq/               # Converted FASTQ from BAM
├── concatenated/         # Merged files (Illumina)
├── concatenated_longreads/ # Merged long read files
└── pipeline_info/        # Nextflow reports
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.html
```

## Profiles

| Profile | Description |
|---------|-------------|
| `standard` | Local execution |
| `cluster` | LSF cluster execution |
| `slurm` | SLURM cluster execution |
| `test` | Test mode with reduced resources |
| `debug` | Debug mode (no cleanup) |

## Requirements

### Software Dependencies

The pipeline requires the following tools (available via Singularity container or Conda):

- **QC**: FastQC, NanoPlot
- **Trimming**: fastp, Trimmomatic, Trim Galore
- **Assembly**: MEGAHIT, MetaSPAdes, Hifiasm, Flye, Unicycler
- **Post-processing**: dnaapler, Circlator, BUSCO
- **Utilities**: samtools, minimap2, pbindex, bam2fastq

### Nextflow Version

Requires Nextflow version **23.10.0** or later.

## Citation

If you use this pipeline, please cite:

- The pipeline itself
- The individual tools used in your analysis

## Author

**Qidong Jia**
