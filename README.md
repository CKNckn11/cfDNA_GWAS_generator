# cfDNA GWAS Pipeline Generator

**cfDNA GWAS Pipeline Generator** is a lightweight, semi-automated pipeline generator for cfDNA-derived germline GWAS analysis. It generates standalone shell scripts from a YAML configuration file and a sample table, supporting both BAM-based and FASTQ-based inputs.

## Introduction

This tool is designed to generate reproducible analysis scripts rather than directly execute all jobs inside Python. After generation, each shell script can be submitted, checked, modified, or rerun independently on a local server or computing cluster.

```bash
$ python generate_pipeline.py -h
usage: cfdna-gwas-generate [-h] [-v] --config CONFIG [--outdir OUTDIR]

cfDNA GWAS Pipeline Generator (Version = 0.1.0): Generate cfDNA germline GWAS shell scripts from a YAML config.

optional arguments:
  -h, --help       show this help message and exit
  -v, --version    show the version of cfdna-gwas-generate and exit.
  --config CONFIG  YAML config file.
  --outdir OUTDIR  Override project.outdir from the YAML config.
```

## Main Processes

The generator supports two cfDNA GWAS input modes.

**1. BAM input mode**

Use this mode when aligned BAM files are already available.

```text
metadata -> bam_qc -> gvcf -> joint -> plink -> pca -> gwas -> plot
```

**2. FASTQ input mode**

Use this mode when raw paired-end FASTQ files need to be aligned first.

```text
metadata -> fastq_qc -> align -> markdup -> bqsr -> bam_qc -> gvcf -> joint -> plink -> pca -> gwas -> plot
```

The generated project contains independent shell scripts for each step, together with a wrapper script:

```bash
bash 00.shell/run_pipeline.sh help
bash 00.shell/run_pipeline.sh all
```

## Installation

Run directly from the source directory:

```bash
cd /400T/ckn/cfDNA_GWAS_generator
python generate_pipeline.py -h
```

Optional local installation:

```bash
python -m pip install -e .
```

After installation, the command-line entry point is:

```bash
cfdna-gwas-generate -h
```

## Quick Start

Generate a BAM-input GWAS project:

```bash
cd /400T/ckn/cfDNA_GWAS_generator
python generate_pipeline.py --config examples/config.bam.example.yaml
```

Generate a FASTQ-input GWAS project:

```bash
python generate_pipeline.py --config examples/config.fastq.example.yaml
```

Or after installation:

```bash
cfdna-gwas-generate --config examples/config.bam.example.yaml
```

Run the generated project:

```bash
cd /400T/ckn/cfDNA_GWAS_generated_bam
bash 00.shell/run_pipeline.sh help
bash 00.shell/run_pipeline.sh all
```

Logs are written to `loginfo/` in the generated project.

## Input Tables

BAM input table:

```text
sample_id    bam    bai    phenotype    age    sex
```

FASTQ input table:

```text
sample_id    fastq1    fastq2    phenotype    age    sex
```

`iid` can be used instead of `sample_id`, and `diagnosis_type` can be used instead of `phenotype`.
