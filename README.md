# cfDNA GWAS Pipeline Generator

This generator builds a standalone shell pipeline from:

- `config.yaml`
- `samples.tsv`

It supports two initial input types:

- `input.type: bam`
- `input.type: fastq`

## BAM Input Table

Required columns:

```text
sample_id	bam
```

Recommended columns:

```text
sample_id	bam	bai	phenotype	age	sex
```

`iid` can be used instead of `sample_id`. `diagnosis_type` can be used instead of `phenotype`.

## FASTQ Input Table

Required columns:

```text
sample_id	fastq1	fastq2
```

Recommended columns:

```text
sample_id	fastq1	fastq2	phenotype	age	sex
```

`iid` can be used instead of `sample_id`. `diagnosis_type` can be used instead of `phenotype`.

## Generate

Install in a virtual environment:

```bash
cd /400T/ckn/cfDNA_GWAS_generator
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Then generate with the installed command:

```bash
cfdna-gwas-generate --config examples/config.bam.example.yaml
```

The old direct-script command is also kept:

```bash
cd /400T/ckn/cfDNA_GWAS_generator
python3 generate_pipeline.py --config config.bam.example.yaml
```

or:

```bash
python3 generate_pipeline.py --config config.fastq.example.yaml
```

## Run Generated Pipeline

```bash
cd /400T/ckn/cfDNA_GWAS_generated_bam
bash 00.shell/run_pipeline.sh metadata
bash 00.shell/run_pipeline.sh bam_qc
```

Run all:

```bash
bash 00.shell/run_pipeline.sh all
```

## Generated Flow

BAM input:

```text
metadata -> bam_qc -> gvcf -> joint -> plink -> pca -> gwas -> plot
```

FASTQ input:

```text
metadata -> fastq_qc -> align -> markdup -> bqsr -> bam_qc -> gvcf -> joint -> plink -> pca -> gwas -> plot
```

Logs are written to `loginfo/` in the generated project.
