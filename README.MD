# Nanoporcini

A pipeline for long-read metabarcoding fungi with Oxford Nanopore Technology.

# Requirements

- [Nextflow](https://www.nextflow.io/) v23.10.0 or later
- A container platform. Either:
  - [Docker](https://www.docker.com/get-started/)
  - or [Singularity](https://docs.sylabs.io/guides/latest/user-guide/)
- A network connection to download dependencies

# Quick start

## Run with default parameters
```bash
nextflow run main.nf \
  --input "data/*.fastq.gz" \
  --outdir output/
```

## Run with custom primers
```bash
nextflow run main.nf \
  --input "data/*.fastq.gz" \
  --primers.fwd "AACTTAAAGGAATTGACGGAAG" \ 
  --primers.rev_rc "GGTAAGCAGAACTGGCG" \
  --outdir output/
```

## Run with custom parameters
```bash
nextflow run main.nf \
  --input "data/*.fastq.gz" \
  -params-file conf/params.yml \
  --outdir output/
```

## Configure resource use (CPU/RAM)
```bash
nextflow run main.nf \
  --input "data/*.fastq.gz" \
  -c conf/envs/local.config \
  --outdir output/
```

# Configuration

See [conf/params.yml](conf/params.yml) as an example of how to configure the parameters of the pipeline.

## Inputs

- FastQ files, one file per sample, file names will be used to identify samples
  - `input`
- Primers 
  - `primers.fwd` - The forward primer sequence 5'-3'.
  - `primers.rev_rc` - The reverse primer sequence (which has been reverse complimented)

## Quality Filtering

- FULL_ITS:
  - minQualityPhred: 20
  - minLength: 300
  - maxLength: 6000
- qc_plot_sample_level
- chimera_filtering:
  - ref_db

## Clustering

- methods
  - VSEARCH
    - min cluster size
  - NanoCLUST
    - min cluster sizes

## Taxonomic Assignments

- dnabarcoder
  - ref_db
  - ref_classifications
  - cutoffs

