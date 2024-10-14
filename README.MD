# Nanoporcini

A pipeline for long-read metabarcoding fungi with Oxford Nanopore Technology.

# Requirements

- [Nextflow](https://www.nextflow.io/) v23.10.0 or later
- A container platform. Either:
  - [Docker](https://www.docker.com/get-started/)
  - or [Singularity](https://docs.sylabs.io/guides/latest/user-guide/)
- A network connection to download dependencies
- A reference database containing ITS sequences
  - [UNITE+INSD 2024](https://zenodo.org/records/12580255)

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

NS5 forward primer and LR6 reverse primer.

# Configuration

Many of the pipeline parameters can be configured to suit your needs.
I recommend using creating a [yaml](https://en.wikipedia.org/wiki/YAML) file to pass parameters to the pipeline.
See [conf/params.yml](conf/params.yml) for an example configuration file that can be given to the run command with the `-params-file` option:

```bash
nextflow run main.nf \
  --input "data/*.fastq.gz" \
  -params-file conf/params.yml \
  --outdir output/
```

## Inputs

### Samples

This pipeline is expecting ONT sequences to be basecalled and demultiplexed already.
I recommend [Dorado](https://github.com/nanoporetech/dorado) for these steps.

- `input` - **Required**
  - type - `string`
  - A set of [fastq](https://en.wikipedia.org/wiki/FASTQ_format) files. One file per sample. File names will be used to identify samples. Uses glob syntax (`*`) to select multiple files.
  - ex) `"path/to/input/data/*.fastq.gz"`

### Primer sequences
- `primers.fwd` 
  - type - `string`
  - The forward primer sequence in 5' to 3' orientation.
  - ex) `"AACTTAAAGGAATTGACGGAAG"` for NS5 primer (see https://unite.ut.ee/primers.php)
- `primers.rev_rc` 
  - type - `string`
  - The reverse primer sequence (which has been reverse complimented)
  - ex) `"GGTAAGCAGAACTGGCG"` for LR6 primer (see https://unite.ut.ee/primers.php)

## Outputs

Specify the output directory with:
- `outdir`
- `qc_plot_sample_level`

## Quality Filtering

- FULL_ITS:
  - minQualityPhred: 20
  - minLength: 300
  - maxLength: 6000
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

## Customising Resource Use (CPU/RAM)

For controlling the number of threads or RAM usage of various tasks see the configuration files in [conf/envs](conf/envs).
These will allow you to use the available resources on your system and improve the pipeline runtime.
Pass this to the run command with the `-c` option:

```bash
nextflow run main.nf \
  --input "data/*.fastq.gz" \
  -c conf/envs/local.config \
  --outdir output/
```


