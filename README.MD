# Nanoporcini

Nanoporcini is an experimental pipeline for long-read [metabarcoding](https://en.wikipedia.org/wiki/Metabarcoding) of fungi with Oxford Nanopore Technologies sequencing.


# Features 

Nanoporcini is written in the [Nextflow](https://www.nextflow.io/) workflow language and features:

- Container support for Docker or Singularity.
- Configurable quality filtering with [cutadapt](https://github.com/marcelm/cutadapt/) and [chopper](https://github.com/wdecoster/chopper).
- Full ITS region extraction with [itsxpress](https://github.com/USDA-ARS-GBRU/itsxpress).
- Chimera detection with [VSEARCH](https://github.com/torognes/vsearch/) 
- A choice of clustering approaches:
  - [VSEARCH](https://github.com/torognes/vsearch/)
  - Custom implementation of [NanoCLUST](https://github.com/genomicsITER/NanoCLUST/)
- Taxonomic classifications using [dnabarcoder](https://github.com/vuthuyduong/dnabarcoder)

# Requirements

- [Nextflow](https://www.nextflow.io/) v23.10.0 or later
- A container platform. Either:
  - [Docker](https://www.docker.com/get-started/)
  - or [Singularity](https://docs.sylabs.io/guides/latest/user-guide/)
- A network connection to download dependencies
- A reference database containing ITS sequences
  - [UNITE+INSD 2024](https://zenodo.org/records/12580255)

# Quick start

## Run with minimal required parameters
```bash
nextflow run https://github.com/aringeri/nanoporcini \
  --input "path/to/data/sample*.fastq.gz" \
  --outdir "nanoporcini_output/" \
  --primers.fwd "AACTTAAAGGAATTGACGGAAG" \ 
  --primers.rev_rc "GGTAAGCAGAACTGGCG" \
  --chimera_filtering.ref_db "path/to/unite/db.fasta" \
  --taxonomic_assignment.dnabarcoder.ref_db "path/to/unite/db.fasta" \
  --taxonomic_assignment.dnabarcoder.ref_classifications "path/to/unite/db.classification" \
  --taxonomic_assignment.dnabarcoder.cutoffs "path/to/unite/db/cutoffs.best.json" 
```
This example uses NS5 forward primer and LR6 reverse primer (see <https://unite.ut.ee/primers.php>).

## Test run with example data

A small example can be run with two samples and 20 reads each. (see [tests/test_data/example/fq](tests/test_data/example/fq)).
This takes ~10 minutes to run on a 4 core macbook pro.

First clone the repository and move into the nanoporcini directory:

```bash
git clone https://github.com/aringeri/nanoporcini
cd nanoporcini 
```

Download the UNITE reference database using the command below or do it through your browser at <https://zenodo.org/records/12580255>. 
This may take a while:

```bash
mkdir -p data/db/unite2024/
curl https://zenodo.org/api/records/12580255/files-archive -o data/db/unite2024/unite2024.zip
unzip data/db/unite2024/unite2024.zip
```

Run the examples:

```bash
nextflow run main.nf \
  -params-file "tests/test_data/example/example-params.yml" \
  --outdir "nanoporcini_output/"
```

# Configuration

Many of the pipeline parameters can be configured to suit your needs.
I recommend using creating a [yaml](https://en.wikipedia.org/wiki/YAML) file to pass parameters to the pipeline.
See [conf/params.yml](conf/params.yml) for an example configuration file that can be given to the run command with the `-params-file` option:

```bash
nextflow run https://github.com/aringeri/nanoporcini \
  --input "path/to/data/sample*.fastq.gz" \
  -params-file conf/params.yml \
  --outdir "nanoporcini_output/"
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
  - The reverse primer sequence (which has been reverse complemented).
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
nextflow run https://github.com/aringeri/nanoporcini \
  --input "path/to/data/sample*.fastq.gz" \
  -c conf/envs/local.config \
  --outdir "nanoporcini_output/"
```

# Further reading

This workflow was developed as part of my Masters degree. 
My thesis can be read online ([https://aringeri.github.io/long-read-ITS-metabarcoding-thesis/](https://aringeri.github.io/long-read-ITS-metabarcoding-thesis/)),
and it outlines the details of the pipeline, reasoning behind design choices and the validation approaches taken.
