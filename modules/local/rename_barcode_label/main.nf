process RENAME_BARCODE_LABEL {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        tuple val(meta), path(zipped_reads)
    
    output:
        tuple val(meta), path( "*.fast*.gz" ), emit: reads


    // relabel each read to have the format @read-id;barcodelabel=$barcode;
    // this is so that vsearch can group clusters by sample id
    script:
    prefix = task.ext.prefix ?: "renamed_reads"
    input_files = "$zipped_reads".tokenize()
    if (input_files.every { it.endsWith("fastq.gz") || it.endsWith("fq.gz") }) {
        ext = "fastq"
        lines_per_record = 4
    } else if (input_files.every { it.endsWith("fasta.gz" || it.endsWith("fa.gz")) }) {
        ext = "fasta"
        lines_per_record = 2
    } else {
        error "input list of zipped reads must have 'fastq.gz', 'fq.gz', 'fasta.gz' or 'fa.gz' ending: $zipped_reads"
    }

    """
    gunzip -c $zipped_reads \\
        | awk 'NR%$lines_per_record==1 {\$1=\$1";barcodelabel="substr(\$8,index(\$8,"=")+1,length(\$8))";"} {print}' \\
        | gzip > ${prefix}.${ext}.gz
    """
}