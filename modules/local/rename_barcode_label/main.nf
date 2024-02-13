process RENAME_BARCODE_LABEL {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(zipped_reads)
    
    output:
    tuple val(meta), path( "*.fast*.gz" ), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    //| awk 'NR%4==1 {sub(/barcode=/,"barcodelabel=")} {print}' \\

    // relabel each read to have the format @read-id;barcodelabel=$barcode;
    // this is so that vsearch can group clusters by sample id
    script:
    prefix = task.ext.prefix ?: "renamed_reads"
    input_files = "$zipped_reads".tokenize()
    if (input_files.every { it.endsWith("fastq.gz") }) {
        ext = "fastq"
    } else if (input_files.every { it.endsWith("fasta.gz") }) {
        ext = "fasta"
    } else {
        error "input list of zipped reads must have 'fastq.gz' or 'fasta.gz' ending: $zipped_reads"
    }

    """
    gunzip -c $zipped_reads \\
        | awk 'NR%4==1 {\$1=\$1";"substr(\$8,0,index(\$8,"=")-1)"label="substr(\$8,index(\$8,"=")+1,length(\$8))";"} {print}' \\
        | gzip > ${prefix}.${ext}.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed -n 1p | sed 's/GNU bash, version //g')
    END_VERSIONS
    """
}