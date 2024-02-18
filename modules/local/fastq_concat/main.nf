process FASTQ_CONCAT {
    label 'process_lo'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(zipped_reads, name: "zipped*.gz", stageAs: "dir*/*")
    
    output:
    tuple val(meta), path( "*.gz" ), emit: merged_reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    prefix = task.ext.prefix ?: "merged_reads"
    input_files = "$zipped_reads".tokenize()
    if (input_files.every { it.endsWith("fastq.gz") }) {
        ext = "fastq"
    } else if (input_files.every { it.endsWith("fasta.gz") }) {
        ext = "fasta"
    } else {
        error "input list of zipped reads must have 'fastq.gz' or 'fasta.gz' ending: $zipped_reads"
    }

    """
    gunzip -c $zipped_reads | gzip > ${prefix}.${ext}.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed -n 1p | sed 's/GNU bash, version //g')
    END_VERSIONS
    """
}