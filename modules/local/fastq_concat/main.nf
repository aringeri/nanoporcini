process FASTQ_CONCAT {
    label 'process_lo'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path("reads*.fastq.gz")
    
    output:
    tuple val(meta), path( "*.fastq.gz" ), emit: merged_reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gunzip -c reads*.fastq.gz | gzip > merged_reads.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed -n 1p | sed 's/GNU bash, version //g')
    END_VERSIONS
    """
}