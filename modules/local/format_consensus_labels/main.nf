process FORMAT_CONSENSUS_LABELS {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(zipped_reads)
    
    output:
    tuple val(meta), path( "*.fasta.gz" ), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gunzip -c $zipped_reads \\
        | awk 'BEGIN{FS=";" } {if (/>/) { print ">"substr(\$1, index(\$1, "=") + 1, length(\$1) )} else {print}}' \\
        | gzip > renamed_reads.fasta.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed -n 1p | sed 's/GNU bash, version //g')
    END_VERSIONS
    """
}