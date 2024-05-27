process FindReadsWithID {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/seqkit:2.6.1--h9ee0642_0':
            'biocontainers/seqkit:2.6.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(centroids_query)
    tuple val(region_meta), path(ref_reads)

    output:
    tuple val(region_meta), path("*.fastq.gz"), emit: matching_reads

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqkit \\
        grep \\
        --threads $task.cpus \\
        -f <(seqkit seq --name --only-id --id-regexp '^[^\\s]+\\s([^\\s,^;]+)' $centroids_query) \\
        $ref_reads \\
        -o ${prefix}.matching.fastq.gz
    """
}
