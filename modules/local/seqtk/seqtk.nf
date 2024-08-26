process seqtk_sample {
    tag "$meta.id (${meta.scenario.count})"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_2':
            'biocontainers/seqtk:1.4--he4a0461_2' }"
    label "small_cpu"
    label "large_mem"

    input:
        tuple val(meta), path(fastx)

    output:
        tuple val(meta), path("${prefix}.*"), emit: subsampled
        path("*.txt"), emit: seed

    script:
    prefix = task.ext.prefix ?: "${fastx.baseName}"
    """
    seqtk sample \\
        -s ${meta.scenario.seed} \\
        $fastx \\
        ${meta.scenario.count} \\
        | gzip > ${prefix}.subsampled.fastq.gz

    echo "${meta.scenario.seed}" > seed.txt

    """
}