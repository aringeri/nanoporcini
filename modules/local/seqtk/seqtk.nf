process seqtk_sample {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_2':
            'biocontainers/seqtk:1.4--he4a0461_2' }"

    input:
        tuple val(meta), path(fastx)

    output:
        tuple val(meta), path("${prefix}.*"), emit: fasta

    script:
    prefix = task.ext.prefix ?: "${fastx.baseName}"
    """
    seqtk sample \\
        -s ${params.sub_sample.seed} \\
        $fastx \\
        ${params.sub_sample.count} \\
        | gzip > ${prefix}.subsampled.fastq.gz

    """
}