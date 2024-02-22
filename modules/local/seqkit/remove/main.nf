process SEQKIT_REMOVE_CHIMERAS {
tag "$meta.id"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.6.1--h9ee0642_0':
        'biocontainers/seqkit:2.6.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq), path(denovo_chimeras), path(ref_chimeras)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: nonchimeras
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqkit \\
        grep \\
        --threads $task.cpus \\
        -v -f <(seqkit seq -n -i $denovo_chimeras $ref_chimeras) \\
        $fastq \\
        -o ${prefix}.nonchimeras.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}