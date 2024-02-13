process SEQKIT_FQ2FA {
    tag "$meta.id"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.6.1--h9ee0642_0':
        'biocontainers/seqkit:2.6.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${prefix}.fasta.gz"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqkit \\
        fq2fa \\
        --threads $task.cpus \\
        $fastq \\
        -o ${prefix}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}