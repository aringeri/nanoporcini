process SEQKIT_FQ2FA {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.6.1--h9ee0642_0':
        'biocontainers/seqkit:2.6.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${prefix}.fasta.gz"), emit: fasta

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqkit \\
        fq2fa \\
        --threads $task.cpus \\
        $fastq \\
        -o ${prefix}.fasta.gz
    """
}