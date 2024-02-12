process ITSXPRESS {
    label 'process_medium'

    conda "bioconda::itsxpress=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/itsxpress:2.0.1--pyhdfd78af_0  ' :
        'biocontainers/itsxpress:2.0.1--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), path("*.fastq.gz"), emit: reads
        tuple val(meta), path("*.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    itsxpress \\
        --fastq $fastq \\
        $args \\
        --threads $task.cpus \\
        --outfile ${prefix}.fastq.gz
    """
}