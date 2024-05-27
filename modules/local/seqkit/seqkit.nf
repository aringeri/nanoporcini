package local.seqkit

process seqkit {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/seqkit:2.6.1--h9ee0642_0':
            'biocontainers/seqkit:2.6.1--h9ee0642_0' }"

    input:
        val(args)
        tuple val(meta), path(fastx)

    output:
        tuple val(meta), path("${prefix}.*"), emit: fasta

    script:
    prefix = task.ext.prefix ?: "${fastx.baseName}"

    def extension = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension       = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    def call_gzip   = extension.endsWith('.gz')
    if ("${prefix}.${extension}" == "$fastx") {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }

    """
    seqkit \\
        ${args} \\
        --threads $task.cpus \\
        $fastx \\
        ${call_gzip ? '| gzip -c' : ''} \\
        > ${prefix}.${extension}
    """
}