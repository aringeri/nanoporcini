process seqkit_add_barcode {
    tag "$meta.id"
    label "small_cpu"
    label "large_mem"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/seqkit:2.6.1--h9ee0642_0':
            'biocontainers/seqkit:2.6.1--h9ee0642_0' }"

    input:
        tuple val(meta), path(fastx)

    output:
        tuple val(meta), path("${prefix}.*"), emit: fasta

    script:
    prefix = task.ext.prefix ?: "${fastx.simpleName}"

    def extension = "fastq"
    if ("$fastx" ==~ /.+\.fasta|.+\.fasta.gz|.+\.fa|.+\.fa.gz|.+\.fas|.+\.fas.gz|.+\.fna|.+\.fna.gz|.+\.fsa|.+\.fsa.gz/ ) {
        extension   = "fasta"
    }
    extension       = fastx.toString().endsWith('.gz') ? "${extension}.gz" : extension
    def call_gzip   = extension.endsWith('.gz')
    if ("${prefix}.with_barcode.${extension}" == "$fastx") {
        error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    }

    """
    seqkit \\
        replace -p '(^\\S*)' -r '\$1;barcodelabel=${meta.id}' \\
        --threads $task.cpus \\
        $fastx \\
        ${call_gzip ? '| gzip -c' : ''} \\
        > ${prefix}.with_barcode.${extension}
    """
}