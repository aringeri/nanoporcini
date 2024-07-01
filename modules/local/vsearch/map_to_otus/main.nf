process VSEARCH_MAP_READS_TO_OTUS {
    tag "$meta.id - $meta.region"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
        tuple val(meta), path(all_reads), path(otus)
 
    output:
        tuple val(meta), path("*.tsv"), emit: otu_tab
        tuple val(meta), path("*.uc"), emit: uc
        tuple val(meta), path("*.log"), emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    vsearch \\
        --usearch_global $all_reads \\
        --db $otus \\
        --id 0.97 \\
        --threads $task.cpus \\
        --strand plus \\
        --sizein \\
        --sizeout \\
        --fasta_width 0 \\
        --qmask none \\
        --dbmask none \\
        --log vsearch.log \\
        --uc ${prefix}.uc \\
        --otutabout ${prefix}.otus.tsv
    """
}
