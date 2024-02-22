process VSEARCH_MAP_READS_TO_OTUS {
    tag "${meta.id}"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
    tuple val(meta), path(all_reads)
    tuple val(meta), path(otus)
 
    output:
    tuple val(meta), path("*.tsv"), emit: otu_tab
    tuple val(meta), path("*.uc"), emit: uc
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    vsearch \\
        --usearch_global $all_reads \\
        --db $otus \\
        --id 0.98 \\
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
        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
