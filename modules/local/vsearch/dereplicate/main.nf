process VSEARCH_DEREPLICATE {
    tag "${meta.id}"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
    tuple val(meta), path(fastx)
 
    output:
    tuple val(meta), path("*.fast*"), emit: reads
    tuple val(meta), path("*.log"), emit: logs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$fastx".endsWith(".fastq.gz")) {
        out_cmd = "--fastqout"
        ext = "fastq"
    } else if ("$fastx".endsWith(".fasta.gz")) {
        out_cmd = "--fastaout"
        ext = "fasta"
    } else {
        error "input file must have ending '.fasta.gz' or 'fastq.gz': $fastx"
    }

    """
    vsearch \\
        --fastx_uniques $fastx \\
        $out_cmd ${prefix}.$ext \\
        --sizein \\
        --sizeout \\
        --log vsearch.log \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
