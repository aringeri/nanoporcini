process VSEARCH_UCHIME_REF {
    tag "$meta.id"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-53dae514294fca7b44842b784ed85a5303ac2d80:7b3365d778c690ca79bc85aaaeb86bb39a2dec69-0':
        'biocontainers/mulled-v2-53dae514294fca7b44842b784ed85a5303ac2d80:7b3365d778c690ca79bc85aaaeb86bb39a2dec69-0' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path('*.nonchimeras.fasta'), emit: nonchimeras
    tuple val(meta), path('*.chimeras.fasta'), emit: chimeras
    tuple val(meta), path('*.log'), optional: true, emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$fasta" == "${prefix}.nonchimeras.fasta") {
        error "input and output names are the same, set prefix in module configuration to disambiguate!"
    }

    """
    vsearch --uchime_ref $fasta \\
        --db $db \\
        --sizein \\
        --sizeout \\
        --fasta_width 0 \\
        --threads $task.cpus \\
        --log vsearch.log \\
        --chimeras ${prefix}.chimeras.fasta \\
        --nonchimeras ${prefix}.nonchimeras.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}