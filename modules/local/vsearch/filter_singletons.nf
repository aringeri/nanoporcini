process filter_singletons {
    tag "$meta.id - $meta.region"
    label "large_cpu"
    label "large_mem"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-53dae514294fca7b44842b784ed85a5303ac2d80:7b3365d778c690ca79bc85aaaeb86bb39a2dec69-0':
        'biocontainers/mulled-v2-53dae514294fca7b44842b784ed85a5303ac2d80:7b3365d778c690ca79bc85aaaeb86bb39a2dec69-0' }"

    input:
        tuple val(meta), path(centroids_fasta_gz)

    output:
        tuple val(meta), path('*.non-singleton.fasta.gz'), emit: non_singleton

    script:
    def centroids_fname="$centroids_fasta_gz"
    if (!("$centroids_fname" ==~ /.+\.fasta.gz/)) {
        error("Expected $centroids_fname to be a '.fasta.gz' file")
    }
    def prefix="$centroids_fname".replaceFirst(/.fasta.gz$/, '') + '.non-singleton'
    """
    vsearch \\
        --sortbysize <(gunzip -c $centroids_fasta_gz) \\
        --threads ${task.cpus} \\
        --sizein \\
        --sizeout \\
        --fasta_width 0 \\
        --minsize 2 \\
        --output ${prefix}.fasta

    gzip -n ${prefix}.fasta
    """
}