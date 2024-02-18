process ITSX {
    label 'process_medium'
    // cpus 4

    conda "bioconda::itsx=1.1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/itsx:1.1.3--hdfd78af_1  ' :
        'biocontainers/itsx:1.1.3--hdfd78af_1' }"

    input:
        tuple val(meta), path(fasta_gz)

    output:
        tuple val(meta), path("*.ITS1.fasta.gz"), optional: true, emit: its1
        tuple val(meta), path("*.ITS2.fasta.gz"), optional: true, emit: its2
        tuple val(meta), path("*.full.fasta.gz"), optional: true, emit: full_its
        tuple val(meta), path("*_no_detections.fasta.gz"), optional: true, emit: no_detections
        tuple val(meta), path("*.chimeric.fasta.gz"), optional: true, emit: chimeric
        tuple val(meta), path("*.positions.txt"), optional: true, emit: positions
        tuple val(meta), path("*.problematic.txt"), optional: true, emit: problematic
        tuple val(meta), path("*.summary.txt"), optional: true, emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ITSx \\
        -i <(gunzip -c $fasta_gz) \\
        -t F \\
        --save_regions ITS1 \\
        --cpu $task.cpus \\
        --detailed_results T \\
        --preserve T \\
        -o ITSx_out
    
    gzip *.fasta
    """
}