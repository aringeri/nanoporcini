process VSEARCH_FILTER_MAX_EE {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
        tuple val(meta), path(fastq)
        val(maxee_rate)
 
    output:
        tuple val(meta), path("output/*.fastq.gz"), emit: filtered_reads
        tuple val(meta), path("*.log"), emit: logs

    script:

    if (!"$fastq".endsWith(".fastq.gz")) {
        error "input file must have '.fastq.gz' ending: $fastq"
    }

    """
    mkdir -p output/

    vsearch --fastx_filter $fastq \\
        --fastq_qmax 90 \\
        --fastq_ascii 33 \\
        --fastq_maxee_rate $maxee_rate \\
        --fastqout output/${fastq.baseName} \\
        --fastq_minlen 500 \\
        --threads $task.cpus \\
        --log output.log

    gzip output/${fastq.baseName}
    """
}
