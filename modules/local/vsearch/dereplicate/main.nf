process VSEARCH_DEREPLICATE {
    tag "$meta.id - $meta.region"
    label 'small_cpu'
    label 'med_mem'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
        tuple val(meta), path(fastx)
 
    output:
        tuple val(meta), path("*.fast*.gz"), emit: reads
        tuple val(meta), path("*.uc"), emit: uc
        tuple val(meta), path("*.log"), emit: logs

    script:

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
    outfile = "${prefix}.$ext"

    """
    vsearch \\
        --fastx_uniques $fastx \\
        $out_cmd $outfile \\
        --strand plus \\
        --sizein \\
        --sizeout \\
        --uc ${prefix}.uc \\
        --log vsearch.log \\
        --threads $task.cpus
    
    gzip $outfile
    """
}
