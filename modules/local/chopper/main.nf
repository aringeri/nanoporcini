process CHOPPER {
    tag "$meta.id $meta.region"

    container "biocontainers/chopper:0.7.0--hdcf5f25_0"

    input:
        tuple val(meta), val(args), path(fastq)
        

    output:
        tuple val(meta), path("*.fastq.gz"), emit: reads
        tuple val(meta), path("*.log")     , optional: true, emit: log

    script:
    def prefix = "$fastq".replaceAll(/.fastq.gz$/, '')
    if ("${prefix}.filtered.fastq.gz" == "$fastq") {
        error "input file matches output file which would causing overwrite: $fastq"
    }
    println(args)

    """
    echo "running chopper with args: $args" > ${prefix}.log

    gunzip -c $fastq | \\
    chopper \\
        --quality $args.minQualityPhred \\
        ${args.minLength != null ? "--minlength $args.minLength" : ""} \\
        --threads $task.cpus \\
        2> >(tee -a ${prefix}.log >&2) \\
        | gzip -n > ${prefix}.filtered.fastq.gz
    """
}
