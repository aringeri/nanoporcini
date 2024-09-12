process dorado_trim_adapters {
    tag "$meta.id - $meta.region"
    label "small_cpu"
    label "med_mem"

    input:
        tuple val(meta), path(untrimmed)

    output:
        tuple val(meta), path("*adapter-trim.fq.gz"), emit: trimmed_reads

    script:
    def prefix = "${meta.id}.adapter-trim.fq.gz"
    if ("$prefix" == "$untrimmed") {
        error("Input and output names are the same")
    }

    """
    dorado trim \\
        --emit-fastq --no-trim-primers \\
        $untrimmed | gzip > $prefix
    """
}