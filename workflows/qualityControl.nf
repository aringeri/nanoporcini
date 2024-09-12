
def collectWithId(id, ch) {
  ch.collect { meta, read -> 
    read
  }.map {
    [ [id: id], it ]
  }
}

def collectByRegionWithId(String id, ch) {
  ch.map { meta, reads -> [ meta.subMap('region'), reads ] }
    .groupTuple()
    .map { meta, reads -> [ meta + [id: id], reads ] }
}


workflow qualityControl {
    take:
        stage_name
        reads // [ Map, FastqFile ]
    
    main:
        combined_reads = collectByRegionWithId(stage_name, reads)

        plotQualityProfile(combined_reads)
        nanoplot_bulk(combined_reads)

        if (params.qc_plot_sample_level) {
            nanoplot(reads.map{ meta, reads -> [ meta + [stage: stage_name], reads ] })
        }
}

process plotQualityProfile {
    tag "$meta.id"
    label "small_cpu"
    label "med_mem"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
            'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"
    when:
        params.qc_quality_profile
    input:
        tuple val(meta), path(fq_files)
    output:
        tuple val(meta), path("*.png"), emit: plots

    script:
    def r_files = "$fq_files".tokenize().collect{"\'$it\'"}.join(',')
    """
    #!/usr/bin/env Rscript
    library(dada2)
    library(ggplot2)

    plot <- plotQualityProfile(c($r_files))
    ggsave('${meta.id}.png')
    """
}

process nanoplot_bulk {
    tag "$meta.id"
    label "small_cpu"
    label "med_mem"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.41.6--pyhdfd78af_0' :
        'biocontainers/nanoplot:1.41.6--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(ontfiles, name: "*.fastq.gz")

    output:
        tuple val(meta), path("*.html")                , emit: html
        tuple val(meta), path("*.png") , optional: true, emit: png
        tuple val(meta), path("*.txt")                 , emit: txt
        tuple val(meta), path("*.log")                 , emit: log

    script:
    def args = task.ext.args ?: ''
    """
    NanoPlot \\
        $args \\
        -t $task.cpus \\
        --fastq $ontfiles
    """
}

process nanoplot {
    tag "$meta.stage $meta.id"
    label "small_cpu"
    label "med_mem"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.41.6--pyhdfd78af_0' :
        'biocontainers/nanoplot:1.41.6--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(ontfile)

    output:
        tuple val(meta), path("*.html")                , emit: html
        tuple val(meta), path("*.png") , optional: true, emit: png
        tuple val(meta), path("*.txt")                 , emit: txt
        tuple val(meta), path("*.log")                 , emit: log

    script:
    def args = task.ext.args ?: ''
    """
    NanoPlot \\
        $args \\
        -t $task.cpus \\
        --fastq $ontfile
    """
}