process NANOPLOT_BULK {
  tag "$meta.id"
  label 'process_low'

  conda "${moduleDir}/environment.yml"
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
    path  "versions.yml"                           , emit: versions

  when:
    task.ext.when == null || task.ext.when

  script:
    def args = task.ext.args ?: ''
    """
    NanoPlot \\
        $args \\
        -t $task.cpus \\
        --fastq $ontfiles

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}