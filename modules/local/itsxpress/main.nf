process ITSXPRESS {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/itsxpress:2.0.1--pyhdfd78af_0  ' :
        'biocontainers/itsxpress:2.0.1--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(fastq)
        val region

    output:
        tuple val(meta), path("*.fastq.gz"), emit: reads
        tuple val(meta), path("*.log"), emit: log


    script:
    def prefix = "FULL_ITS"
    def cmd_region = ""
    if (region == Regions.ITS1) {
        cmd_region = "ITS1"
    } else if (region == Regions.ITS2) {
        cmd_region = "ITS2"
    } else if (region == Regions.FULL_ITS) { 
        cmd_region = "ALL"
    } else {
        error "Region not supported or recognized: $region"
    }
    """
    rm -r tmp/ || true
    mkdir tmp/

    itsxpress \\
        --fastq $fastq \\
        --single_end \\
        --tempdir ./tmp/ \\
        --single_end --taxa Fungi \\
        --region $cmd_region \\
        --threads $task.cpus \\
        --outfile ${fastq.baseName}.${region}.fastq.gz
    """
}