process ITSXPRESS {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/itsxpress:2.0.1--pyhdfd78af_0  ' :
        'biocontainers/itsxpress:2.0.1--pyhdfd78af_0' }"

    input:
        tuple val(meta), path(fastq) // can be fastq or fastq.gz
        val region

    output:
        tuple val(meta), path("*.fastq.gz"), emit: reads
        tuple val(meta), path("*.log"), emit: log


    script:
    if (!(region in ["ITS1", "ITS2", "FULL_ITS"])) {
        error "Region not supported or recognized: $region"
    }

    def cmd_region = ''
    if ("$region" == "FULL_ITS") {
        cmd_region = "ALL"
    } else {
        cmd_region = "$region"
    }
    def outfile = "${fastq.baseName}".replaceAll(/.fastq$/, '')
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
        --outfile ${outfile}.${region}.fastq.gz
    """
}