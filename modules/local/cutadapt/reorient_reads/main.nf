process CUTADAPT_REORIENT_READS {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.6--py39hf95cd2a_1' :
        'biocontainers/cutadapt:4.6--py39hf95cd2a_1' }"

    input:
        tuple val(meta), path(fastq)
    output:
        tuple val(meta), path("*.fastq.gz"), emit: reads
        tuple val(meta), path("cutadapt.log"), emit: log
    

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        fwd_primer = "CTTGGTCATTTAGAGGAAGTAA" // ITSF1 (5'->3')
        rev_primer_rc = "CCCGTCTTGAAACACGG" // LRU3 (reverse complement)
        """
        cutadapt \\
            -a '$fwd_primer;required...$rev_primer_rc;required' \\
            --trimmed-only \\
            --action=none \\
            --rc \\
            --cores=$task.cpus \\
            -o ${prefix}.reoriented.fastq.gz \\
            $fastq \\
            > cutadapt.log
        """
}