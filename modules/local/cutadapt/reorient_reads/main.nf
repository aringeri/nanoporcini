process CUTADAPT_REORIENT_READS {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.6--py39hf95cd2a_1' :
        'biocontainers/cutadapt:4.6--py39hf95cd2a_1' }"

    input:
        tuple val(meta), path(fastq)
    output:
        tuple val(meta), path("*.reoriented.fastq.gz"), emit: reads
        tuple val(meta), path("*.unmatched.fastq.gz"), emit: unmatched
        tuple val(meta), path("*.partial-match.fastq.gz"), emit: partial
        tuple val(meta), path("cutadapt.log"), emit: log
    

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        fwd_primer = "CTTGGTCATTTAGAGGAAGTAA" // ITSF1 (5'->3')
        rev_primer_rc = "CCCGTCTTGAAACACGG" // LRU3 (reverse complement)
        // error = "1"
        action = "retain"
        """
        cutadapt \\
            -a '$fwd_primer;required...$rev_primer_rc;required' \\
            --action=$action \\
            --rc \\
            --cores=$task.cpus \\
            -o ${prefix}.reoriented.fastq.gz \\
            --untrimmed-output ${prefix}.nomatch.fastq.gz \\
            $fastq \\
            > cutadapt.log

        echo 'Searching for partial matches:' >> cutadapt.log

        cutadapt \\
            -g '$fwd_primer' \\
            -a '$rev_primer_rc' \\
            --action=$action \\
            --rc \\
            --cores=$task.cpus \\
            -o ${prefix}.partial-match.fastq.gz \\
            --untrimmed-output ${prefix}.unmatched.fastq.gz \\
            ${prefix}.nomatch.fastq.gz \\
            >> cutadapt.log
        """
}