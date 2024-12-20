workflow nanoclust_consensus {
    take:
        // one file per cluster, per repetition, per sample
        // expecting meta map to contain 'cluster.id' property
        reads_by_cluster_fastq_gz // Channel<(Map<cluster_id>, Reads_FQ_GZ)>

    main:
        draft = draftSelection(reads_by_cluster_fastq_gz)
        aligned = mapReadsToDraft(draft.draft).aligned
        racon = raconConsensus(aligned).racon_output
        fmt_racon = removeReverseComplementTagFromSeqs(racon)
        medaka = medakaConsensus(fmt_racon).consensus

        reps = medaka.map { meta, cluster_consensus -> [ meta.subMap('region', 'scenario', 'cluster_method', 'umap_min_cluster_size'), cluster_consensus ] }
                .groupTuple()
                .map { meta, cluster_consensus -> [ meta + [id: "${meta.scenario.count}/${meta.scenario.rep}"], cluster_consensus ] }
        consensus = concatConsensusSeqs(reps)

    emit:
        consensus = consensus
}

process readCorrection {
    memory { 7.GB }
    container "quay.io/biocontainers/canu:2.2--ha47f30e_0"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path('corrected_reads.correctedReads.fasta'), emit: corrected_reads

    script:
    genomeSize=1000
    minReadLength=400
    stopOnLowCoverage=1
    minInputCoverage=2
    """
    canu -correct -p corrected_reads -nanopore-raw \\
    $reads genomeSize=${genomeSize} \\
    stopOnLowCoverage=${stopOnLowCoverage} minInputCoverage=${minInputCoverage} \\
    minReadLength=${minReadLength} minOverlapLength=200

    gunzip corrected_reads.correctedReads.fasta.gz
    """
}

process draftSelection {
    tag "${meta.scenario.count}/${meta.scenario.rep} - ${meta.umap_min_cluster_size} - Cluster ${meta.cluster.id}"
    container "quay.io/biocontainers/fastani:1.34--h4dfc31f_3"
    label "small_cpu"
    label "med_mem"

    input:
        tuple val(meta), path(reads_by_cluster_fastq_gz)

    output:
        tuple val(meta), path('reads_by_cluster.fastq'), path('*draft_read.fastq', arity: '1'), emit: draft
        tuple val(meta), path('fastani_output.ani'), emit: fastani_output
        tuple val(meta), path('*.log'), emit: logs

    script:
    def num_polishing_reads=params.consensus.num_polishing_reads
    """
    gunzip -c $reads_by_cluster_fastq_gz | head -n \$(($num_polishing_reads*4)) > reads_by_cluster.fastq
    split -l 4 -a 4 reads_by_cluster.fastq split_reads
    find split_reads* > read_list.txt

    fastANI --ql read_list.txt --rl read_list.txt \\
        -o fastani_output.ani \\
        -t ${task.cpus} \\
        -k 16 --fragLen 160 \\
        2> fastani.log

    DRAFT=\$(awk 'NR>0{name[\$1] = \$1; arr[\$1] += \$3; count[\$1] += 1}  END{for (a in arr) {print arr[a] / count[a], name[a] }}' fastani_output.ani | sort -rg | cut -d " " -f2 | head -n1)
    cat \$DRAFT > cluster_${meta.cluster.id}_draft_read.fastq
    """
}

process mapReadsToDraft {
    tag "${meta.scenario.count}/${meta.scenario.rep} - Cluster ${meta.cluster.id}"
    container "quay.io/biocontainers/minimap2:2.28--he4a0461_1"
    label "small_cpu"
    label "med_mem"

    input:
        tuple val(meta), path(reads_by_cluster_fastq), path(draft_read_fastq)
    output:
        tuple val(meta), path(reads_by_cluster_fastq), path(draft_read_fastq), path('*aligned.sam', arity: '1'), emit: aligned
        tuple val(meta), path('*.log'), emit: logs

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont \\
        --no-long-join \\
        -r100 -a $draft_read_fastq $reads_by_cluster_fastq \\
        -o cluster_${meta.cluster.id}_aligned.sam \\
        2> minimap.log
    """
}

process raconConsensus {
    tag "${meta.scenario.count}/${meta.scenario.rep} - Cluster ${meta.cluster.id}"
    container "quay.io/biocontainers/racon:1.5.0--h21ec9f0_2"
    label "large_mem"
    label "small_cpu"

    input:
        tuple val(meta), path(reads_by_cluster_fastq), path(draft_read_fastq), path(aligned_sam)

    output:
        tuple val(meta), path(reads_by_cluster_fastq), path(draft_read_fastq), path('*racon_consensus.fasta', arity: '1'), env(success), emit: racon_output
        tuple val(meta), path('*.log'), emit: logs

    script:
    """
    success=0
    if [[ \$(racon -t ${task.cpus} \\
            --quality-threshold=9 \\
            -w 250 $reads_by_cluster_fastq $aligned_sam $draft_read_fastq \\
            > cluster_${meta.cluster.id}_racon_consensus.fasta \\
            2> racon.log) -eq 0 \\
            && -s cluster_${meta.cluster.id}_racon_consensus.fasta ]]; then
        success=1
    else
        cat $draft_read_fastq \\
            | head -n 2 \\
            | sed '1 s/^@/>/' > cluster_${meta.cluster.id}_racon_consensus.fasta
    fi
    """
}

process medakaConsensus {
    tag "${meta.scenario.count}/${meta.scenario.rep} - Cluster ${meta.cluster.id}"
    container "docker.io/ontresearch/medaka:sha3486abaab0d3b90351617eb8622acf2028edb154"
    label "large_mem"
    label "small_cpu"

    input:
        tuple val(meta), path(reads_by_cluster_fastq), path(draft_read_fastq), path(racon_consensus_fasta), val(success)

    output:
        tuple val(meta), path('*_medaka/consensus.fasta', arity: '1'), emit: consensus
        tuple val(meta), path('*.log'), emit: logs

    script:
    def model="r1041_e82_400bps_sup_g615"
    if (success == "0") {
        log.warn """Racon correction for cluster ${meta.cluster.id} failed due to not enough overlaps. Taking draft read as consensus"""
    }
    """
    medaka_consensus -i $reads_by_cluster_fastq -d $racon_consensus_fasta \\
        -o cluster_${meta.cluster.id}_medaka \\
        -t ${task.cpus} \\
        -m $model \\
        > medaka.log
    """

}

process removeReverseComplementTagFromSeqs {
    tag "${meta.scenario.count}/${meta.scenario.rep} - Cluster ${meta.cluster.id}"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'https://depot.galaxyproject.org/singularity/seqkit:2.6.1--h9ee0642_0':
                'biocontainers/seqkit:2.6.1--h9ee0642_0' }"
    label "med_mem"
    label "small_cpu"

    input:
        tuple val(meta), path(reads_by_cluster_fastq), path(draft_read_fastq), path(racon_consensus_fasta), val(success)

    output:
        tuple val(meta), path('reads_by_cluster_no_rc.fastq'), path(draft_read_fastq), path(racon_consensus_fasta), val(success)

    // Reverse complement ('rc') has been added to the ids of all reads after using cutadapt.
    // This is causing issues in medaka, so trim these off for now
    script:
    """
    seqkit replace -p '\\src' -r '' $reads_by_cluster_fastq -o reads_by_cluster_no_rc.fastq
    """
}

process concatConsensusSeqs {
    tag "${meta.id} - ${meta.region}"
    label 'large_mem'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        tuple val(meta), path(input_files, name: "*.fasta", stageAs: "dir/*.fasta")

    output:
        tuple val(meta), path( "consensus_sequences.fasta" ), emit: merged

    script:

    def input_files = "$input_files".tokenize()

    """
    cat dir/*.fasta > consensus_sequences.fasta
    """
}
