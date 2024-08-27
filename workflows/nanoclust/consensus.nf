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
    tag "${meta.scenario.count}/${meta.scenario.rep} - Cluster ${meta.cluster.id}"
    container "quay.io/biocontainers/fastani:1.34--h4dfc31f_3"
    label "med_cpu"
    label "med_mem"

    input:
        tuple val(meta), path(reads_by_cluster_fasta_gz)

    output:
        tuple val(meta), path('reads_by_cluster.fasta'), path('*draft_read.fasta', arity: '1'), emit: draft
        tuple val(meta), path('fastani_output.ani'), emit: fastani_output
        tuple val(meta), path('*.log'), emit: logs

    script:
    """
    gunzip -c $reads_by_cluster_fasta_gz > reads_by_cluster.fasta
    split -l 2 reads_by_cluster.fasta split_reads
    find split_reads* > read_list.txt

    fastANI --ql read_list.txt --rl read_list.txt \\
        -o fastani_output.ani \\
        -t ${task.cpus} \\
        -k 16 --fragLen 160 \\
        2> fastani.log

    DRAFT=\$(awk 'NR>0{name[\$1] = \$1; arr[\$1] += \$3; count[\$1] += 1}  END{for (a in arr) {print arr[a] / count[a], name[a] }}' fastani_output.ani | sort -rg | cut -d " " -f2 | head -n1)
    cat \$DRAFT > cluster_${meta.cluster.id}_draft_read.fasta
    """
}

process mapReadsToDraft {
    tag "${meta.scenario.count}/${meta.scenario.rep} - Cluster ${meta.cluster.id}"
    container "quay.io/biocontainers/minimap2:2.28--he4a0461_1"
    label "small_cpu"
    label "med_mem"

    input:
        tuple val(meta), path(reads_by_cluster_fasta), path(draft_read_fasta)
    output:
        tuple val(meta), path(reads_by_cluster_fasta), path(draft_read_fasta), path('*aligned.sam', arity: '1'), emit: aligned
        tuple val(meta), path('*.log'), emit: logs

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont \\
        --no-long-join \\
        -r100 -a $draft_read_fasta $reads_by_cluster_fasta \\
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
        tuple val(meta), path(reads_by_cluster_fasta), path(draft_read_fasta), path(aligned_sam)

    output:
        tuple val(meta), path(reads_by_cluster_fasta), path(draft_read_fasta), path('*racon_consensus.fasta', arity: '1'), emit: racon_output
        tuple val(meta), path('*.log'), emit: logs

    script:
    """
    success=1
    racon -t ${task.cpus} \\
        --quality-threshold=9 \\
        -w 250 $reads_by_cluster_fasta $aligned_sam $draft_read_fasta \\
        > cluster_${meta.cluster.id}_racon_consensus.fasta \\
        2> racon.log
    """
}

process medakaConsensus {
    tag "${meta.scenario.count}/${meta.scenario.rep} - Cluster ${meta.cluster.id}"
    container "docker.io/ontresearch/medaka:sha3486abaab0d3b90351617eb8622acf2028edb154"
    label "large_mem"
    label "small_cpu"

    input:
        tuple val(meta), path(reads_by_cluster_fasta), path(draft_read_fasta), path(racon_consensus_fasta)

    output:
        tuple val(meta), path('*_medaka/consensus.fasta', arity: '1'), emit: consensus
        tuple val(meta), path('*.log'), emit: logs

    script:
    def model="r1041_e82_400bps_sup_g615"
    """
    medaka_consensus -i $reads_by_cluster_fasta -d $racon_consensus_fasta \\
        -o cluster_${meta.cluster.id}_medaka \\
        -t ${task.cpus} \\
        -m $model \\
        > medaka.log
    """

}