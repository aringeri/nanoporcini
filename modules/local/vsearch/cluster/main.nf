process VSEARCH_CLUSTER_A {
    tag "$meta.id - $meta.region"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-53dae514294fca7b44842b784ed85a5303ac2d80:7b3365d778c690ca79bc85aaaeb86bb39a2dec69-0':
        'biocontainers/mulled-v2-53dae514294fca7b44842b784ed85a5303ac2d80:7b3365d778c690ca79bc85aaaeb86bb39a2dec69-0' }"

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path('*.aln.gz')                , optional: true, emit: aln
        tuple val(meta), path('*.biom.gz')               , optional: true, emit: biom
        tuple val(meta), path('*.mothur.tsv.gz')         , optional: true, emit: mothur
        tuple val(meta), path('*.otu.tsv.gz')            , optional: true, emit: otu
        tuple val(meta), path('*.out.tsv.gz')            , optional: true, emit: out
        tuple val(meta), path('*.blast.tsv.gz')          , optional: true, emit: blast
        tuple val(meta), path('*.uc.tsv.gz')             , optional: true, emit: uc
        tuple val(meta), path('*.centroids.fasta.gz')    , optional: true, emit: centroids
        tuple val(meta), path('*.consensus.fasta.gz')    , optional: true, emit: consensus
        tuple val(meta), path('*.clusters.fasta*.gz')    , optional: true, emit: clusters
        tuple val(meta), path('*.profile.txt.gz')        , optional: true, emit: profile
        tuple val(meta), path('*.msa.fasta.gz')          , optional: true, emit: msa
        tuple val(meta), path('*.log')                   , optional: true, emit: log


    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (!args2.contains("--cluster_fast") && !args2.contains("--cluster_size") && !args2.contains("--cluster_smallmem") && !args2.contains("--cluster_unoise") ) {
        error "Unknown clustering option provided (${args2})"
    }

    def ext_map = [ "--alnout" : "aln" ,
                    "--biomout" : "biom" ,
                    "--blast6out" : "blast.tsv" ,
                    "--centroids" : "centroids.fasta" ,
                    "--clusters" : "clusters.fasta" ,
                    "--consout" : "consensus.fasta" ,
                    "--mothur_shared_out" : "mothur.tsv" ,
                    "--msaout" : "msa.fasta" ,
                    "--otutabout" : "otu.tsv" ,
                    "--profile" : "profile.txt" ,
                    "--uc" : "uc.tsv" ,
                    "--userout" : "out.tsv"]

    outputs = args3.tokenize().collect { 
            ext = ext_map[(it)]
            if (ext == null) {
                error "Unknown output file format provided (${it})"
            }
            "${it} ${prefix}.${ext}"
        }.join(" ")

    """
    vsearch \\
        $args2 $fasta \\
        $outputs \\
        --threads $task.cpus \\
        $args


    gzip -n ${prefix}.*
    """
}
