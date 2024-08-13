include {
  UnGzip
 } from '../workflows/classify'

workflow nanoclust {
    take:
        sample_reads // Channel<Map, Fastq.GZ> one per repetition, per sample

    main:
        min_cluster_props = Channel.fromList(params.cluster.nanoclust.min_cluster_props)

        umap = UnGzip(sample_reads)
            | kmerFreq
            | umapTransform

        gatherMinClusterSizeStats(umap)

        clusters = hdbscanCluster(umap)

        createOtuTable(clusters)
        splitReadsByCluster(clusters.join(sample_reads)).reads_by_cluster_fastq_gz
            | findMostAbundantSeqsInCluster
}

process kmerFreq {
    tag "${meta.scenario.count}/${meta.scenario.rep}"
    container "docker.io/hecrp/nanoclust-kmer_freqs:latest"

    input:
        tuple val(meta), path(fastq)
    output:
        tuple val(meta), path("*.tsv"), emit: kmer_freqs

    script:
    """
    kmer_freq.py -t $task.cpus \\
        -k 6 \\
        -r $fastq > freqs.tsv
    """
}

process umapTransform {
    tag "${meta.scenario.count}/${meta.scenario.rep}"
    container "docker.io/hecrp/nanoclust-read_clustering:latest"
    containerOptions "${ workflow.containerEngine == 'singularity' ? '--env NUMBA_CACHE_DIR="./tmp/numba_cache"' : '' }"

    input:
        tuple val(meta), path(kmer_freqs)
    output:
        tuple val(meta), path("*output.tsv"), emit: umap_tsv

    script:
    template "umap_transform.py"
}

process hdbscanCluster {
    tag "${meta.scenario.count}/${meta.scenario.rep}"
    container "docker.io/hecrp/nanoclust-read_clustering:latest"

    input:
        tuple val(meta), path(umap_tsv)
    output:
        tuple val(meta), path("*.tsv"), emit: clusters

    script:
    def min_cluster_size_prop=0.005
    """
    #!/usr/bin/env python

    import pandas as pd
    import hdbscan

    umap_out = pd.read_csv("$umap_tsv", delimiter="\t")
    X = umap_out.loc[:, ["D1", "D2"]]

    nseqs = X.shape[0]
    df = pd.DataFrame({})
    min_size = int(max(2, nseqs * $min_cluster_size_prop))
    clusters = hdbscan.HDBSCAN(min_cluster_size=min_size, cluster_selection_epsilon=0.5, core_dist_n_jobs=${task.cpus}).fit_predict(X)
    umap_out['cluster_id'] = clusters

    umap_out.loc[:, ["read", "cluster_id"]].to_csv(f"hdbscan.output.tsv", sep="\t", index=False)
    """
}

process createOtuTable {
    tag "${meta.scenario.count}/${meta.scenario.rep}"
    container "docker.io/hecrp/nanoclust-read_clustering:latest"

    input:
        tuple val(meta), path(clusters_tsv)
    output:
        tuple val(meta), path("*.tsv"), emit: otu_table

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    clusters = pd.read_csv("$clusters_tsv", delimiter="\t")

    # extract barcode to new column
    clusters['barcode'] = clusters['read'].apply(lambda read: read.split(';')[1].split('=')[1])

    # group and count number of reads from each barcode in each cluster
    grouped_counts = clusters.groupby(['cluster_id', 'barcode']).count()

    # place barcode as columns
    otu_table = grouped_counts.unstack('barcode', fill_value=0)
    # drop unused 'read' column level
    otu_table.columns = otu_table.columns.droplevel()

    otu_table.to_csv("otu_table.tsv", sep="\t")
    """
}

process splitReadsByCluster {
    tag "${meta.scenario.count}/${meta.scenario.rep}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/seqkit:2.6.1--h9ee0642_0':
            'biocontainers/seqkit:2.6.1--h9ee0642_0' }"

    input:
        tuple val(meta), path(clusters_tsv), path(reads)
    output:
        tuple val(meta), path("*.fastq.gz"), emit: reads_by_cluster_fastq_gz
        tuple val(meta), path("*.txt"), emit: reads_ids_by_cluster

    script:
    // assuming read id ends with ';'
    // will add ';cluster={cluster_id};' details to reads
    """
    NCLUSTERS=\$(awk 'NR > 1 {print \$2}' $clusters_tsv | sort -nr | uniq | head -n1)

    for ((i = -1 ; i <= \$NCLUSTERS ; i++));
    do
        cluster_id=\$i
        awk -v cluster="\$cluster_id" '(\$2 == cluster) {print \$1}' $clusters_tsv > "cluster_\${cluster_id}_ids.txt"
        seqkit grep \\
            -f "cluster_\${cluster_id}_ids.txt" \\
            $reads \\
            | seqkit replace -p '(.*);' -r "\\\$1;cluster=\${cluster_id};" \\
                -o "cluster_\${cluster_id}_reads.fastq.gz"
    done
    """
}

process findMostAbundantSeqsInCluster {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
        tuple val(meta), path(clusters)

    output:
        tuple val(meta), path("most_abundant.fastq.gz"), emit: most_abundant

    script:
    """
    for cluster in $clusters;
    do
        vsearch \\
            --fastx_uniques "\$cluster" \\
            --fastqout - \\
            --sizeout \\
            --topn 1 \\
            --threads $task.cpus \\
            2>> vsearch.log \\
            >> most_abundant.fastq
    done

    gzip most_abundant.fastq
    """
}

process gatherMinClusterSizeStats {
    container "docker.io/hecrp/nanoclust-read_clustering:latest"
    containerOptions "${ workflow.containerEngine == 'singularity' ? '--env MPLCONFIGDIR="./tmp/mplconfig"' : '' }"

    input:
        tuple val(meta), path(umap_tsv)
    output:
        tuple val(meta), path("*.tsv"), emit: stats

    script:
    """
    #!/usr/bin/env python

    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import hdbscan

    def gather_cluster_stats(umap, props):
        nseqs = umap.shape[0]
        df = pd.DataFrame({})
        for prop in props:
            min_size = int(max(2, nseqs * prop))
            a = hdbscan.HDBSCAN(min_cluster_size=min_size, cluster_selection_epsilon=0.5, core_dist_n_jobs=${task.cpus}).fit_predict(umap)
            counts = np.unique(a, return_counts=True)
            cluster_result = pd.DataFrame({
                'cluster': counts[0],
                'size': counts[1]
            })
            unassigned = sum(cluster_result[cluster_result['cluster'] == -1]['size'])
            num_otus: int = sum(cluster_result['cluster'] != -1)
            df = pd.concat([
                df,
                pd.DataFrame({
                    'thresh': [prop],
                    'min_cluster_size': min_size,
                    'loss': unassigned / nseqs,
                    'unassigned': unassigned,
                    'numOTUs': num_otus
                })
            ])
        return df

    umap_out = pd.read_csv("$umap_tsv", delimiter="\t")
    X = umap_out.loc[:, ["D1", "D2"]]

    props = np.linspace(0, 0.03, 201)
    df = gather_cluster_stats(X, props)
    df.to_csv('min_cluster_size_stats.tsv', sep='\t', index=False)
    """
}