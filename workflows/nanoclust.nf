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
            | gatherMinClusterSizeStats
}

process kmerFreq {
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
    container "docker.io/hecrp/nanoclust-read_clustering:latest"
    containerOptions "${ workflow.containerEngine == 'singularity' ? '--env NUMBA_CACHE_DIR="./tmp/numba_cache"' : '' }"

    input:
        tuple val(meta), path(kmer_freqs)
    output:
        tuple val(meta), path("*output.tsv"), emit: umap_tsv

    script:
    template "umap_transform.py"
}

process gatherMinClusterSizeStats {
    container "docker.io/hecrp/nanoclust-read_clustering:latest"

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
            a = hdbscan.HDBSCAN(min_cluster_size=min_size, cluster_selection_epsilon=0.5).fit_predict(umap)
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