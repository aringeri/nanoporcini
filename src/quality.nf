#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process nanoplot {
    publishDir "$params.outdir/NanoPlot/$outdir/$sampleId", mode: 'copy', overwrite: false
    
    input:
    tuple val(sampleId), path(reads, name: 'reads.fq.gz')
    val outdir
    
    output:
    path "*"

    """
    NanoPlot --fastq_rich reads.fq.gz -o . --verbose
    """
}

// process qualityMulti {
//     publishDir "$params.outdir/NanoPlot/$outdir", mode: 'copy', overwrite: false
    
//     input:
//     val seqs
//     val outdir
    
//     output:
//     path "*"

//     """
//     NanoPlot --fastq_rich $seqs -o . --verbose
//     """
// }