#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.outdir = 'output'

process qualitySingle {
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