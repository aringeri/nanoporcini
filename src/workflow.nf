#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.input = "$baseDir/../data/fq/raw/*.fq.gz"
params.outdir = 'output'

include { nanoplot } from './quality' addParams(outdir: params.outdir)
include { nanoplot as nanoplot2 } from './quality' addParams(outdir: params.outdir)

process nanofilt {
    publishDir "$params.outdir/NanoFilt/q$qThresh/$sampleId", mode: 'copy', overwrite: false

    input:
    tuple val(sampleId), path(reads, name: '*.fq.gz')
    val qThresh

    output:
    // path "*.log" -- ignore log files for now (can be access in the work directory)
    tuple val(sampleId), path("*.fq.gz")

    """
    gunzip -c $reads \
        | NanoFilt -q $qThresh \
        | gzip > trimmed-"$sampleId".fq.gz
    """
}

process porechop {
    input:
    tuple val(sampleId), path('reads.fq.gz')

    output:
    tuple val(sampleId), path("*.fq.gz")

    """
    porechop-runner.py --threads 16 -i reads.fq.gz -o trimmed-reads.fq.gz
    """
}

workflow {
  ch_input = Channel.fromPath( params.input )
	.map{ rawReadPath -> 
        sampleId = rawReadPath.name.replaceAll(".fq.gz\$", "")
        tuple(sampleId, rawReadPath) 
    }
	.view()

  nanoplot(ch_input, Channel.value('raw'))

  without_adapters = porechop(ch_input)

  filteredQ12 = nanofilt(without_adapters, Channel.value(12))
  nanoplot2(filteredQ12, Channel.value('q12'))
}
