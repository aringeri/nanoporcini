#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = "$baseDir/../data/fq/raw/*.fq.gz"
params.outdir = 'output'

include { qualitySingle } from './quality' addParams(outdir: params.outdir)
include { qualitySingle as qualitySingle2 } from './quality' addParams(outdir: params.outdir)

process filtering {
    publishDir "$params.outdir/NanoFilt/q$qThresh/$sampleId", mode: 'copy', overwrite: false

    input:
    tuple val(sampleId), path('reads.fq.gz')
    val qThresh

    output:
    // path "*.log" -- ignore log files for now (can be access in the work directory)
    tuple val(sampleId), path("*.fq.gz")

    """
    gunzip -c reads.fq.gz \
        | NanoFilt -q $qThresh \
        | gzip > trimmed-"$sampleId".fq.gz
    """
}

process test {

    input:
    // path(l, name:'trim.log')
    path(trimmed, name:'trimmed.fq.gz')

    output:
    stdout

    """
    echo $trimmed
    """
}

workflow {
  ch_input = Channel.fromPath( params.input )
	.map{ rawReadPath -> 
        sampleId = rawReadPath.name.replaceAll(".fq.gz\$", "")
        tuple(sampleId, rawReadPath) 
    }
	.view()

  qualitySingle(ch_input, Channel.value('raw'))

  filteredQ12 = filtering(ch_input, Channel.value(12))

  qualitySingle2(filteredQ12, Channel.value('q12'))

  //quality(params.input, 'raw')
  //test()

}
