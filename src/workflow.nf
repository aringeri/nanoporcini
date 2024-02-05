#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.input = "$baseDir/../data/fq/raw/*.fq.gz"
params.outdir = 'output'

include { nanoplot } from './quality' addParams(outdir: params.outdir)
include { nanoplot as nanoplot2 } from './quality' addParams(outdir: params.outdir)
include { nanoplot_summary } from './quality' addParams(outdir: params.outdir)

process nanofilt {
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

process porechop {
    input:
    tuple val(sampleId), path('reads.fq.gz')

    output:
    tuple val(sampleId), path("*.fq.gz")

    """
    porechop-runner.py --threads 16 -i reads.fq.gz | gzip > trimmed-reads.fq.gz
    """
}

process itsx {
  input:
    tuple val(sampleId), path('reads.fq.gz')

  output:
    tuple val(sampleId), path('ITSx_out.*.fasta')
  
  """
  gunzip -c reads.fq.gz > reads.fq
  seqtk seq -A reads.fq > reads.fa
  ITSx -i reads.fa --cpu 8 --save_regions all -t Fungi
  # TODO salvage quality scores
  """
}

workflow {
  ch_input = Channel.fromPath( params.input )
    .map{ rawReadPath -> 
          sampleId = rawReadPath.name.replaceAll(".fq.gz\$", "")
          tuple(sampleId, rawReadPath) 
    }
    .first()
    .view()

  nanoplot(ch_input, Channel.value('raw'))

  without_adapters = porechop(ch_input)

  filteredQ12 = nanofilt(without_adapters, Channel.value(12))
  nanoplot2(filteredQ12, Channel.value('q12'))
}

process mk_subsample {
  input:
    tuple val(sampleId), path('in.fq.gz')

  output:
    tuple val(sampleId), path("*.fq.gz")

  """
  gunzip -c in.fq.gz > in.fq
  seqtk sample -s11 in.fq 100 \
    | gzip > out.fq.gz
  """
}

process nanoplot_conda {
  // conda 'bioconda::nanoplot=1.42.0'
  conda '/opt/miniconda3/miniconda3'

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

process nanoplot_list {
  conda '/opt/miniconda3/miniconda3'
  publishDir "$params.outdir/NanoPlot/$outdir/", mode: 'copy', overwrite: false

  input:
    path('reads*.fq.gz')
    val outdir

  output:
    path "*"

  """
  NanoPlot --fastq_rich reads*.fq.gz -o . --verbose
  """
}

workflow subsample {
  ch_input = Channel.fromPath(params.input)
    .map{ rawReadPath -> 
          sampleId = rawReadPath.name.replaceAll(".fq.gz\$", "")
          tuple(sampleId, rawReadPath) 
    }
    .view()
  
  subs = mk_subsample(ch_input)

  onlyFiles = subs
    .map { tuple -> tuple[1] }
    .collect()
  
  nanoplot_list(onlyFiles, Channel.value('sub100'))
  
  without_adapters = porechop(subs)
  filteredQ15 = nanofilt(without_adapters, Channel.value(15))

  q15Files = filteredQ15
    .map{ tuple -> tuple[1] }
    .collect()

  // Not working -- nanofilt not found - conda?
  nanoplot_summary(q15Files, Channel.value('sub100-Q15'))
}