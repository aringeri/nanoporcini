#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.input = "$baseDir/../data/sub100/*.fastq.gz"
params.unite_db = "$baseDir/../data/db/utax_reference_dataset_all_25.07.2023.fasta.gz"
params.outdir = 'output'

include { NANOPLOT } from '../modules/nf-core/nanoplot/main'
include { NANOPLOT as NANOPLOT_2 } from '../modules/nf-core/nanoplot/main'
include { NANOPLOT_BULK } from '../modules/local/nanoplot_bulk/main'
include { NANOPLOT_BULK as NANOPLOT_BULK_2 } from '../modules/local/nanoplot_bulk/main'
include { PORECHOP_PORECHOP } from '../modules/nf-core/porechop/porechop/main'
include { FILTLONG } from '../modules/nf-core/filtlong'

include { ITSXPRESS } from '../modules/local/itsxpress/main'
include { VSEARCH_CLUSTER } from '../modules/nf-core/vsearch/cluster'
include { VSEARCH_SINTAX } from '../modules/nf-core/vsearch/sintax'

workflow {
    ch_raw_reads = Channel.fromPath(params.input)

    NANOPLOT_BULK(
      ch_raw_reads.collect().map { [ [id: "raw"], it] }
    )

    ch_reads = ch_raw_reads
      .map{ fastq -> 
        [ [id: fastq.name.replaceAll(".fastq.gz\$", "")], fastq]
      }
    // must have ending '.fastq.gz'
    NANOPLOT(ch_reads)

    chopped = PORECHOP_PORECHOP(ch_reads)

    filtered = FILTLONG(
      chopped.reads.map{ ch -> 
        (meta, reads) = ch
        shortreads = []
        [meta, shortreads, reads] 
      }
    )

    NANOPLOT_2(filtered.reads)
    NANOPLOT_BULK_2(
      filtered.reads.collect{ tuple -> 
          (meta, reads) = tuple
          reads
        }.map{ reads -> 
          [ [id: "filtered"], reads ]
        }
    )

    its1 = ITSXPRESS(filtered.reads).reads
    centroids = VSEARCH_CLUSTER(its1).centroids

    VSEARCH_SINTAX(centroids, params.unite_db)
    

}