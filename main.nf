#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.input = "$baseDir/data/sub100/*.fastq.gz"
params.unite_db = "$baseDir/data/db/utax_reference_dataset_all_25.07.2023.fasta.gz"
params.outdir = 'output'

include { NANOPLOT } from './modules/nf-core/nanoplot'
include { NANOPLOT as NANOPLOT_2 } from './modules/nf-core/nanoplot'
include { NANOPLOT_BULK } from './modules/local/nanoplot_bulk'
include { NANOPLOT_BULK as NANOPLOT_BULK_2 } from './modules/local/nanoplot_bulk'
include { PORECHOP_PORECHOP } from './modules/nf-core/porechop/porechop'
include { FILTLONG } from './modules/nf-core/filtlong'

include { ITSXPRESS } from './modules/local/itsxpress'
// include { VSEARCH_RELABEL } from './modules/local/vsearch/relabel'
include { RENAME_BARCODE_LABEL } from './modules/local/rename_barcode_label'
include { FASTQ_CONCAT } from './modules/local/fastq_concat'

include { VSEARCH_DEREPLICATE } from './modules/local/vsearch/dereplicate'

include { VSEARCH_CLUSTER_A } from './modules/local/vsearch/cluster'
include { VSEARCH_UCHIME_DENOVO } from './modules/local/vsearch/uchime_denovo'
include { VSEARCH_UCHIME_REF } from './modules/local/vsearch/uchime_ref'

include { FORMAT_CONSENSUS_LABELS } from './modules/local/format_consensus_labels'
// include { VSEARCH_CLUSTER } from './modules/nf-core/vsearch/cluster'
include { VSEARCH_SINTAX } from './modules/nf-core/vsearch/sintax'

include { PHYLOSEQ } from './modules/local/phyloseq'

workflow {
    ch_raw_reads = Channel.fromPath(params.input, checkIfExists: true)

    NANOPLOT_BULK(
      ch_raw_reads.collect().map { [ [id: "raw"], it] }
    )

    ch_reads = ch_raw_reads
      .map{ fastq -> 
        [ [id: fastq.name.replaceAll(".fastq.gz\$", "")], fastq]
      }
    // must have ending '.fastq.gz'
    // NANOPLOT(ch_reads)

    chopped = PORECHOP_PORECHOP(ch_reads)

    filtered = FILTLONG(
      chopped.reads.map{ ch -> 
        (meta, reads) = ch
        shortreads = []
        [meta, shortreads, reads] 
      }
    )

    // NANOPLOT_2(filtered.reads)
    NANOPLOT_BULK_2(
      filtered.reads.collect{ tuple -> 
          (meta, reads) = tuple
          reads
        }.map{ reads -> 
          [ [id: "filtered"], reads ]
        }
    )

    its1 = ITSXPRESS(filtered.reads).reads
    
    all_reads = its1.collect {
      (meta, read) = it
      read
    }.map {
      [ [id: "all-reads", foo:"bar"], it]
    }
    
    prepped_for_vsearch = FASTQ_CONCAT(all_reads).merged_reads 
      | RENAME_BARCODE_LABEL

    cluster_out = VSEARCH_DEREPLICATE(prepped_for_vsearch.reads).reads
      | VSEARCH_CLUSTER_A
    
    VSEARCH_UCHIME_REF(
      VSEARCH_UCHIME_DENOVO(cluster_out.centroids).nonchimeras, 
      params.unite_db
    )
    
    consensus = FORMAT_CONSENSUS_LABELS(cluster_out.consensus).reads
    tax = VSEARCH_SINTAX(consensus, params.unite_db).tsv

    PHYLOSEQ (
      tax.join(cluster_out.otu).map {
        (meta, tax_tsv, otu_tsv) = it
        [meta["id"], tax_tsv, otu_tsv]
      }
    )
}