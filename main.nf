#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.input = "$baseDir/data/sub100/*.fastq.gz"
// params.unite_db = "$baseDir/data/db/utax_reference_dataset_all_25.07.2023.fasta.gz"
params.unite_db = "$baseDir/data/db/UNITE-full-all-10.15156-BIO-2938070-20230725/sh_general_release_dynamic_s_all_25.07.2023.fasta"
params.sintax_db = "$baseDir/data/db/unite-sintax.fasta"
params.outdir = 'output'

include { NANOPLOT } from './modules/nf-core/nanoplot'
include { NANOPLOT as NANOPLOT_2 } from './modules/nf-core/nanoplot'
include { NANOPLOT_BULK } from './modules/local/nanoplot_bulk'
include { NANOPLOT_BULK as NANOPLOT_BULK_2 } from './modules/local/nanoplot_bulk'
include { NANOPLOT_BULK as NANOPLOT_BULK_3 } from './modules/local/nanoplot_bulk'
include { NANOPLOT_BULK as NANOPLOT_BULK_4 } from './modules/local/nanoplot_bulk'
include { NANOPLOT_SINGLE } from './modules/local/nanoplot_single'
include { PORECHOP_PORECHOP } from './modules/nf-core/porechop/porechop'
include { FILTLONG } from './modules/nf-core/filtlong'

include { ITSXPRESS } from './modules/local/itsxpress'
include { ITSX } from './modules/local/itsx'

include { SEQKIT_FQ2FA } from './modules/local/seqkit/fq2fa'
include { SEQKIT_FQ2FA as SEQKIT_FQ2FA_2 } from './modules/local/seqkit/fq2fa'
include { SEQKIT_REMOVE_CHIMERAS } from './modules/local/seqkit/remove'
include { RENAME_BARCODE_LABEL } from './modules/local/rename_barcode_label'
include { FASTQ_CONCAT } from './modules/local/fastq_concat'

include { VSEARCH_DEREPLICATE } from './modules/local/vsearch/dereplicate'
include { VSEARCH_DEREPLICATE as VSEARCH_DEREPLICATE_2 } from './modules/local/vsearch/dereplicate'

include { VSEARCH_CLUSTER_A } from './modules/local/vsearch/cluster'
include { VSEARCH_UCHIME_DENOVO } from './modules/local/vsearch/uchime_denovo'
include { VSEARCH_UCHIME_REF } from './modules/local/vsearch/uchime_ref'
include { VSEARCH_MAP_READS_TO_OTUS } from './modules/local/vsearch/map_to_otus'

include { CUTADAPT_REORIENT_READS } from './modules/local/cutadapt/reorient_reads'
include { FORMAT_CONSENSUS_LABELS } from './modules/local/format_consensus_labels'
include { VSEARCH_SINTAX } from './modules/nf-core/vsearch/sintax'

include { Classify } from "./workflows/classify"
include { PrepUniteDBForQiime } from "./workflows/qiime_prep_db"
include { LoadTaxTableIntoPhyloseq } from './workflows/phyloseq/import/qiime'

include { PHYLOSEQ } from './modules/local/phyloseq'

include { CLUSTER } from './workflows/vsearch_cluster'

def collectWithId(id, ch) {
  ch.collect {
    (meta, read) = it
    read
  }.map {
    [ [id: id], it]
  }
}

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

    oriented = CUTADAPT_REORIENT_READS(ch_reads)
    NANOPLOT_BULK_3(
      collectWithId("cutadapt-oriented", oriented.reads)
    )

    //chopped = PORECHOP_PORECHOP(ch_reads)
    filtered = FILTLONG(
      oriented.reads.map{ ch -> 
        (meta, reads) = ch
        shortreads = []
        [meta, shortreads, reads] 
      }
    )

    // NANOPLOT_2(filtered.reads)
    NANOPLOT_BULK_2(
      collectWithId("filtered", filtered.reads)
    )

    // all_reads = collectWithId("all-reads-filtered", filtered.reads)
    derep = (ITSXPRESS(filtered.reads).reads 
      | RENAME_BARCODE_LABEL).reads
      | VSEARCH_DEREPLICATE //per sample
    
    // NANOPLOT_SINGLE(derep.reads)

    NANOPLOT_BULK_4(
      collectWithId("full-its-derep", derep.reads)
    )

    uchime_denovo = VSEARCH_UCHIME_DENOVO(derep.reads)
    uchime_ref = VSEARCH_UCHIME_REF(
      SEQKIT_FQ2FA(derep.reads).fasta, 
      params.unite_db)

    nonchimeras = SEQKIT_REMOVE_CHIMERAS(
      derep.reads.join(uchime_denovo.chimeras).join(uchime_ref.chimeras)
    )

    all_reads = FASTQ_CONCAT(collectWithId("all_reads_full_its", nonchimeras.nonchimeras)).merged_reads
    all_reads_derep = VSEARCH_DEREPLICATE_2(all_reads).reads
    
    cluster_out = VSEARCH_CLUSTER_A(all_reads_derep)
    
    all_reads_fa = SEQKIT_FQ2FA_2(all_reads).fasta
    otus = VSEARCH_MAP_READS_TO_OTUS(all_reads_fa, cluster_out.centroids)
    
    uniteDB = PrepUniteDBForQiime(params.unite_db)

    cResults = Classify(
      cluster_out.centroids,
      uniteDB.blastDB,
      uniteDB.taxonomy
    )

    LoadTaxTableIntoPhyloseq(cResults.classifications.map{it[1]})

    tax = VSEARCH_SINTAX(cluster_out.centroids, params.sintax_db).tsv

    PHYLOSEQ (
      tax.join(otus.otu_tab).map {
        (meta, tax_tsv, otu_tsv) = it
        [meta["id"], tax_tsv, otu_tsv]
      }
    )

    // its = ITSXPRESS(derep.reads).reads
    // VSEARCH_DEREPLICATE_2(
    //   its.map { ch -> 
    //     (meta, reads) = ch
    //     meta2 = meta.clone()
    //     meta2.id = "full-its"
    //     [meta2, reads]
    //   }
    // )
    
    /*
    itsx_its1 = ITSX(SEQKIT_FQ2FA(filtered.reads).fasta).its1
    
    all_reads = its1.collect {
      (meta, read) = it
      read
    }.map {
      [ [id: "all-reads", foo:"bar"], it]
    }

    combined = collectWithId("all_reads_itsx_its1", itsx_its1)
    cluster_out = CLUSTER(combined)
      
    // prepped_for_vsearch = FASTQ_CONCAT(all_reads).merged_reads 
    //   | RENAME_BARCODE_LABEL

    // cluster_out = VSEARCH_DEREPLICATE(prepped_for_vsearch.reads).reads
    //   | VSEARCH_CLUSTER_A
    
    // VSEARCH_UCHIME_REF(
    //   VSEARCH_UCHIME_DENOVO(cluster_out.centroids).nonchimeras, 
    //   params.unite_db
    // )
    
    consensus = FORMAT_CONSENSUS_LABELS(cluster_out.consensus).reads
    tax = VSEARCH_SINTAX(consensus, params.unite_db).tsv

    PHYLOSEQ (
      tax.join(cluster_out.otu).map {
        (meta, tax_tsv, otu_tsv) = it
        [meta["id"], tax_tsv, otu_tsv]
      }
    )*/
}