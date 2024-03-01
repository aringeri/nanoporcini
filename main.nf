#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.input = "$baseDir/data/sub100/*.fastq.gz"
// params.unite_db = "$baseDir/data/db/utax_reference_dataset_all_25.07.2023.fasta.gz"
params.unite_db = "$baseDir/data/db/UNITE-full-all-10.15156-BIO-2938070-20230725/sh_general_release_dynamic_s_all_25.07.2023.fasta"
params.sintax_db = "$baseDir/data/db/unite-sintax.fasta"
params.outdir = 'output'
params.classifier = 'blast'

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


include { PrepUniteDBForQiime } from "./workflows/qiime_prep_db"
include { LoadTaxTableIntoPhyloseq } from './workflows/phyloseq/import/qiime'
include { ImportSintaxTaxonomyIntoPhyloseq } from './workflows/phyloseq/import/sintax'
include { ClassifyTaxonomySintax } from './workflows/classify_taxonomy_sintax'
include { ClassifyTaxonomyBlast } from "./workflows/classify"

include { PHYLOSEQ; CreatePhyloseqObject } from './modules/local/phyloseq'

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

    if (params.classifier == "blast") {
      uniteDB = PrepUniteDBForQiime(params.unite_db)

      cResults = ClassifyTaxonomyBlast(
        cluster_out.centroids,
        uniteDB.blastDB,
        uniteDB.taxonomy
      )

      tax_rds = LoadTaxTableIntoPhyloseq(cResults.classifications.map{it[1]})
    } else if (params.classifier == "sintax") {
      tax_rds = VSEARCH_SINTAX(cluster_out.centroids, params.sintax_db)
        .tsv
        .map { 
          (meta, tsv) = it
          tsv
        }
        | ImportSintaxTaxonomyIntoPhyloseq

    } else {
      error "input param 'classifier' not defined or not supported: ${params.classifier}"
    }

    CreatePhyloseqObject(
      // careful with merge when using channels with multiple values.
      otus.otu_tab.merge(tax_rds).map {
        (meta, otu, tax) = it
        // re-order parameters
        [meta, tax, otu]
      }
    )
}