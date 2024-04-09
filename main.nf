#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.input = "$baseDir/data/sub100/*.fastq.gz"
// params.unite_db = "$baseDir/data/db/utax_reference_dataset_all_25.07.2023.fasta.gz"
params.unite_db = "$baseDir/data/db/UNITE-full-all-10.15156-BIO-2938070-20230725/sh_general_release_dynamic_s_all_25.07.2023.fasta"
params.sintax_db = "$baseDir/data/db/unite-sintax.fasta"
params.rdp_lsu_db = "$baseDir/data/db/RDP-LSU/RDP-LSU-training-11-dada2/RDP_LSU_fixed_train_set_v2.fa"
params.rdp_lsu_trained_model_dir = "$baseDir/data/db/RDP-LSU/RDPClassifier_fungiLSU_trainsetNo11_trained/"
params.outdir = 'output'
params.classifier = 'blast'

include { NANOPLOT as NANOPLOT_2 } from './modules/nf-core/nanoplot'
include { 
  NANOPLOT_SINGLE;
  NANOPLOT_SINGLE as NANOPLOT_SINGLE_2;
  NANOPLOT_SINGLE as NANOPLOT_SINGLE_3 
  } from './modules/local/nanoplot_single'
include { PORECHOP_PORECHOP } from './modules/nf-core/porechop/porechop'
include { FILTLONG; FILTLONG as FILTLONG_2 } from './modules/nf-core/filtlong'

include { VSEARCH_FILTER_MAX_EE } from './modules/local/vsearch/filter_ee'
include { VSEARCH_FILTER_MAX_EE as VSEARCH_FILTER_MAX_EE_2 } from './modules/local/vsearch/filter_ee'

include { ITSXPRESS } from './modules/local/itsxpress'
include { ITSX } from './modules/local/itsx'

include { CHOPPER } from './modules/local/chopper'

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
include { VSEARCH_SINTAX } from './modules/nf-core/vsearch/sintax'


include { PrepUniteDBForQiime } from "./workflows/qiime_prep_db"
include { LoadTaxTableIntoPhyloseq } from './workflows/phyloseq/import/qiime'
include { ImportSintaxTaxonomyIntoPhyloseq } from './workflows/phyloseq/import/sintax'
include { ClassifyTaxonomyBlast } from "./workflows/classify"
include { ClassifyRDP } from "./workflows/ClassifyRDP"

include { PHYLOSEQ; CreatePhyloseqObject } from './modules/local/phyloseq'

include { CLUSTER } from './workflows/vsearch_cluster'
include { ExtractRegions } from './workflows/extract_regions'
include { 
  qualityControl; 
  qualityControl as qualityControl_postPrimer;
  qualityControl as qualityControl_itsx;
  qualityControl as qualityControl_q_filter;
  qualityControl as qualityControl_chimera
  } from './workflows/qualityControl'

def collectWithId(id, ch) {
  ch.collect { meta, read -> 
    read
  }.map {
    [ [id: id], it ]
  }
}

def collectByRegionWithId(String id, ch) {
  ch.map { meta, reads -> [ meta.subMap('region'), reads ] }
    .groupTuple()
    .map { meta, reads -> [ meta + [id: id], reads ] }
}

def removeFileEndings(file, String extension, String... rest) {
  for (ext in [extension, *rest]) {
    if (file.endsWith(ext)) {
      return file.replaceAll("$ext\$", "")
    }
  }
  return file
}

workflow {
    ch_raw_reads = Channel.fromPath(params.input, checkIfExists: true)

    ch_reads = ch_raw_reads
      .map{ fastq -> 
        [ [id: removeFileEndings(fastq.name, ".fastq.gz", ".fq.gz")], fastq]
      }
    
    qualityControl('raw_reads_all_samples', ch_reads)

    oriented = CUTADAPT_REORIENT_READS(ch_reads)

    qualityControl_postPrimer('post_primer_trimming', oriented.reads)

    extracted = ExtractRegions(oriented.reads)

    regions = extracted.its1
      .mix(extracted.its2)
      .mix(extracted.full_its)
      .mix(extracted.lsu)

    qualityControl_itsx('post_its_extraction', regions)
    
    filtered_regions = CHOPPER(
      regions.map { meta, reads -> 
        def args = params.qualityFiltering[meta.region]
        [meta, args, reads] 
      }
    )
    qualityControl_q_filter('post_quality_filtering', filtered_regions.reads)

    derep = RENAME_BARCODE_LABEL(filtered_regions.reads).reads
      |  VSEARCH_DEREPLICATE //per sample

    uchime_denovo = VSEARCH_UCHIME_DENOVO(derep.reads)
    fa_for_uchime = SEQKIT_FQ2FA(derep.reads).fasta

    uchime_fa_with_db = fa_for_uchime.map { meta, fa -> 
      db = (meta.region == "LSU") ? params.rdp_lsu_db : params.unite_db
      tuple(meta, fa, db)
    }

    uchime_ref = VSEARCH_UCHIME_REF(uchime_fa_with_db)

    nonchimeras = SEQKIT_REMOVE_CHIMERAS(
      derep.reads.join(uchime_denovo.chimeras).join(uchime_ref.chimeras)
    )
    qualityControl_chimera('post_chimera_filtering', nonchimeras.nonchimeras)

    pooled_reads = FASTQ_CONCAT(collectByRegionWithId("pooled_reads", nonchimeras.nonchimeras)).merged_reads
    NANOPLOT_SINGLE_2(pooled_reads)

    pooled_reads_derep = VSEARCH_DEREPLICATE_2(pooled_reads).reads
    NANOPLOT_SINGLE_3(pooled_reads_derep.map { meta, reads -> 
      [meta + [id: 'post_derep_pooled_reads'], reads]
    })

    cluster_out = VSEARCH_CLUSTER_A(pooled_reads_derep)
    
    pooled_reads_fa = SEQKIT_FQ2FA_2(pooled_reads).fasta
    
    merged_channels = pooled_reads_fa
      .map { meta, reads -> [ meta.subMap(['region', 'id']), reads ] }
      .join( 
        cluster_out.centroids.map{ meta, reads -> [ meta.subMap(['region', 'id']), reads ] }
      )
    otus = VSEARCH_MAP_READS_TO_OTUS(merged_channels)

    centroids = cluster_out.centroids.branch {
      lsu: it[0].region == "LSU"
      its: true
    }

    tax_rds_lsu = ClassifyRDP(centroids.lsu, params.region.LSU.rdp_trained_model_dir).tax_rds

    if (params.classifier == "blast") {
      uniteDB = PrepUniteDBForQiime(params.unite_db)

      cResults = ClassifyTaxonomyBlast(
        centroids.its,
        uniteDB.blastDB,
        uniteDB.taxonomy
      )

      tax_rds_its = LoadTaxTableIntoPhyloseq(cResults.classifications)
    } else if (params.classifier == "sintax") {
      tax_rds_its = VSEARCH_SINTAX(centroids.its, params.sintax_db)
        .tsv
        | ImportSintaxTaxonomyIntoPhyloseq

    } else {
      error "input param 'classifier' not defined or not supported: ${params.classifier}"
    }

    tax_rds = tax_rds_its.mix(tax_rds_lsu)

    CreatePhyloseqObject(
      otus.otu_tab.join(tax_rds).map { meta, otu, tax ->  [meta, tax, otu] } // re-order parameters
    )
}