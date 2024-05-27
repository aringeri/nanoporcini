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

include {
    SEQKIT_FQ2FA
    SEQKIT_FQ2FA as SEQKIT_FQ2FA_2
    SEQKIT_FQ2FA as SEQKIT_FQ2FA_3
} from './modules/local/seqkit/fq2fa'
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
    qualityControl as qualityControl_postPrimer
    qualityControl as qualityControl_itsx
    qualityControl as qualityControl_q_filter
    qualityControl as qualityControl_chimera
} from './workflows/qualityControl'

include {
    FindReadsWithID
    FindReadsWithID as FindReadsWithID_ITS2
    FindReadsWithID as FindReadsWithID_LSU
} from './modules/local/seqkit/FindReadsWithID'

include {
    seqkit
} from './modules/local/seqkit/seqkit'

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

    oriented = CUTADAPT_REORIENT_READS(ch_reads)

    extracted = ExtractRegions(oriented.reads)

    regions = extracted.its1
        .mix(extracted.its2)
        .mix(extracted.full_its)
        .mix(extracted.lsu)

    filtered_reads = CHOPPER(
        regions.map { meta, reads ->
            def args = params.qualityFiltering[meta.region]
            [meta, args, reads]
        }
    ).reads

    filtered_reads_full_its = filtered_reads.filter { meta, reads -> meta.region == "FULL_ITS" }

//    filtered_reads_full_its.view()

    pooled_full_its = FASTQ_CONCAT(collectByRegionWithId("pooled_full_its", filtered_reads_full_its)).merged_reads
        | RENAME_BARCODE_LABEL
    pooled_full_its_derep = VSEARCH_DEREPLICATE(pooled_full_its)

//    pooled_full_its.reads.view()

    uchime_denovo = VSEARCH_UCHIME_DENOVO(pooled_full_its_derep.reads)
//    uchime_denovo.chimeras.view()
    uchime_ref = SEQKIT_FQ2FA(pooled_full_its_derep.reads).fasta
        .map { meta, fasta -> tuple(meta, fasta, params.unite_db) }
        | VSEARCH_UCHIME_REF

//    uchime_ref.chimeras.view()

    full_its_nonchimeras = SEQKIT_REMOVE_CHIMERAS(
        pooled_full_its_derep.reads.join(uchime_denovo.chimeras).join(uchime_ref.chimeras)
    )
//    full_its_nonchimeras.view()

    full_its_centroids = seqkit(
            "seq --only-id --id-regexp '([^\\s,;]+);'",
            VSEARCH_CLUSTER_A(full_its_nonchimeras).centroids
    )
    full_its_otus = VSEARCH_MAP_READS_TO_OTUS(
        SEQKIT_FQ2FA_2(pooled_full_its).join(full_its_centroids)
    )

    region_matches = (
        FindReadsWithID(
            full_its_centroids.first(), // convert to value channel so it can be re-used for each region
            collectByRegionWithId('full_its_centroid_matches', filtered_reads).filter { meta, _ -> meta.region != "FULL_ITS" }
        ) | SEQKIT_FQ2FA_3
    ).branch { meta, reads ->
        lsu: meta.region == "LSU"
        its: true
    }
//    region_matches.its.view()
//    region_matches.lsu.view()

//    its1_centroid_reads = FindReadsWithID(full_its_centroids, filtered_reads.its1.collect{meta, reads -> reads})
//    its2_centroid_reads = FindReadsWithID_ITS2(full_its_centroids, filtered_reads.its2.collect{meta, reads -> reads})
//    lsu_centroid_reads = FindReadsWithID_LSU(full_its_centroids, filtered_reads.lsu.collect{meta, reads -> reads})
//    its1_centroid_reads.view()

    uniteDB = PrepUniteDBForQiime(params.unite_db)

    tax_rds_its = ClassifyTaxonomyBlast(
        full_its_centroids.mix(region_matches.its),
        uniteDB.blastDB,
        uniteDB.taxonomy
    ).classifications | LoadTaxTableIntoPhyloseq

    tax_rds_lsu = ClassifyRDP(region_matches.lsu, params.region.LSU.rdp_trained_model_dir).tax_rds

    tax_rds = tax_rds_its.mix(tax_rds_lsu)
//    tax_rds.view()
//    full_its_otus.otu_tab.view()

    otu_and_tax = tax_rds.combine(full_its_otus.otu_tab.map{ meta, otu -> otu })
            .map { meta, tax, otu ->  [meta, tax, otu] } // re-order parameters

    otu_and_tax.view()
    CreatePhyloseqObject(
        otu_and_tax
    )

    /*

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
    )*/
}