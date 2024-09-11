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
params.trim_adapters = false
params.qc_quality_profile = true

include { dorado_trim_adapters } from './modules/local/dorado/trim'
include { CUTADAPT_REORIENT_READS } from './modules/local/cutadapt/reorient_reads'
include { CHOPPER } from './modules/local/chopper'
include { RENAME_BARCODE_LABEL } from './modules/local/rename_barcode_label'

include { ExtractRegions } from './workflows/extract_regions'
include {
    qualityControl;
    qualityControl as qualityControl_postPrimer
    qualityControl as qualityControl_itsx
    qualityControl as qualityControl_q_filter
    qualityControl as qualityControl_chimera
    qualityControl as qualityControl_adapter
} from './workflows/qualityControl'

include {
    SEQKIT_FQ2FA
    SEQKIT_FQ2FA as SEQKIT_FQ2FA_2
    SEQKIT_FQ2FA as SEQKIT_FQ2FA_3
} from './modules/local/seqkit/fq2fa'

include { SEQKIT_REMOVE_CHIMERAS } from './modules/local/seqkit/remove'

include {
    FASTQ_CONCAT
    FASTQ_CONCAT as FASTQ_CONCAT_2
    FASTQ_CONCAT as FASTQ_CONCAT_3
} from './modules/local/fastq_concat'

include { nanoclust } from './workflows/nanoclust'
include { nanoclust_consensus } from './workflows/nanoclust/consensus'
include { VSEARCH_DEREPLICATE } from './modules/local/vsearch/dereplicate'

include { VSEARCH_CLUSTER } from './modules/local/vsearch/cluster'
include { VSEARCH_UCHIME_DENOVO } from './modules/local/vsearch/uchime_denovo'
include { VSEARCH_UCHIME_REF } from './modules/local/vsearch/uchime_ref'
include { VSEARCH_MAP_READS_TO_OTUS } from './modules/local/vsearch/map_to_otus'

include { seqtk_sample } from './modules/local/seqtk/seqtk'
include { subsample } from './workflows/subsample'
include {
    assignTaxDnabarcoder
    assignTaxDnabarcoder as assignTaxDnabarcoder_vsearch
    assignTaxDnabarcoder as assignTaxDnabarcoder_vsearch_consensus
    assignTaxDnabarcoder as assignTaxDnabarcoder_consensus
} from './workflows/assign_tax_dnabarcoder'

include { PrepUniteDBForQiime } from "./workflows/qiime_prep_db"
include { LoadTaxTableIntoPhyloseq } from './workflows/phyloseq/import/qiime'
include { ClassifyTaxonomyBlast; UnGzip } from "./workflows/classify"
include { ClassifyRDP } from "./workflows/ClassifyRDP"

include { CreatePhyloseqObject; CreatePhyloseqOTUObject } from './modules/local/phyloseq'

include { FindReadsWithID } from './modules/local/seqkit/FindReadsWithID'

include {
    seqkit
    seqkit as seqkit_2
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
        .filter { meta, fastq -> meta.id !in params.exclude }

    qualityControl('01-raw_reads_all_samples', ch_reads)

    renamed = RENAME_BARCODE_LABEL(ch_reads)

    if (params.trim_adapters) {
        renamed = dorado_trim_adapters(renamed)
        qualityControl_adapter('XX-adapter-trimming', renamed)
    }

    oriented = CUTADAPT_REORIENT_READS(renamed)

    qualityControl_postPrimer('02-post_primer_trimming', oriented.reads)

    extracted = ExtractRegions(oriented.reads)

    regions = extracted.its1
        .mix(extracted.its2)
        .mix(extracted.full_its)
        .mix(extracted.lsu)

    qualityControl_itsx('03-post_its_extraction', regions)

    filtered_reads = CHOPPER(
        regions.map { meta, reads ->
            def args = params.qualityFiltering[meta.region]
            [meta, args, reads]
        }
    ).reads

    qualityControl_q_filter('04-post_quality_filtering', filtered_reads)

//     derep = VSEARCH_DEREPLICATE(filtered_reads)

    uchime_denovo = VSEARCH_UCHIME_DENOVO(filtered_reads)

    uchime_ref = SEQKIT_FQ2FA(filtered_reads).fasta
        .map { meta, fasta ->
            def db = (meta.region == 'LSU')
                    ? "$baseDir/data/db/RDP-LSU/rdp_train.LSU.dada2.fasta.gz"
                    : params.chimera_filtering.ref_db
            [meta, fasta, db]
        }
        | VSEARCH_UCHIME_REF

    nonchimeras = SEQKIT_REMOVE_CHIMERAS(
        filtered_reads.join(uchime_denovo.chimeras).join(uchime_ref.chimeras)
    )

    qualityControl_chimera('05-post_chimera_filtering', nonchimeras)

    subsampled = subsample(nonchimeras)

    if ('nanoclust' in params.cluster.methods) {
        pooled_nanoclust = FASTQ_CONCAT_3(
                subsampled.map { meta, reads -> [ meta.subMap('region', 'scenario'), reads ] }
                            .groupTuple()
                            .map { meta, reads -> [ meta + [id: "all_samples"], reads ] }
            ).merged_reads
        nanoclust_result = nanoclust(pooled_nanoclust)

        assignTaxDnabarcoder(
            nanoclust_result.most_abundant_by_cluster
                .map { meta, reads ->
                    [ meta + [cluster_method: "nanoclust_abundant"], reads ]
                }
        )

        assignTaxDnabarcoder_consensus(
            nanoclust_result.consensus_by_cluster
        )
    }

    if ('vsearch' in params.cluster.methods) {
        derep = VSEARCH_DEREPLICATE(subsampled).reads

        pooled = FASTQ_CONCAT(
            derep.map { meta, reads -> [ meta.subMap('region', 'scenario'), reads ] }
                        .groupTuple()
                        .map { meta, reads -> [ meta + [id: "all_samples"], reads ] }
        ).merged_reads

        cluster = VSEARCH_CLUSTER(pooled)

        otus = cluster.otu | UnGzip
        CreatePhyloseqOTUObject(otus)

        cluster.centroids.map { meta, reads ->
            [ meta + [cluster_method: 'vsearch'], reads ]
        } | assignTaxDnabarcoder_vsearch

        vsearch_consensus = cluster.clusters
            .flatMap { meta, clusters ->
                [  clusters.collect { cluster_fname ->
                      (full,id)=(cluster_fname.getName() =~ /cluster_(\-?\d+)/)[0]
                      meta + [ cluster : [ id: id ] ]
                   },
                    clusters
                ].transpose()
            }.map { meta, reads -> [ meta + [cluster_method: 'vsearch_consensus'], reads ] }
            | nanoclust_consensus

        vsearch_consensus | assignTaxDnabarcoder_vsearch_consensus
    }

    if (params.taxonomic_assignment.enabled) {
        centroids = seqkit(
            "seq --only-id --id-regexp '([^\\s,;]+);'",
            cluster.centroids
        )
        all_centroids = FASTQ_CONCAT_2(collectWithId('all_centroids', centroids)).merged_reads

        // Need to remove barcode here so ids match the centroid ids
        filtered_reads_no_barcode = seqkit_2(
                "seq --only-id --id-regexp '([^\\s,;]+);'",
                filtered_reads
        )

        region_matches = (
            FindReadsWithID(
                all_centroids,//.first(), // convert to value channel so it can be re-used for each region
                collectByRegionWithId('centroid_matches', filtered_reads_no_barcode)
            ) | SEQKIT_FQ2FA_3
        ).branch { meta, reads ->
            lsu: meta.region == "LSU"
            its: true
        }

        uniteDB = PrepUniteDBForQiime(params.unite_db)

        tax_rds_its = ClassifyTaxonomyBlast(
            region_matches.its,
            uniteDB.blastDB,
            uniteDB.taxonomy
        ).classifications | LoadTaxTableIntoPhyloseq

        tax_rds_lsu = ClassifyRDP(region_matches.lsu, params.region.LSU.rdp_trained_model_dir).tax_rds

        tax_rds = tax_rds_its.mix(tax_rds_lsu)

        tax_and_otu = tax_rds.map{ meta, tax -> [ meta.subMap('region'), tax ] }
            .join(
                otus.map{ meta, otu -> [ meta.subMap('region'), otu ] }
            )
            .map { meta, tax, otu -> [ meta + [id: "clustered-by-region"], tax, otu ] }

        CreatePhyloseqObject(tax_and_otu)
    }
}