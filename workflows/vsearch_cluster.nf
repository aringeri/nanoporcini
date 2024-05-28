include { FASTQ_CONCAT } from '../modules/local/fastq_concat'
include { RENAME_BARCODE_LABEL } from '../modules/local/rename_barcode_label'
include { VSEARCH_DEREPLICATE } from '../modules/local/vsearch/dereplicate'
include { VSEARCH_CLUSTER } from '../modules/local/vsearch/cluster'

workflow CLUSTER {
    take:
        all_reads
    
    main:
        prepped_for_vsearch = FASTQ_CONCAT(all_reads).merged_reads 
            | RENAME_BARCODE_LABEL

        cluster_out = VSEARCH_DEREPLICATE(prepped_for_vsearch.reads).reads
            | VSEARCH_CLUSTER
    
    emit:
        otu = cluster_out.otu
        centroids = cluster_out.centroids
        consensus = cluster_out.consensus
        
}