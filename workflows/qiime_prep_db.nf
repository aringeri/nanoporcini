include { ImportFastaIntoQiime; ImportTaxonomyFeatureDataIntoQiime } from "./qiime/import" 

process ConvertHeadersToQiimeFmt {
    container "docker.io/biocontainers/biocontainers:v1.2.0_cv1"

    input:
        path(db_fasta)
    
    output:
        path("*.qiime.fasta")
    
    script:
    if (!"$db_fasta".endsWith(".fasta")) {
        error "Expecting input to be a '.fasta' file: $db_fasta"
    }

    """
    awk -vRS=">" -vORS="\\n" -vFS="\\n" -vOFS="\\t" \\
        ' NR>1 {split(\$1, arr, "|"); print(">"arr[3]"_"arr[2]"_"arr[4]"\\n"\$2)} ' \\
        $db_fasta \\
        > ${db_fasta.baseName}.qiime.fasta
    """
}

process MakeBlastDB {
    container "qiime2/amplicon:2023.9"

    input:
        path(db_seqs_qza)
    
    output:
        path("*.blast.qza")

    script:
    """
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime feature-classifier makeblastdb \\
        --i-sequences $db_seqs_qza \\
        --o-database ${db_seqs_qza.baseName}.blast.qza
    """
}

process ExtractTaxonomyAsTSV {
    container "docker.io/biocontainers/biocontainers:v1.2.0_cv1"

    input:
        path(db_fasta)
    
    output:
        path("*.tsv")
    
    script:
    if (!"$db_fasta".endsWith(".fasta")) {
        error "Expecting input to be a '.fasta' file: $db_fasta"
    }

    /* Extract taxonomy from UNITE fasta headers.
    Ex)
        >Glomeraceae|AM076560|SH146432.05FU|refs|k__Fungi;p__Glomeromycota;c__Glomeromycetes;o__Glomerales;f__Glomeraceae;g__;s__uncultured_Glomus
    would become (tsv):
        SH146432.05FU_AM076560_refs  k__Fungi;p__Glomeromycota;c__Glomeromycetes;o__Glomerales;f__Glomeraceae;g__;s__uncultured_Glomus
    */
    """
    awk -vRS=">" -vORS="\\n" -vFS="\\n" -vOFS="\\t" \\
        ' NR>1 {split(\$1, arr, "|"); print(arr[3]"_"arr[2]"_"arr[4]"\\t"arr[5])} ' \\
        $db_fasta \\
        > ${db_fasta.baseName}.tsv
    """
}

workflow PrepUniteDBForQiime {
    take:
        dbFasta

    main:

    blastDB_qza = ConvertHeadersToQiimeFmt(dbFasta)
        | ImportFastaIntoQiime
        | MakeBlastDB

    taxonomy = ExtractTaxonomyAsTSV(dbFasta)
        | ImportTaxonomyFeatureDataIntoQiime

    emit:
        blastDB = blastDB_qza
        taxonomy = taxonomy

}
