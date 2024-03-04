include { ImportFastaIntoQiime } from "./qiime/import" 

process UnGzip {
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        tuple val(meta), path(gzipped)
    
    output:
        tuple val(meta), path("*"), emit: unzipped
    
    script:

    if (!"$gzipped".endsWith(".gz")) {
        error "Expecting input to be a gzipped file ending in .gz: $gzipped"
    }

    """
    gunzip -c $gzipped > ${gzipped.baseName} 
    """
}

process ImportFastaQiime {

    container "qiime2/amplicon:2023.9"

    input:
        tuple val(meta), path(fasta)
    output:
        tuple val(meta), path("*.qza"), emit: sequences
    
    script:
    
    """
    export NUMBA_CACHE_DIR="./numbacache"
    export MPLCONFIGDIR="./mplconfigdir"

    qiime tools import --type 'FeatureData[Sequence]' \\
        --input-path $fasta \\
        --output-path ${fasta.baseName}.qza
    """
}

process ClassifyConsensusBlast {
    cpus 4
    container "qiime2/amplicon:2023.9"

    input:
        tuple val(meta), path(toClassify) // .qza
        path(blastDB)
        path(referenceTax)
    output:
        tuple val(meta), path("classifications.qza"), emit: classifications
        tuple val(meta), path("search-results.qza"), emit: searchResults
    
    script:
    
    """
    export NUMBA_CACHE_DIR="./numbacache"
    export MPLCONFIGDIR="./mplconfigdir"

    qiime feature-classifier classify-consensus-blast \\
        --i-query $toClassify \\
        --i-blastdb $blastDB \\
        --i-reference-taxonomy $referenceTax \\
        --o-classification classifications.qza \\
        --o-search-results search-results.qza \\
        --verbose \\
        --p-num-threads $task.cpus
    """
}

process ExportQiimeData {
    container "qiime2/amplicon:2023.9"

    input:
        tuple val(meta), path(qiimeData) // .qza
    output:
        tuple val(meta), path("output/*"), emit: exported

    script:

    if (!"$qiimeData".endsWith(".qza")) {
        error "Expecting input to be a qiime artifact ending in .qza: $qiimeData"
    }

    """
    export NUMBA_CACHE_DIR="./numbacache"
    export MPLCONFIGDIR="./mplconfigdir"
    mkdir output/

    qiime tools export \\
        --input-path $qiimeData  \\
        --output-path output/
    """
}

workflow ClassifyTaxonomyBlast {
    take:
        reads // reads to classify in fasta format
        qiime_db // BLASTDB -- qza file
        qiime_taxonomy // FeatureData[Taxonomy] -- qza file
    
    main:
        qiime = UnGzip(reads) | ImportFastaQiime

        classifyResults = ClassifyConsensusBlast(
            qiime.sequences,
            qiime_db,
            qiime_taxonomy
        )

        export = ExportQiimeData(classifyResults.classifications)
        // ExportQiimeData(classifyResults.searchResults)
    
    emit:
        classifications = export.exported
}
