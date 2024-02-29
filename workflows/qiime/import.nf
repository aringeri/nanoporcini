process ImportFastaIntoQiime {

    container "qiime2/amplicon:2023.9"

    input:
        path(fasta)
    output:
        path("*.qza"), emit: sequences
    
    script:
    
    """
    qiime tools import --type 'FeatureData[Sequence]' \\
        --input-path $fasta \\
        --output-path ${fasta.baseName}.qza
    """
}

process ImportTaxonomyFeatureDataIntoQiime {

    container "qiime2/amplicon:2023.9"

    input:
        path(taxonomy_tsv) // two columns (ID, Taxonomy), no header
    output:
        path("*.qza")
    
    script:
    
    """
    qiime tools import --type 'FeatureData[Taxonomy]' \\
        --input-format HeaderlessTSVTaxonomyFormat \\
        --input-path $taxonomy_tsv \\
        --output-path ${taxonomy_tsv.baseName}.qza
    """
}
