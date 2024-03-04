process ImportQiime2TaxonomyIntoPhyloseq {
    tag "$prefix"

    container 'ghcr.io/aringeri/phyloseq-speedytax'

    input:
        path(tax_tsv)

    output:
         path("*.rds"), emit: tax_rds

    script:
    if (!"$tax_tsv".endsWith(".tsv")) {
        error "Expecting input to be a '.tsv' file: $tax_tsv"
    }

    prefix  = tax_tsv.baseName
    """
    #!/usr/bin/env Rscript
    library(speedytax)

    taxTable  <- import_qiime2_tax_table("$tax_tsv")
    saveRDS(taxTable, file = paste0("$prefix", ".rds"))
    """
}

process TweakTaxTableForSpeedytax {
    container "docker.io/biocontainers/biocontainers:v1.2.0_cv1"

    input:
        path(tax_tsv)

    output:
        path("output/*.tsv")
    
    script:
    if (!"$tax_tsv".endsWith(".tsv")) {
        error "Expecting input to be a '.tsv' file: $tax_tsv"
    }

    """
    mkdir output/

    # Rename 'Consensus' header to 'Confidence'
    # Replace separator in Taxon column from ';' to '; ' (has a space).
    sed '1s/Consensus/Confidence/' $tax_tsv \\
        | awk -vRS="\n" -vORS="\n" -vFS="\t" -vOFS="\t" 'NR == 1 {print \$0} NR > 1 {gsub(/;/, "; ", \$2); print}' \\
        > output/$tax_tsv
    """
}

workflow LoadTaxTableIntoPhyloseq {
    take:
        tax_table_tsv
    
    main:
    
    tax_rds = TweakTaxTableForSpeedytax(tax_table_tsv)
        | ImportQiime2TaxonomyIntoPhyloseq
    
    emit:
        tax_table_rds = tax_rds
}