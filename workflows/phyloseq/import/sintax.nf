process ImportSintaxTaxonomyIntoPhyloseq {
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

    taxTable <- import_sintax_tax_table("$tax_tsv")
    saveRDS(taxTable, file = paste0("$prefix", ".rds"))
    """
}