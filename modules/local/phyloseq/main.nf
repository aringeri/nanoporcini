process PHYLOSEQ {
    tag "$prefix"
    label 'process_low'

    conda "bioconda::bioconductor-phyloseq=1.44.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/bioconductor-phyloseq:1.44.0--r43hdfd78af_0' :
    //     'biocontainers/bioconductor-phyloseq:1.44.0--r43hdfd78af_0' }"
    container 'ghcr.io/aringeri/phyloseq-speedytax'

    input:
    tuple val(prefix), path(tax_tsv), path(otu_tsv)

    output:
    tuple val(prefix), path("*phyloseq.rds"), emit: rds
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def otu_tsv = "\"${otu_tsv}\""
    def tax_tsv = "\"${tax_tsv}\""
    def prefix  = "\"${prefix}\""
    """
    #!/usr/bin/env Rscript

    #install.packages("speedytax", repos = "https://cran.ms.unimelb.edu.au/")
    library(speedytax)
    suppressPackageStartupMessages(library(phyloseq))

    otu_df  <- read.table($otu_tsv, sep="\\t", header=TRUE, comment.char = "", row.names=1)
    otu_mat <- as.matrix(otu_df)

    OTU     <- otu_table(otu_mat, taxa_are_rows=TRUE)
    TAX  <- import_sintax_tax_table($tax_tsv)
    phy_obj <- phyloseq(OTU, TAX)

    saveRDS(phy_obj, file = paste0($prefix, "_phyloseq.rds"))

    # Version information
    writeLines(c("\\"${task.process}\\":",
        paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
        paste0("    phyloseq: ", packageVersion("phyloseq"))),
        "versions.yml"
    )
    """
}

process CreatePhyloseqObject {
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
         'https://depot.galaxyproject.org/singularity/bioconductor-phyloseq:1.44.0--r43hdfd78af_0' :
         'biocontainers/bioconductor-phyloseq:1.44.0--r43hdfd78af_0' }"

    input:
        tuple val(meta), path(tax_rds), path(otu_tsv)

    output:
        tuple val(meta), path("*.phyloseq.rds"), emit: rds

    script:
    if (!"$tax_rds".endsWith(".rds")) {
        error "Expecting taxonomy table to be an '.rds' file: $tax_rds"
    }
    if (!"$otu_tsv".endsWith(".tsv")) {
        error "Expecting otu table to be an '.tsv' file: $otu_tsv"
    }
    
    
    prefix  = meta.id
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(phyloseq))

    otu_df  <- read.table("$otu_tsv", sep="\\t", header=TRUE, comment.char = "", row.names=1)
    otu_mat <- as.matrix(otu_df)

    OTU     <- otu_table(otu_mat, taxa_are_rows=TRUE)
    TAX     <- readRDS("$tax_rds")
    phy_obj <- phyloseq(OTU, TAX)

    saveRDS(phy_obj, file = paste0("$prefix", ".phyloseq.rds"))
    """
}