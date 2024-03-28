
process UnGzip {
    tag "$meta.id - $meta.region"
    
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

process ConvertRDPResultToTaxTable {
    tag "$meta.id - $meta.region"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
        tuple val(meta), path(rdp_result_tsv)
    
    output:
        tuple val(meta), path("*.tax.tsv"), emit: tax_table_tsv
    
    script:

    if (!"$rdp_result_tsv".endsWith(".tsv")) {
        error "Expecting input to be a RDP result file ending in .tsv: $rdp_result_tsv"
    }

    """
    awk -vRS="\n" -vFS="\t" \\
        'BEGIN {print "id\\tDomain\\tPhylum\\tClass\\tOrder\\tFamily\\tGenus"} {print \$1"\\t"\$6"\\t"\$9"\\t"\$12"\\t"\$15"\\t"\$18"\\t"\$21}' \\
        $rdp_result_tsv \\
        > ${rdp_result_tsv.baseName}.tax.tsv
    """
}


process ClassifyReadsRDP {
    tag "$meta.id - $meta.region"
    container "biocontainers/rdp_classifier:2.14--hdfd78af_0"

    input:
        tuple val(meta), path(reads)
        path(model_dir) // RDP properties file
    output:
        tuple val(meta), path("*.tsv"), emit: classifications

    script:

    if (!"$reads".endsWith(".fasta")) {
        error "Expecting input to be a fasta file: $reads"
    }

    """
    rdp_classifier \\
        -t $model_dir/rRNAClassifier.properties \\
        -o ${reads.baseName}.tsv \\
        $reads
    """
}

process ImportTaxonomyIntoPhyloseq {
    tag "$meta.id - $meta.region"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
         'https://depot.galaxyproject.org/singularity/bioconductor-phyloseq:1.44.0--r43hdfd78af_0' :
         'biocontainers/bioconductor-phyloseq:1.44.0--r43hdfd78af_0' }"

    input:
        tuple val(meta), path(tax_tsv)

    output:
        tuple val(meta), path("*.rds"), emit: tax_rds

    script:
    if (!"$tax_tsv".endsWith(".tsv")) {
        error "Expecting input to be a '.tsv' file: $tax_tsv"
    }

    def prefix  = tax_tsv.baseName
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(phyloseq))

    tax_df  <- read.table("$tax_tsv", sep="\\t", header=TRUE, comment.char = "", row.names=1)
    tax_mat <- as.matrix(tax_df)

    TAX     <- tax_table(tax_mat)
    saveRDS(TAX, file = paste0("$prefix", ".rds"))
    """
}

workflow ClassifyRDP {
    take:
        reads // reads to classify in fasta format
        trained_rdp_model_dir
    
    main:
        plain_reads = UnGzip(reads)

        tax_rds = ClassifyReadsRDP(plain_reads, trained_rdp_model_dir)
            | ConvertRDPResultToTaxTable
            | ImportTaxonomyIntoPhyloseq

        // export = ExportQiimeData(classifyResults.classifications)
        // ExportQiimeData(classifyResults.searchResults)
    
    emit:
         tax_rds = tax_rds
}
