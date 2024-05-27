library(dada2)
library(phyloseq)
library(ggplot2)

setwd("/Users/alex/repos/long-read-ITS-metabarcoding")

qual <- plotQualityProfile(c("data/sub100/20221214.CT.Sample1.BC89.fq.gz", "data/sub1000/20221214.CT.Sample1.BC89.fastq.gz"))

ggsave('scratch/quality%03d.png', qual)

seqs <- getSequences("data/sub100/20221214.CT.Sample1.BC89.fq.gz")

taxa = assignTaxonomy(
  seqs,
  "data/db/RDP-LSU/rdp_train.LSU.dada2.fasta.gz",
  multithread = 2
)

rows <- names(rownames(taxa))
ids <- lapply(ids, FUN = function(x) unlist(strsplit(x, "\\s"))[1])
rownames(taxa) = ids

write.table(taxa, file="scratch/taxa.tsv", sep="\t", row.names=TRUE)

tax_t <- tax_table(taxa)
