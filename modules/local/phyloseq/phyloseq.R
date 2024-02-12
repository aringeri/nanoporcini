
install.packages("speedytax", repos = "https://cran.ms.unimelb.edu.au/")
library(speedytax)
library(phyloseq)
library(ggplot2)


otu <- read.table("data/otu.tsv", header = TRUE, sep = "\t", comment.char = "", row.names=1)
otu_mat <- data.matrix(otu)
#rownames(otu_mat) <- as.vector(otu[,1])

p_otu <- otu_table(otu_mat, taxa_are_rows=TRUE)
p_tax <- import_sintax_tax_table("/data/tax-no-empty.tsv")
taxa_names(p_tax) <- gsub(";.*", "", taxa_names(p_tax))

physeq = phyloseq(p_otu, p_tax)

ggsave("plot.png", device="png", path="/data", plot=plot_bar(physeq, fill="Domain"))

saveRDS(physeq, "data/phyloseq.rds")