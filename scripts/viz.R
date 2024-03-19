library(phyloseq)
library(ggplot2)

setwd('/Users/alex/repos/long-read-ITS-metabarcoding')

#phylo <- readRDS("output/full/phyloseq/all_reads_full_its.phyloseq.rds")
phylo <- readRDS("output/full-itsx-trim-server/phyloseq/all_reads_full_its.phyloseq.rds")
phylo_ns <- filter_taxa(phylo, function (x) {sum(x) > 10}, prune=TRUE)
ntaxa(phylo_ns)
sample_sums(phylo_ns)
taxa_sums(phylo_ns)
sam_data(phylo_ns)
#phylo_ns <- filter_taxa(phylo, function (x) {sum(x > 0) > 1}, prune=TRUE)

phylo_ra <- transform_sample_counts(phylo_ns,function(x) x / sum(x))
ps_genus <- tax_glom(phylo_ra, taxrank = "Genus", NArm = FALSE)
# ps_genusP <- subset_taxa(ps_genus, Genus %in% c("D_5__Flavisolibacter", "D_5__Halomonas",
                                                # "D_5__Thiobacillus", "D_5__Sphingomonas",
                                                # "D_5__Bacillus", "D_5__uncultured Acidobacterium sp.", 
                                                # "D_5__Bradyrhizobium", "D_5__Ohtaekwangia", 
                                                # "D_5__Steroidobacter")


fungi <- subset_taxa(ps_genus, Domain %in% c("k_Fungi"))
agaricomycetes <- subset_taxa(ps_genus, Class %in% c("c_Agaricomycetes"))


plot_heatmap(phylo)
plot_heatmap(phylo_ns)
plot_bar(phylo_ra, fill='Genus', facet_grid='Phylum')

plot_bar(ps_genus, fill='Phylum')
plot_bar(fungi, fill='Phylum')
plot_heatmap(fungi)
plot_heatmap(agaricomycetes)
plot_bar(agaricomycetes, fill="Family", facet_grid="Order")
