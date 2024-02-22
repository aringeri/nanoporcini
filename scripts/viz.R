library(phyloseq)
library(ggplot2)

phylo <- readRDS("output/localtest/phyloseq/its1/all_reads_full_its_phyloseq.rds")
phylo_ns <- filter_taxa(phylo, function (x) {sum(x > 0) > 1}, prune=TRUE)

phylo_ra <- transform_sample_counts(phylo_ns,function(x) x / sum(x))
ps_genus <- tax_glom(phylo_ra, taxrank = "Genus", NArm = FALSE)
# ps_genusP <- subset_taxa(ps_genus, Genus %in% c("D_5__Flavisolibacter", "D_5__Halomonas",
                                                # "D_5__Thiobacillus", "D_5__Sphingomonas",
                                                # "D_5__Bacillus", "D_5__uncultured Acidobacterium sp.", 
                                                # "D_5__Bradyrhizobium", "D_5__Ohtaekwangia", 
                                                # "D_5__Steroidobacter")


agaricomycetes <- subset_taxa(ps_genus, Class %in% c("c_Agaricomycetes"))


plot_heatmap(phylo)
plot_heatmap(phylo_ns)
plot_bar(phylo_ra, fill='Genus', facet_grid='Phylum')
