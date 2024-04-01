library(phyloseq)
library(ggplot2)
library(cowplot)

setwd('/Users/alex/repos/long-read-ITS-metabarcoding')

#phylo <- readRDS("output/full/phyloseq/all_reads_full_its.phyloseq.rds")
its1 <- readRDS("output/full-multi-region/phyloseq/ITS1/pooled_reads/blast/pooled_reads.phyloseq.rds")
its2 <- readRDS("output/full-multi-region/phyloseq/ITS2/pooled_reads/blast/pooled_reads.phyloseq.rds")
full_its <- readRDS("output/full-multi-region/phyloseq/FULL_ITS/pooled_reads/blast/pooled_reads.phyloseq.rds")
lsu <- readRDS("output/full-multi-region/phyloseq/LSU/pooled_reads/blast/pooled_reads.phyloseq.rds")
#lsu <- readRDS("output/sub100-lsu/)
full_its_old <- readRDS("output/full-itsx-trim-server/phyloseq/all_reads_full_its.phyloseq.rds")

filter_at_least <- function(phylo, n) {
  filter_taxa(phylo, function (x) {sum(x) > n}, prune=TRUE)
  #filter_taxa(phylo, function (x) {sum(x > n) > 0}, prune=TRUE)
}
otu_size_gt = 5

phylum_bar_plot <- function(region, level="Phylum", relative=TRUE) {
  filtered <- filter_at_least(region, otu_size_gt)
  print(filtered)
  filtered <- tax_glom(filtered, taxrank = level, NArm = FALSE)
  
  if (relative) {
    relative <- transform_sample_counts(filtered,function(x) x / sum(x))
  } else {
    relative <- filtered
  }

  plot_bar(relative, fill=level)
}

plot_grid(
  phylum_bar_plot(its1, relative=FALSE),
  phylum_bar_plot(its2, relative=FALSE),
  phylum_bar_plot(full_its, relative=FALSE),
  phylum_bar_plot(lsu, relative=FALSE),
  #phylum_bar_plot(full_its_old, relative=FALSE),
  labels = c("its1", "its2", "full its", "lsu"))

its1_fungi <- subset_taxa(its1, Domain %in% c("k_Fungi"))
its2_fungi <- subset_taxa(its2, Domain %in% c("k_Fungi"))
full_its_fungi <- subset_taxa(full_its, Domain %in% c("k_Fungi"))
lsu_fungi <- subset_taxa(lsu, Domain %in% c("Fungi"))
full_its_old_fungi <- subset_taxa(full_its_old, Domain %in% c("k_Fungi"))

plot_grid(
  phylum_bar_plot(its1_fungi),
  phylum_bar_plot(its2_fungi),
  phylum_bar_plot(full_its_fungi),
  phylum_bar_plot(lsu_fungi),
  labels = c("its1", "its2", "full its", "lsu"))

its1_agarico <- subset_taxa(its1, Class %in% c("c_Agaricomycetes"))
its2_agarico <- subset_taxa(its2, Class %in% c("c_Agaricomycetes"))
full_its_agarico <- subset_taxa(full_its, Class %in% c("c_Agaricomycetes"))
lsu_agarico <- subset_taxa(lsu, Class %in% c("Agaricomycetes"))

plot_grid(
  phylum_bar_plot(its1_agarico, "Order", FALSE),
  phylum_bar_plot(its2_agarico, "Order", FALSE),
  phylum_bar_plot(full_its_agarico, "Order", FALSE),
  phylum_bar_plot(lsu_agarico, "Order", FALSE),
  labels = c("its1", "its2", "full its", "lsu"))

#plot_grid(
#  phylum_bar_plot(full_its_agarico, "Genus", FALSE),
#  phylum_bar_plot(lsu_agarico, "Genus", FALSE),
#  labels = c("full its", "lsu"))
  
plot_heatmap(filter_at_least(lsu_agarico, otu_size_gt))
#warntaxa(phylo_ns)
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
