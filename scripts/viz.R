library(phyloseq)
library(ggplot2)
library(cowplot)

setwd('/Users/alex/repos/long-read-ITS-metabarcoding')

#phylo <- readRDS("output/full/phyloseq/all_reads_full_its.phyloseq.rds")
its1 <- tax_table(readRDS("output/full-multi-region/phyloseq/ITS1/pooled_reads/blast/pooled_reads.phyloseq.rds"))
# its1_tax_table <- import_tax_tsv('output/full-multi-region/qiime-export/ITS1/pooled_reads/output/taxonomy.tsv')
# its1_tax_table %>%
#   dplyr::filter(!is.na(Species)) %>%
#   dplyr::mutate(Species = paste0('s_', Species)) %>%
#   merge(its1, by=0, all=TRUE) %>%
#   dplyr::filter(Species.x != Species.y) %>%
#   tibble::column_to_rownames('Row.names') %>%
#   View()


its2 <- readRDS("output/full-multi-region/phyloseq/ITS2/pooled_reads/blast/pooled_reads.phyloseq.rds")
full_its <- readRDS("output/full-multi-region/phyloseq/FULL_ITS/pooled_reads/blast/pooled_reads.phyloseq.rds")
#lsu <- readRDS("output/full-multi-region/phyloseq/LSU/pooled_reads/blast/pooled_reads.phyloseq.rds")
lsu <- readRDS("output/full-multi-region-dada2/phyloseq/LSU/pooled_reads/blast/pooled_reads.phyloseq.rds")
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
lsu_fungi <- subset_taxa(lsu, Kingdom %in% c("Fungi"))
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

glom_and_count <- function(phylo, rank="Species") {
  glom <- tax_glom(phylo, taxrank = rank, NArm = FALSE)
  glom.sums <- data.frame(taxa_sum = taxa_sums(glom))
  taxa  <- merge(tax_table(glom), glom.sums, by=0, all=TRUE)
  # ordered_taxa <- taxa[with(taxa, order(taxa_sum, decreasing = TRUE)), ]
  total <- sum(taxa$taxa_sum)
  taxa$prop <- taxa$taxa_sum / total
  return(taxa)
}

simple_view <- function(df, title) {
  View(df[1:100, c('Row.names', 'Genus', 'Species', 'taxa_sum', 'prop')], title = title)
}

top_species_its1 <- glom_and_count(its1, "Species")
top_species_its2 <- glom_and_count(its2, "Species")
top_species_full_its <- glom_and_count(full_its, "Species")

simple_view(top_species_full_its, "Full ITS")
simple_view(top_species_its2, "ITS2")

merged <- merge(x = top_species_its2, y = top_species_its1, by=c("Family", "Genus", "Species"), all = TRUE)
# d <- data.frame(Species = merged$Species)
d <- merged[, c('Species', 'prop.x', 'prop.y')]
d[is.na(d$prop.x), "prop.x"] <- 0
d[is.na(d$prop.y), "prop.y"] <- 0
t.test(d$prop.x, d$prop.y, paired = TRUE)
# d$prop_diff <- diff$prop.y - diff$prop.x


top_genus_its1 <- glom_and_count(filter_at_least(its1_fungi, otu_size_gt), "Genus")
top_genus_its2 <- glom_and_count(filter_at_least(its2_fungi, otu_size_gt), "Genus")
top_genus_full_its <- glom_and_count(filter_at_least(full_its_fungi, otu_size_gt), "Genus")


d_m <- merge_w_props(top_genus_its1, top_genus_its2)

chisq.test(d_m$taxa_sum.y, p = d_m$prop.x)

merge_w_props <- function(x, y) {
  m <- merge(x = x, y = y, by=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), all = TRUE)
  d_m <- m[, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", 'prop.x', 'prop.y', 'taxa_sum.x', 'taxa_sum.y')]
  d_m[is.na(d_m$prop.x), "prop.x"] <- 0
  d_m[is.na(d_m$prop.y), "prop.y"] <- 0
  d_m[is.na(d_m$taxa_sum.x), "taxa_sum.x"] <- 0
  d_m[is.na(d_m$taxa_sum.y), "taxa_sum.y"] <- 0
  return(d_m)
}

its1.v.its2 <- merge_w_props(top_genus_its1, top_genus_its2)
its1.v.full <- merge_w_props(top_genus_its1, top_genus_full_its)
its2.v.full <- merge_w_props(top_genus_its2, top_genus_full_its)

plot(its1.v.its2$prop.x, its1.v.its2$prop.y)
abline(a=0, b=1)

plot(its1.v.full$prop.x, its1.v.full$prop.y)
abline(a=0, b=1)
plot(its2.v.full$prop.x, its2.v.full$prop.y)
abline(a=0, b=1)

t.test(its1.v.its2$prop.x, its1.v.its2$prop.y, paired = TRUE)
t.test(its1.v.full$prop.x, its1.v.full$prop.y, paired = TRUE)
t.test(its2.v.full$prop.x, its2.v.full$prop.y, paired = TRUE)


ts <- transform_sample_counts(its1, function(x) 100 * x / sum(x))
filter_taxa(ts, function (y) {sum(y) > 0.5 }, prune=TRUE)
transform_sample_counts(full_its, function(x) 100 * x / sum(x))
