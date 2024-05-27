library(phyloseq)
library(ggplot2)
library(magrittr)
library(dplyr)
library(cowplot)

populate_sample_data <- function(phylo, region) {
  s_names <- sample_names(phylo)
  sample_data(
    data.frame(
      SampleName = s_names,
      Region = rep(region, length(s_names)),
      row.names = s_names)
  )
}

phyloseq_dir <- "output/full-multi-region/phyloseq"
its1 <- readRDS("output/full-multi-region/phyloseq/ITS1/pooled_reads/blast/pooled_reads.phyloseq.rds")
sample_data(its1) <- populate_sample_data(its1, "ITS1")

its2 <- readRDS("output/full-multi-region/phyloseq/ITS2/pooled_reads/blast/pooled_reads.phyloseq.rds")
sample_data(its2) <- populate_sample_data(its1, "ITS2")

full_its <- readRDS("output/full-multi-region/phyloseq/FULL_ITS/pooled_reads/blast/pooled_reads.phyloseq.rds")
sample_data(full_its) <- populate_sample_data(full_its, "FULL_ITS")


merge_phylo_region <- function(regions, region_names) {
  merged <- NULL
  for (i in seq_along(regions)) {
    cur <- regions[[i]]
    taxa_names(cur) <- paste0(taxa_names(cur), ".", region_names[i])
    sample_names(cur) <- paste0(sample_names(cur), ".", region_names[i])

    if (is.null(merged)) {
      merged <- cur
    } else {
      merged <- merge_phyloseq(merged, cur)
    }
  }
  return(merged)
}
renamed.full_its <- full_its
taxa_names(renamed.full_its) <- paste0(taxa_names(full_its), ".FULL_ITS")
sample_names(renamed.full_its) <- paste0(sample_names(full_its), ".FULL_ITS")

renamed.its1 <- its1
taxa_names(renamed.its1) <- paste0(taxa_names(its1), ".ITS1")
sample_names(renamed.its1) <- paste0(sample_names(its1), ".ITS1")

merged.its1_and_full_its <- merge_phyloseq(renamed.its1, renamed.full_its)

regions <- list(its1, full_its)

top_n_taxa <- function(phylo, n=20) {
  phylo %>%
    taxa_sums() %>%
    sort(decreasing = TRUE) %>%
    names() %>% .[1:n] %>%
    prune_taxa(phylo)
}

taxa_sums_with_names <- function(phylo) {
  sums <- taxa_sums(phylo) %>% sort(decreasing = TRUE)
  tax_table(phylo)[names(sums), ] %>%
    as.data.frame() %>% as_tibble(rownames = 'OTU') %>%
    mutate(taxa_sum = sums)
}

# top_20_full_its <- top_n_taxa(full_its, 20)
# plot_bar(top_20_full_its, fill='Genus')
# taxa_sums_with_names(top_20_full_its)[, c('Genus', 'taxa_sum')]

taxa_sum <- list()
plots <- list()
for (region in regions) {
  top_20 <- top_n_taxa(region, 20)
  plots[[length(plots) + 1 ]] <- plot_bar(top_20, fill='Genus')
  taxa_sum[[length(taxa_sum) + 1 ]] <- taxa_sums_with_names(top_20)
}

plot_richness(full_its, measures = c("Shannon", "InvSimpson", "Fisher"))
plot_richness(its1, measures = c("Shannon", "InvSimpson", "Fisher"))

subsample.full_its <- prune_samples(
  sapply(sample_names(full_its),
         function(sample) !(sample %in% c('barcode96', 'barcode95')))
  , full_its)
ord.full_its <- ordinate(subsample.full_its, method="NMDS", distance="bray")
plot_ordination(full_its, ord.full_its, type='sample', label='SampleName')

subsample.its1 <- prune_samples(
  sapply(sample_names(its1),
         function(sample) !(sample %in% c('barcode96', 'barcode95')))
  , its1)

ord.its1 <- ordinate(subsample.its1, method="NMDS", distance="bray")
plot_ordination(subsample.its1, ord.its1, type='sample', label='SampleName')

subsample.merged <- merged.its1_and_full_its
# prune_samples(
#   sample_names(merged.its1_and_full_its)[grep('barcode9[5,6]', sample_names(merged.its1_and_full_its), invert = TRUE)],
#   merged.its1_and_full_its
# )
sample_data(subsample.merged)
ord.merged <- ordinate(subsample.merged, method = "NMDS", distance = "bray")
plot_ordination(subsample.merged, ord.merged, type='sample', label='Region', color = "SampleName")


merged.its1_and_its2 <- merge_phylo_region(c(its1, its2, full_its), c("ITS1", "ITS2", "FULL_ITS"))

glom_genus <- merged.its1_and_its2 %>%
  subset_taxa(Domain %in% c("k_Fungi")) %>%
  tax_glom(taxrank = "Genus") %>%
  filter_taxa(function (x) {sum(x) > 1}, prune=TRUE) %>%
  prune_samples(
    sample_names(.)[grep('barcode9[5,6]', sample_names(.), invert = TRUE)],
    .
  )
ord.merged <- ordinate(glom_genus, method = "NMDS", distance = "bray")
plot_ordination(glom_genus, ord.merged, type='samples', label='Region', color = "SampleName")
plot_ordination(glom_genus, ord.merged, type='taxa', color = "Phylum") + facet_wrap(~Phylum, 3)