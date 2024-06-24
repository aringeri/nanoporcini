
library(phyloseq)
library(ggplot2)
library(cowplot)
library(magrittr)
library(tidyr)

setwd("/Users/alex/repos/long-read-ITS-metabarcoding")

import_tax_tsv <- function(path) {
  tsv <- read.table(path, sep='\t', header=TRUE, row.names = 1)
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_df <- as.data.frame(
    tidyr::separate_wider_delim(tsv, delim=';', cols=Taxon, names=ranks, too_few = 'align_start')
  )
  rownames(tax_df) <- rownames(tsv)
  tax_df <- tax_df[, ranks]
  tax_df[, ranks] <- apply(tax_df[, ranks], MARGIN = 2, FUN = function(x) sub('[k,p,c,o,f,g,s]__', '', x))
  return(tax_df)
}

import_tax_tsv_phyloseq <- function(path) {
  tax_df <- import_tax_tsv(path)
  return(
    phyloseq::tax_table(as.matrix(tax_df))
  )
}

import_tax_tsv('output/test-centroid-based/qiime-export/FULL_ITS/centroid_matches/output/taxonomy.tsv')

its1_tax_table <- import_tax_tsv('output/test-centroid-based/qiime-export/ITS1/centroid_matches/output/taxonomy.tsv')
# its1 <- tax_table(readRDS("output/test-centroid-based/phyloseq/ITS1/full_its_centroid_matches/blast/full_its_centroid_matches.phyloseq.rds"))
# its1_tax_table %>% dplyr::mutate(Species = paste0('s_', Species)) %>%
#   merge(its1, by=0, all=TRUE) %>%
#   dplyr::filter(Species.x != Species.y) %>%
#   tibble::column_to_rownames('Row.names') %>%
#   View()
# rownames(its1_tax_table) <- paste0(rownames(its1_tax_table), ".ITS1")
# its1_tax_table$Region <- rep("ITS1", dim(its1_tax_table)[1])
its2_tax_table <- import_tax_tsv('output/test-centroid-based/qiime-export/ITS2/centroid_matches/output/taxonomy.tsv')
#rownames(its2_tax_table) <- paste0(rownames(its2_tax_table), ".ITS2")
# its2_tax_table$Region <- rep("ITS2", dim(its2_tax_table)[1])
full_its_tax_table <- import_tax_tsv('output/test-centroid-based/qiime-export/FULL_ITS/centroid_matches/output/taxonomy.tsv')
# full_its_tax_table$Region <- rep("FULL_ITS", dim(full_its_tax_table)[1])
# rownames(full_its_tax_table) <- paste0(rownames(full_its_tax_table), ".FULL_ITS")
lsu_tax_table <- read.table('output/test-centroid-based/rdp-dada2/LSU/centroid_matches/prefix.tsv', sep='\t', header = TRUE)
lsu_tax_table$Species <- NA
# lsu_tax_table$Region <- rep("LSU", dim(lsu_tax_table)[1])
#rownames(lsu_tax_table) <- paste0(rownames(lsu_tax_table), ".LSU")


named_its1_tax_table <- import_tax_tsv('output/test-centroid-based/qiime-export/ITS1/full_its_centroid_matches/output/taxonomy.tsv')
rownames(named_its1_tax_table) <- paste0(rownames(named_its1_tax_table), ".ITS1")
named_its2_tax_table <- import_tax_tsv('output/test-centroid-based/qiime-export/ITS2/full_its_centroid_matches/output/taxonomy.tsv')
rownames(named_its2_tax_table) <- paste0(rownames(named_its2_tax_table), ".ITS2")
named_full_its_tax_table <- import_tax_tsv('output/test-centroid-based/qiime-export/FULL_ITS/pooled_full_its/output/taxonomy.tsv')
rownames(named_full_its_tax_table) <- paste0(rownames(named_full_its_tax_table), ".FULL_ITS")
named_lsu_tax_table <- read.table('output/test-centroid-based/rdp-dada2/LSU/full_its_centroid_matches/prefix.tsv', sep='\t', header = TRUE)
named_lsu_tax_table$Species <- NA
rownames(named_lsu_tax_table) <- paste0(rownames(named_lsu_tax_table), ".LSU")

combined <- t(rbind(named_its1_tax_table, named_its2_tax_table, named_full_its_tax_table, named_lsu_tax_table))
View(combined[, sort(colnames(combined))])


otu_mat  <- as.matrix(
  read.table("output/test-centroid-based/vsearch-otu-map/FULL_ITS/pooled_full_its/pooled_full_its.otus.tsv", sep="\t", header=TRUE, comment.char = "", row.names=1)
)
otus <- otu_table(otu_mat, taxa_are_rows = TRUE)

glom_and_count <- function(phylo, rank="Genus") {
  glom <- tax_glom(phylo, taxrank = rank, NArm = FALSE)
  glom.sums <- data.frame(taxa_sum = taxa_sums(glom))
  taxa  <- merge(tax_table(glom), glom.sums, by=0, all=TRUE)
  # ordered_taxa <- taxa[with(taxa, order(taxa_sum, decreasing = TRUE)), ]
  total <- sum(taxa$taxa_sum)
  taxa$prop <- taxa$taxa_sum / total
  return(taxa)
}

counts <- list()
i <- 1
for (tax_table in list(full_its_tax_table, its1_tax_table, its2_tax_table, lsu_tax_table)) {
  phylo <- phyloseq::phyloseq(tax_table(as.matrix(tax_table)), otus)
  g <- glom_and_count(phylo)
  counts[[i]] <- g
  i <- i + 1
  # View(g)
}


ps <- list()
i <- 1
for (tax_table in list(its1_tax_table, its2_tax_table, full_its_tax_table, lsu_tax_table)) {
  phylo <- phyloseq::phyloseq(tax_table(as.matrix(tax_table)), otus)
  p <- plot_bar(phylo, fill="Phylum")
  ps[[i]] <- p
  i <- i+1
}

cowplot::plot_grid(plotlist = ps, ncol = 2, labels=c("ITS1", "ITS2", "FULL_ITS", "LSU"))


count_lowest_assignment_rank <- function(tax) {
  ranks <- colnames(tax)
  result <- data.frame(count = rep(0, length(ranks) + 1 ), row.names = c('Unassigned', ranks))
  for (rowname in rownames(tax)) {
    row <- tax[rowname,]
    idx <- which(!is.na(row) & row != 'Unassigned')
    rank <- if (length(idx) == 0) {
      'Unassigned'
    } else {
      ranks[max(idx)]
    }
    result[rank, 'count'] <- result[rank, 'count'] + 1
  }
  return(result)
}
count_lowest_assignment_rank(its1_tax_table)
count_lowest_assignment_rank(its2_tax_table)
count_lowest_assignment_rank(full_its_tax_table)
count_lowest_assignment_rank(lsu_tax_table)

# Find the lowest taxonomic rank where all regions are assigned
find_lowest_assigned_rank_in_OTU <- function(centroid_tax, additional_tax) {
  ranks <- colnames(centroid_tax)
  result <- data.frame(count = rep(0, length(ranks) + 1 ), row.names = c('Unassigned', ranks))

  for (rowname in rownames(centroid_tax)) {
    row <- centroid_tax[rowname,]
    lowest <- !is.na(row) & row != 'Unassigned'
    for (region in additional_tax) {
      if (rowname %in% rownames(region)) {
        # print(lowest)
        r <- region[rowname,]
        lowest <- lowest & !is.na(r) & row != 'Unassigned'
      }
    }
    # print(lowest)
    idx <- which(lowest)
    rank <- if (length(idx) == 0) {
      'Unassigned'
    } else {
      ranks[max(idx)]
    }
    # print(rank)
    result[rank, 'count'] <- result[rank, 'count'] + 1
  }
  return(result)
}
# subset by phylum, class
find_lowest_assigned_rank_in_OTU(full_its_tax_table, list(its1_tax_table, its2_tax_table, lsu_tax_table))
find_lowest_assigned_rank_in_OTU(full_its_tax_table, list(its1_tax_table, its2_tax_table))

find_proportion_of_consistent_assignments_per_OTU <- function(centroid_tax, additional_tax) {
  ranks <- colnames(centroid_tax)
  result <- data.frame(
    count = rep(0, length(ranks) + 1 ),
    prop_nas = rep(0, length(ranks) + 1 ),
    row.names = c(ranks, 'Unassigned'))

  # rowname <- rownames(centroid_tax)[2]
  rnames <- rownames(centroid_tax)
  for (rowname in rnames) {
    row <- centroid_tax[rowname,]
    otu_tax <- as.data.frame(t(sapply(additional_tax, function(df) df[rowname, ]))) %>% rbind(row)

    # print(otu_tax)
    for (rank in ranks) {
      rank_counts <- otu_tax %>% dplyr::count(!!sym(rank))
      max_matching <-
        rank_counts[!is.na(rank_counts[,rank]) & rank_counts[,rank] != 'Unassigned', ] %>%
        dplyr::arrange(desc(n))

      num_nas <- sum(
        rank_counts[is.na(rank_counts[,rank]) | rank_counts[,rank] == 'Unassigned', ]['n']
      )
      # print(max_matching)
      # print(length(max_matching$n) == 0)
      if (length(max_matching$n) == 0) {
        max_matching <- 0
      } else {
        max_matching <- max_matching[1,'n']
      }
      proportion <- max_matching / (length(additional_tax) + 1)
      result[rank, 'count'] <- result[rank, 'count'] + proportion

      prop_nas <- num_nas / (length(additional_tax) + 1)
      result[rank, 'prop_nas'] <- result[rank, 'prop_nas'] + prop_nas
    }
  }
  num_otus <- dim(centroid_tax)[1]
  return(result %>%
           dplyr::mutate(
             proportion = count / num_otus,
             proportion_nas = prop_nas / num_otus,
             dissent = 1 - proportion - proportion_nas
           )
  )
}
props <- find_proportion_of_consistent_assignments_per_OTU(full_its_tax_table, list(its1_tax_table, its2_tax_table, lsu_tax_table))
# props$dissent <- 1 - (props$proportion + props$proportion_nas)
props
# synonym list