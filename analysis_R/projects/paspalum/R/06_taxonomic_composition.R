#####################################################################################
# Project: Paspalum Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 06_taxonomic_composition.R
# Author: Gustavo Goncalves
# Last edit: 19/11/2026
#
# Description:
# Taxonomic summary and relative abundance of bacterial communities at
# Phylum and Family levels, by sample and experimental group.
#####################################################################################

# Load project configuration
source("projects/paspalum/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")

# Load phyloseq object
ps <- readRDS(file.path(objects_path, "phyloseq_object.rds"))
theme_set(theme_bw())

# Create output directory
tax_path <- file.path(results_path, "taxonomic_composition")
dir.create(tax_path, recursive = TRUE, showWarnings = FALSE)

#####################################################################################
# Taxonomic summary (ASV richness, not read abundance)
#####################################################################################

taxa <- as.data.frame(tax_table(ps))
assigned <- colSums(!is.na(taxa) & taxa != "")
assigned.percent <- round(100 * assigned / nrow(taxa), 2)

top_phylum <- head(sort(table(taxa$Phylum), decreasing = TRUE), 15)
top_class  <- head(sort(table(taxa$Class),  decreasing = TRUE), 15)
top_family <- head(sort(table(taxa$Family), decreasing = TRUE), 20)
top_genus  <- head(sort(table(taxa$Genus),  decreasing = TRUE), 20)

capture.output({
  cat("=====================================\n")
  cat("TAXONOMIC COMPOSITION SUMMARY\n")
  cat("Project: Paspalum Metabarcoding 16S\n")
  cat("=====================================\n\n")
  cat("Number of samples:", nsamples(ps), "\n")
  cat("Total ASVs:", ntaxa(ps), "\n")
  cat("Total reads:", sum(sample_sums(ps)), "\n\n")

  cat("TAXONOMIC ASSIGNMENT\n\n")
  print(data.frame(Assigned_ASVs = assigned, Percentage = assigned.percent))
  cat("\nTOP 15 PHYLA BY ASV RICHNESS\n\n")
  print(top_phylum)
  cat("\nTOP 15 CLASSES BY ASV RICHNESS\n\n")
  print(top_class)
  cat("\nTOP 20 FAMILIES BY ASV RICHNESS\n\n")
  print(top_family)
  cat("\nTOP 20 GENERA BY ASV RICHNESS\n\n")
  print(top_genus)
}, file = file.path(tax_path, "taxonomic_summary.txt"))

#####################################################################################
# Prepare relative abundance data
#####################################################################################

ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
metadata <- data.frame(sample_data(ps))
metadata$SampleID <- rownames(metadata)

# Prepare one taxonomic rank at a time. Missing assignments remain visible as
# "Unclassified"; "Other" contains only less abundant classified taxa.
prepare_rank <- function(rank, n_top = 10) {
  ps.rank <- ps.rel
  tax.table <- tax_table(ps.rank)
  values <- as.character(tax.table[, rank])
  values[is.na(values) | values == ""] <- "Unclassified"
  tax.table[, rank] <- values
  tax_table(ps.rank) <- tax.table

  ps.rank <- tax_glom(ps.rank, taxrank = rank)
  otu <- as(otu_table(ps.rank), "matrix")
  if (taxa_are_rows(ps.rank)) otu <- t(otu)

  names.rank <- as.character(tax_table(ps.rank)[, rank])
  abundance <- tapply(colSums(otu), names.rank, sum)
  abundance <- sort(abundance, decreasing = TRUE)

  # Top 10 classified taxa, plus a distinct Unclassified category if present.
  top <- head(names(abundance)[names(abundance) != "Unclassified"], n_top)
  groups <- ifelse(names.rank %in% top, names.rank,
                   ifelse(names.rank == "Unclassified", "Unclassified", "Other"))

  grouped <- t(rowsum(t(otu), group = groups))
  df <- as.data.frame(grouped, check.names = FALSE)
  df$SampleID <- rownames(df)
  df <- df %>%
    pivot_longer(cols = -SampleID, names_to = rank, values_to = "Abundance") %>%
    left_join(metadata, by = "SampleID")

  levels.rank <- c(top, if ("Unclassified" %in% groups) "Unclassified",
                   if ("Other" %in% groups) "Other")
  df[[rank]] <- factor(df[[rank]], levels = levels.rank)

  # The grouped matrix already contains zeroes for absent taxa in each sample.

  mean.acc <- df %>%
    group_by_at(c("Species", "Accession", "Origin", rank)) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop")

  mean.species <- df %>%
    group_by_at(c("Species", "Origin", rank)) %>%
    summarise(MeanAbundance = mean(Abundance), .groups = "drop")

  global <- data.frame(
    Taxon = names(abundance),
    MeanAbundance = as.numeric(abundance) / nsamples(ps),
    row.names = NULL
  )
  global <- global[order(global$MeanAbundance, decreasing = TRUE), ]

  colors <- scales::hue_pal()(length(top))
  names(colors) <- top
  if ("Unclassified" %in% levels.rank) colors <- c(colors, Unclassified = "#888888")
  if ("Other" %in% levels.rank) colors <- c(colors, Other = "#D3D3D3")

  list(df = df, mean.acc = mean.acc, mean.species = mean.species,
       global = global, top = top, colors = colors)
}

phylum <- prepare_rank("Phylum")
family <- prepare_rank("Family")

#####################################################################################
# Save relative abundance tables (semicolon-separated, for Brazilian Excel)
#####################################################################################

write.csv2(phylum$mean.acc,
           file.path(tax_path, "mean_phylum_abundance_by_group.csv"), row.names = FALSE)
write.csv2(phylum$mean.species,
           file.path(tax_path, "mean_phylum_abundance_by_species.csv"), row.names = FALSE)
write.csv2(family$mean.acc,
           file.path(tax_path, "mean_family_abundance_by_group.csv"), row.names = FALSE)
write.csv2(family$mean.species,
           file.path(tax_path, "mean_family_abundance_by_species.csv"), row.names = FALSE)

#####################################################################################
# Relative abundance summary (percentages)
#####################################################################################

capture.output({
  cat("=====================================\n")
  cat("RELATIVE ABUNDANCE SUMMARY\n")
  cat("Project: Paspalum Metabarcoding 16S\n")
  cat("=====================================\n\n")
  cat("Number of samples:", nsamples(ps), "\n")
  cat("Values are mean relative abundance (%).\n")
  cat("'Other' groups lower-abundance classified taxa.\n")
  cat("'Unclassified' groups taxa lacking assignment at the displayed rank.\n\n")

  for (rank in c("Phylum", "Family")) {
    result <- if (rank == "Phylum") phylum else family
    cat("\n=====================================\n")
    cat(toupper(rank), "\n")
    cat("=====================================\n\n")

    cat("TOP 20", toupper(rank), "BY GLOBAL MEAN RELATIVE ABUNDANCE (%)\n\n")
    top.global <- head(result$global, 20)
    top.global$MeanAbundance <- round(100 * top.global$MeanAbundance, 2)
    print(top.global, row.names = FALSE)

    group.summary <- as.data.frame(result$mean.acc)
    group.summary$MeanAbundance <- round(100 * group.summary$MeanAbundance, 2)
    cat("\nMEAN RELATIVE ABUNDANCE BY ACCESSION AND ORIGIN (%)\n\n")
    for (acc in unique(group.summary$Accession)) {
      for (origin in unique(group.summary$Origin)) {
        part <- group.summary[
          group.summary$Accession == acc & group.summary$Origin == origin, ]
        if (nrow(part)) {
          cat("\nAccession:", as.character(acc),
              "| Origin:", as.character(origin), "\n\n")
          print(part[, c(rank, "MeanAbundance")], row.names = FALSE)
        }
      }
    }

    species.summary <- as.data.frame(result$mean.species)
    species.summary$MeanAbundance <- round(100 * species.summary$MeanAbundance, 2)
    cat("\nMEAN RELATIVE ABUNDANCE BY SPECIES AND ORIGIN (%)\n\n")
    for (species in unique(species.summary$Species)) {
      for (origin in unique(species.summary$Origin)) {
        part <- species.summary[
          species.summary$Species == species & species.summary$Origin == origin, ]
        if (nrow(part)) {
          cat("\nSpecies:", as.character(species),
              "| Origin:", as.character(origin), "\n\n")
          print(part[, c(rank, "MeanAbundance")], row.names = FALSE)
        }
      }
    }
  }
}, file = file.path(tax_path, "relative_abundance_summary.txt"))

#####################################################################################
# Plots: mean per accession/origin; mean per species/origin; individual samples
#####################################################################################

plot_rank <- function(result, rank, prefix) {
  rank.label <- rank
  fill.scale <- scale_fill_manual(values = result$colors,
                                  breaks = names(result$colors), drop = FALSE)
  bar.layer <- geom_col(position = position_stack(reverse = TRUE),
                        color = "grey50", linewidth = 0.2)
  y.scale <- scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"))

  p1 <- ggplot(result$mean.acc,
               aes_string(x = "Accession", y = "MeanAbundance", fill = rank)) +
    bar.layer + facet_wrap(~Origin) + fill.scale + y.scale +
    labs(x = "Accession", y = "Mean relative abundance", fill = rank.label)
  print(p1)
  ggsave(file.path(tax_path, paste0(prefix, "_mean_", tolower(rank),
                                   "_accession_origin.png")),
         plot = p1, width = 12, height = 6, dpi = 300)

  p2 <- ggplot(result$mean.species,
               aes_string(x = "Species", y = "MeanAbundance", fill = rank)) +
    bar.layer + facet_wrap(~Origin) + fill.scale + y.scale +
    labs(x = "Species", y = "Mean relative abundance", fill = rank.label)
  print(p2)
  ggsave(file.path(tax_path, paste0(sprintf("%02d", as.integer(prefix) + 1),
                                   "_mean_", tolower(rank), "_species_origin.png")),
         plot = p2, width = 12, height = 6, dpi = 300)

  df <- result$df
  sample.order <- unique(df$SampleID[
    order(df$Accession, df$Origin, df$SampleID)
  ])
  df$SampleID <- factor(df$SampleID, levels = sample.order)
  p3 <- ggplot(df, aes_string(x = "SampleID", y = "Abundance", fill = rank)) +
    bar.layer + fill.scale + y.scale +
    labs(x = "Sample", y = "Relative abundance", fill = rank.label) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  print(p3)
  ggsave(file.path(tax_path, paste0(sprintf("%02d", as.integer(prefix) + 2),
                                   "_", tolower(rank), "_all_samples.png")),
         plot = p3, width = 14, height = 7, dpi = 300)
}

plot_rank(phylum, "Phylum", "01")
plot_rank(family, "Family", "04")
