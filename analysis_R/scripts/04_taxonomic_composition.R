#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 04_taxonomic_composition.R
# Autor: Gustavo Goncalves
# Last edit: 19/08/2026
#
# Description:
# Script to investigate relative abundance of taxas between samples from a phyloseq object.
#
#####################################################################################

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

# Set up initial parameters
path <- "results/objects"; list.files(path) # Set path
ps <- readRDS(file.path(path, "phyloseq_object.rds")) # Load phyloseq object
theme_set(theme_bw()) # Set ggplot theme

# --------------------------------------------------
# Taxonomic summary
# --------------------------------------------------

# Set up taxa object
taxa <- as.data.frame(tax_table(ps))

# ASVs assigned at each taxonomic level
assigned <- colSums(!is.na(taxa));assigned
assigned_percent <- round(100 * assigned / nrow(taxa),2);assigned_percent

# Most abundant taxa by number of ASVs
top_phylum <- head(sort(table(taxa$Phylum), decreasing = TRUE),15);top_phylum
top_class <- head(sort(table(taxa$Class), decreasing = TRUE),15);top_class
top_genus <- head(sort(table(taxa$Genus), decreasing = TRUE),20);top_genus

# --------------------------------------------------
# Taxonomic composition visualization
# --------------------------------------------------

# Obtain relative abundance
ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

# Top phylum
ps.phylum <- tax_glom(ps.rel, taxrank = "Phylum")

top10 <- names(
  sort(
    taxa_sums(ps.phylum),
    decreasing = TRUE
  )[1:10]
)

ps.phylum.top10 <- prune_taxa(
  top10,
  ps.phylum
)

# SampleID
plot_bar(
  ps.phylum.top10,
  x = "SampleID",
  fill = "Phylum"
)

# SampleID per Farm
plot_bar(
  ps.phylum.top10,
  x = "SampleID",
  fill = "Phylum"
) +
  facet_wrap(~Farm, scales = "free_x")

# Facet Farms
plot_bar(
  ps.phylum.top10,
  x = "CropStage",
  fill = "Phylum"
) +
  facet_wrap(~Farm)
