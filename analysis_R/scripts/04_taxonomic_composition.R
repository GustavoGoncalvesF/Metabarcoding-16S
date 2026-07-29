#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 04_taxonomic_composition.R
# Autor: Gustavo Goncalves
# Last edit: 29/07/2026
#
# Description:
# Script to investigate relative abundance of taxas between samples from a phyloseq object.
#
#####################################################################################

# Set up packages
library(phyloseq); packageVersion("phyloseq")

# Set up initial parameters
path <- "results"; list.files(path) # Set path
ps <- readRDS(file.path(path, "phyloseq_object.rds")) # Load phyloseq object
theme_set(theme_bw()) # Set ggplot theme

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
