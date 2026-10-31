#####################################################################################
# Project: Coffee Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 06_taxonomic_composition.R
# Author: Gustavo Goncalves
# Last edit: 07/10/2026
#
# Description:
# Script to investigate taxonomic composition and relative abundance of taxa
# in coffee rhizosphere microbiome samples.
#
#####################################################################################

# Load project configuration
source("projects/coffee/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

# Load phyloseq object
ps <- readRDS(
  file.path(objects_path, "phyloseq_object.rds")
)

# Set ggplot theme
theme_set(theme_bw())


#####################################################################################
# Taxonomic summary
#####################################################################################

# Extract taxonomy table
taxa <- as.data.frame(
  tax_table(ps)
)

# ASVs assigned at each taxonomic level
assigned <- colSums(
  !is.na(taxa)
)

assigned_percent <- round(
  100 * assigned / nrow(taxa),
  2
)

assigned
assigned_percent

# Taxa with highest number of assigned ASVs
top_phylum <- head(
  sort(
    table(taxa$Phylum),
    decreasing = TRUE
  ),
  15
)

top_class <- head(
  sort(
    table(taxa$Class),
    decreasing = TRUE
  ),
  15
)

top_genus <- head(
  sort(
    table(taxa$Genus),
    decreasing = TRUE
  ),
  20
)

top_phylum
top_class
top_genus


#####################################################################################
# Taxonomic composition visualization
#####################################################################################

# Transform counts to relative abundance
ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)


#####################################################################################
# Phylum-level composition
#####################################################################################

# Agglomerate ASVs at Phylum level
ps.phylum <- tax_glom(
  ps.rel,
  taxrank = "Phylum"
)

# Select the 10 most abundant phyla
top10 <- names(
  head(
    sort(
      taxa_sums(ps.phylum),
      decreasing = TRUE
    ),
    10
  )
)

ps.phylum.top10 <- prune_taxa(
  top10,
  ps.phylum
)

# Relative abundance by SampleID
plot_bar(
  ps.phylum.top10,
  x = "SampleID",
  fill = "Phylum"
)

# Relative abundance by SampleID, faceted by Farm
plot_bar(
  ps.phylum.top10,
  x = "SampleID",
  fill = "Phylum"
) +
  facet_wrap(
    ~Farm,
    scales = "free_x"
  )

# Relative abundance by CropStage, faceted by Farm
plot_bar(
  ps.phylum.top10,
  x = "CropStage",
  fill = "Phylum"
) +
  facet_wrap(~Farm)
