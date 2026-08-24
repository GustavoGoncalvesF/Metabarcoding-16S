#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 03_alpha_diversity.R
# Autor: Gustavo Goncalves
# Last edit: 29/07/2026
#
# Description:
# Script to exploratory analysis of alpha diversity from samples of a phyloseq object.
#
#####################################################################################

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

# Set up initial parameters
path <- "results/objects"; list.files(path) # Set path
ps <- readRDS(file.path(path, "phyloseq_object.rds")) # Load phyloseq object
theme_set(theme_bw()) # Set ggplot theme

# Observed + Shannon - Farm x CropStage
plot_richness(
  ps,
  x = "Farm",
  measures = c("Shannon", "Observed"),
  color = "CropStage"
)

# (Observed + Shannon - Farm x CropStage) + facet Year
plot_richness(
  ps,
  x = "Farm",
  measures = c("Shannon", "Observed"),
  color = "CropStage"
) +
  facet_wrap(~Year)

# Observed + Shannon - Farm x pH
plot_richness(
  ps,
  x="Farm",
  measures="Shannon",
  color="pH"
)

