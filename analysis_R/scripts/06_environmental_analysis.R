#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 06_environmental_analysis.R
# Autor: Gustavo Goncalves
# Last edit: 29/07/2026
#
# Description:
# Script to include metadata as pH, altitude, phosphorus and etc in environmental analysis.
#
#####################################################################################

# Set up packages
library(phyloseq); packageVersion("phyloseq")

# Set up initial parameters
path <- "results"; list.files(path) # Set path
ps <- readRDS(file.path(path, "phyloseq_object.rds")) # Load phyloseq object
theme_set(theme_bw()) # Set ggplot theme

# Environment data frame
env <- data.frame(sample_data(ps))[
  ,
  c(
    "pH",
    "Phosphorus",
    "Potassium",
    "CEC",
    "Altitude"
  )
]
pcoa_scores <- as.data.frame(ord.pcoa$vectors) # Extract scores from ord.pcoa

# envfit to environmental factors analysis
ef <- envfit(
  pcoa_scores[,1:2],
  env,
  permutations = 999
); ef
