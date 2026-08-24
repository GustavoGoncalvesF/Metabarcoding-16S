#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 07_core_microbiome.R
# Autor: Gustavo Goncalves
# Last edit: 29/07/2026
#
# Description:
# Script to .
#
#####################################################################################

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

# Set up initial parameters
path <- "results/objects"; list.files(path) # Set path
ps <- readRDS(file.path(path, "phyloseq_object.rds")) # Load phyloseq object
theme_set(theme_bw()) # Set ggplot theme
