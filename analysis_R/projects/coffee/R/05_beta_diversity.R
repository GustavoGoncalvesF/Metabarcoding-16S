#####################################################################################
# Project: Coffee Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 05_beta_diversity.R
# Author: Gustavo Goncalves
# Last edit: 07/10/2026
#
# Description:
# Exploratory analysis of beta diversity in coffee rhizosphere microbiome samples,
# including PCoA ordination, PERMANOVA, and beta-dispersion analysis.
#
#####################################################################################

# Load project configuration
source("projects/coffee/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(vegan); packageVersion("vegan")

# Load phyloseq object
ps <- readRDS(
  file.path(objects_path, "phyloseq_object.rds")
)

# Extract sample metadata
metadata <- data.frame(
  sample_data(ps)
)

# Set ggplot theme
theme_set(theme_bw())


#####################################################################################
# Bray-Curtis dissimilarity
#####################################################################################

bray <- phyloseq::distance(
  physeq = ps,
  method = "bray"
)


#####################################################################################
# PCoA ordination
#####################################################################################

ord.pcoa <- ordinate(
  ps,
  method = "PCoA",
  distance = bray
)

# Plot by CropStage
plot_ordination(
  ps,
  ord.pcoa,
  color = "CropStage"
) +
  geom_point(size = 3)

# Plot by Quality
plot_ordination(
  ps,
  ord.pcoa,
  color = "Quality"
) +
  geom_point(size = 3)

# Plot by Farm
plot_ordination(
  ps,
  ord.pcoa,
  color = "Farm"
) +
  geom_point(size = 3)

# Plot by Farm with group ellipses
plot_ordination(
  ps,
  ord.pcoa,
  color = "Farm"
) +
  geom_point(size = 3) +
  stat_ellipse()


#####################################################################################
# PERMANOVA
#####################################################################################

# Individual effects

adonis2(
  bray ~ Farm,
  data = metadata,
  permutations = 999
)

adonis2(
  bray ~ CropStage,
  data = metadata,
  permutations = 999
)

adonis2(
  bray ~ Quality,
  data = metadata,
  permutations = 999
)

# Marginal effects
adonis2(
  bray ~ Farm + CropStage + Quality,
  data = metadata,
  permutations = 999,
  by = "margin"
)


#####################################################################################
# Beta-dispersion
#####################################################################################

# Test homogeneity of multivariate dispersion among farms
disp.farm <- betadisper(
  bray,
  metadata$Farm
)

permutest(
  disp.farm,
  permutations = 999
)

# Visualize beta-dispersion
boxplot(disp.farm)

# Pairwise comparisons among farms
TukeyHSD(disp.farm)


#####################################################################################
# End
#####################################################################################
