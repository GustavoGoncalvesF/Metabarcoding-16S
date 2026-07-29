#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 05_beta_diversity.R
# Autor: Gustavo Goncalves
# Last edit: 29/07/2026
#
# Description:
# Script to exploratory analysis of beta diversity from samples of a phyloseq object.
#
#####################################################################################

# Set up packages
library(phyloseq); packageVersion("phyloseq")

# Set up initial parameters
path <- "results"; list.files(path) # Set path
ps <- readRDS(file.path(path, "phyloseq_object.rds")) # Load phyloseq object
theme_set(theme_bw()) # Set ggplot theme

# Bray-Curtis (distance matrix)
bray <- phyloseq::distance(
  physeq = ps,
  method = "bray"
)

# PCoA ordinate
ord.pcoa <- ordinate(
  ps,
  method = "PCoA",
  distance = "bray"
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

# Plot by Farm + groups
plot_ordination(
  ps,
  ord.pcoa,
  color = "Farm"
) +
  geom_point(size = 3) +
  stat_ellipse()

# ==================================================================================
# PERMANOVA
# ==================================================================================
adonis2(bray ~ Farm, data = data.frame(sample_data(ps))) # Farm factor
adonis2(bray ~ CropStage, data = data.frame(sample_data(ps))) # CropStage factor
adonis2(bray ~ Quality, data = data.frame(sample_data(ps))) #Quality factor
adonis2(bray ~ Farm + CropStage + Quality, data = data.frame(sample_data(ps)), by = "margin")

# beta-dispersion test
disp.farm <- betadisper(bray, sample_data(ps)$Farm)
permutest(disp.farm)
boxplot(disp.farm) # Visualize beta-dispersion
TukeyHSD(disp.farm) # Show differences between farms
