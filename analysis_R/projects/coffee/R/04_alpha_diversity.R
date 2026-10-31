#####################################################################################
# Project: Coffee Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 04_alpha_diversity.R
# Author: Gustavo Goncalves
# Last edit: 07/10/2026
#
# Description:
# Exploratory analysis of alpha diversity in coffee rhizosphere microbiome samples.
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
# Calculate alpha diversity
#####################################################################################

alpha <- estimate_richness(
  ps,
  measures = c("Observed", "Shannon")
)

# Add sample metadata
metadata <- data.frame(
  sample_data(ps)
)

alpha <- cbind(
  alpha,
  metadata[rownames(alpha), , drop = FALSE]
)


#####################################################################################
# Exploratory plots
#####################################################################################

# Observed richness and Shannon diversity - Farm x CropStage
plot_richness(
  ps,
  x = "Farm",
  measures = c("Observed", "Shannon"),
  color = "CropStage"
)

# Observed richness and Shannon diversity - Farm x CropStage, faceted by Year
plot_richness(
  ps,
  x = "Farm",
  measures = c("Observed", "Shannon"),
  color = "CropStage"
) +
  facet_wrap(~Year)


#####################################################################################
# Alpha diversity x pH
#####################################################################################

# Shannon diversity x pH
ggplot(
  alpha,
  aes(
    x = pH,
    y = Shannon
  )
) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    x = "Soil pH",
    y = "Shannon diversity"
  )

# Observed richness x pH
ggplot(
  alpha,
  aes(
    x = pH,
    y = Observed
  )
) +
  geom_point() +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  labs(
    x = "Soil pH",
    y = "Observed ASV richness"
  )
