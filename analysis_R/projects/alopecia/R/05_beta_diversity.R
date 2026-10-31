#####################################################################################
# Project: Alopecia Microbiome - 16S rRNA Metabarcoding
# Script: 05_beta_diversity.R
# Author: Gustavo Goncalves
# Last edit: 08/10/2026
#
# Description:
# Exploratory beta diversity analysis of scalp microbiome samples using
# Bray-Curtis dissimilarity, PCoA, PERMANOVA, and dispersion tests.
#
#####################################################################################

# Load project configuration
source("projects/alopecia/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(vegan); packageVersion("vegan")

# Load phyloseq object
ps <- readRDS(
  file.path(objects_path, "phyloseq_object.rds")
)

# Set ggplot theme
theme_set(theme_bw())

# Set random seed for reproducibility
set.seed(123)


#####################################################################################
# Prepare data
#####################################################################################

# Extract sample metadata
metadata <- data.frame(
  sample_data(ps)
)

# Transform counts to relative abundance
ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)


#####################################################################################
# Bray-Curtis dissimilarity
#####################################################################################

bray <- phyloseq::distance(
  physeq = ps.rel,
  method = "bray"
)


#####################################################################################
# Principal Coordinates Analysis (PCoA)
#####################################################################################

# Calculate PCoA
ord.pcoa <- ordinate(
  ps.rel,
  method = "PCoA",
  distance = bray
)

# PCoA by Treatment
plot_ordination(
  ps.rel,
  ord.pcoa,
  color = "Treatment"
) +
  geom_point(size = 3)


# PCoA by AlopeciaGrade
plot_ordination(
  ps.rel,
  ord.pcoa,
  color = "AlopeciaGrade"
) +
  geom_point(size = 3)


# PCoA by Age
plot_ordination(
  ps.rel,
  ord.pcoa,
  color = "Age"
) +
  geom_point(size = 3)


#####################################################################################
# PERMANOVA
#####################################################################################

# Treatment effect
adonis.treatment <- adonis2(
  bray ~ Treatment,
  data = metadata,
  permutations = 999
)

adonis.treatment


# Alopecia severity
# Exclude samples without AlopeciaGrade
valid.grade <- !is.na(metadata$AlopeciaGrade)

bray.grade <- as.dist(
  as.matrix(bray)[valid.grade, valid.grade]
)

metadata.grade <- metadata[
  valid.grade,
  ,
  drop = FALSE
]

adonis.grade <- adonis2(
  bray.grade ~ AlopeciaGrade,
  data = metadata.grade,
  permutations = 999
)

adonis.grade


# Age
adonis.age <- adonis2(
  bray ~ Age,
  data = metadata,
  permutations = 999
)

adonis.age


#####################################################################################
# Multivariate dispersion
#####################################################################################

# Test homogeneity of dispersion between treatment groups
disp.treatment <- betadisper(
  bray,
  metadata$Treatment
)

# Permutation test
permutest(
  disp.treatment,
  permutations = 999
)

# Visualize dispersion
boxplot(
  disp.treatment,
  xlab = "Treatment",
  ylab = "Distance to group median"
)
