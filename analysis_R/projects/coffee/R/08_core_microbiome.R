#####################################################################################
# Project: Coffee Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 08_core_microbiome.R
# Author: Gustavo Goncalves
# Last edit: 07/10/2026
#
# Description:
# Script to identify and explore the core microbiome at the Family level
# in coffee rhizosphere microbiome samples.
#
#####################################################################################

# Load project configuration
source("projects/coffee/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(microbiome); packageVersion("microbiome")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")

# Load phyloseq object
ps <- readRDS(
  file.path(objects_path, "phyloseq_object.rds")
)

# Set ggplot theme
theme_set(theme_bw())


#####################################################################################
# Prepare Family-level data
#####################################################################################

# Transform counts to relative abundance
ps.rel <- microbiome::transform(
  ps,
  "compositional"
)

# Agglomerate taxa at Family level
ps.family <- tax_glom(
  ps.rel,
  taxrank = "Family"
)


#####################################################################################
# Prevalence
#####################################################################################

# Calculate prevalence at different detection thresholds
prev.0 <- prevalence(
  ps.family,
  detection = 0,
  sort = TRUE
)

prev.001 <- prevalence(
  ps.family,
  detection = 0.001,
  sort = TRUE
)

# Summarize prevalence
data.frame(
  Taxon = names(prev.001),
  Prevalence = as.numeric(prev.001)
) |>
  head(20)


#####################################################################################
# Core microbiome
#####################################################################################

# Identify core Families
# Detection threshold: 0.1% relative abundance
# Prevalence threshold: 90% of samples
core.family <- core_members(
  ps.family,
  detection = 0.001,
  prevalence = 0.90
)

# Calculate mean relative abundance
mean.abundance <- taxa_sums(ps.family) /
  nsamples(ps.family)

# Create core Family table
core.family.table <- data.frame(
  ASV = core.family,
  Family = as.character(
    tax_table(ps.family)[core.family, "Family"]
  ),
  Prevalence = as.numeric(
    prev.001[core.family]
  ),
  MeanAbundance = mean.abundance[core.family]
)

core.family.table <- core.family.table |>
  arrange(
    desc(MeanAbundance)
  )

core.family.table

# Number of core Families
length(core.family)


#####################################################################################
# Family summary
#####################################################################################

family.summary <- data.frame(
  ASV = names(prev.001),
  Family = as.character(
    tax_table(ps.family)[names(prev.001), "Family"]
  ),
  Prevalence = as.numeric(prev.001),
  MeanAbundance = mean.abundance[names(prev.001)]
)

# Identify actual core members
family.summary$Core <-
  family.summary$ASV %in% core.family


#####################################################################################
# Core microbiome visualization
#####################################################################################

ggplot(
  family.summary,
  aes(
    x = Prevalence,
    y = MeanAbundance,
    color = Core
  )
) +
  geom_point(
    size = 3,
    alpha = 0.7
  ) +
  geom_vline(
    xintercept = 0.90,
    linetype = "dashed"
  ) +
  scale_y_log10() +
  labs(
    x = "Prevalence",
    y = "Mean relative abundance",
    color = "Core microbiome"
  )
