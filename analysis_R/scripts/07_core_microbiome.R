#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 07_core_microbiome.R
# Autor: Gustavo Goncalves
# Last edit: 19/08/2026
#
# Description:
# Script to .
#
#####################################################################################

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(microbiome); packageVersion("microbiome")
library(ggplot2); packageVersion("ggplot2")

# Set up initial parameters
path <- "results/objects"; list.files(path) # Set path
ps <- readRDS(file.path(path, "phyloseq_object.rds")) # Load phyloseq object
theme_set(theme_bw()) # Set ggplot theme

# Transform to compositional abundances
ps.rel <- microbiome::transform(ps,"compositional")

# Working on Family level
ps.family <- tax_glom(ps.rel, taxrank = "Family")

# Calculate prevalence at different detection thresholds
prev.0 <- prevalence(ps.family, detection = 0, sort = TRUE)
prev.001 <- prevalence(ps.family, detection = 0.001, sort = TRUE)

# Summarize prevalence
data.frame(
  Taxon = names(prev.001),
  Prevalence = as.numeric(prev.001)
) |>
  head(20)

# Identify core taxa
core.family <- core_members(
  ps.family,
  detection = 0.001,
  prevalence = 0.90
)

# Extract taxonomy of core taxa
core.family.table <- data.frame(
  ASV = core.family,
  Family = as.character(
    tax_table(ps.family)[core.family, "Family"]
  ),
  Prevalence = as.numeric(
    prev.001[core.family]
  )
)

# Mean relative abundance
mean.abundance <- apply(
  otu_table(ps.family),
  2,
  mean
)

# Core microbiome table
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
  dplyr::arrange(desc(MeanAbundance))

core.family.table

# Number of core Families
length(core.family)

# Plot core Families visualization
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
  geom_hline(
    yintercept = 0.001,
    linetype = "dashed"
  ) +
  scale_y_log10() +
  labs(
    x = "Prevalence",
    y = "Mean relative abundance",
    color = "Core microbiome"
  )
