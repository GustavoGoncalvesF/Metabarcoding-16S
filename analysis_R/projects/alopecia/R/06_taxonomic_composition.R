
#####################################################################################
# Project: Alopecia Microbiome - 16S rRNA Metabarcoding
# Script: 06_taxonomic_composition.R
# Author: Gustavo Goncalves
# Last edit: 08/10/2026
#
# Description:
# Exploratory analysis of taxonomic composition and relative abundance
# of bacterial taxa across scalp microbiome samples and treatment groups.
#
#####################################################################################

# Load project configuration
source("projects/alopecia/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")

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
  !is.na(taxa) & taxa != ""
)

assigned_percent <- round(
  100 * assigned / nrow(taxa),
  2
)

assigned
assigned_percent

# Taxa with highest number of assigned ASVs
top_phylum <- head(
  sort(table(taxa$Phylum), decreasing = TRUE),
  15
)

top_class <- head(
  sort(table(taxa$Class), decreasing = TRUE),
  15
)

top_genus <- head(
  sort(table(taxa$Genus), decreasing = TRUE),
  20
)

top_phylum
top_class
top_genus


#####################################################################################
# Relative abundance
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
# Keep taxa without Phylum assignment
ps.phylum <- tax_glom(
  ps.rel,
  taxrank = "Phylum",
  NArm = FALSE
)

# Extract relative abundance and taxonomy
phylum.data <- psmelt(ps.phylum)

# Assign unclassified taxa to a separate category
phylum.data$Phylum <- as.character(phylum.data$Phylum)

phylum.data$Phylum[
  is.na(phylum.data$Phylum) |
    phylum.data$Phylum == ""
] <- "Unclassified"

# Identify the 10 most abundant phyla
top10 <- phylum.data |>
  group_by(Phylum) |>
  summarise(
    TotalAbundance = sum(Abundance),
    .groups = "drop"
  ) |>
  arrange(desc(TotalAbundance)) |>
  slice_head(n = 10) |>
  pull(Phylum)

# Group remaining phyla as Other
phylum.data <- phylum.data |>
  mutate(
    PhylumGroup = if_else(
      Phylum %in% top10,
      Phylum,
      "Other"
    )
  )

# Aggregate relative abundance by sample
phylum.sample <- phylum.data |>
  group_by(
    Sample,
    Treatment,
    PhylumGroup
  ) |>
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

# Include absent phyla with zero abundance
phylum.sample <- phylum.sample |>
  complete(
    nesting(Sample, Treatment),
    PhylumGroup = unique(c(top10, "Other")),
    fill = list(Abundance = 0)
  )


#####################################################################################
# Taxonomic composition by sample
#####################################################################################

ggplot(
  phylum.sample,
  aes(
    x = Sample,
    y = Abundance,
    fill = PhylumGroup
  )
) +
  geom_col(width = 0.9) +
  facet_wrap(
    ~Treatment,
    scales = "free_x"
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Sample",
    y = "Relative abundance",
    fill = "Phylum"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )


#####################################################################################
# Mean taxonomic composition by Treatment
#####################################################################################

# Calculate mean relative abundance across samples
phylum.mean <- phylum.sample |>
  group_by(
    Treatment,
    PhylumGroup
  ) |>
  summarise(
    MeanAbundance = mean(Abundance),
    .groups = "drop"
  )

# Plot mean relative abundance
ggplot(
  phylum.mean,
  aes(
    x = Treatment,
    y = MeanAbundance,
    fill = PhylumGroup
  )
) +
  geom_col(width = 0.7) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Treatment",
    y = "Mean relative abundance",
    fill = "Phylum"
  )


#####################################################################################
# Genus-level composition
#####################################################################################

# Agglomerate ASVs at Genus level
ps.genus <- tax_glom(
  ps.rel,
  taxrank = "Genus",
  NArm = FALSE
)

# Extract relative abundance and taxonomy
genus.data <- psmelt(ps.genus)

# Assign unclassified taxa to a separate category
genus.data$Genus <- as.character(genus.data$Genus)

genus.data$Genus[
  is.na(genus.data$Genus) |
    genus.data$Genus == ""
] <- "Unclassified"

# Identify the 15 most abundant genera
top15 <- genus.data |>
  group_by(Genus) |>
  summarise(
    TotalAbundance = sum(Abundance),
    .groups = "drop"
  ) |>
  arrange(desc(TotalAbundance)) |>
  slice_head(n = 15) |>
  pull(Genus)

# Group remaining genera as Other
genus.data <- genus.data |>
  mutate(
    GenusGroup = if_else(
      Genus %in% top15,
      Genus,
      "Other"
    )
  )

# Aggregate relative abundance by sample
genus.sample <- genus.data |>
  group_by(
    Sample,
    Treatment,
    GenusGroup
  ) |>
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

# Include absent genera with zero abundance
genus.sample <- genus.sample |>
  complete(
    nesting(Sample, Treatment),
    GenusGroup = unique(c(top15, "Other")),
    fill = list(Abundance = 0)
  )


#####################################################################################
# Genus composition by Treatment
#####################################################################################

# Calculate mean relative abundance across samples
genus.mean <- genus.sample |>
  group_by(
    Treatment,
    GenusGroup
  ) |>
  summarise(
    MeanAbundance = mean(Abundance),
    .groups = "drop"
  )

# Plot mean relative abundance
ggplot(
  genus.mean,
  aes(
    x = Treatment,
    y = MeanAbundance,
    fill = GenusGroup
  )
) +
  geom_col(width = 0.7) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Treatment",
    y = "Mean relative abundance",
    fill = "Genus"
  )


#####################################################################################
# End
#####################################################################################
