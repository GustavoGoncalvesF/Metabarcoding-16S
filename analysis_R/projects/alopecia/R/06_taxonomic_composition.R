#####################################################################################
# Project: Alopecia Microbiome - 16S rRNA Metabarcoding
# Script: 06_taxonomic_composition.R
# Author: Gustavo Goncalves
# Last edit: 11/11/2026
#
# Description:
# Taxonomic summary and relative abundance analysis of bacterial communities
# in scalp microbiome samples.
#
# Sensitivity analysis:
# - All samples
# - Without patient 16 (16A and 16B)
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

# Create output directory
tax_path <- file.path(results_path, "taxonomic_composition")
dir.create(tax_path, recursive = TRUE, showWarnings = FALSE)


#####################################################################################
# Prepare metadata
#####################################################################################

metadata <- data.frame(sample_data(ps))

# Preserve original sample identifiers
metadata$SampleID <- sample_names(ps)

# Identify patients and timepoints
metadata$Patient <- sub("[AB]$", "", metadata$SampleID)

metadata$Timepoint <- factor(
  ifelse(
    grepl("A$", metadata$SampleID),
    "Before",
    "After"
  ),
  levels = c("Before", "After")
)

metadata$Treatment <- factor(
  metadata$Treatment,
  levels = c("Control", "Active")
)

stopifnot(
  nrow(metadata) == 40,
  !anyNA(metadata$Treatment),
  all(table(metadata$Patient, metadata$Timepoint) == 1)
)


#####################################################################################
# Prepare sensitivity datasets
#####################################################################################

ps.no16 <- prune_samples(
  !sample_names(ps) %in% c("16A", "16B"),
  ps
)

ps.no16 <- prune_taxa(
  taxa_sums(ps.no16) > 0,
  ps.no16
)

datasets <- list(
  all_samples = ps,
  without_patient16 = ps.no16
)


#####################################################################################
# Taxonomic summary
#####################################################################################

for (dataset_name in names(datasets)) {

  ps.sub <- datasets[[dataset_name]]

  output_path <- file.path(tax_path, dataset_name)

  dir.create(
    output_path,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Extract taxonomy table
  taxa <- as.data.frame(
    tax_table(ps.sub)
  )

  # ASVs assigned at each taxonomic level
  assigned <- colSums(
    !is.na(taxa) & taxa != ""
  )

  assigned_percent <- round(
    100 * assigned / nrow(taxa),
    2
  )

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

  # Save taxonomic summary
  capture.output({

    cat("=====================================\n")
    cat("TAXONOMIC COMPOSITION SUMMARY\n")
    cat("Project: Alopecia Metabarcoding 16S\n")
    cat("Dataset:", dataset_name, "\n")
    cat("=====================================\n\n")

    cat("Number of samples:", nsamples(ps.sub), "\n")
    cat("Total ASVs:", ntaxa(ps.sub), "\n")
    cat("Total reads:", sum(sample_sums(ps.sub)), "\n\n")

    cat("-------------------------------------\n")
    cat("TAXONOMIC ASSIGNMENT\n")
    cat("-------------------------------------\n\n")

    print(
      data.frame(
        Assigned_ASVs = assigned,
        Percentage = assigned_percent
      )
    )

    cat("\n-------------------------------------\n")
    cat("TOP 15 PHYLA BY ASV RICHNESS\n")
    cat("-------------------------------------\n\n")
    print(top_phylum)

    cat("\n-------------------------------------\n")
    cat("TOP 15 CLASSES BY ASV RICHNESS\n")
    cat("-------------------------------------\n\n")
    print(top_class)

    cat("\n-------------------------------------\n")
    cat("TOP 20 GENERA BY ASV RICHNESS\n")
    cat("-------------------------------------\n\n")
    print(top_genus)

  }, file = file.path(output_path, "taxonomic_summary.txt"))
}


#####################################################################################
# Prepare relative abundance data
#####################################################################################

# Transform counts to relative abundance
ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

# Assign unclassified phyla
ps.tax <- ps.rel

tax.table <- tax_table(ps.tax)

tax.table[, "Phylum"][
  is.na(tax.table[, "Phylum"]) |
    tax.table[, "Phylum"] == ""
] <- "Unclassified"

tax_table(ps.tax) <- tax.table

# Agglomerate ASVs at Phylum level
ps.phylum <- tax_glom(
  ps.tax,
  taxrank = "Phylum"
)

# Extract abundance matrix (samples x taxa)
otu <- as(
  otu_table(ps.phylum),
  "matrix"
)

if (taxa_are_rows(ps.phylum)) {
  otu <- t(otu)
}

# Extract phylum names
phylum.names <- as.character(
  tax_table(ps.phylum)[, "Phylum"]
)

# Calculate global abundance per phylum
phylum.abundance <- tapply(
  colSums(otu),
  phylum.names,
  sum
)

# Rank phyla by global relative abundance
phylum.order <- names(
  sort(phylum.abundance, decreasing = TRUE)
)

# Select top 10 phyla
top10 <- head(phylum.order, 10)

# Group remaining phyla into Other
phylum.group <- ifelse(
  phylum.names %in% top10,
  phylum.names,
  "Other"
)

# Sum abundance by phylum within each sample
otu.grouped <- t(
  rowsum(
    t(otu),
    group = phylum.group
  )
)

# Convert to long-format data frame
phylum.df <- as.data.frame(
  otu.grouped,
  check.names = FALSE
)

phylum.df$SampleID <- rownames(phylum.df)

phylum.df <- phylum.df %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "Phylum",
    values_to = "Abundance"
  )

# Add sample metadata
phylum.df <- phylum.df %>%
  left_join(
    metadata[, c(
      "SampleID",
      "Patient",
      "Treatment",
      "Timepoint"
    )],
    by = "SampleID"
  )

# Check metadata alignment
stopifnot(
  !anyNA(phylum.df$Treatment),
  !anyNA(phylum.df$Timepoint)
)

# Set phylum abundance order
phylum.df$Phylum <- factor(
  phylum.df$Phylum,
  levels = c(top10, "Other")
)

# Order samples numerically by patient, then A/B
sample.order <- metadata$SampleID[
  order(
    as.numeric(metadata$Patient),
    metadata$Timepoint
  )
]

phylum.df$SampleID <- factor(
  phylum.df$SampleID,
  levels = sample.order
)


#####################################################################################
# Plot settings
#####################################################################################

# Define colors for the top 10 phyla
phylum.colors <- scales::hue_pal()(length(top10))

names(phylum.colors) <- top10

# Add light gray for Other
phylum.colors <- c(
  phylum.colors,
  Other = "#D3D3D3"
)


#####################################################################################
# Relative abundance plots and tables
#####################################################################################

for (dataset_name in names(datasets)) {

  output_path <- file.path(tax_path, dataset_name)

  # Select samples
  df <- phylum.df

  if (dataset_name == "without_patient16") {
    df <- df %>%
      filter(Patient != "16")
  }

  df$SampleID <- droplevels(df$SampleID)

  # Save individual relative abundances
  write.csv(
    df,
    file.path(
      output_path,
      "phylum_relative_abundance_by_sample.csv"
    ),
    row.names = FALSE
  )


  ###################################################################################
  # Mean relative abundance by Treatment x Timepoint
  ###################################################################################

  phylum.mean <- df %>%
    group_by(
      Treatment,
      Timepoint,
      Phylum
    ) %>%
    summarise(
      MeanAbundance = mean(Abundance),
      .groups = "drop"
    )

  # Save mean relative abundance table
  write.csv(
    phylum.mean,
    file.path(
      output_path,
      "mean_phylum_abundance_by_group.csv"
    ),
    row.names = FALSE
  )


  ###################################################################################
  # Plot 1 - Mean relative abundance by Treatment
  ###################################################################################

  p1 <- ggplot(
    phylum.mean,
    aes(
      x = Treatment,
      y = MeanAbundance,
      fill = Phylum
    )
  ) +
    geom_col(
      position = position_stack(reverse = TRUE),
      color = "grey50",
      linewidth = 0.2
    ) +
    facet_wrap(~Timepoint) +
    scale_fill_manual(
      values = phylum.colors,
      breaks = c(top10, "Other")
    ) +
    scale_y_continuous(
      labels = scales::percent
    ) +
    labs(
      x = "Treatment",
      y = "Mean relative abundance",
      fill = "Phylum"
    )

  print(p1)

  ggsave(
    file.path(
      output_path,
      "01_mean_phylum_treatment_timepoint.png"
    ),
    plot = p1,
    width = 10,
    height = 6,
    dpi = 300
  )


  ###################################################################################
  # Plot 2 - Mean relative abundance by Timepoint
  ###################################################################################

  p2 <- ggplot(
    phylum.mean,
    aes(
      x = Timepoint,
      y = MeanAbundance,
      fill = Phylum
    )
  ) +
    geom_col(
      position = position_stack(reverse = TRUE),
      color = "grey50",
      linewidth = 0.2
    ) +
    facet_wrap(~Treatment) +
    scale_fill_manual(
      values = phylum.colors,
      breaks = c(top10, "Other")
    ) +
    scale_y_continuous(
      labels = scales::percent
    ) +
    labs(
      x = "Timepoint",
      y = "Mean relative abundance",
      fill = "Phylum"
    )

  print(p2)

  ggsave(
    file.path(
      output_path,
      "02_mean_phylum_timepoint_treatment.png"
    ),
    plot = p2,
    width = 10,
    height = 6,
    dpi = 300
  )


  ###################################################################################
  # Plot 3 - Relative abundance of individual samples
  ###################################################################################

  p3 <- ggplot(
    df,
    aes(
      x = SampleID,
      y = Abundance,
      fill = Phylum
    )
  ) +
    geom_col(
      position = position_stack(reverse = TRUE),
      color = "grey50",
      linewidth = 0.2
    ) +
    scale_fill_manual(
      values = phylum.colors,
      breaks = c(top10, "Other")
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

  print(p3)

  ggsave(
    file.path(
      output_path,
      "03_phylum_all_samples.png"
    ),
    plot = p3,
    width = 16,
    height = 7,
    dpi = 300
  )


  ###################################################################################
  # Plot 4 - Individual samples grouped by Treatment
  ###################################################################################

  p4 <- p3 +
    facet_grid(
      ~Treatment,
      scales = "free_x",
      space = "free_x"
    )

  print(p4)

  ggsave(
    file.path(
      output_path,
      "04_phylum_samples_by_treatment.png"
    ),
    plot = p4,
    width = 16,
    height = 7,
    dpi = 300
  )

}

#####################################################################################
# End
#####################################################################################
