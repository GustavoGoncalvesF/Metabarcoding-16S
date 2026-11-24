#####################################################################################
# Project: Alopecia Microbiome - 16S rRNA Metabarcoding
# Script: 05_beta_diversity.R
# Author: Gustavo Goncalves
# Last edit: 11/11/2026
#
# Description:
# Exploratory beta diversity analysis using Bray-Curtis dissimilarity,
# PCoA, paired PERMANOVA, and multivariate dispersion tests.
#
# Sensitivity analysis:
# - All samples
# - Without patient 16 (samples 16A and 16B)
#
#####################################################################################

# Load project configuration
source("projects/alopecia/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(vegan); packageVersion("vegan")
library(permute)

# Load phyloseq object
ps <- readRDS(
  file.path(objects_path, "phyloseq_object.rds")
)

# Set ggplot theme
theme_set(theme_bw())

# Set random seed for reproducibility
set.seed(123)

# Create output directory
beta_path <- file.path(results_path, "beta_diversity")
dir.create(beta_path, recursive = TRUE, showWarnings = FALSE)


#####################################################################################
# Prepare data
#####################################################################################

# Extract sample metadata
metadata <- data.frame(sample_data(ps))

# Experimental variables
metadata$Patient <- factor(metadata$Patient)

metadata$Timepoint <- factor(
  metadata$Timepoint,
  levels = c("Before", "After")
)

# Check sample alignment
stopifnot(
  identical(rownames(metadata), sample_names(ps)),
  !anyNA(metadata$Patient),
  !anyNA(metadata$Timepoint)
)

# Confirm paired design
print(table(metadata$Timepoint))
print(table(metadata$Patient, metadata$Timepoint))


#####################################################################################
# Beta diversity analysis function
#####################################################################################

run_beta_diversity <- function(ps, output_path) {

  dir.create(
    output_path,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Remove ASVs absent from this subset
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)

  # Extract metadata
  metadata <- data.frame(sample_data(ps))

  metadata$Patient <- factor(metadata$Patient)

  metadata$Timepoint <- factor(
    metadata$Timepoint,
    levels = c("Before", "After")
  )

  # Confirm complete patient pairs
  pair_table <- table(
    metadata$Patient,
    metadata$Timepoint
  )

  stopifnot(all(pair_table == 1))

  # Transform counts to relative abundance
  ps.rel <- transform_sample_counts(
    ps,
    function(x) x / sum(x)
  )


  ###################################################################################
  # Bray-Curtis dissimilarity
  ###################################################################################

  bray <- phyloseq::distance(
    physeq = ps.rel,
    method = "bray"
  )

  # Align metadata with distance matrix
  metadata <- metadata[
    attr(bray, "Labels"),
    ,
    drop = FALSE
  ]

  saveRDS(
    bray,
    file.path(output_path, "bray_curtis.rds")
  )


  ###################################################################################
  # Principal Coordinates Analysis (PCoA)
  ###################################################################################

  ord.pcoa <- ordinate(
    ps.rel,
    method = "PCoA",
    distance = bray
  )

  # PCoA by Timepoint
  p1 <- plot_ordination(
    ps.rel,
    ord.pcoa,
    color = "Timepoint"
  ) +
    geom_point(size = 3) +
    labs(
      title = "PCoA - Bray-Curtis",
      color = "Timepoint"
    )

  print(p1)

  ggsave(
    file.path(output_path, "01_pcoa_timepoint.png"),
    plot = p1,
    width = 8,
    height = 6,
    dpi = 300
  )


  # PCoA with paired samples connected
  p2 <- plot_ordination(
    ps.rel,
    ord.pcoa,
    color = "Timepoint"
  ) +
    geom_line(
      aes(group = Patient),
      color = "grey60",
      alpha = 0.6
    ) +
    geom_point(size = 3) +
    labs(
      title = "PCoA - Paired samples",
      color = "Timepoint"
    )

  print(p2)

  ggsave(
    file.path(output_path, "02_pcoa_paired.png"),
    plot = p2,
    width = 8,
    height = 6,
    dpi = 300
  )


  ###################################################################################
  # PERMANOVA - Paired design
  ###################################################################################

  # Restrict permutations to within-patient swaps
  permutation_design <- permute::how(
    nperm = 999,
    blocks = metadata$Patient
  )

  # Adjust for patient identity
  adonis.result <- adonis2(
    bray ~ Patient + Timepoint,
    data = metadata,
    permutations = permutation_design,
    by = "margin"
  )

  # Extract Timepoint effect
  adonis.timepoint <- adonis.result[
    "Timepoint",
    ,
    drop = FALSE
  ]

  print(adonis.timepoint)

  write.csv(
    as.data.frame(adonis.timepoint),
    file.path(output_path, "03_permanova_timepoint.csv")
  )


  ###################################################################################
  # Multivariate dispersion
  ###################################################################################

  disp.timepoint <- betadisper(
    bray,
    metadata$Timepoint
  )

  dispersion.test <- permutest(
    disp.timepoint,
    permutations = 999
  )

  print(dispersion.test)

  capture.output(
    dispersion.test,
    file = file.path(
      output_path,
      "04_dispersion_timepoint.txt"
    )
  )

  png(
    file.path(output_path, "05_dispersion_timepoint.png"),
    width = 2400,
    height = 1800,
    res = 300
  )

  boxplot(
    disp.timepoint,
    xlab = "Timepoint",
    ylab = "Distance to group median",
    main = "Multivariate dispersion"
  )

  dev.off()


  ###################################################################################
  # Return results
  ###################################################################################

  invisible(list(
    bray = bray,
    ordination = ord.pcoa,
    permanova = adonis.timepoint,
    dispersion = dispersion.test
  ))
}


#####################################################################################
# Main analysis - All samples
#####################################################################################

beta_all <- run_beta_diversity(
  ps,
  file.path(beta_path, "all_samples")
)


#####################################################################################
# Sensitivity analysis - Without patient 16
#####################################################################################

# Remove the complete patient pair
ps.no16 <- subset_samples(
  ps,
  !(Patient %in% c("16"))
)

ps.no16 <- prune_taxa(
  taxa_sums(ps.no16) > 0,
  ps.no16
)

beta_no16 <- run_beta_diversity(
  ps.no16,
  file.path(beta_path, "without_patient16")
)


#####################################################################################
# Save analysis objects
#####################################################################################

saveRDS(
  beta_all,
  file.path(beta_path, "beta_all_results.rds")
)

saveRDS(
  beta_no16,
  file.path(beta_path, "beta_without_patient16_results.rds")
)

#####################################################################################
# End
#####################################################################################
