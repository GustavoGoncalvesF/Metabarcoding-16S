#####################################################################################
# Project: Paspalum Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 04_alpha_diversity.R
# Author: Gustavo Goncalves
# Last edit: 09/11/2026
#
# Description:
# Exploratory analysis of alpha diversity in Paspalum rhizosphere samples.
#
#####################################################################################

# Load project configuration
source("projects/paspalum/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")

# Load phyloseq object
ps <- readRDS(
  file.path(objects_path, "phyloseq_object.rds")
)

# Set ggplot theme
theme_set(theme_bw())

# Create output directory
alpha_path <- file.path(results_path, "alpha_diversity")
dir.create(alpha_path, recursive = TRUE, showWarnings = FALSE)

#####################################################################################
# Calculate alpha diversity
#####################################################################################

# Calculate alpha diversity
alpha <- estimate_richness(
  ps,
  measures = c("Observed", "Shannon")
)

# Restore original sample identifiers
rownames(alpha) <- sample_names(ps)

# Extract sample metadata
metadata <- data.frame(
  sample_data(ps)
)

# Combine diversity metrics and metadata
alpha <- cbind(
  alpha,
  metadata[rownames(alpha), , drop = FALSE]
)

# Add sample identifiers
alpha$SampleID <- rownames(alpha)

# Check metadata
head(alpha[, c(
  "SampleID", "Species", "Accession",
  "Origin", "Observed", "Shannon"
)])

#####################################################################################
# Alpha diversity summary
#####################################################################################

# Save individual alpha diversity metrics
write.csv2(
  alpha,
  file.path(alpha_path, "alpha_diversity_metrics.csv"),
  row.names = FALSE
)

# Function to summarize alpha diversity by group
alpha_summary <- function(data, group) {

  aggregate(
    data[, c("Observed", "Shannon")],
    by = list(Group = data[[group]]),
    FUN = function(x) c(
      N = length(x),
      Mean = mean(x),
      SD = sd(x),
      Median = median(x),
      Min = min(x),
      Max = max(x)
    )
  )
}

# Save summary report
capture.output({

  cat("=====================================\n")
  cat("ALPHA DIVERSITY SUMMARY\n")
  cat("Project: Paspalum Metabarcoding 16S\n")
  cat("=====================================\n\n")

  cat("Number of samples:", nrow(alpha), "\n\n")

  cat("-------------------------------------\n")
  cat("GENERAL SUMMARY\n")
  cat("-------------------------------------\n\n")

  print(
    summary(alpha[, c("Observed", "Shannon")])
  )

  cat("\n-------------------------------------\n")
  cat("ALPHA DIVERSITY BY ORIGIN\n")
  cat("-------------------------------------\n\n")

  print(alpha_summary(alpha, "Origin"))

  cat("\n-------------------------------------\n")
  cat("ALPHA DIVERSITY BY SPECIES\n")
  cat("-------------------------------------\n\n")

  print(alpha_summary(alpha, "Species"))

  cat("\n-------------------------------------\n")
  cat("ALPHA DIVERSITY BY ACCESSION\n")
  cat("-------------------------------------\n\n")

  print(alpha_summary(alpha, "Accession"))

  cat("\n-------------------------------------\n")
  cat("ALPHA DIVERSITY BY ACCESSION AND ORIGIN\n")
  cat("-------------------------------------\n\n")

  for (acc in unique(alpha$Accession)) {

    cat("\nAccession:", as.character(acc), "\n\n")

    print(
      alpha_summary(
        alpha[alpha$Accession == acc, , drop = FALSE],
        "Origin"
      )
    )
  }

}, file = file.path(alpha_path, "alpha_diversity_summary.txt"))


#####################################################################################
# Exploratory plots
#####################################################################################

# Observed richness and Shannon diversity - Origin x Species
p1 <- plot_richness(
  ps,
  x = "Origin",
  measures = c("Observed", "Shannon"),
  color = "Species"
)

print(p1)

ggsave(
  file.path(alpha_path, "01_origin_species.png"),
  plot = p1,
  width = 10,
  height = 6,
  dpi = 300
)


# Observed richness and Shannon diversity - Origin x Accession
p2 <- plot_richness(
  ps,
  x = "Origin",
  measures = c("Observed", "Shannon"),
  color = "Accession"
)

print(p2)

ggsave(
  file.path(alpha_path, "02_origin_accession.png"),
  plot = p2,
  width = 10,
  height = 6,
  dpi = 300
)


# Observed richness and Shannon diversity - Accession x Origin
p3 <- plot_richness(
  ps,
  x = "Accession",
  measures = c("Observed", "Shannon"),
  color = "Origin"
)

print(p3)

ggsave(
  file.path(alpha_path, "03_accession_origin.png"),
  plot = p3,
  width = 10,
  height = 6,
  dpi = 300
)


# Observed richness and Shannon diversity - Accession x Origin, faceted by Species
p4 <- p3 + facet_wrap(~Species)

print(p4)

ggsave(
  file.path(alpha_path, "04_accession_origin_by_species.png"),
  plot = p4,
  width = 12,
  height = 8,
  dpi = 300
)
