
#####################################################################################
# Project: Alopecia Microbiome - 16S rRNA Metabarcoding
# Script: 04_alpha_diversity.R
# Author: Gustavo Goncalves
# Last edit: 11/11/2026
#
# Description:
# Exploratory alpha diversity analysis using Observed richness and Shannon.
# Analyses by Treatment, AlopeciaGrade, Age, and Timepoint.
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

# Recover original sample identifiers
sample_ids <- sample_names(ps)

idx <- match(
  rownames(alpha),
  make.names(sample_ids)
)

stopifnot(
  !anyNA(idx),
  !anyDuplicated(idx)
)

alpha$SampleID <- sample_ids[idx]

# Extract sample metadata
metadata <- data.frame(sample_data(ps))

# Align metadata with alpha diversity
metadata <- metadata[
  match(alpha$SampleID, metadata$SampleID),
  ,
  drop = FALSE
]

# Validate alignment
stopifnot(
  identical(alpha$SampleID, metadata$SampleID)
)

# Combine diversity metrics and metadata
alpha <- cbind(
  alpha,
  metadata[
    , !colnames(metadata) %in% "SampleID",
    drop = FALSE
  ]
)

# Patient identification
alpha$Patient <- factor(
  sub("[AB]$", "", alpha$SampleID)
)

# Sampling timepoint
alpha$Timepoint <- factor(
  ifelse(
    grepl("A$", alpha$SampleID),
    "Before",
    "After"
  ),
  levels = c("Before", "After")
)

# Validate paired design
stopifnot(
  nrow(alpha) == 40,
  nlevels(alpha$Patient) == 20,
  all(table(alpha$Patient, alpha$Timepoint) == 1)
)

# Summary statistics
print(summary(alpha[, c("Observed", "Shannon")]))

# Check missing metadata
print(colSums(is.na(alpha)))


#####################################################################################
# Prepare sensitivity datasets
#####################################################################################

# All samples
alpha_all <- alpha

# Without patient 16
alpha_no16 <- droplevels(
  alpha[as.character(alpha$Patient) != "16", , drop = FALSE]
)

# Check sample counts
stopifnot(
  nrow(alpha_all) == 40,
  nrow(alpha_no16) == 38
)

print(table(alpha_all$Timepoint))
print(table(alpha_no16$Timepoint))


#####################################################################################
# Generate plots and statistics
#####################################################################################

datasets <- list(
  all_samples = alpha_all,
  without_patient16 = alpha_no16
)

for (dataset_name in names(datasets)) {

  df <- datasets[[dataset_name]]

  output_path <- file.path(alpha_path, dataset_name)

  dir.create(
    output_path,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Save alpha diversity table
  write.csv(
    df,
    file.path(output_path, "alpha_diversity.csv"),
    row.names = FALSE
  )


  ###################################################################################
  # Alpha diversity by Treatment
  ###################################################################################

  # Observed richness
  p1 <- ggplot(
    df,
    aes(x = Treatment, y = Observed, fill = Treatment)
  ) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2.5) +
    labs(
      x = "Treatment",
      y = "Observed ASV richness"
    ) +
    guides(fill = "none")

  print(p1)

  ggsave(
    file.path(output_path, "01_observed_treatment.png"),
    plot = p1,
    width = 8,
    height = 6,
    dpi = 300
  )


  # Shannon diversity
  p2 <- ggplot(
    df,
    aes(x = Treatment, y = Shannon, fill = Treatment)
  ) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2.5) +
    labs(
      x = "Treatment",
      y = "Shannon diversity"
    ) +
    guides(fill = "none")

  print(p2)

  ggsave(
    file.path(output_path, "02_shannon_treatment.png"),
    plot = p2,
    width = 8,
    height = 6,
    dpi = 300
  )


  ###################################################################################
  # Alpha diversity by AlopeciaGrade
  ###################################################################################

  # Observed richness
  p3 <- ggplot(
    df,
    aes(x = AlopeciaGrade, y = Observed, fill = AlopeciaGrade)
  ) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2.5) +
    labs(
      x = "Alopecia grade",
      y = "Observed ASV richness"
    ) +
    guides(fill = "none")

  print(p3)

  ggsave(
    file.path(output_path, "03_observed_alopecia_grade.png"),
    plot = p3,
    width = 8,
    height = 6,
    dpi = 300
  )


  # Shannon diversity
  p4 <- ggplot(
    df,
    aes(x = AlopeciaGrade, y = Shannon, fill = AlopeciaGrade)
  ) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2.5) +
    labs(
      x = "Alopecia grade",
      y = "Shannon diversity"
    ) +
    guides(fill = "none")

  print(p4)

  ggsave(
    file.path(output_path, "04_shannon_alopecia_grade.png"),
    plot = p4,
    width = 8,
    height = 6,
    dpi = 300
  )


  ###################################################################################
  # Alpha diversity and Age
  ###################################################################################

  # Shannon diversity vs Age
  p5 <- ggplot(
    df,
    aes(x = Age, y = Shannon, color = Treatment)
  ) +
    geom_point(size = 3) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black"
    ) +
    labs(
      x = "Age (years)",
      y = "Shannon diversity"
    )

  print(p5)

  ggsave(
    file.path(output_path, "05_shannon_age.png"),
    plot = p5,
    width = 8,
    height = 6,
    dpi = 300
  )


  # Observed richness vs Age
  p6 <- ggplot(
    df,
    aes(x = Age, y = Observed, color = Treatment)
  ) +
    geom_point(size = 3) +
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "black"
    ) +
    labs(
      x = "Age (years)",
      y = "Observed ASV richness"
    )

  print(p6)

  ggsave(
    file.path(output_path, "06_observed_age.png"),
    plot = p6,
    width = 8,
    height = 6,
    dpi = 300
  )


  ###################################################################################
  # Paired alpha diversity by Timepoint
  ###################################################################################

  # Observed richness
  p7 <- ggplot(
    df,
    aes(x = Timepoint, y = Observed, group = Patient)
  ) +
    geom_line(color = "grey60", alpha = 0.5) +
    geom_point(aes(color = Timepoint), size = 2.5) +
    labs(
      title = "Observed richness before and after treatment",
      x = "Timepoint",
      y = "Observed ASVs"
    ) +
    guides(color = "none")

  print(p7)

  ggsave(
    file.path(output_path, "07_observed_paired.png"),
    plot = p7,
    width = 8,
    height = 6,
    dpi = 300
  )


  # Shannon diversity
  p8 <- ggplot(
    df,
    aes(x = Timepoint, y = Shannon, group = Patient)
  ) +
    geom_line(color = "grey60", alpha = 0.5) +
    geom_point(aes(color = Timepoint), size = 2.5) +
    labs(
      title = "Shannon diversity before and after treatment",
      x = "Timepoint",
      y = "Shannon diversity"
    ) +
    guides(color = "none")

  print(p8)

  ggsave(
    file.path(output_path, "08_shannon_paired.png"),
    plot = p8,
    width = 8,
    height = 6,
    dpi = 300
  )


  ###################################################################################
  # Paired Wilcoxon tests
  ###################################################################################

  # Extract Before samples
  before <- df[
    df$Timepoint == "Before",
    ,
    drop = FALSE
  ]

  # Extract After samples
  after <- df[
    df$Timepoint == "After",
    ,
    drop = FALSE
  ]

  # Match patients
  after <- after[
    match(
      as.character(before$Patient),
      as.character(after$Patient)
    ),
    ,
    drop = FALSE
  ]

  # Validate pairing
  stopifnot(
    nrow(before) > 1,
    nrow(before) == nrow(after),
    identical(
      as.character(before$Patient),
      as.character(after$Patient)
    )
  )

  # Observed richness
  test_observed <- wilcox.test(
    before$Observed,
    after$Observed,
    paired = TRUE,
    exact = FALSE
  )

  # Shannon diversity
  test_shannon <- wilcox.test(
    before$Shannon,
    after$Shannon,
    paired = TRUE,
    exact = FALSE
  )

  print(test_observed)
  print(test_shannon)

  capture.output(
    test_observed,
    test_shannon,
    file = file.path(
      output_path,
      "09_wilcoxon_paired.txt"
    )
  )
}


#####################################################################################
# Save analysis objects
#####################################################################################

saveRDS(
  alpha_all,
  file.path(alpha_path, "alpha_all_results.rds")
)

saveRDS(
  alpha_no16,
  file.path(alpha_path, "alpha_without_patient16_results.rds")
)

#####################################################################################
# End
#####################################################################################
