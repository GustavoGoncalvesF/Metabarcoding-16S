#####################################################################################
# Project: Alopecia Microbiome - 16S rRNA Metabarcoding
# Script: 04_alpha_diversity.R
# Author: Gustavo Goncalves
# Last edit: 08/10/2026
#
# Description:
# Exploratory analysis of alpha diversity in scalp microbiome samples,
# including treatment groups, alopecia severity, and patient age.
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


#####################################################################################
# Alpha diversity metrics
#####################################################################################

# Calculate alpha diversity
alpha <- estimate_richness(
  ps,
  measures = c("Observed", "Shannon")
)

# Extract sample metadata
metadata <- data.frame(
  sample_data(ps)
)

# Combine diversity metrics and metadata
alpha <- cbind(
  alpha,
  metadata[rownames(alpha), , drop = FALSE]
)

# Summary statistics
summary(
  alpha[, c("Observed", "Shannon")]
)


#####################################################################################
# Alpha diversity by Treatment
#####################################################################################

# Observed richness
ggplot(
  alpha,
  aes(
    x = Treatment,
    y = Observed,
    fill = Treatment
  )
) +
  geom_boxplot(
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.5
  ) +
  labs(
    x = "Treatment",
    y = "Observed ASV richness"
  ) +
  guides(fill = "none")


# Shannon diversity
ggplot(
  alpha,
  aes(
    x = Treatment,
    y = Shannon,
    fill = Treatment
  )
) +
  geom_boxplot(
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.5
  ) +
  labs(
    x = "Treatment",
    y = "Shannon diversity"
  ) +
  guides(fill = "none")


#####################################################################################
# Alpha diversity by AlopeciaGrade
#####################################################################################

# Observed richness
ggplot(
  alpha,
  aes(
    x = AlopeciaGrade,
    y = Observed,
    fill = AlopeciaGrade
  )
) +
  geom_boxplot(
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.5
  ) +
  labs(
    x = "Alopecia grade",
    y = "Observed ASV richness"
  ) +
  guides(fill = "none")


# Shannon diversity
ggplot(
  alpha,
  aes(
    x = AlopeciaGrade,
    y = Shannon,
    fill = AlopeciaGrade
  )
) +
  geom_boxplot(
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2.5
  ) +
  labs(
    x = "Alopecia grade",
    y = "Shannon diversity"
  ) +
  guides(fill = "none")


#####################################################################################
# Alpha diversity and Age
#####################################################################################

# Shannon diversity vs Age
ggplot(
  alpha,
  aes(
    x = Age,
    y = Shannon,
    color = Treatment
  )
) +
  geom_point(
    size = 3
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  labs(
    x = "Age (years)",
    y = "Shannon diversity"
  )


# Observed richness vs Age
ggplot(
  alpha,
  aes(
    x = Age,
    y = Observed,
    color = Treatment
  )
) +
  geom_point(
    size = 3
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  labs(
    x = "Age (years)",
    y = "Observed ASV richness"
  )
