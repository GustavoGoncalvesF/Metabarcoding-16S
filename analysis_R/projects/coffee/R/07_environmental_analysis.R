#####################################################################################
# Project: Coffee Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 07_environmental_analysis.R
# Author: Gustavo Goncalves
# Last edit: 07/10/2026
#
# Description:
# Script to investigate associations between environmental variables and
# microbial community composition.
#
#####################################################################################

# Load project configuration
source("projects/coffee/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")

# Load phyloseq object
ps <- readRDS(
  file.path(objects_path, "phyloseq_object.rds")
)

# Set ggplot theme
theme_set(theme_bw())


#####################################################################################
# Environmental data
#####################################################################################

env <- data.frame(sample_data(ps))[
  ,
  c(
    "pH",
    "Phosphorus",
    "Potassium",
    "CEC",
    "Altitude"
  )
]


#####################################################################################
# Bray-Curtis dissimilarity
#####################################################################################

bray <- phyloseq::distance(
  physeq = ps,
  method = "bray"
)


#####################################################################################
# PCoA + envfit
#####################################################################################

# PCoA ordination
ord.pcoa <- ordinate(
  ps,
  method = "PCoA",
  distance = bray
)

# Extract PCoA scores
pcoa_scores <- as.data.frame(
  ord.pcoa$vectors
)

# Fit environmental variables onto the PCoA ordination
ef <- envfit(
  pcoa_scores[, 1:2],
  env,
  permutations = 999
)

ef


#####################################################################################
# Distance-based redundancy analysis (dbRDA)
#####################################################################################

dbrda <- capscale(
  bray ~ pH + Phosphorus + Potassium + CEC + Altitude,
  data = env
)

dbrda

# Global significance
anova.global <- anova(
  dbrda,
  permutations = 999
)

# Marginal effects
anova.margin <- anova(
  dbrda,
  by = "margin",
  permutations = 999
)

# Variance inflation factors
vif <- vif.cca(dbrda)

# R-squared
r2 <- RsquareAdj(dbrda)


#####################################################################################
# Environmental analysis report
#####################################################################################

sink(
  file.path(
    results_path,
    "environmental_analysis_report.txt"
  )
)

cat("\n")
cat("=====================================\n")
cat("ENVIRONMENTAL ANALYSIS REPORT\n")
cat("Project:", project_name, "\n")
cat("=====================================\n\n")

# Environmental variables
cat("Environmental variables:\n")
cat(
  paste(
    "-",
    colnames(env),
    collapse = "\n"
  )
)
cat("\n\n")


# envfit
cat("ENVFIT\n")
cat("-------------------------------------\n")

ef.table <- data.frame(
  Variable = rownames(ef$vectors$arrows),
  R2 = round(ef$vectors$r, 3),
  Pvalue = ef$vectors$pvals
)

print(ef.table)

cat("\n")


# dbRDA
cat("dbRDA\n")
cat("-------------------------------------\n")

cat(
  "Variance explained:",
  round(100 * r2$r.squared, 2),
  "%\n"
)

cat(
  "Adjusted R²:",
  round(100 * r2$adj.r.squared, 2),
  "%\n\n"
)


# Global ANOVA
cat("GLOBAL MODEL\n")
cat("-------------------------------------\n")

print(anova.global)

cat("\n")


# Marginal effects
cat("MARGINAL EFFECTS\n")
cat("-------------------------------------\n")

print(anova.margin)

cat("\n")


# VIF
cat("MULTICOLLINEARITY (VIF)\n")
cat("-------------------------------------\n")

print(
  round(vif, 3)
)

cat(
  "\nMaximum VIF:",
  round(max(vif), 3),
  "\n\n"
)


# Automatic interpretation
cat("INTERPRETATION\n")
cat("-------------------------------------\n")

if (anova.global$`Pr(>F)`[1] < 0.05) {
  cat("Global dbRDA model is significant.\n")
} else {
  cat("Global dbRDA model is not significant.\n")
}

if (max(vif) < 5) {
  cat("No evidence of problematic multicollinearity.\n")
} else {
  cat("High multicollinearity detected.\n")
}

cat(
  "Environmental variables explain ",
  round(100 * r2$r.squared, 2),
  "% of the community variation (",
  round(100 * r2$adj.r.squared, 2),
  "% adjusted R²).\n",
  sep = ""
)

sink()


#####################################################################################
# End
#####################################################################################
