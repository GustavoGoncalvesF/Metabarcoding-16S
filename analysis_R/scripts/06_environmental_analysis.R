#####################################################################################
# Project: Coffe Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 06_environmental_analysis.R
# Autor: Gustavo Goncalves
# Last edit: 29/07/2026
#
# Description:
# Script to include metadata as pH, altitude, phosphorus and etc in environmental analysis.
#
#####################################################################################

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")

# Set up initial parameters
path <- "results/objects"; list.files(path) # Set path
ps <- readRDS(file.path(path, "phyloseq_object.rds")) # Load phyloseq object
theme_set(theme_bw()) # Set ggplot theme

# Environment data frame
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

# Bray-Curtis (distance matrix)
bray <- phyloseq::distance(
  physeq = ps,
  method = "bray"
)

# PCoA ordinate
ord.pcoa <- ordinate(
  ps,
  method = "PCoA",
  distance = "bray"
)
pcoa_scores <- as.data.frame(ord.pcoa$vectors) # Extract scores from ord.pcoa

# envfit to environmental factors analysis
ef <- envfit(
  pcoa_scores[,1:2],
  env,
  permutations = 999
); ef

# distance-based redundancy analysis
dbrda <- capscale(
  bray ~ pH + Phosphorus + Potassium + CEC + Altitude,
  data = env
); dbrda

anova.global <- anova(dbrda, permutations = 999) # Global ANOVA
anova.margin <- anova(dbrda, by = "margin") # Marginal ANOVA

# VIF
vif <- vif.cca(dbrda)

# R^2 adjust
r2 <- RsquareAdj(dbrda)

# ==================================================
# Environmental analysis report
# ==================================================
sink(file.path(dirname(path), "environmental_analysis_report.txt"))

cat("\n")
cat("=====================================\n")
cat("ENVIRONMENTAL ANALYSIS REPORT\n")
cat("=====================================\n\n")

# Variables
cat("Environmental variables:\n")
cat(paste("-", colnames(env), collapse = "\n"))
cat("\n\n")

# --------------------------------------------------
# envfit
# --------------------------------------------------

cat("ENVFIT\n")
cat("-------------------------------------\n")

ef.table <- data.frame(
  Variable = rownames(ef$vectors$arrows),
  R2 = round(ef$vectors$r, 3),
  Pvalue = ef$vectors$pvals
)

print(ef.table)

cat("\n")

# --------------------------------------------------
# dbRDA
# --------------------------------------------------

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

# --------------------------------------------------
# Global ANOVA
# --------------------------------------------------

cat("GLOBAL MODEL\n")
cat("-------------------------------------\n")

print(anova.global)

cat("\n")

# --------------------------------------------------
# Marginal effects
# --------------------------------------------------

cat("MARGINAL EFFECTS\n")
cat("-------------------------------------\n")

print(anova.margin)

cat("\n")

# --------------------------------------------------
# VIF
# --------------------------------------------------

cat("MULTICOLLINEARITY (VIF)\n")
cat("-------------------------------------\n")

print(round(vif, 3))

cat(
  "\nMaximum VIF:",
  round(max(vif), 3),
  "\n\n"
)

# --------------------------------------------------
# Automatic interpretation
# --------------------------------------------------

cat("INTERPRETATION\n")
cat("-------------------------------------\n")

if(anova.global$`Pr(>F)`[1] < 0.05){
  cat("✓ Global dbRDA model is significant.\n")
}else{
  cat("✗ Global dbRDA model is not significant.\n")
}

if(max(vif) < 5){
  cat("✓ No evidence of problematic multicollinearity.\n")
}else{
  cat("⚠ High multicollinearity detected.\n")
}

cat(
  "✓ Environmental variables explain ",
  round(100 * r2$r.squared, 2),
  "% of the community variation (",
  round(100 * r2$adj.r.squared, 2),
  "% adjusted R²).\n",
  sep = ""
)
sink()
