#####################################################################################
# Project: Paspalum Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 05_beta_diversity.R
# Author: Gustavo Goncalves
# Last edit: 19/11/2026
#
# Description:
# Exploratory beta diversity analysis using Bray-Curtis, PCoA, PERMANOVA,
# and multivariate dispersion tests. Saves figures and a statistical TXT report.
#####################################################################################

# Load project configuration
source("projects/paspalum/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(vegan); packageVersion("vegan")

# Load phyloseq object
ps <- readRDS(file.path(objects_path, "phyloseq_object.rds"))

# Set ggplot theme and random seed
theme_set(theme_bw())
set.seed(123)

# Create output directory
beta_path <- file.path(results_path, "beta_diversity")
dir.create(beta_path, recursive = TRUE, showWarnings = FALSE)

#####################################################################################
# Prepare data and Bray-Curtis dissimilarity
#####################################################################################

metadata <- data.frame(sample_data(ps))
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
bray <- phyloseq::distance(ps.rel, method = "bray")
metadata <- metadata[attr(bray, "Labels"), , drop = FALSE]

#####################################################################################
# Principal Coordinates Analysis (PCoA)
#####################################################################################

ord.pcoa <- ordinate(ps.rel, method = "PCoA", distance = bray)

p1 <- plot_ordination(ps.rel, ord.pcoa, color = "Origin") +
  geom_point(size = 3)
print(p1)
ggsave(file.path(beta_path, "01_pcoa_origin.png"), p1,
       width = 8, height = 6, dpi = 300)

p2 <- plot_ordination(ps.rel, ord.pcoa, color = "Species") +
  geom_point(size = 3)
print(p2)
ggsave(file.path(beta_path, "02_pcoa_species.png"), p2,
       width = 8, height = 6, dpi = 300)

p3 <- plot_ordination(ps.rel, ord.pcoa, color = "Accession") +
  geom_point(size = 3)
print(p3)
ggsave(file.path(beta_path, "03_pcoa_accession.png"), p3,
       width = 8, height = 6, dpi = 300)

p4 <- plot_ordination(ps.rel, ord.pcoa, color = "Origin") +
  geom_point(size = 3) + facet_wrap(~Species)
print(p4)
ggsave(file.path(beta_path, "04_pcoa_origin_by_species.png"), p4,
       width = 10, height = 6, dpi = 300)

p5 <- plot_ordination(ps.rel, ord.pcoa, color = "Species") +
  geom_point(size = 3) + facet_wrap(~Origin)
print(p5)
ggsave(file.path(beta_path, "05_pcoa_species_by_origin.png"), p5,
       width = 10, height = 6, dpi = 300)

#####################################################################################
# General PERMANOVA
#####################################################################################

adonis.origin <- adonis2(bray ~ Origin, data = metadata, permutations = 999)
adonis.species <- adonis2(bray ~ Species, data = metadata, permutations = 999)
adonis.accession <- adonis2(bray ~ Accession, data = metadata, permutations = 999)

#####################################################################################
# General multivariate dispersion
#####################################################################################

disp.origin <- betadisper(bray, metadata$Origin)
test.disp.origin <- permutest(disp.origin, permutations = 999)

disp.species <- betadisper(bray, metadata$Species)
test.disp.species <- permutest(disp.species, permutations = 999)

disp.accession <- betadisper(bray, metadata$Accession)
test.disp.accession <- permutest(disp.accession, permutations = 999)
tukey.accession <- TukeyHSD(disp.accession)

png(file.path(beta_path, "06_dispersion_origin.png"),
    width = 2400, height = 1800, res = 300)
boxplot(disp.origin, xlab = "Origin", ylab = "Distance to group median")
dev.off()

png(file.path(beta_path, "07_dispersion_species.png"),
    width = 2400, height = 1800, res = 300)
boxplot(disp.species, xlab = "Species", ylab = "Distance to group median")
dev.off()

png(file.path(beta_path, "08_dispersion_accession.png"),
    width = 2400, height = 1800, res = 300)
boxplot(disp.accession, xlab = "Accession", ylab = "Distance to group median")
dev.off()

#####################################################################################
# BAG1: dispersion among accessions and species PERMANOVA
#####################################################################################

ps.bag1 <- subset_samples(ps.rel, Origin == "BAG1")
bray.bag1 <- phyloseq::distance(ps.bag1, method = "bray")
metadata.bag1 <- data.frame(sample_data(ps.bag1))
metadata.bag1 <- metadata.bag1[attr(bray.bag1, "Labels"), , drop = FALSE]

disp.bag1 <- betadisper(bray.bag1, metadata.bag1$Accession)
test.disp.bag1 <- permutest(disp.bag1, permutations = 999)
tukey.bag1 <- TukeyHSD(disp.bag1)

png(file.path(beta_path, "10_dispersion_accession_bag1.png"),
    width = 2400, height = 1800, res = 300)
boxplot(disp.bag1, xlab = "Accession", ylab = "Distance to group median")
dev.off()

ord.bag1 <- ordinate(ps.bag1, method = "PCoA", distance = bray.bag1)
p.species <- plot_ordination(ps.bag1, ord.bag1,
                             color = "Species", shape = "Accession") +
  geom_point(size = 3)
print(p.species)
ggsave(file.path(beta_path, "09_pcoa_species_bag1.png"), p.species,
       width = 8, height = 6, dpi = 300)

adonis.species.bag1 <- adonis2(
  bray.bag1 ~ Species, data = metadata.bag1, permutations = 999
)

#####################################################################################
# Origin adjusted for Accession (balanced subset)
#####################################################################################

# BGP-337 is present only in BAG1, so exclude it for the BAG1 x BAG2 comparison.
ps.balanced <- subset_samples(ps.rel, Accession != "BGP-337")
sample_data(ps.balanced)$Accession <- droplevels(
  sample_data(ps.balanced)$Accession
)

bray.balanced <- phyloseq::distance(ps.balanced, method = "bray")
metadata.balanced <- data.frame(sample_data(ps.balanced))
metadata.balanced <- metadata.balanced[
  attr(bray.balanced, "Labels"), , drop = FALSE
]

# Permutations restricted to within-accession exchanges.
# Only the Origin row is interpreted for this constrained permutation design.
adonis.origin.adjusted <- adonis2(
  bray.balanced ~ Accession + Origin,
  data = metadata.balanced,
  permutations = 999,
  strata = metadata.balanced$Accession,
  by = "margin"
)

disp.origin.balanced <- betadisper(
  bray.balanced, metadata.balanced$Origin
)
test.disp.origin.balanced <- permutest(
  disp.origin.balanced,
  permutations = 999,
  strata = metadata.balanced$Accession
)

png(file.path(beta_path, "11_dispersion_origin_balanced.png"),
    width = 2400, height = 1800, res = 300)
boxplot(disp.origin.balanced, xlab = "Origin",
        ylab = "Distance to group median")
dev.off()

#####################################################################################
# Origin effect within each Accession
#####################################################################################

accessions <- c("BGP-210", "BGP-293", "BGP-486")
adonis.by.accession <- list()

for (acc in accessions) {
  ps.acc <- subset_samples(ps.rel, Accession == acc)
  bray.acc <- phyloseq::distance(ps.acc, method = "bray")
  metadata.acc <- data.frame(sample_data(ps.acc))
  metadata.acc <- metadata.acc[attr(bray.acc, "Labels"), , drop = FALSE]

  adonis.by.accession[[acc]] <- adonis2(
    bray.acc ~ Origin,
    data = metadata.acc,
    permutations = 999
  )
}

#####################################################################################
# Save statistical summary
#####################################################################################

capture.output({
  cat("=====================================\n")
  cat("BETA DIVERSITY SUMMARY\n")
  cat("Project: Paspalum Metabarcoding 16S\n")
  cat("=====================================\n\n")
  cat("Number of samples:", nsamples(ps), "\n")
  cat("Distance: Bray-Curtis (sample relative abundances)\n")
  cat("Permutations requested: 999\n\n")

  cat("-------------------------------------\n")
  cat("GENERAL PERMANOVA\n")
  cat("-------------------------------------\n\n")
  cat("Origin:\n")
  print(adonis.origin)
  cat("\nSpecies:\n")
  print(adonis.species)
  cat("\nAccession:\n")
  print(adonis.accession)

  cat("\n-------------------------------------\n")
  cat("GENERAL MULTIVARIATE DISPERSION\n")
  cat("-------------------------------------\n\n")
  cat("Origin:\n")
  print(test.disp.origin)
  cat("\nSpecies:\n")
  print(test.disp.species)
  cat("\nAccession:\n")
  print(test.disp.accession)
  cat("\nTukey HSD - Accession dispersion:\n")
  print(tukey.accession)

  cat("\n-------------------------------------\n")
  cat("BAG1 ANALYSES\n")
  cat("-------------------------------------\n\n")
  cat("PERMANOVA - Species within BAG1:\n")
  print(adonis.species.bag1)
  cat("\nDispersion - Accessions within BAG1:\n")
  print(test.disp.bag1)
  cat("\nTukey HSD - BAG1 accession dispersion:\n")
  print(tukey.bag1)

  cat("\n-------------------------------------\n")
  cat("ORIGIN ADJUSTED FOR ACCESSION\n")
  cat("-------------------------------------\n\n")
  cat("Balanced subset: BGP-210, BGP-293, BGP-486 (18 samples)\n")
  cat("Permutations restricted within accession.\n")
  cat("Interpret Origin row only for constrained test.\n\n")
  print(adonis.origin.adjusted)

  cat("\nDispersion - Origin (balanced subset):\n")
  print(test.disp.origin.balanced)

  cat("\n-------------------------------------\n")
  cat("ORIGIN EFFECT WITHIN EACH ACCESSION\n")
  cat("-------------------------------------\n\n")
  for (acc in accessions) {
    cat("Accession:", acc, "\n")
    print(adonis.by.accession[[acc]])
    cat("\n")
  }

  cat("-------------------------------------\n")
  cat("INTERPRETATION NOTES\n")
  cat("-------------------------------------\n\n")
  cat("General tests are exploratory and do not mutually adjust for covariates.\n")
  cat("Species and Accession are confounded/nested in this sampling design.\n")
  cat("BAG1 species test is descriptive/exploratory, not independent accession-level replication.\n")
  cat("PERMANOVA is sensitive to group dispersion; inspect betadisper results.\n")
  cat("Within-accession tests have few samples and limited possible permutations.\n")

}, file = file.path(beta_path, "beta_diversity_summary.txt"))

# Also display the key adjusted result in the console
print(adonis.origin.adjusted)
cat("\nSaved statistical summary to:",
    file.path(beta_path, "beta_diversity_summary.txt"), "\n")
