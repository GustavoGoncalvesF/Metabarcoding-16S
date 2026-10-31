#####################################################################################
# Project: Alopecia Microbiome - 16S rRNA Metabarcoding
# Script: 03_phyloseq_object.R
# Author: Gustavo Goncalves
# Last edit: 08/10/2026
#
# Description:
# Script to integrate ASV counts, taxonomic assignments, and sample metadata
# into a phyloseq object for the alopecia microbiome project.
#
#####################################################################################

# Load project configuration
source("projects/alopecia/config.R")

# Set up packages
library(phyloseq); packageVersion("phyloseq")

# --------------------------------------------------
# Load DADA2 outputs
# --------------------------------------------------

seqtab.nochim <- readRDS(
  file.path(objects_path, "seqtab-nochim.rds")
)

taxa <- readRDS(
  file.path(objects_path, "taxonomy.rds")
)

# Load sample metadata
metadata <- read.csv(
  file.path(project_path, "metadata", "sample_data.csv"),
  row.names = 1,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

# --------------------------------------------------
# Validate sample and ASV identifiers
# --------------------------------------------------

# Check whether sample IDs match
if (!setequal(rownames(metadata), rownames(seqtab.nochim))) {
  stop("Sample IDs in metadata do not match those in the ASV table.")
}

# Reorder metadata to match ASV table
metadata <- metadata[
  rownames(seqtab.nochim),
  ,
  drop = FALSE
]

# Check whether ASV sequences match taxonomy
if (!identical(colnames(seqtab.nochim), rownames(taxa))) {
  stop("ASV sequences in the ASV table do not match the taxonomy table.")
}

# --------------------------------------------------
# Prepare sample metadata
# --------------------------------------------------

# Numeric variables
metadata$Age <- as.numeric(metadata$Age)

# Experimental treatment
metadata$Treatment <- factor(
  metadata$Treatment,
  levels = c("Control", "Active")
)

# Alopecia severity
metadata$AlopeciaGrade <- factor(
  metadata$AlopeciaGrade,
  levels = c(1, 2, 3),
  ordered = TRUE
)

# Gray hair level
metadata$GrayHairLevel <- factor(
  metadata$GrayHairLevel,
  levels = c("None", "Low", "Medium", "High"),
  ordered = TRUE
)

# Hair loss duration
metadata$HairLossDurationMonths <- factor(
  metadata$HairLossDurationMonths,
  levels = c("LessThan3", "From3To12", "MoreThan12"),
  ordered = TRUE
)

# Quality of life
metadata$QualityOfLife <- factor(
  metadata$QualityOfLife,
  levels = c("Poor", "Medium", "Good", "VeryGood"),
  ordered = TRUE
)

# Diet quality
metadata$DietQuality <- factor(
  metadata$DietQuality,
  levels = c("Poor", "Fair", "Good_and_balanced"),
  ordered = TRUE
)

# Binary variables
binary.variables <- c(
  "WeakOrSheddingHair",
  "OilyOrDandruff",
  "Coloring",
  "Bleaching",
  "Straightening",
  "BlowDryerWithBrushing",
  "BlowDryerWithoutBrushing",
  "FlatIron",
  "Stressed",
  "Anxious"
)

for (variable in binary.variables) {
  metadata[[variable]] <- factor(
    metadata[[variable]],
    levels = c("No", "Yes")
  )
}

# Additional categorical variables
metadata$HairEthnicity <- factor(metadata$HairEthnicity)

metadata$CurrentHairCondition <- factor(
  metadata$CurrentHairCondition,
  levels = c("Healthy", "Unhealthy")
)

# --------------------------------------------------
# Validate metadata
# --------------------------------------------------

# Check sample alignment
stopifnot(
  identical(
    rownames(metadata),
    rownames(seqtab.nochim)
  )
)

# Check treatment groups
print(table(metadata$Treatment, useNA = "ifany"))

# Check missing values
print(colSums(is.na(metadata)))

# --------------------------------------------------
# Construct phyloseq object
# --------------------------------------------------

ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(metadata)
)

# Add sample identifiers as metadata
sample_data(ps)$SampleID <- sample_names(ps)

# Add DNA sequences
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

ps <- merge_phyloseq(ps, dna)

# Rename ASVs
taxa_names(ps) <- paste0(
  "ASV",
  seq_len(ntaxa(ps))
)

# --------------------------------------------------
# Save phyloseq object
# --------------------------------------------------

saveRDS(
  ps,
  file.path(objects_path, "phyloseq_object.rds")
)

# Final summary
ps
nsamples(ps)
ntaxa(ps)

