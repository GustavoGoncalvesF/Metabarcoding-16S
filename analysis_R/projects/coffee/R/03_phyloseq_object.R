#####################################################################################
# Project: Coffee Rhizosphere Microbiome - 16S rRNA Metabarcoding
# Script: 03_phyloseq_object.R
# Author: Gustavo Goncalves
# Last edit: 07/10/2026
#
# Description:
# Script to create and adjust the phyloseq object for the coffee rhizosphere
# microbiome project.
#
#####################################################################################

# Load project configuration
source("projects/coffee/config.R")

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

# Load metadata
metadata <- read.csv(
  file.path(project_path, "metadata", "sample_data.csv"),
  sep = ";",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

# Verify sample names
if (!setequal(rownames(metadata), rownames(seqtab.nochim))) {
  stop(
    "Sample IDs in metadata do not match those in the ASV table."
  )
}

# Reorder metadata to match the ASV table
metadata <- metadata[
  rownames(seqtab.nochim),
  ,
  drop = FALSE
]

# Verify ASV sequences
if (!identical(colnames(seqtab.nochim), rownames(taxa))) {
  stop(
    "ASV sequences in the ASV table do not match the taxonomy table."
  )
}

# Convert variable types
metadata$Replicate <- factor(metadata$Replicate)
metadata$Treatment <- factor(metadata$Treatment)
metadata$Farm <- factor(metadata$Farm)
metadata$CoffeeSpecies <- factor(metadata$CoffeeSpecies)
metadata$Cultivar <- factor(metadata$Cultivar)
metadata$Quality <- factor(metadata$Quality)
metadata$CropStage <- factor(metadata$CropStage)
metadata$SampleType <- factor(metadata$SampleType)
metadata$Year <- factor(metadata$Year)

metadata$Altitude <- as.numeric(metadata$Altitude)
metadata$Potassium <- as.numeric(metadata$Potassium)
metadata$Phosphorus <- as.numeric(metadata$Phosphorus)
metadata$pH <- as.numeric(metadata$pH)
metadata$CEC <- as.numeric(metadata$CEC)

# Construct phyloseq object
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxa),
  sample_data(metadata)
)

# Add SampleID column
sample_data(ps)$SampleID <- sample_names(ps)

# Store DNA sequences in the refseq slot
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

ps <- merge_phyloseq(
  ps,
  dna
)

# Rename taxa using short ASV identifiers
taxa_names(ps) <- paste0(
  "ASV",
  seq_len(ntaxa(ps))
)

# Save phyloseq object
saveRDS(
  ps,
  file.path(objects_path, "phyloseq_object.rds")
)
